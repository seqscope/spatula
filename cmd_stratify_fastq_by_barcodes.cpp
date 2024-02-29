#include "spatula.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "htslib/bgzf.h"
#include "seq_utils.h"
#include "fastq_utils.h"
#include <ctime>
//#include <set>
#include <map>
//#include <unordered_set>
#include <sys/stat.h>
#include <sys/types.h>

/////////////////////////////////////////////////////////////////////////////////
// stratify-fastq-by-barcodes : Stratify FASTQ files by spatial barcodes and UMIs
/////////////////////////////////////////////////////////////////////////////////
int32_t cmdStratifyFASTQByBarcodes(int32_t argc, char **argv)
{
    std::string fq1f;
    std::string fq2f;
    std::string outprefix;
    int32_t startSBCD = 1; // start position of spatial barcode
    int32_t lenSBCD = 27;  // length of spatial barcode to read
    bool fq2SBCD = false;  // spatial barcode is in Read 2
    int32_t startUMI = 31; // start position of UMI
    int32_t lenUMI = 9;    // length of UMI to read
    bool fq2UMI = false;   // UMI is in Read 2
    int32_t nbatch = 512;  // Number of batches to process in memory
    int32_t min_freq = 1;  // Minimum frequency of SBCD-UMI pairs to be considered as a valid pair
    int32_t max_freq = INT32_MAX; // Maximum frequency of SBCD-UMI pairs to be considered as a valid pair

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("fq1", &fq1f, "FASTQ file for read 1")
    LONG_STRING_PARAM("fq2", &fq2f, "FASTQ file for read 2")

    LONG_PARAM_GROUP("Barcode and UMI options", NULL)
    LONG_INT_PARAM("start-sbcd", &startSBCD, "1-based start position of spatial barcode")
    LONG_INT_PARAM("len-sbcd", &lenSBCD, "Length of spatial barcode to read")
    LONG_PARAM("fq2-sbcd", &fq2SBCD, "Spatial barcode is in Read 2")
    LONG_INT_PARAM("start-umi", &startUMI, "1-based start position of UMI")
    LONG_INT_PARAM("len-umi", &lenUMI, "Length of UMI to read")
    LONG_PARAM("fq2-umi", &fq2UMI, "UMI is in Read 2")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &outprefix, "Output prefix")
    LONG_INT_PARAM("nbatch", &nbatch, "Number of batches to process")
    LONG_INT_PARAM("min-freq", &min_freq, "Minimum frequency of SBCD-UMI pairs to be considered as a valid pair")
    LONG_INT_PARAM("max-freq", &max_freq, "Maximum frequency of SBCD-UMI pairs to be considered as a valid pair")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    notice("Analysis started");

    if ( fq1f.empty() || fq2f.empty() )
        error("Missing required option --fq1 and --fq2");

    // read the FASTQ files
    htsFile *hf1 = hts_open(fq1f.c_str(), "r");
    if (hf1 == NULL) 
        error("Cannot read %s", fq1f.c_str());

    htsFile *hf2 = hts_open(fq2f.c_str(), "r");
    if (hf2 == NULL) 
        error("Cannot read %s", fq2f.c_str());

    notice("Reading FASTQ files %s %s", fq1f.c_str(), fq2f.c_str());

    // split the FASTQ files based on spatial barcodes
    fq_rec_t rec1, rec2;
    char buf1[65536], buf2[65536];
    // std::unordered_map<uint128_t, uint32_t, uint128_hash> sbcd_umi_hist;
    uint64_t nrecs = 0;

    // create batch files, which splits the FASTQ files based on the hash of spatial barcodes
    // std::vector<htsFile*> batch_files;
    std::vector<FILE*> batch_files;
    std::vector<std::string> batch_filenames;
    for(int32_t i=0; i < nbatch; ++i) {
        //std::string batch_filename = outprefix + ".tmp.batch" + std::to_string(i) + ".tsv.gz";
        std::string batch_filename = outprefix + ".tmp.batch" + std::to_string(i) + ".bin";
        batch_filenames.push_back(batch_filename);
        //batch_files.push_back(hts_open(batch_filename.c_str(), "wz"));
        batch_files.push_back(fopen(batch_filename.c_str(), "wb"));
        if ( batch_files.back() == NULL )
            error("Cannot open %s for writing", batch_filename.c_str());
    }

    // STEP 1 : read the FASTQ files and split to the batch files
    uint64_t u64buf1[255], u64buf2[255];
    uint64_t verbose_unit_read = 1000000;
    uint64_t verbose_unit_write = 100000;
    while (true)
    {
        bool ret1 = hf1 == NULL ? true : read_fastq_record(hf1, rec1);
        bool ret2 = hf2 == NULL ? true : read_fastq_record(hf2, rec2);

        if ( (ret1 != ret2) && ( hf1 != NULL ) && ( hf2 != NULL ) )
            error("FASTQ files %s and %s have different number of records", fq1f.c_str(), fq2f.c_str());

        if (ret1 == false || ret2 == false) // EOF reached
            break;

        // extract SBCD sequence to buf1
        fq_rec_t *prec = fq2SBCD ? &rec2 : &rec1;
        if ( prec->seq.l < startSBCD + lenSBCD - 1)
            error("Read is too short (length = %d) to extract the SBCD from %d to %d", prec->seq.l, startSBCD, startSBCD + lenSBCD);
        strncpy(buf1, prec->seq.s + startSBCD - 1, lenSBCD);

        // extract UMI sequence to buf2
        prec = fq2UMI ? &rec2 : &rec1;
        if ( prec->seq.l < startUMI + lenUMI - 1 )
            error("Read is too short (length = %d) to extract the UMI from %d to %d", prec->seq.l, startUMI, startUMI + lenUMI);
        strncpy(buf2, prec->seq.s + startUMI - 1, lenUMI);

        uint64_t n1 = seq2nt5(buf1, lenSBCD);
        uint64_t n2 = seq2nt5(buf2, lenUMI);

        // create a batch file based on the hash of SBCD
        int32_t ibatch = (std::hash<uint64_t>()(n1) % nbatch);
        //htsFile* wh = batch_files[ibatch];
        FILE* wh = batch_files[ibatch];

        // write the FASTQ record to the batch file, ignoring the base quality scores
        //hprintf(wh, "%llu\t%llu\t%s\t%s\n", n1, n2, rec1.seq.s, rec2.seq.s);
        int32_t nchunks1 = seq2nt5multi(rec1.seq.s, rec1.seq.l, u64buf1, MAX_NT5_UNIT_64);
        int32_t nchunks2 = seq2nt5multi(rec2.seq.s, rec2.seq.l, u64buf2, MAX_NT5_UNIT_64);
        uint16_t l1 = (uint16_t)rec1.seq.l;
        uint16_t l2 = (uint16_t)rec2.seq.l;

        fwrite(&n1, sizeof(uint64_t), 1, wh);
        fwrite(&n2, sizeof(uint64_t), 1, wh);
        fwrite(&l1, sizeof(uint16_t), 1, wh);
        fwrite(&l2, sizeof(uint16_t), 1, wh);
        fwrite(u64buf1, sizeof(uint64_t), nchunks1, wh);
        fwrite(u64buf2, sizeof(uint64_t), nchunks2, wh);

        ++nrecs;

        if ( nrecs % verbose_unit_read == 0 )
            notice("Processing %llu FASTQ records...", nrecs);
    }
    notice("Finished processing %llu FASTQ records...", nrecs);
    // write UINT64_MAX to the end of the batch files
    for(int32_t i=0; i < nbatch; ++i) {
        FILE* wh = batch_files[i];
        uint64_t u64max = UINT64_MAX;
        fwrite(&u64max, sizeof(uint64_t), 1, wh);
        fclose(wh);
    }

    if ( hf1 != NULL ) hts_close(hf1);
    if ( hf2 != NULL ) hts_close(hf2);

    // close batch file handle
    //for(int32_t i=0; i < nbatch; ++i) {
    //    hts_close(batch_files[i]);
    //}

    notice("Finished reading %llu records and writing temporary files", nrecs);

    // STEP 2 : Read the batch files and store the stratified SBCD-UMI pairs

    uint64_t nrec_pass = 0, nrec_fail = 0, nuniq_pass = 0, nuniq_fail = 0;

    // process each batch file separately
    htsFile* wf = hts_open((outprefix + ".merged.tsv.gz").c_str(), "wz");
    for(int32_t i=0; i < nbatch; ++i) {
        //tsv_reader* tr = new tsv_reader(batch_filenames[i].c_str());
        //if ( tr == NULL )
        //    error("Cannot open %s for reading", batch_filenames[i].c_str());
        FILE* rh = fopen(batch_filenames[i].c_str(), "rb");
        if ( rh == NULL )
            error("Cannot open %s for reading", batch_filenames[i].c_str());

        // stratify the SBCD-UMI pairs based on the hash of SBCD
        std::map<uint64_t, std::map<uint64_t, std::vector<std::string> > > sbcd_umi_str;
        uint64_t n1, n2;
        uint64_t u64buf1[255], u64buf2[255];
        uint16_t l1, l2;

        if ( fread(&n1, sizeof(uint64_t), 1, rh) != 1 )
            error("Cannot read the first record from %s", batch_filenames[i].c_str());
        while (n1 != UINT64_MAX) {
            if ( fread(&n2, sizeof(uint64_t), 1, rh) != 1 )
                error("Cannot read the second record from %s", batch_filenames[i].c_str()); 
            if ( fread(&l1, sizeof(uint16_t), 1, rh) != 1 )
                error("Cannot read the length of SBCD from %s", batch_filenames[i].c_str());
            if ( fread(&l2, sizeof(uint16_t), 1, rh) != 1 )
                error("Cannot read the length of UMI from %s", batch_filenames[i].c_str());

            // read the encoded sequences
            int32_t nchunks1 = (l1 + MAX_NT5_UNIT_64 - 1) / MAX_NT5_UNIT_64;
            int32_t nchunks2 = (l2 + MAX_NT5_UNIT_64 - 1) / MAX_NT5_UNIT_64;
            if ( fread(u64buf1, sizeof(uint64_t), nchunks1, rh) != nchunks1 )
                error("Cannot read the encoded SBCD sequence from %s", batch_filenames[i].c_str());
            if ( fread(u64buf2, sizeof(uint64_t), nchunks2, rh) != nchunks2 )
                error("Cannot read the encoded UMI sequence from %s", batch_filenames[i].c_str());
            
            // convert the encoded sequences to strings
            if ( nt5multi2seq(u64buf1, l1, buf1, MAX_NT5_UNIT_64) == false )
                error("Cannot convert the SBCD to sequence");
            if ( nt5multi2seq(u64buf2, l2, buf2, MAX_NT5_UNIT_64) == false )
                error("Cannot convert the UMI to sequence");
            buf1[l1] = '\t';
            strncpy(buf1 + l1 + 1, buf2, l2);
            buf1[l1 + 1 + l2] = '\0';
            sbcd_umi_str[n1][n2].emplace_back(buf1);

            // read the next record
            fread(&n1, sizeof(uint64_t), 1, rh);
        }

        // while ( tr->read_line() ) {
        //     n1 = tr->uint64_field_at(0);
        //     n2 = tr->uint64_field_at(1);
        //     const char* s1 = tr->str_field_at(2);
        //     const char* s2 = tr->str_field_at(3);
        //     sbcd_umi_str[n1][n2].push_back(std::string(s1) + "\t" + std::string(s2));
        // }

        // write the stratified SBCD-UMI pairs to the output file
        std::map<uint64_t, std::map<uint64_t, std::vector<std::string> > >::iterator it;
        for(it = sbcd_umi_str.begin(); it != sbcd_umi_str.end(); ++it) {
            std::map<uint64_t, std::vector<std::string> >::iterator it2;
            for(it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
                if ( ( it2->second.size() >= min_freq ) && ( it2->second.size() <= max_freq ) ) {
                    std::vector<std::string>& v = it2->second;
                    for(int32_t i=0; i < (int32_t)v.size(); ++i) {
                        hprintf(wf, "%llu\t%llu\t%zu\t%s\n", it->first, it2->first, v.size(), v[i].c_str());
                    }
                    ++nuniq_pass;
                    nrec_pass += v.size();
                }
                else {
                    ++nuniq_fail;
                    nrec_fail += it2->second.size();
                }
                if ( (nuniq_pass + nuniq_fail) % verbose_unit_write == 0 )
                    notice("Writing %llu (%llu pass, %llu failed) SBCD-UMI pairs... in batch %d/%d", nuniq_pass + nuniq_fail, nuniq_pass, nuniq_fail, i+1, nbatch);
            }
        }
        //tr->close();
        //delete tr;
        fclose(rh);
        remove(batch_filenames[i].c_str());
    }
    hts_close(wf);

    notice("Wrote %llu (%llu unique) records and skipped %llu (%llu unique) records", nrec_pass, nuniq_pass, nrec_fail, nuniq_fail);

    notice("Analysis finished");

    return 0;
}
