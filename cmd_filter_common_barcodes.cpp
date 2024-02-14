#include "spatula.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "htslib/bgzf.h"
#include "seq_utils.h"
#include "fastq_utils.h"
#include <ctime>
#include <set>
#include <unordered_set>
#include <sys/stat.h>
#include <sys/types.h>

///////////////////////////////////////////////////////////////////////////////////
// filter-common-barcodes : Filter FASTQ file to filter out rare barcode+UMI pairs
///////////////////////////////////////////////////////////////////////////////////
int32_t cmdFilterCommonBarcodes(int32_t argc, char **argv)
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
    int32_t min_freq = 2;  // Minimum frequency of SBCD-UMI pairs to be considered as a valid pair

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
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    notice("Analysis started");

    if ( fq1f.empty() || fq2f.empty() )
        error("Missing required option --fq1 and --fq2");

    htsFile *hf1 = hts_open(fq1f.c_str(), "r");
    if (hf1 == NULL) 
        error("Cannot read %s", fq1f.c_str());

    htsFile *hf2 = hts_open(fq2f.c_str(), "r");
    if (hf2 == NULL) 
        error("Cannot read %s", fq2f.c_str());

    notice("Reading FASTQ files %s %s", fq1f.c_str(), fq2f.c_str());

    fq_rec_t rec1, rec2;
    char buf1[255], buf2[255];
    // std::unordered_map<uint128_t, uint32_t, uint128_hash> sbcd_umi_hist;
    uint64_t nrecs = 0;

    // keep track of the total UMI counts in memory
    std::map<uint64_t, uint64_t> umi_hist;

    // create batch files
    std::vector<BGZF*> batch_files;
    std::vector<std::string> batch_filenames;
    for(int32_t i=0; i < nbatch; ++i) {
        std::string batch_filename = outprefix + ".tmp.batch" + std::to_string(i) + ".gz";
        batch_filenames.push_back(batch_filename);
        batch_files.push_back(bgzf_open(batch_filename.c_str(), "w"));
        if ( batch_files.back() == NULL )
            error("Cannot open %s for writing", batch_filename.c_str());
    }

    // pre-calculate log factorials
    double logfacs[65536];
    logfacs[0] = 0;
    for(int32_t i=1; i < 65536; ++i) {
        logfacs[i] = logfacs[i-1] + log(double(i));
    }

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
        ++umi_hist[n2];

        // create a batch file based on the hash of SBCD
        int32_t ibatch = (std::hash<uint64_t>()(n1) % nbatch);
        if ( bgzf_write(batch_files[ibatch], &n1, sizeof(uint64_t)) != sizeof(uint64_t) )
            error("Cannot write to the batch file %s", batch_filenames[ibatch].c_str());
        if ( bgzf_write(batch_files[ibatch], &n2, sizeof(uint64_t)) != sizeof(uint64_t) )
            error("Cannot write to the batch file %s", batch_filenames[ibatch].c_str());

        //uint128_t key = static_cast<uint128_t>(n1) << 64 | n2;
        //sbcd_umi_hist[key]++;
        ++nrecs;
    }
    if ( hf1 != NULL ) hts_close(hf1);
    if ( hf2 != NULL ) hts_close(hf2);

    // close batch file handle
    for(int32_t i=0; i < nbatch; ++i) {
        bgzf_close(batch_files[i]);
    }

    notice("Finished reading %llu records and write temporary files", nrecs);

    // maintain valid SBCD+UMI pairs
    std::unordered_set<uint128_t, uint128_hash> valid_sbcd_umi;
    uint64_t npass = 0, nfail = 0;

    // process each batch file separately
    for(int32_t i=0; i < nbatch; ++i) {
        BGZF* bf = bgzf_open(batch_filenames[i].c_str(), "r");
        if ( bf == NULL )
            error("Cannot open %s for reading", batch_filenames[i].c_str());

        // calculate the histogram of SBCD+UMI pairs
        std::map<uint64_t, std::map<uint64_t, uint64_t> > sbcd_umi_hist;
        uint64_t n1, n2;
        while ( bgzf_read(bf, &n1, sizeof(uint64_t)) == sizeof(uint64_t) ) {
            if ( bgzf_read(bf, &n2, sizeof(uint64_t)) != sizeof(uint64_t) )
                error("Cannot read from the batch file %s", batch_filenames[i].c_str());
            sbcd_umi_hist[n1][n2]++;
        }

        // calculate the probability of observing the specific counts
        std::map<uint64_t, std::map<uint64_t, uint64_t> >::iterator it;
        for(it = sbcd_umi_hist.begin(); it != sbcd_umi_hist.end(); ++it) {
            std::map<uint64_t, uint64_t>::iterator it2;
            for(it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
                if ( it2->second >= min_freq ) {
                    // valid barcode
                    uint128_t key = static_cast<uint128_t>(it->first) << 64 | it2->first;
                    valid_sbcd_umi.insert(key);
                    ++npass;
                }
                else {
                    ++nfail;
                }
            }

            // uint64_t n_sbcd = 0;
            // std::map<uint64_t, uint64_t>::iterator it2;
            // for(it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
            //     n_sbcd += it2->second;
            // }
            // for(it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
            //     double p_umi = umi_hist[it2->first] / (double)nrecs;
            //     double p_sbcd = n_sbcd / (double) nrecs;

            //     // calculate the probability under the null
            //     double rate = p_umi * p_sbcd;
            //     double lograte = log(rate);
            //     double pdups = 1.0 - rate * exp(-rate) / (1.0-exp(-rate));
            //     double dpois = pow(rate, it2->second) * exp(-rate) / exp(logfacs[it2->second]) / ( 1 - exp(-rate) );
            //     // printf("%llu\t%llu\t%llu\t%.5lg\t%5lg\t%.5lg\t%.5lg\t%.5lg\n", it->first, it2->first, it2->second, p_umi, p_sbcd, rate, pdups, dpois);
            // }
        }
        bgzf_close(bf);
    }

    notice("Identified %llu valid SBCD-UMI pairs and filtered out %llu pairs", npass, nfail);

    notice("Removing temporary files");
    // close batch file handle
    for(int32_t i=0; i < nbatch; ++i) {
        // remove the temporary batch files
        remove(batch_filenames[i].c_str());
    }

    hf1 = hts_open(fq1f.c_str(), "r");
    if (hf1 == NULL) 
        error("Cannot read %s", fq1f.c_str());

    hf2 = hts_open(fq2f.c_str(), "r");
    if (hf2 == NULL) 
        error("Cannot read %s", fq2f.c_str());

    notice("Reading FASTQ files %s %s again and write filtered FASTQs", fq1f.c_str(), fq2f.c_str());

    std::string out1f = outprefix + ".filt.R1.fastq.gz";
    std::string out2f = outprefix + ".filt.R2.fastq.gz";
    htsFile* wf1 = hts_open(out1f.c_str(), "wz");
    if (wf1 == NULL) 
        error("Cannot read %s", fq1f.c_str());

    htsFile* wf2 = hts_open(out2f.c_str(), "wz");
    if (wf2 == NULL) 
        error("Cannot read %s", fq2f.c_str());

    uint64_t nwritten = 0;
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

        uint128_t key = static_cast<uint128_t>(n1) << 64 | n2;

        if ( valid_sbcd_umi.find(key) != valid_sbcd_umi.end() ) {
            // write the record to the output
            write_fastq_record(wf1, rec1);
            write_fastq_record(wf2, rec2);
            ++nwritten;
        }
    }
    notice("Finished writing %llu FASTQ records", nwritten);
    if ( hf1 != NULL ) hts_close(hf1);
    if ( hf2 != NULL ) hts_close(hf2);
    if ( wf1 != NULL ) hts_close(wf1);
    if ( wf2 != NULL ) hts_close(wf2);


    // // calculate non-unique SBCD+UMI pairs
    // std::unordered_map<uint128_t, uint32_t, uint128_hash>::iterator it;
    // uint64_t ndups = 0, sumdups = 0;
    // for(it = sbcd_umi_hist.begin(); it != sbcd_umi_hist.end(); ++it) {
    //     if ( it->second > 1 ) {
    //         ++ndups;
    //         sumdups += it->second;
    //     }
    // }

    // notice("Read %llu records\nFound %zu unique (%.5lf) SBCD-UMI pairs\n%llu (%.5lf) of unique pairs were observed 2+ times,\ncovering %llu (%.5lf) reads", 
    //     nrecs, sbcd_umi_hist.size(), sbcd_umi_hist.size()/(double)nrecs, ndups, ndups/(double)sbcd_umi_hist.size(), sumdups, sumdups/(double)nrecs);

    // Read the FASTQ files once again and write the valid SBCD-UMI pairs


    notice("Analysis finished");

    return 0;
}
