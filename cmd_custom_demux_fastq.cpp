#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "htslib/bgzf.h"
#include "seq_utils.h"
#include "fastq_utils.h"
#include "multiproc_compressor.h"
#include <ctime>
#include <map>
#include <sys/stat.h>
#include <sys/types.h>

int32_t hamming_dist_Nspecial(const char *s1, const char *s2, int32_t len, bool count_N = false)
{
    int32_t dist = 0;
    for (int32_t i = 0; i < len; ++i)
    {
        if (s1[i] != s2[i])
        {
            if (s1[i] == 'N')
            {
                if (count_N)
                    ++dist;
            }
            else
            {
                ++dist;
            }
        }
    }
    return dist;
}

/////////////////////////////////////////////////////////////////////////////////
// custom-demux-fastq : Demultiplex FASTQ files based on custom index file
/////////////////////////////////////////////////////////////////////////////////
int32_t cmdCustomDemuxFASTQ(int32_t argc, char **argv)
{
    std::string fqI1f;
    std::string fqI2f;
    std::string fqR1f;
    std::string fqR2f;
    std::string samplef;
    std::string outprefix;
    std::string outsuffixR1(".R1.fastq.gz");
    std::string outsuffixR2(".R2.fastq.gz");
    std::string outsuffixI1(".I1.fastq.gz");
    std::string outsuffixI2(".I2.fastq.gz");
    std::string keyword_ambiguous("ambiguous");
    std::string cmd_compress("gzip -c");
    int32_t verbose_unit_read = 1000000;
    bool consider_N_as_mismatch = false; // consider N as mismatch
    int32_t max_mismatch = 2;            // maximum mismatch allowed
    int32_t min_diff = 2;                // minimum difference between the best and second best match

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("R1", &fqR1f, "FASTQ file for read 1")
    LONG_STRING_PARAM("R2", &fqR2f, "FASTQ file for read 2")
    LONG_STRING_PARAM("I1", &fqI1f, "FASTQ file for index 1")
    LONG_STRING_PARAM("I2", &fqI2f, "FASTQ file for index 2")
    LONG_STRING_PARAM("sample", &samplef, "Sample sheet file containing [ID] [I1] [I2]")

    LONG_PARAM_GROUP("Settings", NULL)
    LONG_STRING_PARAM("cmd", &cmd_compress, "Command to compress the output files")
    LONG_PARAM("consider-N-as-mismatch", &consider_N_as_mismatch, "Consider N as mismatch (default: false)")
    LONG_INT_PARAM("max-mismatch", &max_mismatch, "Maximum number of mismatch allowed")
    LONG_INT_PARAM("min-diff", &min_diff, "Minimum difference between the best and second best match")
    LONG_INT_PARAM("verbose-chunk", &verbose_unit_read, "Number of records to print output messages (default: 1000000)")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &outprefix, "Output prefix")
    LONG_STRING_PARAM("suffix-R1", &outsuffixR1, "Output suffix for read 1")
    LONG_STRING_PARAM("suffix-R2", &outsuffixR2, "Output suffix for read 2")
    LONG_STRING_PARAM("suffix-I1", &outsuffixI1, "Output suffix for index 1")
    LONG_STRING_PARAM("suffix-I2", &outsuffixI2, "Output suffix for index 2")
    LONG_STRING_PARAM("ambiguous", &keyword_ambiguous, "The keyword for ambiguous samples (default: ambiguous)")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    notice("Analysis started");

    // At the minimum R1 and I1 must be provided
    if (fqR1f.empty() || fqI1f.empty() || outprefix.empty() || samplef.empty() )
        error("Missing required option --R1, --I1, --sample, --out");

    // read the FASTQ files
    htsFile *hR1 = hts_open(fqR1f.c_str(), "r");
    if (hR1 == NULL)
        error("Cannot read %s", fqR1f.c_str());
    notice("Successfully opened %s", fqR1f.c_str());

    htsFile *hR2 = NULL;
    if (!fqR2f.empty())
    {
        hR2 = hts_open(fqR2f.c_str(), "r");
        if (hR2 == NULL)
            error("Cannot read %s", fqR2f.c_str());
        notice("Successfully opened %s", fqR2f.c_str());
    }

    htsFile *hI1 = hts_open(fqI1f.c_str(), "r");
    if (hI1 == NULL)
        error("Cannot read %s", fqI1f.c_str());
    notice("Successfully opened %s", fqI1f.c_str());

    htsFile *hI2 = NULL;
    if (!fqI2f.empty())
    {
        hI2 = hts_open(fqI2f.c_str(), "r");
        if (hI2 == NULL)
            error("Cannot read %s", fqI2f.c_str());
        notice("Successfully opened %s", fqI2f.c_str());
    }

    notice("Reading the input FASTQ files");

    // read the sample sheet
    dataframe_t sample_df(samplef.c_str());
    std::vector<std::string> sample_ids;
    std::vector<std::string> sample_i1s;
    std::vector<std::string> sample_i2s;

    // check if the requird columns exist
    if ( !sample_df.has_column("ID") )
        error("Cannot find ID column in %s", samplef.c_str());
    if ( !sample_df.has_column("I1") )
        error("Cannot find I1 column in %s", samplef.c_str());
    for (int32_t i = 0; i < sample_df.nrows; ++i)
    {
        sample_ids.push_back(sample_df.get_str_elem(i, "ID"));
        sample_i1s.push_back(sample_df.get_str_elem(i, "I1"));
        if (sample_df.has_column("I2"))
        {
            sample_i2s.push_back(sample_df.get_str_elem(i, "I2"));
        }
    }
    std::vector<std::string> sample_prefixes = sample_ids;
    sample_prefixes.push_back(keyword_ambiguous);

    fq_rec_t recR1, recR2, recI1, recI2;
    uint64_t nrecs = 0;

    multiproc_compressor_t out_pipe(cmd_compress.c_str());
    std::vector<int32_t> idx_R1s;
    std::vector<int32_t> idx_R2s;
    std::vector<int32_t> idx_I1s;
    std::vector<int32_t> idx_I2s;

    std::vector<int32_t> mismatches(sample_ids.size(), 0);
    std::vector<uint64_t> nwritten(sample_prefixes.size(), 0);

    notice("Opening the output files for writing");
    // open output files for writing
    for (int32_t i = 0; i < (int32_t)sample_prefixes.size(); ++i)
    {
        std::string outR1f = outprefix + "." + sample_prefixes[i] + outsuffixR1;
        std::string outR2f = outprefix + "." + sample_prefixes[i] + outsuffixR2;
        std::string outI1f = outprefix + "." + sample_prefixes[i] + outsuffixI1;
        std::string outI2f = outprefix + "." + sample_prefixes[i] + outsuffixI2;
        
        out_pipe.add_filename(outR1f.c_str());
        idx_R1s.push_back(out_pipe.size() - 1);
        if ( hR2 != NULL )
        {
            out_pipe.add_filename(outR2f.c_str());
            idx_R2s.push_back(out_pipe.size() - 1);
        }
        out_pipe.add_filename(outI1f.c_str());
        idx_I1s.push_back(out_pipe.size() - 1);
        if ( hI2 != NULL )
        {
            out_pipe.add_filename(outI2f.c_str());
            idx_I2s.push_back(out_pipe.size() - 1);
        }
    }

    // open the output files using pipes
    if ( !out_pipe.open_pipes() )
        error("Cannot open the output files");

    notice("Started reading the input FASTQ files...");
    while (true)
    {
        bool retR1 = hR1 == NULL ? true : read_fastq_record(hR1, recR1);
        bool retR2 = hR2 == NULL ? true : read_fastq_record(hR2, recR2);
        bool retI1 = hI1 == NULL ? true : read_fastq_record(hI1, recI1);
        bool retI2 = hI2 == NULL ? true : read_fastq_record(hI2, recI2);

        if (!(retR1 && retR2 && retI1 && retI2))
        { // EOF reached in one of the files
            // make sure that all non-null files reached EOF
            // each file handle should be either NULL of reached EOF
            if ( !((retR1 == false) && (hR2 == NULL || retR2 == false) && (retI1 == false) && (hI2 == NULL || retI2 == false)) )
            {
                notice("R1 : %s, %s", hR1 == NULL ? "Empty" : "Not Empty", retR1 ? "Not EOF" : "EOF");
                notice("R2 : %s, %s", hR2 == NULL ? "Empty" : "Not Empty", retR2 ? "Not EOF" : "EOF");
                notice("I1 : %s, %s", hI1 == NULL ? "Empty" : "Not Empty", retI1 ? "Not EOF" : "EOF");
                notice("I2 : %s, %s", hI2 == NULL ? "Empty" : "Not Empty", retI2 ? "Not EOF" : "EOF");
                error("FASTQ files %s, %s, %s, and %s have different number of records", fqR1f.c_str(), fqR2f.c_str(), fqI1f.c_str(), fqI2f.c_str());
            }
            break;
        }

        // calculate the hamming distance for each sample
        for (int32_t i = 0; i < (int32_t)sample_ids.size(); ++i)
        {
            int32_t dist1 = hamming_dist_Nspecial(recI1.seq.s, sample_i1s[i].c_str(), recI1.seq.l, consider_N_as_mismatch);
            int32_t dist2 = hI2 == NULL ? 0 : hamming_dist_Nspecial(recI2.seq.s, sample_i2s[i].c_str(), recI2.seq.l, consider_N_as_mismatch);
            mismatches[i] = dist1 + dist2;
        }

        // find the best match
        int32_t imin = 0, inext = -1;
        for (int32_t i = 1; i < sample_ids.size(); ++i)
        {
            if (mismatches[i] < mismatches[imin])   // update imin
            {
                inext = imin;
                imin = i;
            }
            else if (inext == -1 || mismatches[i] < mismatches[inext]) // update inext
            {
                inext = i;
            }
        }

        if ( ( mismatches[imin] > max_mismatch ) || (inext != -1 && (mismatches[inext] - mismatches[imin]) < min_diff) ) {
            imin = sample_prefixes.size() - 1;
        }
        ++nwritten[imin];

        write_fastq_record_fd( out_pipe.get_fd(idx_R1s[imin]), recR1 );
        if (hR2 != NULL)
            write_fastq_record_fd( out_pipe.get_fd(idx_R2s[imin]), recR2 );
        write_fastq_record_fd( out_pipe.get_fd(idx_I1s[imin]), recI1 );
        if (hI2 != NULL)
            write_fastq_record_fd( out_pipe.get_fd(idx_R2s[imin]), recI2 );

        ++nrecs;

        if (nrecs % verbose_unit_read == 0)
            notice("Processing %llu FASTQ records...", nrecs);
    }
    notice("Finished processing %llu FASTQ records...", nrecs);

    // report the number of records per sample
    for (int32_t i = 0; i < (int32_t)sample_prefixes.size(); ++i)
    {
        notice("Sample %s:\t%llu\t(%.3lf %%) records", sample_prefixes[i].c_str(), nwritten[i], nwritten[i]*100.0/nrecs);
    }

    // close the input files
    hts_close(hR1);
    if (hR2 != NULL)
        hts_close(hR2);
    hts_close(hI1);
    if (hI2 != NULL)
        hts_close(hI2);

    // close pipe outputs
    out_pipe.close_pipes();

    notice("Analysis finished");

    return 0;
}
