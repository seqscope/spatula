#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "seq_utils.h"
#include "generic_utils.h"
#include "sge.h"
#include <ctime>
#include <map>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
// merge-pseudobulk : Merge multiple pseudobulk matrices into one
////////////////////////////////////////////////////////////////////////////////
int32_t cmdMergePseudobulk(int32_t argc, char **argv)
{
    std::string listf;
    std::vector<std::string> v_tsvf;
    std::string outf;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input/Output options", NULL)
    LONG_STRING_PARAM("list", &listf, "File containing list of pseudobulk tsv files to merge")
    LONG_MULTI_STRING_PARAM("tsv", &v_tsvf, "Input pseudobulk tsv files to merge (can be specified multiple times)")
    LONG_STRING_PARAM("out", &outf, "Output filename")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    if ( v_tsvf.empty() && listf.empty() )
        error("--tsv or --list must be specified");
    if ( outf.empty() )
        error("--out must be specified");
    if ( ( !v_tsvf.empty()) && ( !listf.empty() ) )
        error("Cannot specify both --tsv and --list at the same time");
    
    notice("Analysis started");

    if ( !listf.empty() ) {
        // read the list of TSV files
        tsv_reader tf(listf.c_str());
        tf.delimiter = '\t';
        while ( tf.read_line() ) {
            if ( tf.nfields < 1 ) {
                error("Empty line in the list file %s", listf.c_str());
            }
            const char* tsvf = tf.str_field_at(0);
            v_tsvf.push_back(tsvf);
        }
    }
    notice("Found %zu pseudobulk files to merge", v_tsvf.size());
    if ( v_tsvf.size() < 2 ) {
        error("At least two pseudobulk files are required to merge");
    }

    notice("Loading pseudobulk matrix: ", v_tsvf[0].c_str());
    pseudobulk_matrix pbm(v_tsvf[0].c_str());
    for (size_t i = 1; i < v_tsvf.size(); ++i) {
        const std::string& tsvf = v_tsvf[i];
        notice("Merging pseudobulk matrix from %s", tsvf.c_str());
        pseudobulk_matrix pbm2(tsvf.c_str());
        if ( !pbm.add(pbm2) ) {
            error("Cannot merge pseudobulk matrix from %s", tsvf.c_str());
        }
    }

    notice("Writing merged pseudobulk matrix to %s", outf.c_str());
    if ( !pbm.write(outf.c_str()) ) {
        error("Cannot write merged pseudobulk matrix to %s", outf.c_str());
    }
    notice("Analysis finished");
    return 0;
}
