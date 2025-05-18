#include "spatula.h"
#include "qgenlib/dataframe.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"
#include "sge2.h"
#include "file_utils.h"
#include <cmath>
#include <ctime>
#include <regex>
#include <cstring>
#include <set>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>
#include <deque>
#include <memory>
#include <fstream>
#include "nanoflann/nanoflann.hpp"

template<typename T>
struct PointCloudXtra {
    struct PointXtra {
        T x;
        T y;
        std::vector<std::string> vals;
        PointXtra() : x(0), y(0) {}
        PointXtra(T x, T y) : x(x), y(y) {}
    };

    using coord_t = T;  //!< The type of each coordinate

    std::vector<PointXtra> pts;

    inline size_t kdtree_get_point_count() const { return pts.size(); }
    inline T kdtree_get_pt(const size_t idx, const size_t dim) const
    {
        if (dim == 0)
            return pts[idx].x;
        else 
            return pts[idx].y;
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const
    {
        return false;
    }
};

using kd_tree_f2_xtra_t = nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<float, PointCloudXtra<float>>,
        PointCloudXtra<float>, 2>;

class DecodeJoiner {
public:
    std::string in_mol_tsv;   // TSV file containing individual molecules
    std::string out_prefix;   // Output Prefix
    std::string tmp_dir;      // Temporary directory for intermediate files

    std::vector<std::string> decode_prefixes;
    std::vector<std::string> decode_tsvs;
    std::vector<int32_t> decode_nfields;

    std::string out_suffix_tsv = ".tsv";            // suffix for the output TSV file
    //std::string out_suffix_hist = ".dist.hist.tsv";    // suffix for the histogram of match distance
    //std::string out_suffix_summary = ".summary.tsv";   // suffix for the summary file

    bool mol_colnames_exclude = true;
    std::set<std::string> mol_colnames_set;

    // output column names for tsv files
    std::string colname_mol_X = "X";
    std::string colname_mol_Y = "Y";
    //std::string colname_mol_feature("gene");
    std::string colname_decode_X = "x";
    std::string colname_decode_Y = "y";
    //std::string colname_decode_feature("feature");

    double tile_size = 500;   // tile size in um
    //double bin_size  = 1.0;   // unit to group temporary output
    double max_dist  = 0.5;   // maximum distance between pixel-level output and the molecules 
    int32_t precision = 3;    // Output precision below the decimal point

    int32_t out_max_k = 1;            // maximum num of pixel-level factors to include in the joined output 
    int32_t out_max_p = 1;            // maximum num of pixel-level probs to include in the joined output
    double mu_scale = 1.0;            // scale factor for the mu values

    int32_t n_threads = 1;      // number of threads to use for processing

    void init_decode_prefix_tsvs(std::vector<std::string>& decode_prefix_tsvs) {
        // parse pix_prefix_tsvs
        for(int32_t i=0; i < (int32_t)decode_prefix_tsvs.size(); ++i) {
            std::vector<std::string> toks;
            split(toks, ",", decode_prefix_tsvs[i].c_str());
            if ( toks.size() != 2 ) {
                error("Cannot parse %s", decode_prefix_tsvs[i].c_str());
            }
            decode_prefixes.push_back(toks[0]);
            decode_tsvs.push_back(toks[1]);
            decode_nfields.push_back(-1);
        }
    }

    void init_mol_colnames(std::string& csv_colnames_include, std::string& csv_colnames_exclude) {
        if ( csv_colnames_exclude.empty() ) {
            if ( csv_colnames_include.empty() ) {
                mol_colnames_exclude = true;
            }
            else {
                mol_colnames_exclude = false;
                std::vector<std::string> colnames;
                split(colnames, ",", csv_colnames_include);
                for(int32_t i=0; i < (int32_t)colnames.size(); ++i) {
                    mol_colnames_set.insert(colnames[i]);
                }
            }
        }
        else if ( csv_colnames_include.empty() ) {
            mol_colnames_exclude = true;
            std::vector<std::string> colnames;
            split(colnames, ",", csv_colnames_exclude);
            for(int32_t i=0; i < (int32_t)colnames.size(); ++i) {
                mol_colnames_set.insert(colnames[i]);
            }
        }
        else {
            error("Cannot specify both --csv-colnames-exclude and --csv-colnames-include");
        }    
    }

protected:
    std::vector<std::thread> split_threads;
    std::vector<std::thread> merge_threads;
    std::mutex global_mutex;
    std::map<uint64_t, std::set<int32_t> > tile2idxs; 

    virtual void merge_tile_worker(int32_t threadId,
        uint64_t tile_id,
        const std::set<int32_t>& idxs
    )
    {
        notice("merge_tile_worker (thread %d, tile %llu) started", threadId, tile_id);
        // make sure that index 0 is always the anchor file
        std::string anchor_file;
        std::vector<std::string> tile_files;
        std::vector<bool> tile_exists;
        if ( idxs.find(0) == idxs.end() ) {
            error("Anchor file must appear first but absent in the tile %llu", tile_id);
        }
        anchor_file = tmp_dir + "/files/" + std::to_string(0) + "/" + std::to_string(tile_id) + ".tsv";
        if ( !fileExists(anchor_file) ) {
            error("Cannot find anchor file %s", anchor_file.c_str());
        }
        for(int32_t i=0; i < decode_prefixes.size(); ++i) {
            std::string tile_file = tmp_dir + "/files/" + std::to_string(i+1) + "/" + std::to_string(tile_id) + ".tsv";
            tile_files.push_back(tile_file);
            if ( idxs.find(i+1) == idxs.end() ) {
                tile_exists.push_back(false);
                notice("WARNING: Tile %llu not found in file %d", tile_id, i+1);
                continue;
            }
            else {
                if ( !fileExists(tile_file) ) {
                    error("Cannot find tile file %s", tile_file.c_str());
                }
                tile_exists.push_back(true);
            }
        }

        std::vector<PointCloudXtra<float> > clouds(decode_prefixes.size());
        for(int32_t i=0; i < decode_prefixes.size(); ++i) {
            if ( !tile_exists[i] ) continue;
            PointCloudXtra<float>& cloud = clouds[i];
            std::string tile_file = tile_files[i];
            tsv_reader tr(tile_file.c_str());
            while( tr.read_line() ) {
                float x = (float)tr.double_field_at(0);
                float y = (float)tr.double_field_at(1);
                cloud.pts.emplace_back(x, y);
                PointCloudXtra<float>::PointXtra& pt = cloud.pts.back();
                for(int32_t j=2; j < tr.nfields; ++j) {
                    pt.vals.push_back(tr.str_field_at(j));
                }
            }
            notice("Loaded %zu points from file %d", cloud.pts.size(), cloud.pts.size());
        }

        std::vector<std::unique_ptr<kd_tree_f2_xtra_t> > kd_trees;
        for(int32_t i=0; i < decode_prefixes.size(); ++i) {
            if ( !tile_exists[i] ) {
                kd_trees.push_back(nullptr);
                continue;
            }
            notice("Building kd-tree for file %d", i+1);
            notice("clouds[%d].pts.size() = %zu", i, clouds[i].pts.size());
            kd_trees.push_back(std::unique_ptr<kd_tree_f2_xtra_t>(new kd_tree_f2_xtra_t(2, clouds[i], {10})));
            kd_trees.back()->buildIndex();
            notice("Finished building kd-tree for file %d", i+1);
        }

        // read the anchor file and find the nearest neighbors
        // write the output line for the tile
        std::string out_file = tmp_dir + "/tiles/" + std::to_string(tile_id) + ".tsv";
        FILE* wf = fopen(out_file.c_str(), "w");
        tsv_reader tr(anchor_file.c_str());
        float query_pt[2];
        const size_t num_neighbors = 1;
        uint64_t nlines = 0;
        while( tr.read_line() ) {
            float x = (float)tr.double_field_at(0);
            float y = (float)tr.double_field_at(1);
            query_pt[0] = x;
            query_pt[1] = y;
            //notice("Anchor: (%f, %f), precision = %d", x, y, precision);
            // write the output line to the tile
            fprintf(wf, "%.*f\t%.*f", precision, x, precision, y);
            // write additional columns
            for(int32_t i=2; i < tr.nfields; ++i) {
                fprintf(wf, "\t%s", tr.str_field_at(i));
            }
            for(int32_t i=0; i < decode_prefixes.size(); ++i) {
                if ( !tile_exists[i] ) {
                    for(int32_t j=2; j < decode_nfields[i]; ++j) {
                        fprintf(wf, "\tNA");
                    }
                    continue;
                }
                std::vector<uint32_t> indices(1);
                std::vector<float> dists(1);
                //notice("foo");
                size_t num_found = kd_trees[i]->knnSearch(query_pt, num_neighbors, &indices[0], &dists[0]);
                //notice("bar");

                if ( num_found > 0 ) {
                    //notice("Found %zu points (%u) in file %d with dist = %f", num_found, indices[0], i+1, dists[0]);
                    PointCloudXtra<float>::PointXtra& pt = clouds[i].pts[indices[0]];
                    // write the output line to the tile
                    if ( dists[0] <= (float)(max_dist * max_dist) ) { // closest point within the max distance
                        // write the output line to the tile
                        if ( pt.vals.size() != decode_nfields[i] - 2 ) {
                            error("Mismatch in the number of fields in file %d", i+1);
                        }
                        for(int32_t j=0; j < (int32_t)pt.vals.size(); ++j) {
                            fprintf(wf, "\t%s", pt.vals[j].c_str());
                        }
                    }
                    else {
                        for(int32_t j=2; j < decode_nfields[i]; ++j) {
                            fprintf(wf, "\tNA");
                        }
                    }
                    //notice("File %d: (%f, %f)", i+1, pt.x, pt.y);
                }
                else {
                    for(int32_t j=2; j < decode_nfields[i]; ++j) {
                        fprintf(wf, "\tNA");
                    }
                }
            }
            fprintf(wf, "\n");
            ++nlines;
            //if ( nlines > 10 ) break; // for testing
            if ( nlines % 1000000 == 0 ) {
                notice("Thread %d, Tile %llu, Processed %llu lines", threadId, tile_id, nlines);
            }
        }
        fclose(wf);

        // remove the kd-trees
        for(int32_t i=0; i < decode_prefixes.size(); ++i) {
            if ( !tile_exists[i] ) continue;
            kd_trees[i].reset();
        }
        // remove the tile files
        for(int32_t i=0; i < decode_prefixes.size(); ++i) {
            if ( !tile_exists[i] ) continue;
            std::string tile_file = tile_files[i];
            if ( std::remove(tile_file.c_str()) != 0 ) {
                error("Cannot remove %s", tile_file.c_str());
            }
        }
        // remove the anchor file
        if ( std::remove(anchor_file.c_str()) != 0 ) {
            error("Cannot remove %s", anchor_file.c_str());
        }
        notice("merge_tile_worker (thread %d, tile %llu) finished processing %llu lines", threadId, tile_id, nlines);
    }

    virtual void split_file_worker(int32_t threadId,  
        int32_t idx_file, 
        const std::string& in_file, 
        const std::string& colname_x, 
        const std::string& colname_y, 
        const std::set<std::string>& colnames_set,
        bool colnames_exclude,
        const std::string& decode_prefix
    ) 
    {
        notice("split_file_worker (thread %d, file %d) started", threadId, idx_file);
        // read the file
        tsv_reader tr(in_file.c_str());

        // create a new path for the input file        
        std::string tmpDir = tmp_dir + "/files/" + std::to_string(idx_file);
        if ( !makePath(tmpDir) ) {
            error("Cannot create directory %s", tmpDir.c_str());
        }

        // read the header line and identify column indices to extract
        int32_t idx_x = -1, idx_y = -1;
        std::vector<int32_t> idxs_include;
        std::vector<std::string> colnames_include;

        // read the first line
        if ( !tr.read_line() ) {
            error("Cannot read the first line of %s", in_file.c_str());
        }
        // parse the header line
        std::map<std::string, int32_t> col2idx;
        for(int32_t i=0; i < tr.nfields; ++i) {
            const char* s = tr.str_field_at(i);
            if ( s[0] == '#' ) ++s;
            if ( col2idx.find(s) != col2idx.end() ) {
                error("Duplicate column name %s in %s", s, in_file.c_str());
            }
            col2idx[s] = i;
        }

        // determine the column indices to include
        if ( col2idx.find(colname_x) == col2idx.end() ) {
            error("Cannot find column %s in the header of %s", colname_x.c_str(), in_file.c_str());
        }
        idx_x = col2idx[colname_x];
        if ( col2idx.find(colname_y) == col2idx.end() ) {
            error("Cannot find column %s in the header of %s", colname_y.c_str(), in_file.c_str());
        }
        idx_y = col2idx[colname_y];

        for(int32_t i=0; i < tr.nfields; ++i) {
            if ( i == idx_x || i == idx_y ) continue; // skip the x and y columns
            const char* s = tr.str_field_at(i);
            if ( s[0] == '#' ) ++s;
            if ( colnames_set.find(s) != colnames_set.end() ) { // found the column name from the set
                // add the column if colnames_exclude is false
                if ( !colnames_exclude ) {
                    idxs_include.push_back(i);
                    colnames_include.push_back(s);
                }
            }
            else if ( colnames_exclude )  {
                idxs_include.push_back(i);
                colnames_include.push_back(s);
            }
        } 

        if ( idx_file > 0 ) {
            decode_nfields[idx_file-1] = colnames_include.size() + 2;
        }
        // write the header file
        std::string hdr_file = tmpDir + "/hdr.tsv";
        FILE* fp_hdr = fopen(hdr_file.c_str(), "w");
        if ( fp_hdr == NULL ) {
            error("Cannot open %s for writing", hdr_file.c_str());
        }
        fprintf(fp_hdr, "%s\t%s", colname_mol_X.c_str(), colname_mol_Y.c_str());
        for(int32_t i=0; i < (int32_t)colnames_include.size(); ++i) {
            fprintf(fp_hdr, "\t%s%s", decode_prefix.c_str(), colnames_include[i].c_str());
        }
        fprintf(fp_hdr, "\n");
        fclose(fp_hdr);

        std::map<uint64_t, FILE*> tile2files;
        std::map<uint64_t, FILE*>::iterator it_tile2files;
        // read the data lines
        int32_t nlines = 0;
        while( tr.read_line() ) {
            double x = tr.double_field_at(idx_x);
            double y = tr.double_field_at(idx_y);
            int32_t xbin = (int32_t)(x / tile_size);
            int32_t ybin = (int32_t)(y / tile_size);
            uint64_t tile_id = ((uint64_t)xbin << 32) | (uint64_t)ybin;

            it_tile2files = tile2files.find(tile_id);
            FILE* fp = NULL;
            if ( it_tile2files == tile2files.end() ) { // create a new file
                std::string tmpFile = tmpDir + "/" + std::to_string(tile_id) + ".tsv";
                fp = fopen(tmpFile.c_str(), "w");
                if ( fp == NULL ) {
                    error("Cannot open %s for writing", tmpFile.c_str());
                }
                tile2files[tile_id] = fp;
            }
            else { // use the existing file
                fp = it_tile2files->second;
            }

            // write the output line to the tile
            fprintf(fp, "%s\t%s", tr.str_field_at(idx_x), tr.str_field_at(idx_y));
            // write additional columns
            for(int32_t i=0; i < (int32_t)idxs_include.size(); ++i) {
                fprintf(fp, "\t%s", tr.str_field_at(idxs_include[i]));
            }
            fprintf(fp, "\n");
            ++nlines;
            if ( nlines % 1000000 == 0 ) {
                notice("split_file_worker thread %d: %d lines processed", threadId, nlines);
            }
        }
        // close the files
        for(it_tile2files = tile2files.begin(); it_tile2files != tile2files.end(); ++it_tile2files) {
            fclose(it_tile2files->second);
        }

        // update the global tile2idxs map
        {
            std::lock_guard<std::mutex> lock(global_mutex);
            for(it_tile2files = tile2files.begin(); it_tile2files != tile2files.end(); ++it_tile2files) {
                uint64_t tile_id = it_tile2files->first;
                tile2idxs[tile_id].insert(idx_file);
            }
        }

        notice("SplitThread %d: Finished processing %d lines", threadId, nlines);
        notice("split_file_worker (thread %d, file %d) started", threadId, idx_file);
    }

    bool merge_temp_files() {
        notice("Merging temporary files to %s", (out_prefix + out_suffix_tsv).c_str());
        std::ofstream outfile(out_prefix + out_suffix_tsv);
        if (!outfile.is_open()) {
            error("Cannot open %s for writing", (out_prefix + out_suffix_tsv).c_str());
        }
        // merge the headers first
        for(int32_t i=0; i <= (int32_t)decode_prefixes.size(); ++i) {
            std::string hdr_file = tmp_dir + "/files/" + std::to_string(i) + "/hdr.tsv";
            tsv_reader tr(hdr_file.c_str());
            if ( !tr.read_line() ) {
                error("Cannot read the header line of %s", (tmp_dir + "/files/" + std::to_string(i) + "/hdr.tsv").c_str());
            }
            if ( i == 0 ) {
                outfile << tr.str_field_at(0) << "\t" << tr.str_field_at(1);
            }
            for(int32_t j=2; j < tr.nfields; ++j) {
                outfile << "\t" << tr.str_field_at(j);
            }
            std::remove(hdr_file.c_str());
        }
        outfile << std::endl;

        // merge each tile
        std::map<uint64_t, std::set<int32_t> >::iterator it_tile2idxs;
        for(it_tile2idxs = tile2idxs.begin(); it_tile2idxs != tile2idxs.end(); ++it_tile2idxs) {
            uint64_t tile_id = it_tile2idxs->first;
            notice("Merging tile %llu", tile_id);
            std::string tile_file = tmp_dir + "/tiles/" + std::to_string(tile_id) + ".tsv";
            std::ifstream tmpFile(tile_file, std::ios::binary);
            if ( tmpFile ) {
                outfile << tmpFile.rdbuf();
                tmpFile.close();
                std::remove(tile_file.c_str());
            }
            else {
                error("Cannot open %s for reading", tile_file.c_str());
            }
        }
        notice("Finished merging temporary files to %s", (out_prefix + out_suffix_tsv).c_str());
        return true;
    }

    bool launch_merge_worker_threads() {
        std::map<uint64_t, std::set<int32_t> >::iterator it_tile2idxs;

        merge_threads.clear();
        merge_threads.reserve(tile2idxs.size());

        std::deque<size_t> running_thread_indices;

        notice("Launching up to %d merge worker threads concurrently.", n_threads);

        int32_t thread_id = 0;
        for(it_tile2idxs = tile2idxs.begin(); it_tile2idxs != tile2idxs.end(); ++it_tile2idxs) {
            if ( running_thread_indices.size() >= n_threads ) {
                size_t oldest_thread_vector_idx = running_thread_indices.front();

                if (merge_threads[oldest_thread_vector_idx].joinable()) {
                    //notice("Reached thread limit (%d). Waiting for oldest thread (vector index %zu) to finish...", n_threads, oldest_thread_vector_idx);
                    merge_threads[oldest_thread_vector_idx].join(); // Wait for completion
                    //notice("Thread (vector index %zu) finished.", oldest_thread_vector_idx);
                }
                // Remove the index of the completed thread from the front of the deque
                running_thread_indices.pop_front();
            }
            uint64_t tile_id = it_tile2idxs->first;
            std::set<int32_t> idx_files = it_tile2idxs->second;

            merge_threads.emplace_back([this, thread_id, tile_id, idx_files]() {
                this->merge_tile_worker(thread_id, tile_id, idx_files);
            });
            running_thread_indices.push_back(merge_threads.size() - 1);
            ++thread_id;
        }
        return true;
    }

    bool launch_split_worker_threads() {
        // set the parameters for the split worker threads
        std::vector<std::string> in_files;
        std::vector<std::string> colname_xs;
        std::vector<std::string> colname_ys;
        std::vector<std::set<std::string> > colnames_sets;
        std::vector<bool> colnames_excludes;
        std::vector<std::string> append_prefixes;
        in_files.push_back(in_mol_tsv);
        colname_xs.push_back(colname_mol_X);
        colname_ys.push_back(colname_mol_Y);
        colnames_sets.push_back(mol_colnames_set);
        colnames_excludes.push_back(mol_colnames_exclude);
        append_prefixes.push_back("");

        for(int32_t i=0; i < (int32_t)decode_prefixes.size(); ++i) {
            in_files.push_back(decode_tsvs[i]);
            colname_xs.push_back(colname_decode_X);
            colname_ys.push_back(colname_decode_Y);
            std::set<std::string> colnames_set;
            for(int32_t i=0; i < out_max_k; ++i) {
                std::string colname = "K" + std::to_string(i+1);
                colnames_set.insert(colname);
            }
            for(int32_t i=0; i < out_max_p; ++i) {
                std::string colname = "P" + std::to_string(i+1);
                colnames_set.insert(colname);
            }
            colnames_sets.push_back(colnames_set);
            colnames_excludes.push_back(false);
            append_prefixes.push_back(decode_prefixes[i]);
        }

        split_threads.clear();
        split_threads.reserve(in_files.size());

        std::deque<size_t> running_thread_indices;

        notice("Launching up to %d split worker threads concurrently.", n_threads);
    
        for(int32_t i=0; i < (int32_t)in_files.size(); ++i) {
            if (running_thread_indices.size() >= n_threads) {
                size_t oldest_thread_vector_idx = running_thread_indices.front();

                if (split_threads[oldest_thread_vector_idx].joinable()) {
                    //notice("Reached thread limit (%d). Waiting for oldest thread (vector index %zu) to finish...", n_threads, oldest_thread_vector_idx);
                    split_threads[oldest_thread_vector_idx].join(); // Wait for completion
                    //notice("Thread (vector index %zu) finished.", oldest_thread_vector_idx);
                }
                // Remove the index of the completed thread from the front of the deque
                running_thread_indices.pop_front();
            }
            int32_t file_idx = static_cast<int32_t>(i); // Match function signature if needed, or keep as size_t
            std::string file = in_files[i];
            std::string col_x = colname_xs[i];
            std::string col_y = colname_ys[i];
            std::set<std::string> colnames = colnames_sets[i]; // Explicit copy
            bool exclude = colnames_excludes[i];     
            std::string append_prefix = append_prefixes[i];      

            split_threads.emplace_back([this, file_idx, file, col_x, col_y, colnames, exclude, append_prefix]() {
                // Call the member function using the captured values
                // Note: threadId and idx_file in the original worker might be the same, adjust if needed.
                // Assuming threadId can just be file_idx here.
                this->split_file_worker(file_idx, file_idx, file, col_x, col_y, colnames, exclude, append_prefix);
            });

            running_thread_indices.push_back(split_threads.size() - 1);
        }

        notice("All %zu split worker tasks have been launched.", in_files.size());
        return true;
    }

    bool join_split_worker_threads() {
        notice("Waiting for %zu launched split worker threads to complete...", split_threads.size());
        size_t joined_count = 0;
        for (size_t i = 0; i < split_threads.size(); ++i) {
            if (split_threads[i].joinable()) {
                split_threads[i].join();
                joined_count++;
            }
        }
        notice("All %zu split worker threads finished.", joined_count);
        split_threads.clear(); // Optional: Clear the vector after joining
        return true;
    }

    bool join_merge_worker_threads() {
        notice("Waiting for %zu launched merge worker threads to complete...", merge_threads.size());
        size_t joined_count = 0;
        for (size_t i = 0; i < merge_threads.size(); ++i) {
            if (merge_threads[i].joinable()) {
                merge_threads[i].join();
                joined_count++;
            }
        }
        notice("All %zu merge worker threads finished.", joined_count);
        merge_threads.clear(); // Optional: Clear the vector after joining
        return true;
    }
public:
    bool run() {
        notice("Launching %d threads to split the input files", n_threads);
        if ( !launch_split_worker_threads() ) {
            error("Error in launching split worker threads");
            return false;
        }
        if ( !join_split_worker_threads() ) {
            error("Error in joining split worker threads");
            return false;
        }
        notice("Merging individual tiles with spatial joining");
        if ( !launch_merge_worker_threads() ) {
            error("Error in launching merge worker threads");
            return false;
        }
        if ( !join_merge_worker_threads() ) {
            error("Error in joining merge worker threads");
            return false;
        }
        notice("Merging temporary files and writing the index");
        if ( !merge_temp_files() ) {
            error("Error in merging temporary files");
            return false;
        }
        return true;
    }
};


/////////////////////////////////////////////////////////////////////////////////////////
// join-pixel-decode : Join FICTURE's pixel-level output with raw transcript-level TSV files
/////////////////////////////////////////////////////////////////////////////////////////
int32_t cmdJoinPixelDecode(int32_t argc, char **argv)
{
    DecodeJoiner dj;
 
    std::vector<std::string> decode_prefix_tsvs; // vector of "[prefix],[tsv-path]" pairs for pixel-level projections
    std::string csv_colnames_include; // comma-separated column names to include in the output TSV file
    std::string csv_colnames_exclude; // comma-separated column names to exclude in the output TSV file

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Key Input/Output Options", NULL)
    LONG_STRING_PARAM("mol-tsv", &dj.in_mol_tsv, "TSV file containing individual molecules")
    LONG_MULTI_STRING_PARAM("decode-prefix-tsv", &decode_prefix_tsvs, "TSV file containing pixel-level factors in [prefix],[tsv-path] format")
    LONG_STRING_PARAM("out-prefix", &dj.out_prefix, "Output prefix for the joined TSV files")
    LONG_STRING_PARAM("tmp-dir", &dj.tmp_dir, "Temporary directory for intermediate files")

    LONG_PARAM_GROUP("Key Parameters", NULL)
    LONG_DOUBLE_PARAM("tile-size", &dj.tile_size, "Tile size to create temporary files for binning")
    //LONG_DOUBLE_PARAM("bin-size", &dj.bin_size, "Bin size for grouping the pixel-level output for indexing")
    LONG_DOUBLE_PARAM("max-dist", &dj.max_dist, "Maximum distance in um to consider a match")
    LONG_DOUBLE_PARAM("mu-scale", &dj.mu_scale, "Scale factor for the resolution - divide by mu_scale in the output")
    LONG_DOUBLE_PARAM("precision", &dj.precision, "Output precision below the decimal point")
    LONG_INT_PARAM("threads", &dj.n_threads, "Number of threads to use for processing")

    LONG_PARAM_GROUP("Expected columns in input and output", NULL)
    LONG_STRING_PARAM("colname-mol-x", &dj.colname_mol_X, "Column name for X-axis for molecular TSV")
    LONG_STRING_PARAM("colname-mol-y", &dj.colname_mol_Y, "Column name for Y-axis for molecular TSV")
    LONG_STRING_PARAM("colname-decode-x", &dj.colname_mol_X, "Column name for X-axis for decoded TSV")
    LONG_STRING_PARAM("colname-decode-y", &dj.colname_mol_Y, "Column name for Y-axis for decoded TSV")
    LONG_STRING_PARAM("colnames-include", &csv_colnames_include, "Comma-separated column names to include in the output TSV file")
    LONG_STRING_PARAM("colnames-exclude", &csv_colnames_exclude, "Comma-separated column names to exclude in the output TSV file")
    LONG_INT_PARAM("out-max-k", &dj.out_max_k, "Maximum number of pixel-level factors to include in the joined output. (Default : 1)")
    LONG_INT_PARAM("out-max-p", &dj.out_max_p, "Maximum number of pixel-level posterior probabilities to include in the joined output. (Default : 1)")

    LONG_PARAM_GROUP("Output File suffixes", NULL)
    LONG_STRING_PARAM("out-suffix-tsv", &dj.out_suffix_tsv, "Suffix for the output TSV file")
    //LONG_STRING_PARAM("out-suffix-hist", &dj.out_suffix_hist, "Suffix for the histogram of match distance")
    //LONG_STRING_PARAM("out-suffix-summary", &dj.out_suffix_summary, "Suffix for the summary file")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    // check the input files
    if ( dj.in_mol_tsv.empty() || decode_prefix_tsvs.empty() || dj.out_prefix.empty() )
        error("--in-mol-tsv, --pix-prefix-tsv, and --out-prefix must be specified");

    dj.init_decode_prefix_tsvs(decode_prefix_tsvs);
    dj.init_mol_colnames(csv_colnames_include, csv_colnames_exclude);

    if ( dj.tmp_dir.empty() ) {
        dj.tmp_dir = dj.out_prefix + ".tmp";
    }
    makePath(dj.tmp_dir);
    makePath(dj.tmp_dir + "/files");
    makePath(dj.tmp_dir + "/tiles");

    uint32_t max_threads = std::thread::hardware_concurrency();
    if ( dj.n_threads > max_threads )
        dj.n_threads = max_threads;
    notice("Using %d threads for processing", dj.n_threads);

    notice("Analysis Started");

    if ( !dj.run() ) {
        error("Failed to run the join-pixel-decode command");
        exit(1);
    }

    // remove the temporary directory
    for(int32_t i=0; i <= (int32_t)decode_prefix_tsvs.size(); ++i) {
        std::string tmpDir = dj.tmp_dir + "/files/" + std::to_string(i);
        if ( !removeDir(tmpDir) ) {
            error("Cannot remove %s", tmpDir.c_str());
        }
    }
    if ( !removeDir(dj.tmp_dir + "/files") ) {
        error("Cannot remove %s", (dj.tmp_dir + "/files").c_str());
    }
    if ( ! removeDir(dj.tmp_dir + "/tiles") ) {
        error("Cannot remove %s", (dj.tmp_dir + "/tiles").c_str());
    }
    if ( !removeDir(dj.tmp_dir) ) {
        error("Cannot remove %s", dj.tmp_dir.c_str());
    }

    notice("Analysis Finished");

    return 0;
}
