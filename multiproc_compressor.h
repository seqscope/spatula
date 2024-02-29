#include <iostream>
#include <vector>
#include <string>
#include <unistd.h>
#include <fcntl.h>
#include <sys/wait.h>
#include "qgenlib/qgen_error.h"
#include "qgenlib/qgen_utils.h"


class multiproc_compressor_t {
private:
    //std::string cmd;
    std::vector<std::string> filenames;
    std::vector<int> pipes;
    char** arglist;
    int32_t argc;

public:
    multiproc_compressor_t(const char* command) {
        std::vector<std::string> toks;
        split(toks, " ", command);
        int32_t ntoks = (int32_t)toks.size();
        arglist = (char **)malloc((ntoks + 1) * sizeof(char *));
        for (int32_t i = 0; i < ntoks; ++i)
        {
            arglist[i] = strdup(toks[i].c_str());
        }
        arglist[ntoks] = NULL;
        argc = ntoks;
    }

    ~multiproc_compressor_t() {
        for (int32_t i = 0; i < argc; ++i)
        {
            free(arglist[i]);
        }
        free(arglist);
    }

    int32_t size() {
        return (int32_t)filenames.size();
    }

    void add_filename(const char* filename) {
        filenames.push_back(filename);
    }

    int32_t get_fd(int32_t i) {
        if (i < 0 || i >= pipes.size()) {
            error("Invalid index %d", i);
            return -1;
        }
        return pipes[i];
    }

    bool open_pipes() {
        for (int32_t i=0; i < filenames.size(); ++i) {
            std::string& filename = filenames[i];
            int pipefd[2];
            if (pipe(pipefd) == -1) {
                error("Cannot create pipe for output file %s", filename.c_str());
                return false;
            }

            pid_t pid = fork();
            if (pid == -1) {
                error("Cannot fork for output file %s", filename.c_str());
                exit(EXIT_FAILURE);
            } else if (pid == 0) { // Child process
                close(pipefd[1]); // Close unused write end

                // Redirect stdout to write end of the pipe
                if (dup2(pipefd[0], STDIN_FILENO) == -1) {
                    error("Cannot redirect stdout to write end of the pipe for output file %s", filename.c_str());
                    return false;
                }

                // Open output file for redirection
                int out_fd = open(filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
                if (out_fd == -1) {
                    error("Cannot open the output file %s for writing", filename.c_str());
                    return false;
                }

                // Redirect standard output to the output file
                if (dup2(out_fd, STDOUT_FILENO) == -1) {
                    error("Cannot redirect stdout to the output file %s", filename.c_str());
                    return false;
                }

                // Execute the command
                //execlp("sh", "sh", "-c", cmd.c_str(), (char *)NULL);
                execvp(arglist[0], arglist);
                error("Cannot execute the command %s", arglist[0]);
                return false;
            } else { // Parent process
                close(pipefd[0]); // Close unused read end
                pipes.push_back(pipefd[1]); // Save write end of the pipe
            }
        }
        notice("Successfully opened %d pipes", (int32_t)pipes.size());
        return true;
    }

    // void write_to_pipe(int idx, const char* s) {
    //     if (idx < 0 || idx >= pipes.size()) {
    //         std::cerr << "Invalid index" << std::endl;
    //         return;
    //     }
    //     write(pipes[idx], s, strlen(s));
    // }

    void close_pipes() {
        for (int32_t i=0; i < pipes.size(); ++i) {
            close(pipes[i]);
        }

        // Wait for all child processes to finish
        while (wait(NULL) > 0);
    }
};
