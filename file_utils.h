#ifndef __FILE_UTILS_H
#define __FILE_UTILS_H

// copied from https://stackoverflow.com/questions/675039/how-can-i-create-directory-tree-in-c-linux

#include <iostream>
#include <string>
#include <sys/stat.h> // stat
#include <errno.h>    // errno, ENOENT, EEXIST
#if defined(_WIN32)
#include <direct.h>   // _mkdir
#endif

bool isDirExist(const std::string& path);
bool dirExists(const std::string& path);
bool fileExists(const std::string& path);
bool makePath(const std::string& path);
bool removeDir(const std::string& path);

/*
int main(int argc, char* ARGV[])
{
    for (int i=1; i<argc; i++)
    {
        std::cout << "creating " << ARGV[i] << " ... " << (makePath(ARGV[i]) ? "OK" : "failed") << std::endl;
    }
    return 0;
}
*/

#endif
