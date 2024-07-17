# Installing spatula

## Installing spatula

Spatula contains `htslib` and `qgenlib` as submodules, so you need to clone the repository recursively to install the required libraries, and build the submodules before building the `spatula` package. An example instruction is given below.

```sh
## STEP 1 : CLONE THE REPOSITORY
## clone the current snapshot of this repository
git clone --recursive https://github.com/seqscope/spatula.git

## move to the spatula directory
cd spatula

## STEP 2 : BUILD THE SUBMODULES
## move to the submodules directory
cd submodules

## build the submodules using build.sh script
sh -x build.sh

## move to the spatula directory
cd ..

## STEP 3 : BUILD SPATULA
## create a build directory
mkdir build
cd build

## Run cmake to configure the build
cmake ..

## Build the spatula package
make
```

If `cmake` is not found, you need to install [cmake](https://cmake.org/) in your system.

## (Optional) Customized specification of the library path

In case any required libraries is missing in `cmake`, you may specify customized installing path by replacing "cmake .." with the following options:

```sh
## If qgenlib is missing or installed in a different directory
$ cmake -DQGEN_INCLUDE_DIRS=/qgenlib_absolute_path/include
        -DHTS_LIBRARIES=/qgenlib_absolute_path/lib/libqgen.a ..

## If htslib is missing or installed in a different directory
$ cmake -DHTS_INCLUDE_DIRS=/htslib_absolute_path/include/  \
        -DHTS_LIBRARIES=/htslib_absolute_path/libhts.a ..

## If bzip2 is missing or installed in a different directory
$ cmake -DBZIP2_INCLUDE_DIRS=/bzip2_absolute_path/include/ \
        -DBZIP2_LIBRARIES=/bzip2_absolute_path/lib/libbz2.a ..

## If lzma is missing or installed in a different directory
$ cmake -DLZMA_INCLUDE_DIRS=/lzma_absolute_path/include/ \
        -DLZMA_LIBRARIES=/lzma_absolute_path/lib/liblzma.a ..

## You may combine the multiples options above if needed.
## Other missing libraries can be handled in a similar way.
```

## Testing spatula

To test whether build was successful, you can run the following command:

```sh
## Current directory: /path/to/install/spatula/build
$ ../bin/spatula --help
```
