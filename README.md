# spatula : A C++ toolkit for spatial transcriptomics

## Overview 

`spatula` is a collection of C++ tools that are designed to help analyze 
sub-micron resolution
spatial transcriptomics data such as Seq-Scope. 
These tools are under active development, so they may change frequently. 

## Documentation 

A detailed documentation, including installation guide and description of individual tools, can be found at 
[https://seqscope.github.io/spatula/](https://seqscope.github.io/spatula/)

## Installation

You can install `spatula` by following the instructions below:

```bash
## clone the repository
git clone --recursive https://github.com/seqscope/spatula.git
cd spatula

## build the submodules
cd submodules
sh -x build.sh
cd ..

## build spatula
mkdir build
cd build
cmake ..
make

## list available package
../bin/spatula --help
```

## References

* Kang, H. M. (2026) spatula : A C++ toolkit for spatial transcriptomics. Zenodo https://doi.org/10.5281/zenodo.18327577
