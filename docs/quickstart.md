# Quickstart for spatula

## Installing spatula

Before installing `spatula`, you need to install
[htslib](https://github.com/samtools/htslib) and 
[qgenlib](https://github.com/hyunminkang/qgenlib) 
in the **same directory** you
want to install `spatula` (i.e. `spatula`, `qgenlib`, and `htslib` should be
siblings directories). You also need [cmake](https://cmake.org/) installed in your system.

After installing `htslib` and `qgenlib`, you can clone the current snapshot of this repository to install as well

```sh
$ git clone https://github.com/seqscope/spatula.git
$ cd spatula
$ mkdir build
$ cd build
$ cmake ..
$ make

## List the available tools
$ ../bin/spatula --help
```

If you encounter any difficulties, see [Install](install.md) for more details.