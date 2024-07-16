# Quickstart for spatula

## Installing spatula

Please follow the instruction below to install `spatula`

```sh
git clone --recursive https://github.com/seqscope/spatula.git
cd spatula
cd submodules
sh -x build.sh
cd ..
mkdir build
cmake ..
make
```

## List available tools

To list the available tools, run the following command:

```sh
../bin/spatula --help
```

If you encounter any difficulties, see [Install](install.md) for more details.