## build htslib
cd htslib
autoreconf -i
./configure
make
cd ..

## build qgenlib
cd qgenlib
mkdir -p build
cd build
cmake ..
make
cd ../../
