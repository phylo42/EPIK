#! /bin/bash

err_report() {
    echo "Error on line $1"
    exit
}

trap 'err_report $LINENO' ERR


cd xpas 

# Compile xpas
mkdir -p bin && cd bin
cmake -DHASH_MAP="USE_TSL_HOPSCOTCH_MAP" -DCMAKE_CXX_FLAGS="-O3" -DCMAKE_BUILD_TYPE=Release ..
make -j4 xpas_dna xpas_aa

cd ../..

# Compile rappas2
mkdir -p bin && cd bin
cmake -DHASH_MAP="USE_TSL_HOPSCOTCH_MAP" -DCMAKE_CXX_FLAGS="-O3" -DCMAKE_BUILD_TYPE=Release ..
make -j4 rappas2-dna rappas2-aa


