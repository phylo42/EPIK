#! /bin/bash

err_report() {
    echo "Error on line $1"
    exit
}

trap 'err_report $LINENO' ERR


# Build RAPPAS.jar
cd rappas
ant -f build-cli.xml dist

cd ../

# Compile rappas2
mkdir -p bin && cd bin
cmake -DHASH_MAP="USE_TSL_HOPSCOTCH_MAP" -DCMAKE_CXX_FLAGS="-O3" ..
make -j4