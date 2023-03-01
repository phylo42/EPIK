#! /bin/bash

err_report() {
    echo "Error on line $1"
    exit
}

trap 'err_report $LINENO' ERR

mkdir -p bin && cd bin
cmake -DCMAKE_CXX_FLAGS="-O3" -DCMAKE_BUILD_TYPE=Release ..
make -j4 


