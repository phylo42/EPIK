# EPIK: Evolutionary Placement with Informative K-mers
EPIK is a program for rapid alignment-free phylogenetic placement, the successor of [RAPPAS](https://github.com/phylo42/RAPPAS).

## Installation

### Prerequisites

- Boost Libraries >=1.6
- CMake >= 3.10
- GCC compiler must support c++17
- zlib
- rapidjson
- click

On Debian-like systems they can be installed with:
```
sudo apt install build-essential cmake libboost-dev libboost-serialization-dev libboost-filesystem-dev libboost-iostreams-dev libboost-program-options-dev zlib1g-dev rapidjson-dev libquadmath0 python3-pip
pip install click
```

### Clone and build
```
git clone --recursive https://github.com/phylo42/EPIK epik
cd epik && mkdir -p bin && cd bin
cmake ..
make -j4
```

## Usage


### Phylogenetic placement
```
python epik.py place -i DATABASE -s [nucl|amino] -o OUTPUT_DIR INPUT_FASTA
```
See `python epik.py place --help` for more information.

### Building databases

To compute phylo-k-mer databases, use [IPK](https://github.com/phylo42/IPK).


## Other

### Code quality

EPIK scores 7.6/10 in the code quality evaluation according to [softwipe](https://github.com/adrianzap/softwipe) [1]. To run it:

```
softwipe --cmake --cpp -x third-party,i2l/third-party,i2l/tests/catch2,i2l/examples --no-execution .
```


## References
[1] A. Zapletal, D. Höhler, C. Sinz, A. Stamatakis (2021) The SoftWipe tool and benchmark for assessing coding standards adherence of scientific software Sci Rep 11, 10015 (2021). https://doi.org/10.1038/s41598-021-89495-8
