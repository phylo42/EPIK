# EPIK: Evolutionary Placement with Informative K-mers


[![build](https://github.com/phylo42/EPIK/actions/workflows/build.yml/badge.svg)](https://github.com/phylo42/EPIK/actions/workflows/build.yml)
<a>
<img src="https://img.shields.io/badge/softwipe-7.6-green" />
</a>
        
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

### Install
You can use `epik.py` from the directory where it was built or install it system-wide or for a single user to make `epik.py` visible from any directory.

For a system-wide installation (requires elevated permissions):
```
sudo cmake --install .
```

Alternatively, to install for the current user, choose a directory where you want to install the tool. For instance, you might choose `/home/$USER/opt` or any other directory that you prefer. Replace `DIRECTORY` in the commands below with your chosen directory path:

```
cmake --install . --prefix DIRECTORY
export PATH=DIRECTORY/bin:$PATH
```
Remember to export the `DIRECTORY/bin` to your `PATH`. You can do this manually each time or add the export command to your shell initialization scripts (e.g., `.bashrc`).


## Usage


### Phylogenetic placement
```
epik.py place -i DATABASE -s [nucl|amino] -o OUTPUT_DIR INPUT_FASTA
```
See `epik.py place --help` for more information.

### Building databases

To compute phylo-k-mer databases, use [IPK](https://github.com/phylo42/IPK).


## Other

### Code quality

Code quality evaluation with [softwipe](https://github.com/adrianzap/softwipe) [1]:
```
softwipe --cmake --cpp -x third-party,i2l/third-party,i2l/tests/catch2,i2l/examples --no-execution .
```


## References
[1] A. Zapletal, D. Höhler, C. Sinz, A. Stamatakis (2021) The SoftWipe tool and benchmark for assessing coding standards adherence of scientific software Sci Rep 11, 10015 (2021). https://doi.org/10.1038/s41598-021-89495-8
