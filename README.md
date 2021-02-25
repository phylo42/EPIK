# RAPPAS2
This is the new version of [RAPPAS](https://github.com/phylo42/RAPPAS), currently under development.

## Installation

### Prerequisites

- Boost Libraries >=1.6
- CMake >= 3.10
- GCC compiler must support c++17
- zlib
- rapidjson

In debian, these can be installed with:
```
sudo apt install build-essential cmake libboost-dev libboost-serialization-dev libboost-filesystem-dev libboost-iostreams-dev libboost-program-options-dev zlib1g-dbg rapidjson-dev libquadmath0
```

### Clone and build
```
git clone --recursive https://github.com/phylo42/rappas2.git
cd rappas2
./build.sh
```

## Usage

### Building databases

The functionality of constructing new databases of phylo k-mers has been moved to [xpas](https://github.com/phylo42/xpas/tree/master), a standalone tool for phylo k-mer database construction, which is a submodule of this repository.

Instead of running commands
```
python rappas2.py build OPTIONS...
```
which were supported in earlier versions of RAPPAS2, run:

```
python xpas.py OPTIONS...
```

### Phylogenetic placement
```
python rappas2.py place -i DATABASE -s [nucl|amino] -o OUTPUT_DIR INPUT_FASTA
```
See `python rappas2.py place --help` for more information.
