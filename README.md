# RAPPAS2
This is the new version of [RAPPAS](https://github.com/phylo42/RAPPAS).

## Installation

### Prerequisites

- Java v11 and Ant
- Boost Libraries >=1.65
- CMake >= 3.10
- GCC compiler must support c++17

In debian, these can be installed with:
```
sudo apt install build-essential cmake libboost-dev libboost-serialization-dev libboost-filesystem-dev libboost-iostreams-dev openjdk-11-jdk ant
```

### Clone and build
```
git clone --recursive https://github.com/phylo42/rappas2.git
cd rappas2
./build.sh
```

## Usage
```
python rappas2.py build -b `which phyml` -r reference.fasta -t tree.newick -m GTR -k 8 -w work_dir
```

See `python rappas2.py --help` and `python rappas2.py build --help` for more information.
