# EPIK: Evolutionary Placement with Informative K-mers


[![build](https://github.com/phylo42/EPIK/actions/workflows/build.yml/badge.svg)](https://github.com/phylo42/EPIK/actions/workflows/build.yml)
<a>
<img src="https://img.shields.io/badge/softwipe-7.6-green" />
</a>

**Please cite:**  [![doi](https://img.shields.io/static/v1?label=doi&message=10.1093/bioinformatics/btad692&color=blue)](https://doi.org/10.1093/bioinformatics/btad692) [1]
        
EPIK is a program for rapid alignment-free phylogenetic placement, the successor of [RAPPAS](https://github.com/phylo42/RAPPAS).

## Installation via Bioconda

It is advised to install the package in a new environment, because our C++ dependencies are strict and may clash with other packages (requiring libboost in particular).
We also recommend to use `mamba, which is faster in solving environment dependencies.
```
conda create -n epik
conda activate epik
conda config set channel_priority strict

# note that we install both ipk (database creation) and epik (phylogenetic placement)
mamba install ipk epik
```

Rapid test:
```
# get some test alignment and tree
wget https://github.com/phylo42/IPK/blob/main/tests/data/D652/reference.fasta
wget https://github.com/phylo42/IPK/blob/main/tests/data/D652/tree.rooted.newick

# activate conda environment
conda activate epik

# build database with IPK : using 1 CPU and default phylogenetic model parameters
# a better approach would be to use appropriate parameters, see documentation
ipk.py build --refalign reference.fasta --reftree tree.rooted.newick --states nucl --workdir . --model GTR

# place with EPIK
epik.py place -i DB.ipk -s nucl -o . reference.fasta

# jplace results
cat placements_reference.fasta.jplace

# you can do post-analyses with the excellent 'gappa' package
# (available in bioconda too, see https://github.com/lczech/gappa)
```


## Installation via compilation

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
pip3 install click
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
To place queries to a phylogenetic tree, you need to first preprocess it with IPK and make a phylo-k-mer database (see [here](https://github.com/phylo42/IPK) for detail). Queries should be in non-compressed fasta format. An example of placement command (see below for possible parameters values):
```
epik.py place -i DATABASE -s [nucl|amino] -o OUTPUT_DIR INPUT_FASTA
```
If EPIK is not installed, run `./epik.py` from the EPIK directory instead. 

### Parameters

| Option    | Meaning                                                                                                                                                                 | Default |
|-----------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------|
| -i        | The path to the phylo-k-mer database to use for placement.                                                                                                              |         |
| -s        | States, `nucl` for DNA and `amino` for proteins                                                                                                                         | nucl    |
| --omega   | The user-defined threshold. Can be set higher than the one used when database was created. (If you are not sure, ignore this parameter.)                                | 1.5     |
| --mu      | The proportion of the database to keep when filtering. Mutually exclusive with `--max-ram`. Should be a value in (0.0, 1.0]                                             | 1.0     |
| --max-ram | The maximum amount of memory used to keep the database content. Mutually exclusive with `--mu`. Sets an approximate limit to EPIK's RAM consumption (i.e. the given limit might be exceeded but EPIK will consider it). Examples: 512, 256K, 42M, 4.2G.                    |         |
| --threads | Number of parallel threads used for placement. EPIK should be compiled with OpenMP support enabled, i.e. `EPIK_OMP=ON`. (If you compile as we recommend, it is enabled) | 1       |

Also, see `epik.py place --help` for information.


## Other

### Code quality

Code quality evaluation with [softwipe](https://github.com/adrianzap/softwipe) [2]:
```
softwipe --cmake --cpp -x third-party,i2l/third-party,i2l/tests/catch2,i2l/examples --no-execution .
```


## References
[1] Romashchenko, N., Linard, B., Pardi, F., & Rivals, E. (2023). EPIK: precise and scalable evolutionary placement with informative k-mers. Bioinformatics, 39(12), btad692. https://doi.org/10.1093/bioinformatics/btad692

[2] Zapletal, A., Höhler, D., Sinz, C., & Stamatakis, A. (2021). The SoftWipe tool and benchmark for assessing coding standards adherence of scientific software. Scientific reports, 11(1), 10015. https://doi.org/10.1038/s41598-021-89495-8
