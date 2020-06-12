# RAPPAS2
This is the new version of [RAPPAS](https://github.com/phylo42/RAPPAS).

## Installation
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
