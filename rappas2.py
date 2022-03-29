#!/usr/bin/env python3
"""
RAPPAS2 wrapper script.

Since not all the functionality of RAPPAS is yet implemented in RAPPAS2,
one has to call RAPPAS to produce intermediate files needed to run RAPPAS2.
This script makes usage of RAPPAS2 transparent by calling RAPPAS where needed,
thus providing a clean CLI interface which is very similar to the one from RAPPAS.
"""

__author__ = "Nikolai Romashchenko"
__license__ = "MIT"


import os
from os import path
from pathlib import Path
import click
import subprocess
from pathlib import Path


@click.group()
def rappas():
    """
    RAPPAS2

    N. Romashchenko, B. Linard, F. Pardi, E. Rivals
    """
    pass


@rappas.command()
@click.option('-i', '--database',
              required=True,
              type=click.Path(dir_okay=False, file_okay=True, exists=True),
              help="Input database.")
@click.option('-s', '--states',
              type=click.Choice(['nucl', 'amino']),
              default='nucl', show_default=True,
              required=True,
              help="States used in analysis.")
@click.option('-o', '--outputdir',
              required=True,
              type=click.Path(dir_okay=True, file_okay=False),
              help="Output directory.")
@click.option('--threads',
             type=int,
             default=4, show_default=True,
             help="Number of threads used.")
@click.argument('input_files', type=click.Path(exists=True), nargs=-1)
def place(database, states, outputdir, threads, input_files):
    """
    Places .fasta files using the input RAPPAS2 database.

    \tpython rappas2.py place -s [nucl|amino] -i db.rps -o output file.fasta [file2.fasta ...]

    """
    current_dir = os.path.dirname(os.path.realpath(__file__))

    if states == 'nucl':
        rappas_bin = f"{current_dir}/bin/rappas/rappas2-dna"
    else:
        rappas_bin = f"{current_dir}/bin/rappas/rappas2-aa"
    

    command = [
        rappas_bin,
        str(database),
        str(outputdir),
        str(threads)
    ]
    command.extend(input_files)
    subprocess.call(command)


if __name__ == "__main__":
    rappas()
