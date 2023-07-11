#!/usr/bin/env python3
"""
EPIK: Evolutionary Placement with Informative K-mers
"""

__author__ = "Nikolai Romashchenko"
__license__ = "MIT"


import os
import click
import subprocess


@click.group()
def epik():
    """
    EPIK: Evolutionary Placement with Informative K-mers

    N. Romashchenko, B. Linard, F. Pardi, E. Rivals
    """
    pass


@epik.command()
@click.option('-i', '--database',
              required=True,
              type=click.Path(dir_okay=False, file_okay=True, exists=True),
              help="Input database.")
@click.option('-s', '--states',
              type=click.Choice(['nucl', 'amino']),
              default='nucl', show_default=True,
              required=True,
              help="States used in analysis.")
@click.option('--omega',
              type=float,
              default=1.5,
              help="User omega value, determines the score threhold.")
@click.option('--mu',
              type=float,
              default=1.0,
              help="The proportion of the database to keep.")
@click.option('-o', '--outputdir',
              required=True,
              type=click.Path(dir_okay=True, file_okay=False),
              help="Output directory.")
@click.option('--threads',
             type=int,
             default=4, show_default=True,
             help="Number of threads used.")
@click.argument('input_file', type=click.Path(exists=True))
def place(database, states, omega, mu, outputdir, threads, input_file):
    """
    Places .fasta files using the input IPK database.

    \tpython epik.py place -s [nucl|amino] -i db.rps -o output file.fasta [file2.fasta ...]

    """
    place_queries(database, states, omega, mu, outputdir, threads, input_file)


def place_queries(database, states, omega, mu, outputdir, threads, input_file):
    current_dir = os.path.dirname(os.path.realpath(__file__))

    if states == 'nucl':
        epik_bin = f"{current_dir}/bin/epik/epik-dna"
    else:
        epik_bin = f"{current_dir}/bin/epik/epik-aa"

    command = [
        epik_bin, 
        "-d", str(database),
        "-q", str(input_file),
        "-j", str(threads),
        "--omega", str(omega),
        "--mu", str(mu),
        "-o", str(outputdir),
    ]
    command.extend(input_file)
    print(" ".join(s for s in command))
    return subprocess.call(command)


if __name__ == "__main__":
    epik()
