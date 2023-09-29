#!/usr/bin/env python3
"""
EPIK: Evolutionary Placement with Informative K-mers
"""

__author__ = "Nikolai Romashchenko"
__license__ = "MIT"


import os
import click
import subprocess


__version__ = "0.2.0"


@click.group()
@click.version_option(__version__)
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
              type=click.Path(dir_okay=True, file_okay=False, exists=True),
              help="Output directory.")
@click.option('--threads',
             type=int,
             default=1, show_default=True,
             help="Number of threads used.")
@click.option('--max-ram',
             type=str,
             default="", show_default=True,
             help="Approximate RAM limit to use. Database may not be fully loaded")
@click.argument('input_file', type=click.Path(exists=True))
def place(database, states, omega, mu, outputdir, threads, max_ram, input_file):
    """
    Places .fasta files using the input IPK database.

    epik.py place -s [nucl|amino] -i DB.ipk -o output file.fasta [file2.fasta ...]

    Examples:
    \tepik.py place -i DB.ipk -o temp --max-ram 4G --threads 8 query.fasta

    """
    place_queries(database, states, omega, mu, outputdir, threads, max_ram, input_file)


def place_queries(database, states, omega, mu, outputdir, threads, max_ram, input_file):
    current_dir = os.path.dirname(os.path.realpath(__file__))
    
    # If EPIK is installed, look for the binary in the installed location,
    # otherwise it is run from sources
    epik_bin_dir = f"{current_dir}" if os.path.exists(f"{current_dir}/epik-dna") else f"{current_dir}/bin/epik"

    if states == 'nucl':
        epik_bin = f"{epik_bin_dir}/epik-dna"
    else:
        epik_bin = f"{epik_bin_dir}/epik-aa"

    command = [
        epik_bin, 
        "-d", str(database),
        "-q", str(input_file),
        "-j", str(threads),
        "--omega", str(omega),
        "--mu", str(mu),
        "-o", str(outputdir),
    ]
    if max_ram:
        command.extend(["--max-ram", max_ram])
    command.append(input_file)
    print(" ".join(s for s in command))
    return subprocess.call(command)


if __name__ == "__main__":
    epik()
