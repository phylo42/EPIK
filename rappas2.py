#!/usr/bin/env python3
"""
RAPPAS / RAPPAS2 wrapper script.

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


NUCL_MODELS = ['JC69', 'HKY85', 'K80', 'F81', 'TN93', 'GTR']
AMINO_MODELS = ['LG', 'WAG', 'JTT', 'Dayhoff', 'DCMut', 'CpREV', 'mMtREV', 'MtMam', 'MtArt']
ALL_MODELS = NUCL_MODELS + AMINO_MODELS


KMER_FILTERS = ["no-filter", "entropy", "max-deviation", "max-difference", "random"]



def validate_filter(ctx, param, value):
    value = value.lower()
    if value not in KMER_FILTERS:
        valid_values = ', '.join(v for v in KMER_FILTERS)
        raise click.BadParameter('Filter must be one of: ' + valid_values)
    return value


@rappas.command()
@click.option('-b', '--arbinary',
              type=click.Path(exists=True),
              required=True,
              help="Binary file for marginal AR, currently 'phyml' and "
                   "'baseml' (from PAML) are supported.")
@click.option('-r', '--refalign', 
              type=click.Path(exists=True),
              required=True,
              help="""Reference alignment in fasta format.
                  It must be the multiple alignment used to build the reference tree.""")
@click.option('-t', '--reftree',
              type=click.Path(exists=True),
              required=True,
              help=" Reference tree in the newick format")
@click.option('-s', '--states',
              type=click.Choice(['nucl', 'amino']),
              default='nucl', show_default=True,
              required=True,
              help="States used in analysis.")
@click.option('-v', '--verbosity',
              type=int,
              default=0, show_default=True,
              help="Verbosity level: -1=none ; 0=default ; 1=high")
@click.option('-w', '--workdir',
              required=True,
              type=click.Path(dir_okay=True, file_okay=False),
              help="Working directory for temp files.")
@click.option('--write-reduction',
              type=click.Path(file_okay=True, dir_okay=False),
              help=" Write reduced alignment to file.")
#@click.option('--dbfilename',
#              type=str,
#              default="db.rps", show_default=True, 
#              help="Output database filename.")
@click.option('-a', '--alpha',
              type=float,
              default=1.0, show_default=True,
              help="Gammma shape parameter used in ancestral reconstruction.")
@click.option('-c', '--categories',
              type=int,
              default=4, show_default=True,
              help="Number of categories used in ancestral reconstruction.")
#@click.option('-g', '--ghosts',
#              type=int,
#              default=1, show_default=True,
#              help="Number of ghost nodes injected per branch.")
@click.option('-k', '--k',
             type=int,
             default=8, show_default=True,
             help="k-mer length used at DB build.")
@click.option('-m', '--model',
             type=click.Choice(ALL_MODELS),
             required=True,
             help="Model used in AR, one of the following:\n"
                  f"nucl: {', '.join(x for x in NUCL_MODELS)}\n"
                  f"amino: {', '.join(x for x in AMINO_MODELS)}")
@click.option('--arparameters',
             type=str,
             help="""Parameters passed to the software used for
                  anc. seq. reconstuct. Overrides -a,-c,-m options.
                  Value must be quoted by ' or ". Do not set options
                  -i,-u,--ancestral (managed by RAPPAS).""")
@click.option('--convert-uo',
              is_flag=True,
              help="U, O amino acids are converted to C, L.")
@click.option('--force-root',
              is_flag=True,
              help="""Root input tree (if unrooted) by adding a root
                  node on righmost branch of the trifurcation.""")
#@click.option('--gap-jump-thresh',
#              type=float,
#              deafult=0.3, show_default=True,
#              help="Gap ratio above which gap jumps are activated.")
@click.option('--no-reduction',
              is_flag=True,
              help="""Do not operate alignment reduction. This will 
                  keep all sites of input reference alignment and 
                  may produce erroneous ancestral k-mers.""")
@click.option('--ratio-reduction',
              type=float,
              default=0.99, show_default=True,
              help="""Ratio for alignment reduction, e.g. sites 
                holding >99% gaps are ignored.""")
@click.option('--omega',
              type=float,
              default=1.5, show_default=True,
              help="""Modifier levelling the threshold used during
                  phylo-kmer filtering, T=(omega/#states)^k""")
@click.option('--filter',
              callback=validate_filter,
              default="no-filter", show_default=True)
@click.option('-u', '--mu',
              type=float,
              default=0.5, show_default=True,
              help="""K-mer filter threshold""")
@click.option('--use-unrooted',
              is_flag=True,
              help="""Confirms you accept to use an unrooted reference
                  tree (option -t). The trifurcation described by the
                  newick file will be considered as root. Be aware that
                  meaningless roots may impact accuracy.""")
@click.option('--ardir',
             type=click.Path(exists=True, dir_okay=True, file_okay=False),
             help="""Skip ancestral sequence reconstruction, and 
                  uses outputs from the specified directory.""")
@click.option('--aronly',
             is_flag=True,
             default=False, show_default=True,
             help="Dev option. Run only ancestral reconstruction and tree extension. No database will be built")
@click.option('--threads',
             type=int,
             default=4, show_default=True,
             help="Number of threads used.")
def build(arbinary, #database,
          refalign, reftree, states, verbosity,
          workdir, write_reduction, #dbfilename,
          alpha, categories, #ghosts,
          k, model, arparameters, convert_uo, force_root, #gap_jump_thresh,
          no_reduction, ratio_reduction, omega,
          filter, mu, use_unrooted, ardir,
          threads, aronly):
    """
    Builds a database of phylo k-mers.

    Minimum usage:

    \tpython rappas2.py build -s [nucl|amino] -b ARbinary -w workdir -r alignment.fasta -t tree.newick

    """
    # create working directory
    Path(workdir).mkdir(parents=True, exist_ok=True)

    rappas1_results_dir = workdir
    if ardir:
        rappas1_results_dir = Path(ardir).parent

    # auxilary files produced by RAPPAS
    extended_tree = f"{rappas1_results_dir}/extended_trees/extended_tree_withBL.tree"
    ar_seq_txt = f"{rappas1_results_dir}/AR/extended_align.phylip_phyml_ancestral_seq.txt"
    extended_tree_node_mapping = f"{rappas1_results_dir}/extended_trees/extended_tree_node_mapping.tsv"
    artree_id_mapping = f"{rappas1_results_dir}/AR/ARtree_id_mapping.tsv"

    # check if auxilary files are alredy present
    files_required = [extended_tree, ar_seq_txt, extended_tree_node_mapping, artree_id_mapping]
    current_dir = os.path.dirname(os.path.realpath(__file__))

    # if not, run RAPPAS1
    return_code = 0
    if not all(path.exists(f) for f in files_required):
        rappas_jar = f"{current_dir}/rappas/dist/RAPPAS.jar"

        command = [
            "java", #"-Xmx1G",
            "-jar", rappas_jar,
            "--phase", "b",
            "--arbinary", str(arbinary),
            #"--database", database,
            "--refalign", str(refalign),
            "--reftree", str(reftree),
            "--states", states,
            "--verbosity", str(verbosity),
            "--workdir", str(workdir),
            "--alpha", str(alpha),
            "--categories", str(categories),
            "--k", str(k),
            "--model", model,
            "--ratio-reduction", str(ratio_reduction),
            "--omega", str(omega),
            "--threads", str(threads),
            "--aronly"
        ]
        if ardir:
            command.extend(["--ardir", ardir])
        if write_reduction:
            command.append("--write-reduction")
        if arparameters:
            command.extend(["--arparameters", arparameters])
        if convert_uo:
            command.append("--convertUO")
        if force_root:
            command.append("--force-root")
        if no_reduction:
            command.append("--no-reduction")
        if use_unrooted:
            command.append("--use_unrooted")
        return_code = subprocess.call(command)

    # exit for --aronly
    if aronly:
        print("AR only requested. No database was built")
        return

    # exit if RAPPAS failed
    if return_code != 0:
        raise RuntimeError(f"RAPPAS failed with exit code {return_code}.")

    # run rappas2 if RAPPAS succeed
    if return_code == 0:
        if states == 'nucl':
            rappas_bin = f"{current_dir}/bin/build/rappas-buildn"
        else:
            raise RuntimeError("Proteins are not supported yet.")

        command = [
            rappas_bin,
            "-t", str(reftree),
            "-x", extended_tree,
            "-a", ar_seq_txt,
            "-e", extended_tree_node_mapping,
            "-m", artree_id_mapping,
            "-w", str(workdir),
            "-k", str(k),
            "-o", str(omega),
            "--" + filter.lower(),
            "-u", str(mu),
            "-j", str(threads)
        ]
        subprocess.call(command)

        hashmaps_dir = f"{workdir}/hashmaps"
        subprocess.call(["rm", "-r", hashmaps_dir])


@rappas.command()
@click.option('-i', '--database',
              required=True,
              type=click.Path(dir_okay=False, file_okay=True, exists=True),
              help="Input database.")
@click.option('-o', '--outputdir',
              required=True,
              type=click.Path(dir_okay=True, file_okay=False),
              help="Output directory.")
@click.option('--threads',
             type=int,
             default=4, show_default=True,
             help="Number of threads used.")
@click.argument('input_files', type=click.Path(exists=True), nargs=-1)
def place(database, outputdir, threads, input_files):
    """
    Places .fasta files using the input RAPPAS2 database.

    \tpython rappas2.py place -i db.rps -o output file.fasta [file2.fasta ...]

    """
    current_dir = os.path.dirname(os.path.realpath(__file__))
    rappas_bin = f"{current_dir}/bin/place/rappas-placen"

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
