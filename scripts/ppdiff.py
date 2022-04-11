#!/usr/bin/env python3

# This script runs RAPPAS and XPAS+RAPPAS2 for a given dataset
# (reference tree, alignemnt, queries) with given parameters
# and compares if resulting placements are equivalent.
# See config.json for an examples of a configuration file.
#
#   Usage: python place_compare.py CONFIG

__author__ = "Nikolai Romashchenko"
__license__ = "MIT"

import click
import os
import json
import distutils
import subprocess
from pathlib import Path
from click.testing import CliRunner
import jplace_diff as jd
from xpas import xpas
from rappas2 import rappas2


# abstract class for RAPPAS and XPAS+RAPPAS2
# has properties for run parameters
class RappasBase:
    def __init__(self, config):
        self.config = config

    # property for workdir
    @property
    def workdir(self):
        return self.config["args"]["workdir"]

    # property for k
    @property
    def k(self):
        return self.config["args"]["k"]

    # property for omega
    @property
    def omega(self):
        return self.config["args"]["omega"]

    # property for queries
    @property
    def queries(self):
        return self.config["args"]["queries"]

    @property
    def placement_file(self):
        return os.path.join(self.workdir, f"placements_{os.path.basename(self.queries)}.jplace")

    @property
    def database_file(self):
        # the subclasses should override this
        raise NotImplementedError

    @property
    def renamed_placement_file(self):
        raise NotImplementedError

    # property that checks if the database file exists
    @property
    def database_exists(self):
        return os.path.isfile(self.database_file)

    # runs RAPPAS / XPAS+RAPPAS2
    def run(self) -> str:
        # create the workdir if not exists
        Path(self.workdir).mkdir(parents=True, exist_ok=True)

        # build the database if not already done
        if not self.database_exists:
            self.build_database()
        else:
            print(f"Database file {self.database_file} already exists. Skipping build")

        # place the queries if not already done
        if not os.path.isfile(self.renamed_placement_file):
            self.place_queries()
            # rename placement file to avoid confusion between RAPPAS and RAPPAS2
            print("Rename: " + self.placement_file + " -> " + self.renamed_placement_file)
            os.rename(self.placement_file, self.renamed_placement_file)
        else:
            print(f"Placement file {self.renamed_placement_file} already exists. Skipping placement")

        return self.renamed_placement_file


class Rappas(RappasBase):
    @property
    def database_file(self):
        return os.path.join(self.workdir, f"DB_session_k{self.k}_o{self.omega}.union")

    @property
    def renamed_placement_file(self):
        # we rename placement files to make the difference between the results of RAPPAS and RAPPAS2
        return f"{self.placement_file}.rappas"

    def build_database(self):
        print(f"Building database {self.database_file}")

        # the command to run RAPPAS
        command = [
            "java",
            "-jar", self.config["soft"]["rappas"],
            "-b", self.config["soft"]["arbinary"],
            "-p", "b",
            "-s", self.config["args"]["states"],
            "-k", self.k,
            "--omega", self.omega,
            "-t", self.config["args"]["tree"],
            "-r", self.config["args"]["alignment"],
            "-w", self.workdir,
            "--ratio-reduction", self.config["args"]["reduction-ratio"],
            "--ardir", self.config["args"]["ardir"],
            "--gap-jump-thresh", "1.0"
        ]

        use_unrooted = "use-unrooted" in self.config["args"] and \
                       self.config["args"]["use-unrooted"]

        if use_unrooted:
            command.append("--use_unrooted")

        print(" ".join(s for s in command))
        subprocess.call(command)

        if not self.database_exists:
            raise Exception(f"Error! Database file {self.database_file} not found after building")

    # method that places queries with RAPPAS
    def place_queries(self):
        command = [
            "java",
            "-jar", self.config["soft"]["rappas"],
            "-b", self.config["soft"]["arbinary"],
            "-p", "p",
            "-s", self.config["args"]["states"],
            "-k", self.k,
            "--omega", self.omega,
            "-d", self.database_file,
            "-w", self.workdir,
            "-q", self.queries
        ]

        # run RAPPAS
        print(" ".join(s for s in command))
        return_code = subprocess.call(command)

        # check the return code of the RAPPAS call
        if return_code != 0:
            raise Exception(f"Error! RAPPAS returned {return_code}")

        # check that the placement file exists
        if not os.path.isfile(self.placement_file):
            raise Exception(f"Error! Placement file {self.placement_file} not found after running RAPPAS")


class Rappas2(RappasBase):
    @property
    def database_file(self):
        return os.path.join(self.workdir, f"DB_k{self.k}_o{self.omega}.rps")

    @property
    def renamed_placement_file(self):
        # we rename placement files to make the difference between the results of RAPPAS and RAPPAS2
        return f"{self.placement_file}.rappas2"

    # builds the database with XPAS
    def build_database(self):
        print(f"Building database {self.database_file}")

        verbosity = 0
        write_reduction = False
        alpha = 1.0
        categories = 4
        model = "GTR"
        arparameters = ""
        convert_uo = True
        no_reduction = False
        filter = "no-filter"
        f = 1.0
        mu = 1.0
        use_unrooted = "use-unrooted" in self.config["args"] and self.config["args"]["use-unrooted"]
        merge_branches = False
        aronly = False
        keep_positions = False
        uncompressed = False
        threads = 2
        #score_model = self.config["args"]["model"]

        xpas.build_database(self.config["soft"]["arbinary"],  # database,
                            self.config["args"]["alignment"], self.config["args"]["tree"],
                            self.config["args"]["states"], verbosity,
                            self.workdir, write_reduction,  # dbfilename,
                            alpha, categories,  # ghosts,
                            self.k, model, arparameters, convert_uo,  # gap_jump_thresh,
                            no_reduction, self.config["args"]["reduction-ratio"], self.omega,
                            filter, #score_model,
                            f, mu, use_unrooted, merge_branches,
                            self.config["args"]["ardir"], aronly,
                            keep_positions, uncompressed,
                            threads)

        if not self.database_exists:
            raise Exception(f"Error! Database file {self.database_file} not found after building")

    # method that places queries with RAPPAS2
    def place_queries(self):
        threads = 2
        return_code = rappas2.place_queries(self.database_file, self.config["args"]["states"],
                                            self.workdir, threads, [self.queries])

        # check the return code of the RAPPAS call
        if return_code != 0:
            raise Exception(f"Error! RAPPAS2 returned {return_code}")

        # check that the placement file exists
        if not os.path.isfile(self.placement_file):
            raise Exception(f"Error! Placement file {self.placement_file} not found after running RAPPAS2")


# checks that file exists
def check_file_exists(file):
    if not os.path.isfile(file):
        raise Exception(f"Error! File {file} not found")


@click.command()
@click.argument('config_file', type=click.Path(exists=True))
def run(config_file: str) -> None:
    with open(config_file, "r") as jsonfile:
        config = json.load(jsonfile)

        # create RAPPAS and RAPPAS2 objects
        rappas = Rappas(config)
        rappas2 = Rappas2(config)

        # run software
        rappas_placement = rappas.run()
        rappas2_placement = rappas2.run()

        # check that the placement files exist
        check_file_exists(rappas_placement)
        check_file_exists(rappas2_placement)

        # compare the two placement files
        jd.jplace_diff(rappas_placement, rappas2_placement)


if __name__ == "__main__":
    run()
