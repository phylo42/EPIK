#!/usr/bin/env python3

# This script compares two .jplace-formatted files to understand how different placements are.
# They are supposed to be two placement results for the same set of queries.
# Usage: python jplace.py FILE1 FILE2

__author__ = "Nikolai Romashchenko"
__license__ = "MIT"


import click
import json
import tqdm
from copy import deepcopy
from typing import Dict, List, Union, Mapping


Number = Union[int, float]
SeqID = str

EPSILON = 1e-4


class PlacementRecord:
    """
    A container for a placement record, e.g.
    Example: [1, -0.1, 0.9, 0.1, 0.0]
    """
    def __init__(self, values: List[Number], fields: List[str]) -> None:
        self._values = values
        self._fields = fields

    def __getattr__(self, item: str):
        """
        'fields' indicated which fields are reported in the .jplace,
        so the list of values can be ordered differently. This method returns
        a value from the list of values by its name.
        Examples of : "edge_num"
        """
        if item not in self._fields:
            raise RuntimeError(f"Wrong field: {item}. "
                               f"Fields listed in the file: {self._fields}")

        return self._values[self._fields.index(item)]


class PlacedSeq:
    """
    A container for a placed sequence.
    """
    def __init__(self,
                 # can be one placement (list) or list of placements
                 placements: Union[PlacementRecord, List[PlacementRecord]],
                 # can be a list with one name or a list of lists with name multiplicity
                 names: Union[List[str], List[List[Union[str, int]]]]) -> None:
        self._placements = placements
        self._names = names

    @staticmethod
    def from_dict(placement_dict: Dict, fields: List[str]) -> "PlacedSeq":
        """
        Creates a PlacedSeq from a placement dictionary.
        Example:
            {
                "p": [...]
                "n": [...]
            }
        """
        placements = [PlacementRecord(p, fields) for p in placement_dict["p"]]

        # Sequence name can be in the field "n" or "nm"
        names_key = "n" if "n" in placement_dict else "nm"
        assert names_key in placement_dict
        names = placement_dict[names_key]

        return PlacedSeq(placements, names)

    @property
    def placements(self):
        return self._placements

    @property
    def names(self):
        return self._names

    @property
    def sequence_name(self):
        # if "name multiplicity"
        if type(self._names[0]) == list:
            return self._names[0][0]
        # if just a name
        else:
            return self._names[0]


class JplaceParser:
    """
    Parses .jplace file, creating a DOM-like data structure using
    PlacedSeq and PlacementRecord.
    """
    def __init__(self, input_file: str) -> None:
        self._input_file = input_file
        self._placements = {}

    def parse(self) -> None:
        """
        Parser the input file, creating a map SeqID -> PlacedSeq in self._placements.
        WARNING: It parses only "edge_num" and "likelihood" fields.
        """

        self._placements = {}
        with open(self._input_file) as jplace_file:
            content = json.load(jplace_file)

            # .jplace file has to have "fields" that determines the order of
            # output fields for each placement. Make sure it is there
            assert "fields" in content, f'{self._input_file} must contain "fields"'
            fields = content["fields"]

            # Make sure the most important two fields are present
            required_fields = ["edge_num", "likelihood"]
            assert all(field in fields for field in required_fields), "Error while parsing " \
                f"{self._input_file}: fields must declare {required_fields}"

            # check if .jplace has at least one placement
            assert "placements" in content,  "Error while parsing " \
                f'{self.input_file}: input file must have the "placements" section.'

            for placement_dict in content["placements"]:
                placed_seq = PlacedSeq.from_dict(placement_dict, fields)

                for name, multiplicity in placed_seq.names:
                    self._placements[name] = placed_seq

    @property
    def placements(self) -> Mapping[SeqID, PlacedSeq]:
        return self._placements


def conditional_print(string, boolean):
    if boolean:
        print(string)


def set_almost_equals(a, b, epsilon=EPSILON):
    """
    Compares two sets for almost equality.
    """
    if len(a) != len(b):
        return False

    la, lb = sorted(list(a)), sorted(list(b))

    for i in range(len(la)):
        if abs(la[i] - lb[i]) > epsilon:
            return False

    return True

def check_first(result1, result2) -> int:
    # if both empty, return 1 matching record
    if len(result1.placements) == len(result2.placements) == 0:
        return 1
    
    # if only one is empty, return a non-match
    if (len(result1.placements) == 0) or (len(result2.placements) == 0):
        return 0
    
    placement1 = result1.placements[0]
    placement2 = result2.placements[0]

    return placement1.edge_num == placement2.edge_num


@click.command()
@click.argument('jplace1', type=click.Path(exists=True))
@click.argument('jplace2', type=click.Path(exists=True))
@click.option('--only-first', is_flag=True, default=False)
def jplace_diff(jplace1: str, jplace2: str, only_first: bool) -> None:
    # parse the input files
    parser1 = JplaceParser(jplace1)
    parser1.parse()

    parser2 = JplaceParser(jplace2)
    parser2.parse()

    num_seqs = len(parser1.placements)
    num_matches = 0
    for name, result1 in tqdm.tqdm(parser1.placements.items()):

        # get the placements for the same sequence from the second file
        result2 = parser2.placements[name]

        if only_first:
            num_matches += check_first(result1, result2)
        else:
            records1 = dict((rec.edge_num, rec.likelihood) for rec in result1.placements)
            records2 = dict((rec.edge_num, rec.likelihood) for rec in result2.placements)
            edge_union = set(records1.keys()).union(set(records2.keys()))

            # check if scores of placements are all the same.
            # If so, we can skip the sequence (no matter what branches are reported)
            scores1 = set(rec.likelihood for rec in result1.placements)
            scores2 = set(rec.likelihood for rec in result2.placements)

            if set_almost_equals(scores1, scores2):
                num_matches += 1
                continue

            # compare placements one by one
            found_mismatch = False
            for edge in edge_union:
                if edge not in records1:
                    conditional_print(f'\n{name}:', not found_mismatch)
                    found_mismatch = True
                    print(f"\t{edge} is not in the first file")
                elif edge not in records2:
                    conditional_print(f'\n{name}:', not found_mismatch)
                    found_mismatch = True
                    print(f"\t{edge} is not in the second file")
                # if found in both, check likelihoods
                elif abs(10**records1[edge] - 10**records2[edge]) > EPSILON:
                    conditional_print(f'\n{name}:', not found_mismatch)
                    found_mismatch = True
                    print(f'\t[{edge}] {records1[edge]} != {records2[edge]}')

            if not found_mismatch:
                num_matches += 1

    print(f"\n{num_matches}/{num_seqs} placements match.")


if __name__ == "__main__":
    jplace_diff()
