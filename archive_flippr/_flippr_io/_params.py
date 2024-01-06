# ---------------------------------------------------------------------
# FragPipe Limited-Proteolysis Processor by Edgar Manriquez-Sandoval
# Stephen D. Fried Lab @ The Johns Hopkins University
# ---------------------------------------------------------------------

import itertools

from os import mkdir
from csv import reader
from pathlib import Path
from datetime import datetime
from argparse import Namespace

from . import _msg

CUTSITE_IDS = [
    "Protein ID",
    "Protein",
    "Gene",
    "Half Tryptic",
    "Cut Site",
    "Cut Site Residue",
    "Cut Site Position",
]

PEPTIDE_IDS = [
    "Peptide Sequence",
    "Modified Sequence",
    "Assigned Modifications",
    "Protein ID",
    "Protein",
    "Gene",
    "Prev AA",
    "Next AA",
    "Start",
    "End",
    "Peptide Length",
    "Half Tryptic",
    "Cut Site",
    "Cut Site Residue",
    "Cut Site Position",
]


# not a good way of doing this,
# 'do as I say not as I do' sort of thing
# error handling should not be THIS integrated with the class
class Params:
    def __init__(self, args: Namespace) -> None:
        self.args = args
        self.time = datetime.now()

    def validate_args(self) -> None:
        self.__validate_replicate_value()
        self.__validate_combined_ion_file()
        self.__validate_combined_protein_file()
        self.__generate_out_paths()

    def fetch_condition_parameters(self, cond: str) -> list:
        return self.condition_replicates[cond], self.cond_out_paths[cond]

    def __validate_replicate_value(self):
        if self.args.replicates is None:
            _msg.error("Missing value for : -n/--replicates")
        else:
            replicates = self.args.replicates[0]

        if replicates <= 1:
            _msg.error(
                f"The value `{replicates:d}` does not satisfy the arguments : -n/--replicates",
                "Replicate value must be greater than one (1)",
            )

        if replicates < 3:
            _msg.notice(
                f"The replicate value `{replicates:d}` is less than `3`",
                "Analysis confidence could be affected",
            )

        self.replicates = replicates

    def __validate_combined_ion_file(self):
        ion_file = Path(self.args.filename)

        if not ion_file.exists():
            _msg.error(f"The file `{ion_file.name:s}` does not exist")

        if ion_file.suffix != ".tsv":
            _msg.error("Expected file `combined_ions.tsv`")

        control_annotation = self.args.control[0]

        control_replicates = tuple(
            [
                f"{c}_{i+1}"
                for c, i in itertools.product(
                    [control_annotation], range(self.replicates)
                )
            ]
        )

        condition_annotations = [d for cond in self.args.conditions for d in cond]

        condition_replicates = [
            f"{d}_{i+1}"
            for d, i in itertools.product(condition_annotations, range(self.replicates))
        ]

        with open(ion_file, "r") as open_file:
            main_file_header = next(reader(open_file, delimiter="\t"))

        if not all(
            True if f"{c} Intensity" in main_file_header else False
            for c in control_replicates
        ):
            _msg.error(f"Control experimental annotation not found in `{ion_file}`")

        if not all(
            True if f"{d} Intensity" in main_file_header else False
            for d in condition_replicates
        ):
            _msg.error(
                f"Condition experimental annotation(s) not found in `{ion_file}`"
            )

        condition_replicates = {
            k: v
            for k, v in zip(
                condition_annotations,
                list(
                    itertools.zip_longest(
                        *[iter(condition_replicates)] * self.replicates
                    )
                ),
            )
        }

        self.ion_file = ion_file
        self.control_annotation = control_annotation
        self.control_replicates = control_replicates
        self.condition_annotations = condition_annotations
        self.condition_replicates = condition_replicates

    def __validate_combined_protein_file(self):
        norm_file = self.args.norm_filename

        if norm_file and (
            self.args.norm_control is None or self.args.norm_conditions is None
        ):
            _msg.error(
                "Normalized file was pass but does not satisfy the arguments : -zc/--norm_control and -zd/--norm_conditions"
            )

        if norm_file is not None:
            norm_file = Path(norm_file)

            norm_control_annotation = self.args.norm_control[0]

            norm_control_replicates = tuple(
                [
                    f"{c}_{i+1}"
                    for c, i in itertools.product(
                        [norm_control_annotation], range(self.replicates)
                    )
                ]
            )

            norm_condition_annotations = [
                d for cond in self.args.norm_conditions for d in cond
            ]

            norm_condition_replicates = [
                f"{d}_{i+1}"
                for d, i in itertools.product(
                    norm_condition_annotations, range(self.replicates)
                )
            ]

            if not norm_file.exists():
                _msg.error(f"The file `{norm_file.name:s}` does not exist")

            if norm_file.suffix != ".tsv":
                _msg.error("Expected file `combined_protein.tsv` or similar")

            with open(norm_file, "r") as open_file:
                norm_file_header = next(reader(open_file, delimiter="\t"))

            if not all(
                True if f"{c} MaxLFQ Intensity" in norm_file_header else False
                for c in norm_control_replicates
            ):
                _msg.error(
                    f"Normalized control experimental annotation not found in `{norm_file}`"
                )

            if not all(
                True if f"{d} MaxLFQ Intensity" in norm_file_header else False
                for d in norm_condition_replicates
            ):
                _msg.error(
                    f"Normalized condition experimental annotation(s) not found in `{norm_file}`"
                )

            if len(norm_condition_annotations) == 1:
                norm_condition_replicates = {
                    k: v
                    for k in self.condition_annotations
                    for v in list(
                        itertools.zip_longest(
                            *[iter(norm_condition_replicates)] * self.replicates
                        )
                    )
                }

            else:
                if len(norm_condition_annotations) == len(self.condition_annotations):
                    norm_condition_replicates = {
                        k: v
                        for k, v in zip(
                            self.condition_annotations,
                            list(
                                itertools.zip_longest(
                                    *[iter(norm_condition_replicates)] * self.replicates
                                )
                            ),
                        )
                    }

                else:
                    _msg.error(
                        f"Mismatching numeber of Normalized condition experimental annotation(s) : {len(norm_condition_annotations)} and Condition experimental annotation(s) : {len(self.condition_annotations)}",
                        "The number of Normalized condition experimental annotation(s) must be 1 or the same as Condition experimental annotation(s)",
                    )
        else:
            (
                norm_control_annotation,
                norm_control_replicates,
                norm_condition_annotations,
                norm_condition_replicates,
            ) = [None] * 4

        self.norm_file = norm_file
        self.norm_control_annotation = norm_control_annotation
        self.norm_control_replicates = norm_control_replicates
        self.norm_condition_annotations = norm_condition_annotations
        self.norm_condition_replicates = norm_condition_replicates

    def __generate_out_paths(self):
        if self.args.prefix:
            self.prefix = f"{self.args.prefix}_"
        else:
            self.prefix = ""

        t = self.time
        self.out_path = Path(
            self.ion_file.parent,
            f"{self.prefix}FLIPPR_OUTPUT_{t.year:04d}{t.month:02d}{t.day:02d}{t.hour:02d}{t.minute:02d}{t.second:02d}",
        )
        self.cond_out_paths = {
            d: Path(self.out_path, f"{self.prefix}{d}_v_{self.control_annotation}")
            for d in self.condition_annotations
        }

        mkdir(self.out_path)
        [mkdir(p) for p in self.cond_out_paths.values()]
