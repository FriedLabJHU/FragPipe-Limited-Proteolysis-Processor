from functools import cached_property
from pathlib import Path

import polars as pl

from . import spr as _spr
from . import validate as _validate
from .parameters import rcParams
from .uniprot import fasta as _fasta


class Study:
    """
    A FLiPPR Study is used to coordinate Limited Proteolysis (LiP) data processing for a single experiment analyzed with FragPipe.
    Studies can contain LiP data alone or LiP and TrP data; with the inclusion of a FASTA file for automatic metadata generation.
    TrP data should result from the same experimental conditions as the LiP data, but without including a broad-spectrum protease like Proteinase K (PK).

    Args:
        lip (str | Path): Directory or path to the FragPipe LFQ output data for Limited Proteolysis (LiP) experiment data.
        trp (str | Path, optional): Directory or path to the FragPipe LFQ output data for Trypsin Only (TrP) experiment data. Defaults to None.
        fasta (str | Path, optional): Path to the FASTA file containing the same protein sequences used in the LiP/TrP experiments. Defaults to None.
        method (str): Method used to quantify ion intensities in FragPipe. `lfq` - Label-Free Quantification, `silac` - Stable Isotope Labelling (currently not incorporated). Defaults to `lfq`.

    Examples:
        Analysis of a single dataset
        >>> flippr.Study(lip = "PXD025926/FragPipe/LiP_LFQ")

        Analysis of a single dataset with metadata integration
        >>> flippr.Study(lip = "PXD025926/FragPipe/LiP_LFQ", fasta = "UP000000625.fasta")

        Including protein-level normalization from trypsin-only data
        >>> flippr.Study(lip = "PXD025926/FragPipe/LiP_LFQ", trp = "PXD025926/FragPipe/TrP_LFQ")

    """

    def __init__(self, lip: str | Path, trp: str | Path | None = None, fasta: str | Path | None = None, method: str = "lfq") -> None:

        _lip, _trp, _fasta_path, _method = _validate._validate_study(lip, trp, fasta, method)

        self.lip: Path = _lip
        self.trp: Path | None = _trp
        self._fasta_path: Path | None = _fasta_path
        self.method: str = _method

        self._pid: int = 0
        self.processes: dict = dict()
        self.results: dict = dict()

    @property
    def samples(self) -> dict:
        """
        Returns the names of the samples or experimental annotations found in the FragPipe outputs for LiP and TrP (if include) datasets.

        Examples:
            Including only LiP data
            >>> e_coli_refolding = flippr.Study(lip = "PXD025926/FragPipe/LiP_LFQ")
            >>> print(e_coli_refolding.samples)
            {"LiP": ("Native", "Refolded_01_min", "Refolded_05_min", "Refolded_30_min")}

            Including LiP and TrP data
            >>> e_coli_refolding = flippr.Study(lip = "PXD025926/FragPipe/LiP_LFQ", trp = "PXD025926/FragPipe/TrP_LFQ")
            >>> print(e_coli_refolding.samples)
            {"LiP": ("Native", "Refolded_01_min", "Refolded_05_min", "Refolded_30_min"),
             "TrP": ("Native", "Refolded")}

        """

        if self.trp is not None:
            return {"LiP": self.lip_samples, "TrP": self.trp_samples}

        return {"LiP": self.lip_samples}

    @property
    def lip_samples(self) -> set:
        """
        Returns the names of the samples or experimental annotations found in the FragPipe outputs for LiP datasets.

        """

        return self.__read_annotation(self.lip)

    @property
    def trp_samples(self) -> set:
        """
        Returns the names of the samples or experimental annotations found in the FragPipe outputs for TrP datasets.

        """

        return self.__read_annotation(self.trp)

    @cached_property
    def fasta(self) -> pl.DataFrame | None:
        """
        Returns a three column dataframe containing Protein IDs, Sequences, and Valid ID if the Protein ID is interpretable as an official Uniprot ID.

        """
        if self._fasta_path is not None:
            return _fasta._read_fasta(self._fasta_path)

        return None

    def add_process(
        self, pid: object | None, lip_ctrl: str, lip_test: str, n_rep: int, trp_ctrl: str | None = None, trp_test: str | None = None, trp_n_rep: int | None = None) -> None:
        """
        Adding a process to calculate the fold-change between test and control conditions.
        Normalization is optional and must include a TrP experiment in the `Study` set-up.

        Args:
            pid (Any): Process ID used to index the `Study().results` dictionary after running `Study().run()`.
            lip_ctrl (str): Control condition sample name from the LiP experiment.
            lip_test (str): Test condition sample name from the LiP experiment.
            n_rep (int): Number of replicates in the LiP experiment.
            trp_ctrl (str, optional): Control condition sample name from the TrP experiment.
            trp_test (str, optional): Test condition sample name from the TrP experiment.
            trp_n_rep (int, optional): Number of replicates in the TrP experiment.

        Examples:
            Simple fold-change calculation between two conditions
            >>> study.add_process("sample", "WT", "Drug", 3)

            Fold-change calculation with normalization
            >>> study.add_process("sample", "WT", "Drug", 3, "WT_TrP", "DMSO_TrP", 3)

            Add multiple processes to the same study
            >>> study.add_process("cond_1", "WT", "Cond1", 5)
            >>> study.add_process("cond_2", "WT", "Cond2", 5, "WT_TrP", "DMSO_TrP", 5)
        """

        if pid is None:
            pid = self._pid
            self._pid += 1

        self.processes.update(
            {
                pid: _spr.Process(
                    self.fasta,
                    self.lip,
                    self.trp,
                    pid,
                    lip_ctrl,
                    lip_test,
                    n_rep,
                    trp_ctrl,
                    trp_test,
                    trp_n_rep,
                )
            }
        )

    def run(self) -> dict[object: _spr.Result]:
        """
        Run the processes added to the study.
        Study parameters can be changed by editing the `flippr.rcParams` dictionary.
        """

        self.results = {pid: proc.run(rcParams) for pid, proc in self.processes.items()}

        return self.results

    def __read_annotation(self, liptrp: Path) -> set:
        def __strip_fragpipe_repilicate(sample: str) -> str:
            return "_".join([_ for _ in sample.split("_")[:-1]])

        annotation = liptrp.joinpath("experiment_annotation.tsv")

        annotation_df = pl.read_csv(annotation, separator="\t")

        annotation_dict = annotation_df.select(pl.col("sample")).to_dict(
            as_series=False
        )

        return set(__strip_fragpipe_repilicate(sample) for sample in annotation_dict["sample"])
