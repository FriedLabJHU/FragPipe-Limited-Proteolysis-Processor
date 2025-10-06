from . import __about__

from pathlib import Path
from typing import Optional

import polars as pl

from . import datatypes as _types
from . import validate as _validate
from . import reader as _reader
from .parameters import rcParams

__version__ = __about__.__version__

class Study:
    """
    A FLiPPR Study is used to coordinate Limited Proteolysis (LiP) data processing for a single experiment analyzed with FragPipe.
    Studies can contain LiP data alone or LiP and TrP data; with the inclusion of a FASTA file for automatic metadata generation.
    TrP data should result from the same experimental conditions as the LiP data, but without including a broad-spectrum protease like Proteinase K (PK).

    Args:
        lip (str | Path): Directory or path to the FragPipe LFQ output data for Limited Proteolysis (LiP) experiment data.
        trp (str | Path, optional): Directory or path to the FragPipe LFQ output data for Trypsin Only (TrP) experiment data. Defaults to None.
        method (str): Acquisition and search method used to quantify ions intensities in FragPipe. `dda` - Data-dependent; `dia` - Data-independent. Defaults to `dda`.

    Examples:
        Analysis of a single dataset
        >>> flippr.Study(lip = "PXD025926/FragPipe/LiP_LFQ")

        Include protein-level normalization from another dataset
        >>> flippr.Study(lip = "PXD025926/FragPipe/LiP_LFQ", trp = "PXD025926/FragPipe/TrP_LFQ")

    """

    def __init__(self, 
                 lip: str | Path, 
                 trp: Optional[str | Path] = None, 
                 method: str = "dda"
    ) -> None:
        """
        docstring
        """

        _lip, _trp, _method = _validate._validate_study(lip, trp, method)

        self.lip: Path = _lip
        self.trp: Optional[Path] = _trp
        self.method: str = _method
        self.processes: dict[str, _types.Process] = dict()
        self.results: dict[str, _types.Result] = dict()

    @property
    def samples(self) -> dict[str, set[str]]:
        """
        Returns the names of the samples or experimental annotations found in the FragPipe outputs for LiP and TrP (if include) datasets.
        
        """
        
        lip_samples = self._get_samples(self.lip)

        if self.trp is not None:
            trp_samples = self._get_samples(self.trp)

            return {"LiP": lip_samples, "TrP": trp_samples}

        return {"LiP": lip_samples}

    def _get_samples(self, path: Path) -> set[str]:
        """
        Returns the names of the samples from experimental annotations.

        """

        annot = _reader._read_experiment_annotation(path)

        samples = [info.get("Sample Name") for info in annot.values()]

        samples = list(filter(None, samples))

        samples = ["_".join(name.split("_")[:-1]) for name in samples]

        return set(samples)


    def add_process(
            self, 
            pid: str, 
            lip_ctrl: str, 
            lip_test: str, 
            n_rep: _types.replicate,
            trp_ctrl: Optional[str] = None, 
            trp_test: Optional[str] = None, 
            trp_n_rep: Optional[_types.replicate] = None
    ) -> None:
        """
        Adding a process to calculate the fold-change between test and control conditions.
        Normalization is optional and must include a TrP experiment in the `Study` set-up.

        Args:
            pid (str): Process ID used to index the `Study().results` dictionary after running `Study().run()`.
            lip_ctrl (str): Control condition sample name from the LiP experiment.
            lip_test (str): Test condition sample name from the LiP experiment.
            n_rep (int | tuple(int, int), tuple(tuple(int, ...), tuple(int, ...))): Number of replicates in the LiP experiment.
            trp_ctrl (str, optional): Control condition sample name from the TrP experiment.
            trp_test (str, optional): Test condition sample name from the TrP experiment.
            trp_n_rep (int, optional): Number of replicates in the TrP experiment.

        Examples:
            Simple fold-change calculation between two conditions
            >>> study.add_process("sample", "WT", "Drug", 3)

            Fold-change calculation with normalization
            >>> study.add_process("sample", "WT", "Drug", 3, "WT_TrP", "DMSO_TrP", 3)

            Add multiple processes to the same study
            >>> study.add_process("Lo_Dose", "WT", "Drug_Lo", 3, "WT_TrP", "DMSO_TrP", 3)
            >>> study.add_process("Hi_Dose", "WT", "Drug_Hi", 3, "WT_TrP", "DMSO_TrP", 3)
        
        """

        self.processes.update(
            {
                pid: _types.Process(
                    rcParams,
                    self.lip,
                    self.trp,
                    self.method,
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

    def run(self) -> dict[str, _types.Result]:
        """
        Run the processes added to the study.
        Global `Study` parameters can be changed by editing the `flippr.rcParams` dictionary.
        
        """

        self.results = {pid: proc.run() for pid, proc in self.processes.items()}

        return self.results
