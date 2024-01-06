import polars as pl

from pathlib import Path
from functools import cached_property

from . import flippr as _flippr
from . import validate as val
from .uniprot import fasta as faa


# TODO TOC
# Metadata integrations
# Saving feature
# Adjustable parameters for protein summary
# Imputation schemes (much later)

class Study:
    """Organizes a FLiPPR Study"""

    def __init__(
        self,
        lip: str | Path,
        trp: str | Path | None = None,
        fasta: str | Path | None = None,
        method: str = "lfq",
    ):
        """ docstring """
        
        _lip, _trp, _fasta, _method = \
            val._validate_study(lip, trp, fasta, method)

        self.lip: Path = _lip
        self.trp: Path | None = _trp
        self._fasta: Path | None = _fasta
        self.method: str = _method
        
        self._pid: int = 0
        self.processes: dict = dict()
        self.results: dict = dict()


    @property
    def samples(self) -> dict:
        if self.trp_samples is not None:
            return {"LiP" : self.lip_samples, "TrP" : self.trp_samples}
        
        return {"LiP" : self.lip_samples}

    @property
    def lip_samples(self) -> set:
        return self._read_annotation(self.lip, "LiP")

    @property
    def trp_samples(self) -> set | None:
        if self.trp is not None:
            return self._read_annotation(self.trp, "TrP")
        
        return None
    
    @cached_property
    def fasta(self) -> pl.DataFrame | None:
        if self._fasta is not None:
            return faa._read_fasta(self._fasta)
        
        return None
    
    def add_process(
        self,
        pid: object | None,
        lip_ctrl: str,
        lip_test: str,
        n_rep: int,
        trp_ctrl: str | None = None,
        trp_test: str | None = None,
        trp_n_rep: int | None = None,
    ) -> None:
        if pid is None:
            pid = self._pid
            self._pid += 1

        self.processes.update({pid:_flippr.Process(self.lip, self.trp, pid, lip_ctrl, lip_test, n_rep, trp_ctrl, trp_test, trp_n_rep)})

    def run(self) -> dict[str, _flippr.Result]:
        
        self.results = {pid: proc.run() for pid, proc in self.processes.items()}

        return self.results


    def _read_annotation(self, liptrp: Path, label: str) -> set:

        def __strip_fragpipe_repilicate(sample: str) -> str:
            return "_".join([_ for _ in sample.split("_")[:-1]])
        
        annotation = liptrp.joinpath("experiment_annotation.tsv")

        annotation_df = pl.read_csv(annotation, separator="\t")

        annotation_dict = annotation_df.select(pl.col("sample")).to_dict(as_series = False)
        
        return set(__strip_fragpipe_repilicate(sample) for sample in annotation_dict["sample"])
    