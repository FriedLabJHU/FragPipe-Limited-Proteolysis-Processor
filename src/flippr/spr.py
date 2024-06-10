from __future__ import annotations

from functools import cached_property
from pathlib import Path

import polars as pl

from . import combine as _combine
from . import functions as _functions
from .parameters import (
    _FLIPPR_PROTEIN_SUMMARY_COLUMNS,
    _LFQ_FP_CONSTANT_ION_COLUMNS,
    _LFQ_FP_CONSTANT_PROTEIN_COLUMNS,
    _LFQ_FP_VARIABLE_ION_COLUMNS,
    _LFQ_FP_VARIABLE_PROTEIN_COLUMNS,
)
from .uniprot import metadata as _metadata


class Process:
    """Organizes individual process runs"""

    def __init__(
        self,
        fasta: pl.DataFrame,
        lip_path: Path,
        trp_path: Path | None,
        pid: object | None,
        lip_ctrl: str,
        lip_test: str,
        n_rep: int,
        trp_ctrl: str | None = None,
        trp_test: str | None = None,
        trp_n_rep: int | None = None,
    ) -> None:
        """docstring"""

        self._is_trp_norm = all(
            trp is not None for trp in [trp_path, trp_ctrl, trp_test, trp_n_rep]
        )

        self._fasta: pl.DataFrame | None = fasta
        self._pid: object = pid

        self._lip_path: Path = lip_path
        self._lip_ctrl_name: str = lip_ctrl
        self._lip_test_name: str = lip_test
        self._lip_ctrl: list[str] = [f"{lip_ctrl}_{n+1}" for n in range(n_rep)]
        self._lip_test: list[str] = [f"{lip_test}_{n+1}" for n in range(n_rep)]
        self._n_rep: int = n_rep
        self._trp_path: Path | None = None
        self._trp_ctrl_name: str | None = None
        self._trp_test_name: str | None = None
        self._trp_ctrl: list[str] | None = None
        self._trp_test: list[str] | None = None
        self._trp_n_rep: int | None = None

        if self._is_trp_norm:
            self._trp_path = trp_path
            self._trp_ctrl_name = trp_ctrl
            self._trp_test_name = trp_test
            self._trp_ctrl = [f"{trp_ctrl}_{n+1}" for n in range(trp_n_rep)]
            self._trp_test = [f"{trp_test}_{n+1}" for n in range(trp_n_rep)]
            self._trp_n_rep = trp_n_rep

    @property
    def name(self) -> str:
        return f"{self._lip_ctrl_name}_v_{self._lip_test_name}"

    def run(self, rcParams: dict):
        result = Result(self, rcParams)
        return result


class Result:
    """Organizes a FLiPPR Result"""

    def __init__(self, cls: Process, rcParams: dict) -> None:
        """doctstring"""

        self.rcParams: dict = rcParams
        self._fasta = cls._fasta

        self.name: str = cls.name
        self._is_trp_norm: bool = cls._is_trp_norm

        self._fc = "FC"
        self._n_rep = cls._n_rep
        self._lip_path: Path = cls._lip_path
        self._trp_path: Path = cls._trp_path
        self._lip_ctrl_name: str = cls._lip_ctrl_name
        self._lip_test_name: str = cls._lip_test_name
        self._lip_ctrl_ints: list[str] = [f"{rep} Intensity" for rep in cls._lip_ctrl]
        self._lip_test_ints: list[str] = [f"{rep} Intensity" for rep in cls._lip_test]
        self._trp_n_rep: int | None = None
        self._trp_ctrl_name: str | None = None
        self._trp_test_name: str | None = None
        self._trp_ctrl_ints: list[str] | None = None
        self._trp_test_ints: list[str] | None = None
        self._trp_ctrl_cols: list[str] | None = None
        self._trp_test_cols: list[str] | None = None
        self._trp_process_cols: list[str] | None = None
        self._norm_factors: pl.DataFrame | None = None

        self._lip_ctrl_cols: list[str] = [
            f"{rep} {var}"
            for var in _LFQ_FP_VARIABLE_ION_COLUMNS
            for rep in cls._lip_ctrl
        ]

        self._lip_test_cols: list[str]  = [
            f"{rep} {var}"
            for var in _LFQ_FP_VARIABLE_ION_COLUMNS
            for rep in cls._lip_test
        ]

        self._lip_process_cols: list[str]  = (
            _LFQ_FP_CONSTANT_ION_COLUMNS
            + self._lip_ctrl_cols
            + self._lip_test_cols
        )

        args = {"ctrl_name":self._lip_ctrl_name, "test_name":self._lip_test_name,
                "ctrl_ints":self._lip_ctrl_ints, "test_ints": self._lip_test_ints,
                "n_rep": self._n_rep,
                "rcParams": self.rcParams}

        self._ion = self._read_fragpipe_tsv(self._lip_path, self._lip_process_cols, "combined_ion.tsv")
        self._ion = self.run(self._ion, args)
        self._ion = self.clean_up(self._ion, args)

        if self._is_trp_norm:
            trp_prot_int_val = self.rcParams.get("trp_protein.intensity_value", "MaxLFQ Intensity")

            self._trp_n_rep = cls._trp_n_rep
            self._trp_ctrl_name = cls._trp_ctrl_name
            self._trp_test_name = cls._trp_test_name
            self._trp_ctrl_ints = [f"{rep} {trp_prot_int_val}" for rep in cls._trp_ctrl]
            self._trp_test_ints = [f"{rep} {trp_prot_int_val}" for rep in cls._trp_test]

            self._trp_ctrl_cols = [
                f"{rep} {var}"
                for var in _LFQ_FP_VARIABLE_PROTEIN_COLUMNS
                for rep in cls._trp_ctrl
            ]

            self._trp_test_cols = [
                f"{rep} {var}"
                for var in _LFQ_FP_VARIABLE_PROTEIN_COLUMNS
                for rep in cls._trp_test
            ]

            self._trp_process_cols = (
                _LFQ_FP_CONSTANT_PROTEIN_COLUMNS
                + self._trp_ctrl_cols
                + self._trp_test_cols
            )

            trp_args = {"ctrl_name":self._trp_ctrl_name, "test_name":self._trp_test_name,
                    "ctrl_ints":self._trp_ctrl_ints, "test_ints": self._trp_test_ints,
                    "n_rep": self._trp_n_rep,
                    "rcParams": self.rcParams}

            self._trp_norm = self._read_fragpipe_tsv(self._trp_path, self._trp_process_cols, "combined_protein.tsv")
            self._trp_norm = self.run(self._trp_norm, trp_args)
            self._ion = _functions._normalize_ratios(self._ion, self._trp_norm, rcParams)
            self._fc = "Normalized FC" # Generated after running `._normalize_ratios()`


    def run(self, df: pl.DataFrame, args: dict) -> pl.DataFrame:
        # Can be performed on ion, mod_pep, pep, or protein
        df = _functions._cull_intensities(df, **args)
        df = _functions._add_alt_hypothesis(df, **args)
        df = _functions._impute_aon_intensities(df, **args)
        df = _functions._add_ttest(df, **args)
        df = _functions._add_fdr(df, **args)
        df = _functions._add_ratio(df, **args)
        df = _functions._log2(df, self._fc)
        df = _functions._log10(df, "P-value")
        df = _functions._log10(df, "Adj. P-value")

        return df

    def clean_up(self, df: pl.DataFrame, args: dict) -> pl.DataFrame:
        # Only meant to be performed on the lip ions
        df = _functions._add_start_end_aa(df, **args)
        df = _functions._add_half_trpytic(df, **args)
        df = _functions._add_cut_sites(df, **args)

        return df

    def _read_fragpipe_tsv(self, path: Path, cols: list[str], filename: str):
        df: pl.DataFrame = pl.read_csv(path.joinpath(filename), separator="\t").select(pl.col(cols))
        return df

    def _add_metadata(self, df: pl.DataFrame) -> pl.DataFrame:
        if self._fasta is not None:
            md = _metadata.bio_metadata(self._fasta)
            return df.join(md, on="Protein ID")


    @property
    def ion(self) -> pl.DataFrame:
        return self._ion
    
    @property
    def trp_protein(self) -> pl.DataFrame:
        return self._trp_norm

    @cached_property
    def modified_peptide(self) -> pl.DataFrame:
        return _combine.combine_by(self.ion, by="MODIFIED PEPTIDE", fc=self._fc)

    @cached_property
    def peptide(self) -> pl.DataFrame:
        return _combine.combine_by(self.ion, by="PEPTIDE", fc=self._fc)

    @cached_property
    def cut_site(self) -> pl.DataFrame:
        return _combine.combine_by(self.ion, by="CUT SITE", fc=self._fc)

    @cached_property
    def protein_summary(self) -> pl.DataFrame:
        self._proteins = self.ion.group_by(by="Protein ID", maintain_order=True).agg(
            pl.col(_FLIPPR_PROTEIN_SUMMARY_COLUMNS).first()
        )

        __mod = _combine.summary_by(self.modified_peptide, by="Modified Peptides", fc=self._fc, rcParams=self.rcParams)

        __pep = _combine.summary_by(self.peptide, by="Peptides", fc=self._fc, rcParams=self.rcParams)

        __cut = _combine.summary_by(self.cut_site, by="Cut Sites", fc=self._fc, rcParams=self.rcParams)

        self._proteins = (
            self._proteins
            .join(__mod, on="Protein ID")
            .join(__pep, on="Protein ID")
            .join(__cut, on="Protein ID")
        )

        if self._fasta is not None:
            return self._add_metadata(self._proteins)
        else:
            return self._proteins
