from __future__ import annotations

import polars as pl
from pathlib import Path
from typing import Optional, Any, cast
from functools import cached_property

from . import combine as _combine
from . import functions as _functions
from . import validate as _validate
from . import reader as _reader
from .parameters import (
    _FLIPPR_ION_COLUMNS,
    _FLIPPR_PROTEIN_COLUMNS,
    _FLIPPR_PROTEIN_SUMMARY_COLUMNS,
)

type replicate = int | tuple[int, int] | tuple[tuple[int, ...], tuple[int, ...]]

class Process:
    """Organizes individual process runs"""

    def __init__(
        self,
        rcParams: dict[str, Any],
        lip_path: Path,
        trp_path: Optional[Path],
        method: str,
        pid: str,
        lip_ctrl: str,
        lip_test: str,
        n_rep: replicate,
        trp_ctrl: Optional[str] = None,
        trp_test: Optional[str] = None,
        trp_n_rep: Optional[replicate] = None,
    ) -> None:
        """docstring"""

        self._rcParams: dict[str, Any] = rcParams
        self._method: str = method
        self._pid: str = pid

        self._lip_path: Path = lip_path
        self._lip_ctrl_name: str = lip_ctrl
        self._lip_test_name: str = lip_test
        self._lip_ctrl_rep_list: list[str]
        self._lip_test_rep_list: list[str]
        self._lip_ctrl_n_rep: int
        self._lip_test_n_rep: int
        (
            self._lip_ctrl_rep_list, 
            self._lip_test_rep_list, 
            self._lip_ctrl_n_rep, 
            self._lip_test_n_rep
        ) = self._create_replicate_variables(lip_ctrl,  lip_test, n_rep)
        
        self._trp_path: Optional[Path] = None
        self._trp_ctrl_name: Optional[str] = None
        self._trp_test_name: Optional[str] = None
        self._trp_ctrl_rep_list: Optional[list[str]] = None
        self._trp_test_rep_list: Optional[list[str]] = None
        self._trp_ctrl_n_rep: Optional[int] = None
        self._trp_test_n_rep: Optional[int] = None

        self._is_trp_norm = all(trp is not None for trp in [trp_path, trp_ctrl, trp_test, trp_n_rep])
        if self._is_trp_norm:
            assert trp_path is not None
            assert trp_ctrl is not None
            assert trp_test is not None
            assert trp_n_rep is not None

            self._trp_path = trp_path
            self._trp_ctrl_name = trp_ctrl
            self._trp_test_name = trp_test
            (
                self._trp_ctrl_rep_list, 
                self._trp_test_rep_list, 
                self._trp_ctrl_n_rep, 
                self._trp_test_n_rep
            ) = self._create_replicate_variables(trp_ctrl,  trp_test, trp_n_rep)

    def run(self):
        return Result(self)
    
    def _create_replicate_variables(self, ctrl: str, test: str, rep: replicate) -> tuple[list[str], list[str], int, int]:
        ctrl_rep_list: list[str]
        test_rep_list: list[str]
        ctrl_n_rep: int
        test_n_rep: int

        match _validate._validate_replicate(rep):
            case "int":
                rep = cast(int, rep)

                ctrl_rep_list = [f"{ctrl}_{i+1}" for i in range(rep)]
                test_rep_list = [f"{test}_{i+1}" for i in range(rep)]
                ctrl_n_rep = rep
                test_n_rep = rep

            case "tuple":
                rep = cast(tuple[int, int], rep)

                ctrl_rep_list = [f"{ctrl}_{i+1}" for i in range(rep[0])]
                test_rep_list = [f"{test}_{i+1}" for i in range(rep[1])]
                ctrl_n_rep = rep[0]
                test_n_rep = rep[1]

            case "tuple_tuple":
                rep = cast(tuple[tuple[int, ...], tuple[int, ...]], rep)

                ctrl_rep_list = [f"{ctrl}_{i}" for i in rep[0]]
                test_rep_list = [f"{test}_{i}" for i in rep[1]]
                ctrl_n_rep = len(rep[0])
                test_n_rep = len(rep[1])

            case _:
                raise ValueError("Input error.")
            
        return ctrl_rep_list, test_rep_list, ctrl_n_rep, test_n_rep
    
    def _ion_intensity_cols(self, rep_list: list[str]) -> list[str]:
        return [f"{rep} Intensity" for rep in rep_list]
    
    def _trp_intensity_cols(self, rep_list: list[str]) -> list[str]:
        match self._method:
            case "dda":
                trp_prot_int_val = self._rcParams.get(
                    "trp_protein.intensity_value",
                    "MaxLFQ Intensity"
                )
            case "dia":
                trp_prot_int_val = "Intensity"
            case _:
                raise ValueError("Input error.")

        return [f"{rep} {trp_prot_int_val}" for rep in rep_list]
    
    @cached_property
    def _ctrl_ion_int_cols(self) -> list[str]:
        return self._ion_intensity_cols(self._lip_ctrl_rep_list)
    
    @cached_property
    def _test_ion_int_cols(self) -> list[str]:
        return self._ion_intensity_cols(self._lip_test_rep_list)
    
    @cached_property
    def _ctrl_trp_int_cols(self) -> list[str]:
        assert self._trp_ctrl_rep_list is not None
        return self._trp_intensity_cols(self._trp_ctrl_rep_list)
    
    @cached_property
    def _test_trp_int_cols(self) -> list[str]:
        assert self._trp_test_rep_list is not None
        return self._trp_intensity_cols(self._trp_test_rep_list)

    @cached_property
    def _ion_columns(self) -> list[str]:
        return _FLIPPR_ION_COLUMNS + self._ctrl_ion_int_cols + self._test_ion_int_cols

    @cached_property
    def _trp_columns(self) -> list[str]:
        return _FLIPPR_PROTEIN_COLUMNS + self._ctrl_trp_int_cols + self._test_trp_int_cols

class Result:
    """Organizes a FLiPPR Result"""

    def __init__(self, cls: Process,) -> None:
        """doctstring"""

        self._fc: str = "FC"
        self._rcParams: dict[str, Any] = cls._rcParams
        
        self.trp_args: Optional[dict[str, Any]] = None
        self._trp_norm: Optional[pl.DataFrame] = None

        self.args: dict[str, Any] = {
            "ctrl_name":    cls._lip_ctrl_name, 
            "test_name":    cls._lip_test_name,
            "ctrl_ints":    cls._ctrl_ion_int_cols, 
            "test_ints":    cls._test_ion_int_cols,
            "ctrl_n_rep":   cls._lip_ctrl_n_rep,
            "test_n_rep":   cls._lip_test_n_rep,
            "rcParams":     cls._rcParams
        }
        
        self._ion = _reader._read_ion(cls._lip_path, cls._method)
        self._ion = self._ion.select(cls._ion_columns)
        self._ion = self.run(self._ion, self.args)
        self._ion = self.clean_up(self._ion, self.args)

        if cls._is_trp_norm:
            assert cls._trp_path is not None

            self.trp_args = {
                "ctrl_name":    cls._trp_ctrl_name,
                "test_name":    cls._trp_test_name,
                "ctrl_ints":    cls._ctrl_trp_int_cols,
                "test_ints":    cls._test_trp_int_cols,
                "ctrl_n_rep":   cls._trp_ctrl_n_rep,
                "test_n_rep":   cls._trp_test_n_rep,
                "rcParams":     cls._rcParams
            }

            self._trp_norm = _reader._read_trp(cls._trp_path, cls._method)
            self._trp_norm = self.run(self._trp_norm, self.trp_args)
            self._ion = _functions._normalize_ratios(self._ion, self._trp_norm, cls._rcParams)
            self._fc = "Normalized FC" # Generated after running `._normalize_ratios()`
            self._ion = _functions._log2(self._ion, self._fc)
            


    def run(self, df: pl.DataFrame, args: dict) -> pl.DataFrame:
        # Can be performed on ion, mod_pep, pep, or protein
        df = _functions._cull_intensities(df, **args)
        df = _functions._add_alt_hypothesis(df, **args)
        df = _functions._impute_aon_intensities(df, **args)
        df = _functions._add_ttest(df, **args)
        df = _functions._add_fdr(df, **args)
        df = _functions._add_ratio(df, **args)
        df = _functions._log2(df, self._fc)
        df = _functions._neg_log10(df, "P-value")
        df = _functions._neg_log10(df, "Adj. P-value")

        return df

    def clean_up(self, df: pl.DataFrame, args: dict) -> pl.DataFrame:
        # Only meant to be performed on the lip ions
        df = _functions._add_start_end_aa(df, **args)
        df = _functions._add_half_trpytic(df, **args)
        df = _functions._add_cut_sites(df, **args)

        return df

    
    @property
    def name(self) -> str:
        """
            Returns the name of the process from the given LiP control and test conditions.
            
        """
        return self.args.get("ctrl_name", "ctrl") + "_v_" + self.args.get("test_name", "test")
    
    @property
    def ion(self) -> pl.DataFrame:
        """
            ion dataframe
        """
        return self._ion
    
    @ion.setter
    def ion(self, ion: pl.DataFrame) -> None:
        """
            ion setter
        """
        self._ion = ion
    
    @property
    def trp_protein(self) -> pl.DataFrame | None:
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
        self._proteins = self.ion.group_by("Protein ID", maintain_order=True).agg(
            pl.col(_FLIPPR_PROTEIN_SUMMARY_COLUMNS).first()
        )

        __mod = _combine.summary_by(self.modified_peptide, by="Modified Peptides", fc=self._fc, rcParams=self._rcParams)

        __pep = _combine.summary_by(self.peptide, by="Peptides", fc=self._fc, rcParams=self._rcParams)

        __cut = _combine.summary_by(self.cut_site, by="Cut Sites", fc=self._fc, rcParams=self._rcParams)

        self._proteins = (
            self._proteins
            .join(__mod, on="Protein ID")
            .join(__pep, on="Protein ID")
            .join(__cut, on="Protein ID")
        )

        return self._proteins
