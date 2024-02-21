import polars as pl

from pathlib import Path
from functools import cached_property

from . import validate as _validate
from . import combine as _combine
from . import functions as _functions
from .uniprot import metadata as _metadata
from .parameters import LFQ_FP_CONSTANT_ION_COLUMNS, LFQ_FP_VARIABLE_ION_COLUMNS
from .parameters import LFQ_FP_CONSTANT_PROTEIN_COLUMNS, LFQ_FP_VARIABLE_PROTEIN_COLUMNS
from .parameters import FLIPPR_PROTEIN_SUMMARY_COLUMNS


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

        _validate._validate_replicate_value(n_rep, "n_rep")

        if trp_n_rep is None:
            trp_n_rep = n_rep
        else:
            _validate._validate_replicate_value(trp_n_rep, "trp_n_rep")

        self._is_trp_norm = all(
            trp is not None for trp in [trp_path, trp_ctrl, trp_test]
        )

        if trp_path is not None:
            if trp_ctrl is not None:
                _validate._validate_sample_annotation(
                    trp_path, trp_ctrl, trp_n_rep, "trp_ctrl"
                )
            if trp_test is not None:
                _validate._validate_sample_annotation(
                    trp_path, trp_test, trp_n_rep, "trp_test"
                )

        _validate._validate_sample_annotation(lip_path, lip_ctrl, n_rep, "lip_ctrl")
        _validate._validate_sample_annotation(lip_path, lip_test, n_rep, "lip_test")

        self._fasta: pl.DataFrame | None = fasta
        self._pid: object = pid

        self._lip_path: Path = lip_path
        self._lip_ctrl_name: str = lip_ctrl
        self._lip_test_name: str = lip_test
        self._lip_ctrl: list[str] = [f"{lip_ctrl}_{n+1}" for n in range(n_rep)]
        self._lip_test: list[str] = [f"{lip_test}_{n+1}" for n in range(n_rep)]
        self._n_rep: int = n_rep

        if self._is_trp_norm:
            # TODO: Allow TrP normalization to occur with "MaxLFQ Intensity" or "Intensity" from `combined_protein.tsv`
            self._trp_path: Path | None = trp_path
            self._trp_ctrl: list[str] = [f"{trp_ctrl}_{n+1}" for n in range(trp_n_rep)]
            self._trp_test: list[str] = [f"{trp_test}_{n+1}" for n in range(trp_n_rep)]
            self._trp_n_rep: int = trp_n_rep

    @property
    def name(self) -> str:
        return f"{self._lip_ctrl_name}_v_{self._lip_test_name}"

    def run(self, **kwargs):
        result = Result(self, **kwargs)
        return result


class Result:
    """Organizes a FLiPPR Result"""

    def __init__(self, cls: Process, **kwargs) -> None:
        """doctstring"""

        self.name = cls.name
        # TODO: Make this able to switch between LFQ and SILAC
        
        self._lip_ctrl_cols = [
            f"{rep} {var}"
            for var in LFQ_FP_VARIABLE_ION_COLUMNS
            for rep in cls._lip_ctrl
        ]

        self._lip_ctrl_ints = [f"{rep} Intensity" for rep in cls._lip_ctrl]

        self._lip_test_cols = [
            f"{rep} {var}"
            for var in LFQ_FP_VARIABLE_ION_COLUMNS
            for rep in cls._lip_test
        ]

        self._lip_test_ints = [f"{rep} Intensity" for rep in cls._lip_test]
        self._lip_process_cols = (
            LFQ_FP_CONSTANT_ION_COLUMNS 
            + self._lip_ctrl_cols 
            + self._lip_test_cols
        )
        self._ions: pl.DataFrame = pl.read_csv(
            cls._lip_path.joinpath("combined_ion.tsv"), 
            separator="\t"
        )
        self._ions = self._ions.select(self._lip_process_cols)
        self._ions = _functions._cull_intensities(
            self._ions, 
            self._lip_ctrl_ints, 
            self._lip_test_ints, 
            cls._n_rep,
            kwargs.get("max_missing_values", 1)
        )
        self._ions = _functions._add_alt_hypothesis(
            self._ions, 
            self._lip_ctrl_ints, 
            self._lip_test_ints, 
            cls._n_rep
        )
        self._ions = _functions._impute_aon_intensities(
            self._ions,
            self._lip_ctrl_ints, 
            self._lip_test_ints,
            kwargs.get("aon_mean", 1e4),
            kwargs.get("aon_std", 1e3)
        )
        self._ions = _functions._add_start_end_aa(self._ions)
        self._ions = _functions._add_half_trpytic(self._ions)
        self._ions = _functions._add_cut_sites(self._ions)
        self._ions = _functions._add_ttest(
            self._ions, self._lip_ctrl_ints, self._lip_test_ints
        )
        self._ions = _functions._add_fdr(self._ions)
        self._ions = _functions._add_ratio(
            self._ions, self._lip_ctrl_ints, self._lip_test_ints
        )

        if (
            cls._is_trp_norm
            and cls._trp_path is not None
            and cls._trp_ctrl is not None
            and cls._trp_test is not None
            and cls._trp_n_rep is not None
        ):
            
            self._trp_ctrl_cols = [
                f"{rep} {var}"
                for var in LFQ_FP_VARIABLE_PROTEIN_COLUMNS
                for rep in cls._trp_ctrl
            ]

            self._trp_ctrl_ints = [f"{rep} MaxLFQ Intensity" for rep in cls._trp_ctrl]

            self._trp_test_cols = [
                f"{rep} {var}"
                for var in LFQ_FP_VARIABLE_PROTEIN_COLUMNS
                for rep in cls._trp_test
            ]

            self._trp_test_ints = [f"{rep} MaxLFQ Intensity" for rep in cls._trp_test]

            self._trp_process_cols = (
                LFQ_FP_CONSTANT_PROTEIN_COLUMNS
                + self._trp_ctrl_cols
                + self._trp_test_cols
            )

            self._norm_factors = pl.read_csv(
                cls._trp_path.joinpath("combined_protein.tsv"),
                separator="\t"
            )
            self._norm_factors.select(self._trp_process_cols)
            self._norm_factors = _functions._cull_intensities(
                self._norm_factors,
                self._trp_ctrl_ints,
                self._trp_test_ints,
                cls._trp_n_rep,
            )
            self._norm_factors = _functions._add_alt_hypothesis(
                self._norm_factors,
                self._trp_ctrl_ints,
                self._trp_test_ints,
                cls._trp_n_rep,
            )
            self._norm_factors = _functions._impute_aon_intensities(
                self._norm_factors,
                self._trp_ctrl_ints,
                self._trp_test_ints,
                kwargs.get("aon_mean", 1e4),
                kwargs.get("aon_std", 1e3)
            )
            self._norm_factors = _functions._add_ttest(
                self._norm_factors,
                self._trp_ctrl_ints,
                self._trp_test_ints
            )
            self._norm_factors = _functions._add_ratio(
                self._norm_factors,
                self._trp_ctrl_ints,
                self._trp_test_ints
            )
            self._ions = _functions._normalize_ratios(
                self._ions,
                self._norm_factors
            )

        self._ions = _functions._log2(
            self._ions, 
            "FC"
        )

        if cls._is_trp_norm:
            self._ions = _functions._log2(self._ions, "Normalized FC")

        self._ions = _functions._log10(self._ions, "P-value")
        self._ions = _functions._log10(self._ions, "Adj. P-value")
        self.ion = self._ions
        self._is_trp_norm = cls._is_trp_norm
        self._fc = "FC"

        if self._is_trp_norm:
            self._fc = "Normalized FC"

        self._fasta = cls._fasta

    def _add_metadata(self, df: pl.DataFrame):
        if self._fasta is not None:
            md = _metadata.bio_metadata(self._fasta)
            return df.join(md, on="Protein ID")
    
    @cached_property
    def modified_peptide(self):
        return _combine.combine_by(self.ion, by="MODIFIED PEPTIDE", fc=self._fc)

    @cached_property
    def peptide(self):
        return _combine.combine_by(self.ion, by="PEPTIDE", fc=self._fc)

    @cached_property
    def cut_site(self):
        return _combine.combine_by(self.ion, by="CUT SITE", fc=self._fc)

    @cached_property
    def protein_summary(self):
        self._proteins = self.ion.group_by(by="Protein ID", maintain_order=True).agg(
            pl.col(FLIPPR_PROTEIN_SUMMARY_COLUMNS).first()
        )

        __mod = _combine.summary_by(
            self.modified_peptide, by="Modified Peptides", fc=self._fc
        )

        __pep = _combine.summary_by(self.peptide, by="Peptides", fc=self._fc)

        __cut = _combine.summary_by(self.cut_site, by="Cut Sites", fc=self._fc)

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
