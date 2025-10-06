import polars as pl
import polars.selectors as cs
from typing import Any
from pathlib import Path

from .parameters import (
    _DIA_FP_CONSTANT_ION_COLUMNS,
    _DIA_RENAME_FP_ION,
    _DIA_RENAME_DIANN_ION,
    _DIA_RENAME_DIANN_PROTEIN
)


def _read_ion(path: Path, method: str) -> pl.DataFrame:
    match method:
        case "dda":
            dda_ion_df = pl.read_csv(path.joinpath("combined_ion.tsv"), separator="\t")
            dda_ion_df = dda_ion_df.with_columns(cs.contains("Intensity").fill_null(strategy="zero"))

            return dda_ion_df
        
        case "dia":
            annot = _read_experiment_annotation(path)

            dia_ion_df = pl.read_csv(path.joinpath("dia-quant-output/report.pr_matrix.tsv"), separator="\t")
            fp_ion_df = pl.read_csv(path.joinpath("ion.tsv"), separator="\t").select(_DIA_FP_CONSTANT_ION_COLUMNS)

            dia_ion_df = _rename_dia_columns(dia_ion_df, annot)
            dia_ion_df = _add_dia_ion_data(dia_ion_df, fp_ion_df)
            dia_ion_df = dia_ion_df.with_columns(cs.contains("Intensity").fill_null(strategy="zero"))

            return dia_ion_df
        
        case _:
            raise ValueError("Input error.")
        
def _read_trp(path: Path, method: str) -> pl.DataFrame:
    match method:
        case "dda":
            dda_trp_df = pl.read_csv(path.joinpath("combined_protein.tsv"), separator="\t")
            dda_trp_df = dda_trp_df.with_columns(cs.contains("Intensity").fill_null(strategy="zero"))

            return dda_trp_df
        
        case "dia":
            annot = _read_experiment_annotation(path)

            dia_trp_df = pl.read_csv(path.joinpath("dia-quant-output/report.pg_matrix.tsv"), separator="\t")

            dia_trp_df = _rename_dia_columns(dia_trp_df, annot, "trp")
            dia_trp_df = dia_trp_df.with_columns(cs.contains("Intensity").fill_null(strategy="zero"))

            return dia_trp_df
        
        case _:
            raise ValueError("Input error.")

def _read_experiment_annotation(path: Path) -> dict[str, dict[str, str]]:
    df = pl.read_csv(
        path.joinpath("experiment_annotation.tsv"), 
        separator="\t"
    )

    annot: dict[str, dict[str, str]] = {}
    for row in df.iter_rows():
        file, sample, sample_name, condition, replicate = row
        annot.update({
            file.split(".").pop(0): # will fail if someone tries to use the same file as multiple replicates
            {
                "Sample": sample,
                "Sample Name": sample_name,
                "Condition": condition,
                "Replicate": replicate
            },
        })

    return annot

def _file_to_sample_name(annot: dict[str, dict[str, str]]) -> dict[str, str]:
    return {file: data.get("Sample Name", "") for file, data in annot.items()}

def _rename_dia_columns(df: pl.DataFrame, annot: dict[str, dict[str, str]], data_type: str = "ion") -> pl.DataFrame:
    other_cols = _DIA_RENAME_DIANN_ION
    if data_type != "ion":
        other_cols =_DIA_RENAME_DIANN_PROTEIN
    
    cols = df.columns

    rename = {}
    for file, sample in _file_to_sample_name(annot).items():
        for col in cols:
            if file in col:
                rename.update({col: f"{sample} Intensity"})

    return df.rename(rename).rename(other_cols)

def _add_dia_ion_data(dia_df: pl.DataFrame, ion_df: pl.DataFrame) -> pl.DataFrame:
        dia_df = dia_df.with_columns(
            (pl.col("Protein ID")+pl.col("Peptide Sequence")).alias("Unique ID")
        )

        ion_df = ion_df.with_columns(
            (pl.col("Protein ID")+pl.col("Peptide Sequence")).alias("Unique ID")
        ).drop("Protein ID", "Peptide Sequence").unique("Unique ID")

        ion_df = ion_df.rename(_DIA_RENAME_FP_ION)

        return dia_df.join(ion_df, on="Unique ID", how="left")
