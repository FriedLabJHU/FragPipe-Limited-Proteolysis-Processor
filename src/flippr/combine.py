import numpy as np
import polars as pl
import scipy as sp

from .functions import _log2, _log10
from .parameters import _FLIPPR_COMBINE_KEY

COMB_NAME_COLUMN = {
    "CUT SITE": "Cut Site ID",
    "PEPTIDE": "Peptide Sequence",
    "MODIFIED PEPTIDE": "Modified Sequence",
}

def combine_by(df: pl.DataFrame, by: str, fc: str) -> pl.DataFrame:

    combined = (
        df.group_by(by=["Protein ID", COMB_NAME_COLUMN[by]], maintain_order=True)
        .agg(
            pl.col(_FLIPPR_COMBINE_KEY[by]).first(),
            pl.col(["P-value", "Adj. P-value", "CV", fc])
            .filter(
                pl.col("T-test").sign() == pl.col("T-test").sign().sum().sign()
            )
        )
        .with_columns(
            pl.when(pl.col("P-value").list.len() > 0)
            .then(pl.col("P-value"))
            .otherwise([1.0])
            .alias("P-value"),
            pl.when(pl.col("Adj. P-value").list.len() > 0)
            .then(pl.col("Adj. P-value"))
            .otherwise([1.0])
            .alias("Adj. P-value"),
            pl.when(pl.col("CV").list.len() > 0)
            .then(pl.col("CV"))
            .otherwise([0.0])
            .alias("CV"),
            pl.when(pl.col(fc).list.len() > 0)
            .then(pl.col(fc))
            .otherwise([0.0])
            .alias(fc)
        )
        .with_columns(
            pl.col("P-value")
            .map_elements(lambda x: sp.stats.combine_pvalues(x)[1], return_dtype=pl.Float64)
            .alias("P-value"),
            pl.col("Adj. P-value")
            .map_elements(lambda x: sp.stats.combine_pvalues(x)[1], return_dtype=pl.Float64)
            .alias("Adj. P-value"),
            pl.col("CV")
            .map_elements(lambda x: np.nanmax(x), return_dtype=pl.Float64)
            .alias("CV"),
            pl.col(fc)
            .map_elements(lambda x: np.nanmedian(x), return_dtype=pl.Float64)
            .alias(fc)
        )
    )

    combined = _log2(combined, fc)
    combined = _log10(combined, "P-value")
    combined = _log10(combined, "Adj. P-value")

    return combined


def summary_by(df: pl.DataFrame, by: str, fc: str) -> pl.DataFrame:
    val = (
        df.group_by(by="Protein ID", maintain_order=True)
        .agg(pl.col("CV").filter(pl.col("-Log10 P-value") > 0).len())
        .rename({"CV": f"No. of Valid {by}"})
    )

    sig = (
        df.group_by(by="Protein ID", maintain_order=True)
        .agg(
            pl.col("CV")
            .filter(
                (
                    (pl.col(f"Log2 {fc}").abs() >= 1)
                    & (pl.col("-Log10 P-value") >= -np.log10(0.01))
                )
                | (
                    (pl.col(f"Log2 {fc}").abs() >= 6)
                    & (pl.col("-Log10 P-value") >= 1.8)
                )
            )
            .len()
        )
        .rename({"CV": f"No. of Significant {by} (P-value)"})
    )

    sigsig = (
        df.group_by(by="Protein ID", maintain_order=True)
        .agg(
            pl.col("CV")
            .filter(

                    (pl.col(f"Log2 {fc}").abs() >= 1)
                    & (pl.col("-Log10 Adj. P-value") >= -np.log10(0.05))

            )
            .len()
        )
        .rename({"CV": f"No. of Significant {by} (Adj. P-value)"})
    )

    return val.join(sig, on="Protein ID").join(sigsig, on="Protein ID")
