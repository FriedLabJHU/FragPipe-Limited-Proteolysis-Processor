import polars as pl
import scipy as sp
import numpy as np

from .functions import _log2, _log10
from .parameters import FLIPPR_COMBINE_KEY

COMB_NAME_COLUMN = {
    "CUT SITE": "Cut Site ID",
    "PEPTIDE": "Peptide Sequence",
    "MODIFIED PEPTIDE": "Modified Sequence"
}

def combine_by(df: pl.DataFrame, by: str, fc: str) -> pl.DataFrame:
    meta = \
    df.group_by(by = "Protein ID", maintain_order=True).agg(
    pl.col(FLIPPR_COMBINE_KEY[by]).first())

    data = \
    df.group_by(by = ["Protein ID", COMB_NAME_COLUMN[by]], maintain_order=True).agg(
        pl.col(["P-value", "Adj. P-value", "CV", fc])
        .filter(
            pl.col("T-test").sign() == pl.col("T-test").sign().sum().sign()
        )
    ).with_columns(
        pl.when(pl.col("P-value").list.len() > 0)
        .then(pl.col("P-value"))
        .otherwise(pl.lit([1.0]))
        .alias("P-value"),

        pl.when(pl.col("Adj. P-value").list.len() > 0)
        .then(pl.col("Adj. P-value"))
        .otherwise(pl.lit([1.0]))
        .alias("Adj. P-value"),

        pl.when(pl.col("CV").list.len() > 0)
        .then(pl.col("CV"))
        .otherwise(pl.lit([0.0]))
        .alias("CV"),

        pl.when(pl.col(fc).list.len() > 0)
        .then(pl.col(fc))
        .otherwise(pl.lit([0.0]))
        .alias(fc),
    ).select(
        pl.col("Protein ID"),

        pl.col("P-value").map_elements(lambda x: sp.stats.combine_pvalues(x)[1])
        .alias("P-value"),

        pl.col("Adj. P-value").map_elements(lambda x: sp.stats.combine_pvalues(x)[1])
        .alias("Adj. P-value"),

        pl.col("CV").map_elements(lambda x: np.nanmax(x))
        .alias("CV"),

        pl.col(fc).map_elements(lambda x: np.nanmedian(x))
        .alias(fc)
    )

    combined = meta.join(data, left_on = "Protein ID", right_on = "Protein ID")
    
    combined = _log2(combined, fc)
    combined = _log10(combined, "P-value")
    combined = _log10(combined, "Adj. P-value")

    return combined

def summary_by(df: pl.DataFrame, by: str, fc: str) -> pl.DataFrame:
    val = df.group_by(by = "Protein ID", maintain_order=True).agg(pl.col("CV").filter(
        (pl.col("-Log10 P-value") > 0)
    ).len()).rename({"CV": f"No. of Valid {by}"})

    sig = df.group_by(by = "Protein ID", maintain_order=True).agg(pl.col("CV").filter(
        ((pl.col(f"Log2 {fc}").abs() > 1) & (pl.col("-Log10 P-value") > -np.log10(0.01))) |
        ((pl.col(f"Log2 {fc}").abs() > 6) & (pl.col("-Log10 P-value") > 1.8))
    ).len()).rename({"CV": f"No. of Significant {by} (P-value)"})

    sigsig = df.group_by(by = "Protein ID", maintain_order=True).agg(pl.col("CV").filter(
        ((pl.col(f"Log2 {fc}").abs() > 1) & (pl.col("-Log10 Adj. P-value") > -np.log10(0.05)))
    ).len()).rename({"CV": f"No. of Significant {by} (Adj. P-value)"})

    return val.join(sig, on = "Protein ID").join(sigsig, on = "Protein ID")
