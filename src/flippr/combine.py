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
            pl.col("CV").list.max()
            .alias("CV"),
            pl.col(fc).list.median()
            .alias(fc)
        )
    )

    combined = _log2(combined, fc)
    combined = _log10(combined, "P-value")
    combined = _log10(combined, "Adj. P-value")

    return combined


def summary_by(df: pl.DataFrame, by: str, fc: str, rcParams: dict) -> pl.DataFrame:

    prot_fc_sig = rcParams.get("protein.fc_sig_thresh", 1.0)
    prot_pv_sig = rcParams.get("protein.pval_sig_thresh", 0.01)
    prot_apv_sig = rcParams.get("protein.adj_pval_sig_thresh", 0.05)
    prot_hfc_sig = rcParams.get("protein.high_fc_sig_thresh", 6.0)
    prot_hpv_sig = rcParams.get("protein.high_pval_sig_thresh", 0.016)

    val = (
        df.group_by(by="Protein ID", maintain_order=True)
        .agg(pl.col("CV").filter(pl.col("-Log10 P-value").gt(0)).len())
        .rename({"CV": f"No. of Valid {by}"})
    )

    sig = (
        df.group_by(by="Protein ID", maintain_order=True)
        .agg(
            pl.col("CV")
            .filter(
                (
                    (pl.col(f"Log2 {fc}").abs().ge(prot_fc_sig))
                    & (pl.col("P-value").le(prot_pv_sig))
                )
                | (
                    (pl.col(f"Log2 {fc}").abs().ge(prot_hfc_sig))
                    & (pl.col("P-value").le(prot_hpv_sig))
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

                    (pl.col(f"Log2 {fc}").abs().ge(prot_fc_sig))
                    & (pl.col("Adj. P-value").le(prot_apv_sig))

            )
            .len()
        )
        .rename({"CV": f"No. of Significant {by} (Adj. P-value)"})
    )

    return val.join(sig, on="Protein ID").join(sigsig, on="Protein ID")
