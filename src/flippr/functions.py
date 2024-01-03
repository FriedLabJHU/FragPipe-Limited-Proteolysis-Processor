import polars as pl
import numpy as np
import scipy as sp

def _cull_intensities(df: pl.DataFrame, ctrl: list[str], test: list[str], n_rep: int) -> pl.DataFrame:
    c_nz = df.select(ctrl).map_rows(np.count_nonzero).to_series()
    t_nz = df.select(test).map_rows(np.count_nonzero).to_series()

    # TODO: Allow the tolerance of `-1` to be changed
    df = \
    df.filter(
        ((c_nz == n_rep) & (t_nz >= n_rep-1)) | \
        ((c_nz >= n_rep-1) & (t_nz == n_rep)) | \
        ((c_nz == n_rep) & (t_nz == 0)) | \
        ((c_nz == 0) & (t_nz == n_rep))
    )

    df = __null_to_nan(df, ctrl, test)

    return df

def __null_to_nan(df: pl.DataFrame, ctrl: list[str], test: list[str]) -> pl.DataFrame:
    df = \
    df.with_columns(
        [
            pl.col(ctrl)
            .map_elements(lambda x: x if x else np.nan),
            pl.col(test)
            .map_elements(lambda x: x if x else np.nan)
        ]
    )
    
    return df


def _add_alt_hypothesis(df: pl.DataFrame, ctrl: list[str], test: list[str], n_rep: int) -> pl.DataFrame:
    c_nz = df.select(ctrl).map_rows(np.count_nonzero).to_series()
    t_nz = df.select(test).map_rows(np.count_nonzero).to_series()

    df = \
    df.with_columns(
        pl.when((c_nz == n_rep) & (t_nz == 0))
        .then(pl.lit("greater"))
        .when((c_nz == 0) & (t_nz == n_rep))
        .then(pl.lit("less"))
        .otherwise(pl.lit("two-sided"))
        .alias("Alternative Hypothesis")
    )

    return df

def _impute_aon_intensities(df: pl.DataFrame, ctrl: list[str], test: list[str]) -> pl.DataFrame:
    df = \
    df.with_columns(
        [
            pl.when(
                pl.all_horizontal(pl.col(c_i).is_nan() for c_i in ctrl) & 
                pl.all_horizontal(pl.col(t_i).is_not_nan() for t_i in test)
            )
            .then(pl.lit(sp.stats.norm.rvs(1e4, 1e3, 1)))
            .otherwise(pl.col(C_I))
            .alias(C_I) 
        for C_I in ctrl
        ]
    ).with_columns(
        [
            pl.when(
                pl.all_horizontal(pl.col(c_i).is_not_nan() for c_i in ctrl) & 
                pl.all_horizontal(pl.col(t_i).is_nan() for t_i in test)
            )
            .then(pl.lit(sp.stats.norm.rvs(1e4, 1e3, 1)))
            .otherwise(pl.col(T_I))
            .alias(T_I) 
        for T_I in test
        ]
    )

    return df

def _add_ttest(df: pl.DataFrame, ctrl: list[str], test: list[str]) -> pl.DataFrame:
    
    def drop_nans(x: list) -> np.ndarray:
        arr_x = np.array(x)
        return arr_x[~np.isnan(arr_x)]
    
    ttest = \
    df.select(
        pl.concat_list(pl.col(ctrl)),
        pl.concat_list(pl.col(test)),
        pl.col("Alternative Hypothesis")
    ).map_rows(lambda x:
        sp.stats.ttest_ind(
            drop_nans(x[0]),
            drop_nans(x[1]),
            alternative = x[2],
            equal_var = False
        ),
    ).rename({"column_0": "T-test", "column_1": "P-value"})

    return df.hstack(ttest, in_place = True)

def _add_fdr(df: pl.DataFrame) -> pl.DataFrame:
    df = df.sort(by = ["Protein ID", "P-value"], descending=[False, False])

    adj_pval_df = \
    df.group_by(by = "Protein ID", maintain_order = True).agg(
        pl.map_groups(
        exprs = "P-value",
        function = lambda x: sp.stats.false_discovery_control(x[0])
        ).alias("Adj. P-value")
    ).explode("Adj. P-value").select(pl.col("Adj. P-value"))

    return df.hstack(adj_pval_df, in_place = True)

def _add_ratio(df: pl.DataFrame, ctrl: list[str], test: list[str]) -> pl.DataFrame:
    log2fc_df = \
    df.select(
        pl.concat_list(pl.col(ctrl)),
        pl.concat_list(pl.col(test)),
    ).map_rows(lambda x:
        (
            np.nanstd(x[1]) / np.nanmean(x[1]),
            np.nanmean(x[1]) / np.nanmean(x[0])
        )
    ).rename({"column_0": "CV", "column_1": "FC"})

    return df.hstack(log2fc_df, in_place = True)

def _normalize_ratios(lip_df: pl.DataFrame, trp_df: pl.DataFrame) -> pl.DataFrame:
    norm = trp_df.select(pl.col("Protein ID", "P-value", "FC"))

    norm = norm.select(
        pl.col("Protein ID"),
        pl.when(
            (-pl.col("P-value").log10() > 2) & (pl.col("FC").log(base = 2).abs() > 1)
        )
        .then(
            pl.col("FC")
        )
        .otherwise(pl.lit(0.0))
        .alias("Normalization Factor")
    )

    lip_df = lip_df.join(norm, on = "Protein ID")

    lip_df = \
    lip_df.with_columns(
        pl.when(pl.col("Normalization Factor") > 0)
        .then(pl.col("FC") / pl.col("Normalization Factor"))
        .otherwise(pl.col("FC"))
        .alias("Normalized FC")
    )

    return lip_df


def _add_start_end_aa(df: pl.DataFrame) -> pl.DataFrame:
    df = \
    df.with_columns(
    [
        pl.col("Peptide Sequence").str.slice(0,1).alias("Start AA"),
        pl.col("Peptide Sequence").str.slice(-1).alias("End AA")
    ]
    )

    return df

def _add_half_trpytic(df: pl.DataFrame) -> pl.DataFrame:
    df = \
    df.with_columns(
    pl.when(
        (pl.col("Prev AA").is_in(["K", "R", "-"]) & ~pl.col("End AA").is_in(["K", "R"]) & ~pl.col("Next AA").is_in(["-"]))
    )
    .then(pl.lit("C_SEMI"))
    .when(
        ((~pl.col("Prev AA").is_in(["K", "R"]) & pl.col("End AA").is_in(["K", "R"]))) & ~pl.col("Prev AA").is_in(["-"]) |
        ((~pl.col("Prev AA").is_in(["K", "R"]) & pl.col("Next AA").is_in(["-"])))
    )
    .then(pl.lit("N_SEMI"))
    .otherwise(pl.lit("FULL_TRP"))
    .alias("Cleavage Type")
    ).with_columns(
        pl.when(pl.col("Cleavage Type") == "FULL_TRP")
        .then(pl.lit(False))
        .otherwise(pl.lit(True))
        .alias("Half Tryptic")
    )

    return df

def _add_cut_sites(df: pl.DataFrame) -> pl.DataFrame:
    df = \
    df.with_columns(
        pl.when(
            pl.col("Cleavage Type") == "C_SEMI"
        )
        .then(pl.format("{}{}", pl.col("Next AA"), pl.col("End") + 1))
        .when(
            pl.col("Cleavage Type") == "N_SEMI"
        )
        .then(pl.format("{}{}", pl.col("Start AA"), pl.col("Start")))
        .when(
            pl.col("Cleavage Type") == "FULL_TRP"
        )
        .then(
            pl.format("{}{}-{}{}", pl.col("Prev AA"), pl.col("Start") - 1, pl.col("End AA"), pl.col("End"))
        )
        .alias("Cut Site")
    ).with_columns(
    pl.format("{}_{}", pl.col("Protein ID"), pl.col("Cut Site"))
    .alias("Cut Site ID")
    )

    return df

def _log2(df: pl.DataFrame, col: str) -> pl.DataFrame:

    df = \
    df.with_columns(
        pl.col(col).log(base = 2)
        .alias(f"Log2 {col}")
    )
    return df

def _log10(df: pl.DataFrame, col: str) -> pl.DataFrame:

    df = \
    df.with_columns(
        -pl.col(col).log10()
        .alias(f"-Log10 {col}")
    )
    return df