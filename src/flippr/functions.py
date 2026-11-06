from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import polars as pl
import polars.selectors as cs
import scipy as sp


type RcParams = Mapping[str, Any]


def _zero_count_col(condition: str) -> str:
    return f"{condition} ZC"


def _intensity_list_col(condition: str) -> str:
    return f"{condition} Intensity"


def _mean_col(condition: str) -> str:
    return f"{condition} Mean"


def _std_col(condition: str) -> str:
    return f"{condition} Std"


def _ttest_ind_from_stats(stats: dict[str, Any]) -> list[float]:
    result = sp.stats.ttest_ind_from_stats(**stats, equal_var=False)
    return [float(result.statistic), float(result.pvalue)]


def _false_discovery_control(p_values: Sequence[float]) -> list[float]:
    return sp.stats.false_discovery_control(p_values).tolist()


def _add_missing_intensity_counts(
    df: pl.DataFrame,
    ctrl_name: str,
    test_name: str,
    ctrl_ints: list[str],
    test_ints: list[str],
) -> pl.DataFrame:
    return df.with_columns(
        pl.concat_list(ctrl_ints).list.count_matches(0).alias(_zero_count_col(ctrl_name)),
        pl.concat_list(test_ints).list.count_matches(0).alias(_zero_count_col(test_name)),
    )


def _filter_quantifiable_rows(
    df: pl.DataFrame,
    ctrl_name: str,
    test_name: str,
    ctrl_n_rep: int,
    test_n_rep: int,
    max_missing: int,
) -> pl.DataFrame:
    ctrl_zc = pl.col(_zero_count_col(ctrl_name))
    test_zc = pl.col(_zero_count_col(test_name))

    return df.filter(
        (ctrl_zc.le(max_missing) & test_zc.eq(0))
        | (ctrl_zc.eq(0) & test_zc.le(max_missing))
        | (ctrl_zc.eq(ctrl_n_rep) & test_zc.eq(0))
        | (ctrl_zc.eq(0) & test_zc.eq(test_n_rep))
    )


def _mark_remaining_zeroes_missing(
    df: pl.DataFrame,
    ctrl_ints: list[str],
    test_ints: list[str],
) -> pl.DataFrame:
    return df.with_columns(cs.by_name(ctrl_ints + test_ints).replace(0.0, None))


def _filter_intensities_by_missingness(
    df: pl.DataFrame,
    ctrl_name: str,
    test_name: str,
    ctrl_ints: list[str],
    test_ints: list[str],
    ctrl_n_rep: int,
    test_n_rep: int,
    rcParams: RcParams,
    **_: Any,
) -> pl.DataFrame:
    max_missing = int(rcParams.get("ion.missing_intensity_thresh", 1))

    df = _add_missing_intensity_counts(df, ctrl_name, test_name, ctrl_ints, test_ints)
    df = _filter_quantifiable_rows(df, ctrl_name, test_name, ctrl_n_rep, test_n_rep, max_missing)
    return _mark_remaining_zeroes_missing(df, ctrl_ints, test_ints)


def _add_ttest_alternative(
    df: pl.DataFrame,
    ctrl_name: str,
    test_name: str,
    ctrl_n_rep: int,
    test_n_rep: int,
    **_: Any,
) -> pl.DataFrame:
    ctrl_zc = pl.col(_zero_count_col(ctrl_name))
    test_zc = pl.col(_zero_count_col(test_name))

    return df.with_columns(
        pl.when(ctrl_zc.eq(ctrl_n_rep) & test_zc.eq(0))
        .then(pl.lit("less"))
        .when(ctrl_zc.eq(0) & test_zc.eq(test_n_rep))
        .then(pl.lit("greater"))
        .otherwise(pl.lit("two-sided"))
        .alias("Alternative Hypothesis")
    )


def _make_imputation_draws(
    height: int,
    ctrl_n_rep: int,
    test_n_rep: int,
    rcParams: RcParams,
) -> pl.DataFrame:
    loc = rcParams.get("ion.aon_impute_loc", 1e4)
    scale = rcParams.get("ion.aon_impute_scale", 1e3)
    seed = rcParams.get("ion.aon_impute_seed", 0)
    rng = np.random.default_rng(seed)

    return pl.DataFrame(
        {
            "CTRL_IMP": rng.normal(loc=loc, scale=scale, size=(height, ctrl_n_rep)).tolist(),
            "TEST_IMP": rng.normal(loc=loc, scale=scale, size=(height, test_n_rep)).tolist(),
        }
    )


def _add_condition_intensity_lists(
    df: pl.DataFrame,
    ctrl_name: str,
    test_name: str,
    ctrl_ints: list[str],
    test_ints: list[str],
    ctrl_n_rep: int,
    test_n_rep: int,
) -> pl.DataFrame:
    ctrl_zc = pl.col(_zero_count_col(ctrl_name))
    test_zc = pl.col(_zero_count_col(test_name))

    return df.with_columns(
        pl.when(ctrl_zc.eq(ctrl_n_rep) & test_zc.eq(0))
        .then(pl.col("CTRL_IMP"))
        .otherwise(pl.concat_list(ctrl_ints))
        .alias(_intensity_list_col(ctrl_name)),
        pl.when(ctrl_zc.eq(0) & test_zc.eq(test_n_rep))
        .then(pl.col("TEST_IMP"))
        .otherwise(pl.concat_list(test_ints))
        .alias(_intensity_list_col(test_name)),
    )


def _write_intensity_lists_to_replicates(
    df: pl.DataFrame,
    ctrl_name: str,
    test_name: str,
    ctrl_ints: list[str],
    test_ints: list[str],
) -> pl.DataFrame:
    df = df.with_columns(
        pl.col(_intensity_list_col(ctrl_name)).list.get(i).alias(col)
        for i, col in enumerate(ctrl_ints)
    )
    return df.with_columns(
        pl.col(_intensity_list_col(test_name)).list.get(i).alias(col)
        for i, col in enumerate(test_ints)
    )


def _drop_nulls_from_condition_lists(
    df: pl.DataFrame,
    ctrl_name: str,
    test_name: str,
) -> pl.DataFrame:
    ctrl_intensity = _intensity_list_col(ctrl_name)
    test_intensity = _intensity_list_col(test_name)

    return df.with_columns(
        pl.col(ctrl_intensity).list.drop_nulls().alias(ctrl_intensity),
        pl.col(test_intensity).list.drop_nulls().alias(test_intensity),
    )


def _drop_imputation_columns(
    df: pl.DataFrame,
    ctrl_name: str,
    test_name: str,
) -> pl.DataFrame:
    return df.drop("CTRL_IMP", "TEST_IMP", _zero_count_col(ctrl_name), _zero_count_col(test_name))


def _impute_all_or_none_intensities(
    df: pl.DataFrame,
    ctrl_name: str,
    test_name: str,
    ctrl_ints: list[str],
    test_ints: list[str],
    ctrl_n_rep: int,
    test_n_rep: int,
    rcParams: RcParams,
    **_: Any,
) -> pl.DataFrame:
    imputed = _make_imputation_draws(df.height, ctrl_n_rep, test_n_rep, rcParams)

    df = df.hstack(imputed)
    df = _add_condition_intensity_lists(
        df,
        ctrl_name,
        test_name,
        ctrl_ints,
        test_ints,
        ctrl_n_rep,
        test_n_rep,
    )
    df = _write_intensity_lists_to_replicates(df, ctrl_name, test_name, ctrl_ints, test_ints)
    df = _drop_nulls_from_condition_lists(df, ctrl_name, test_name)
    return _drop_imputation_columns(df, ctrl_name, test_name)


def _add_intensity_summary_stats(
    df: pl.DataFrame,
    ctrl_name: str,
    test_name: str,
) -> pl.DataFrame:
    return df.with_columns(
        pl.col(_intensity_list_col(ctrl_name)).list.mean().alias(_mean_col(ctrl_name)),
        pl.col(_intensity_list_col(ctrl_name)).list.std().alias(_std_col(ctrl_name)),
        pl.col(_intensity_list_col(test_name)).list.mean().alias(_mean_col(test_name)),
        pl.col(_intensity_list_col(test_name)).list.std().alias(_std_col(test_name)),
    )


def _add_welch_ttest_stats(
    df: pl.DataFrame,
    ctrl_name: str,
    test_name: str,
) -> pl.DataFrame:
    return df.with_columns(
        pl.struct(
            mean1=pl.col(_mean_col(ctrl_name)),
            std1=pl.col(_std_col(ctrl_name)),
            nobs1=pl.col(_intensity_list_col(ctrl_name)).list.len(),
            mean2=pl.col(_mean_col(test_name)),
            std2=pl.col(_std_col(test_name)),
            nobs2=pl.col(_intensity_list_col(test_name)).list.len(),
            alternative=pl.col("Alternative Hypothesis"),
        )
        .map_elements(_ttest_ind_from_stats, return_dtype=pl.List(pl.Float64))
        .alias("Stats")
    )


def _unpack_ttest_stats(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns(
        pl.col("Stats").list.first().alias("T-test"),
        pl.col("Stats").list.last().alias("P-value"),
    )


def _add_welch_ttest(
    df: pl.DataFrame,
    ctrl_name: str,
    test_name: str,
    **_: Any,
) -> pl.DataFrame:
    df = _add_intensity_summary_stats(df, ctrl_name, test_name)
    df = _add_welch_ttest_stats(df, ctrl_name, test_name)
    return _unpack_ttest_stats(df)


def _adjust_p_values_by_protein(df: pl.DataFrame, **_: Any) -> pl.DataFrame:
    df = df.sort(by=["Protein ID", "P-value"], descending=[False, False])

    adjusted_p_values = (
        df.group_by(by="Protein ID", maintain_order=True)
        .agg(pl.col("P-value"))
        .with_columns(
            pl.col("P-value")
            .map_elements(_false_discovery_control, return_dtype=pl.List(pl.Float64))
            .alias("Adj. P-value")
        )
        .explode("Adj. P-value")
        .select("Adj. P-value")
    )

    return df.hstack(adjusted_p_values)


def _drop_statistics_intermediates(
    df: pl.DataFrame,
    ctrl_name: str,
    test_name: str,
) -> pl.DataFrame:
    return df.drop(
        _mean_col(ctrl_name),
        _std_col(ctrl_name),
        _mean_col(test_name),
        _std_col(test_name),
        "Stats",
        _intensity_list_col(ctrl_name),
        _intensity_list_col(test_name),
    )


def _add_fold_change_and_cv(
    df: pl.DataFrame,
    ctrl_name: str,
    test_name: str,
    **_: Any,
) -> pl.DataFrame:
    df = df.with_columns(
        (pl.col(_mean_col(test_name)) / pl.col(_mean_col(ctrl_name))).alias("FC"),
        (pl.col(_std_col(test_name)) / pl.col(_mean_col(test_name))).alias("CV"),
    )
    return _drop_statistics_intermediates(df, ctrl_name, test_name)


def _select_trp_normalization_factors(
    trp_norm: pl.DataFrame,
    rcParams: RcParams,
) -> pl.DataFrame:
    trp_prot_fc = rcParams.get("trp_protein.fc_sig_thresh", 1.0)
    trp_prot_pval = rcParams.get("trp_protein.pval_sig_thresh", 0.01)

    return trp_norm.select(cs.by_name("Protein ID", "P-value", "FC")).select(
        pl.col("Protein ID"),
        pl.when(
            pl.col("P-value").le(trp_prot_pval)
            & pl.col("FC").log(base=2).abs().ge(trp_prot_fc)
        )
        .then(pl.col("FC"))
        .otherwise(pl.lit(0.0))
        .alias("Normalization Factor")
    )


def _apply_trp_normalization(
    df: pl.DataFrame,
    normalization_factors: pl.DataFrame,
) -> pl.DataFrame:
    df = df.join(normalization_factors, on="Protein ID", how="left")

    return df.with_columns(
        pl.col("Normalization Factor").fill_null(0.0)
    ).with_columns(
        pl.when(pl.col("Normalization Factor").gt(0))
        .then(pl.col("FC") / pl.col("Normalization Factor"))
        .otherwise(pl.col("FC"))
        .alias("Normalized FC")
    )


def _normalize_fold_change_by_trp(
    df: pl.DataFrame,
    trp_norm: pl.DataFrame,
    rcParams: RcParams,
    **_: Any,
) -> pl.DataFrame:
    normalization_factors = _select_trp_normalization_factors(trp_norm, rcParams)
    return _apply_trp_normalization(df, normalization_factors)


def _add_peptide_terminal_residues(df: pl.DataFrame, **_: Any) -> pl.DataFrame:
    return df.with_columns(
        pl.col("Peptide Sequence").str.head(1).alias("Start AA"),
        pl.col("Peptide Sequence").str.tail(1).alias("End AA"),
    )


def _cleavage_type_expr() -> pl.Expr:
    return (
        pl.when(
            pl.col("Prev AA").is_in(["K", "R", "-"])
            & pl.col("End AA").is_in(["K", "R"])
            & ~pl.col("Next AA").is_in(["-"])
        )
        .then(pl.lit("FULL_TRP"))
        .when(
            pl.col("Prev AA").is_in(["K", "R"])
            & pl.col("Next AA").is_in(["-"])
        )
        .then(pl.lit("FULL_TRP"))
        .when(
            pl.col("Prev AA").is_in(["M"])
            & pl.col("Start").eq(2)
            & pl.col("End AA").is_in(["K", "R"])
        )
        .then(pl.lit("FULL_TRP"))
        .when(
            pl.col("Prev AA").is_in(["K", "R", "-"])
            & ~pl.col("End AA").is_in(["K", "R"])
            & ~pl.col("Next AA").is_in(["-"])
        )
        .then(pl.lit("C_SEMI"))
        .when(
            pl.col("Prev AA").is_in(["M"])
            & pl.col("Start").eq(2)
            & ~pl.col("End AA").is_in(["K", "R"])
        )
        .then(pl.lit("C_SEMI"))
        .when(
            ~pl.col("Prev AA").is_in(["K", "R", "-"])
            & pl.col("End AA").is_in(["K", "R"])
            & ~pl.col("Next AA").is_in(["-"])
        )
        .then(pl.lit("N_SEMI"))
        .when(
            ~pl.col("Prev AA").is_in(["K", "R", "-"])
            & pl.col("Next AA").is_in(["-"])
        )
        .then(pl.lit("N_SEMI"))
        .otherwise(pl.lit(None))
        .alias("Cleavage Type")
    )


def _add_tryptic_cleavage_annotations(df: pl.DataFrame, **_: Any) -> pl.DataFrame:
    return (
        df.with_columns(_cleavage_type_expr())
        .filter(~pl.col("Cleavage Type").is_null())
        .with_columns(pl.col("Cleavage Type").ne("FULL_TRP").alias("Half Tryptic"))
    )


def _cut_site_expr() -> pl.Expr:
    return (
        pl.when(pl.col("Cleavage Type").eq("C_SEMI"))
        .then(pl.format("{}{}", pl.col("Next AA"), pl.col("End") + 1))
        .when(pl.col("Cleavage Type").eq("N_SEMI"))
        .then(pl.format("{}{}", pl.col("Start AA"), pl.col("Start")))
        .when((pl.col("Cleavage Type").eq("FULL_TRP")) & (pl.col("Prev AA").ne("-")))
        .then(
            pl.format(
                "{}{}-{}{}",
                pl.col("Prev AA"),
                pl.col("Start") - 1,
                pl.col("End AA"),
                pl.col("End"),
            )
        )
        .when((pl.col("Cleavage Type").eq("FULL_TRP")) & (pl.col("Prev AA").eq("-")))
        .then(
            pl.format(
                "{}{}-{}{}",
                pl.col("Start AA"),
                pl.col("Start"),
                pl.col("End AA"),
                pl.col("End"),
            )
        )
        .alias("Cut Site")
    )


def _add_cut_site_ids(df: pl.DataFrame, **_: Any) -> pl.DataFrame:
    return (
        df.with_columns(_cut_site_expr())
        .with_columns(
            pl.format("{}_{}", pl.col("Protein ID"), pl.col("Cut Site")).alias("Cut Site ID")
        )
    )


def _add_log2_column(df: pl.DataFrame, col: str) -> pl.DataFrame:
    return df.with_columns(
        pl.when(pl.col(col).eq(0))
        .then(0.0)
        .otherwise(pl.col(col).log(base=2))
        .alias(f"Log2 {col}")
    )


def _add_neg_log10_column(df: pl.DataFrame, col: str) -> pl.DataFrame:
    return df.with_columns((-pl.col(col).log10()).alias(f"-Log10 {col}"))
