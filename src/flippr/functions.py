import numpy as np
import polars as pl
import scipy as sp


def _cull_intensities(df: pl.DataFrame,
                      ctrl_name: str,  test_name: str,
                      ctrl_ints: list[str],  test_ints: list[str],
                      n_rep: int,
                      rcParams: dict, **kwargs) -> pl.DataFrame:

    max_missing = rcParams.get("ion.missing_intensity_thresh", 1)

    df = \
    df.with_columns(# Count the number of zeros intensity replicates
        pl.concat_list(ctrl_ints).list.count_matches(0).alias(f"{ctrl_name} ZC"),
        pl.concat_list(test_ints).list.count_matches(0).alias(f"{test_name} ZC"),
    ).filter(# Cull based on zero count (ZC)
        (pl.col(f"{ctrl_name} ZC").le(max_missing) & pl.col(f"{test_name} ZC").eq(0))
        | (pl.col(f"{ctrl_name} ZC").eq(0)           & pl.col(f"{test_name} ZC").le(max_missing))
        | (pl.col(f"{ctrl_name} ZC").eq(n_rep)       & pl.col(f"{test_name} ZC").eq(0))
        | (pl.col(f"{ctrl_name} ZC").eq(0)           & pl.col(f"{test_name} ZC").eq(n_rep))
    ).with_columns(# Replace zero values remaining with NaN
        pl.col(ctrl_ints).replace(0.0, float("nan")),
        pl.col(test_ints).replace(0.0, float("nan"))
    )

    return df


def _add_alt_hypothesis(df: pl.DataFrame,
                        ctrl_name: str,  test_name: str,
                        n_rep: int,
                        rcParams: dict, **kwargs) -> pl.DataFrame:

    df = \
    df.with_columns(# Add alternative hypothesis
        pl.when(pl.col(f"{ctrl_name} ZC").eq(n_rep) & pl.col(f"{test_name} ZC").eq(0))
        .then(pl.lit("greater"))
        .when(pl.col(f"{ctrl_name} ZC").eq(0) & pl.col(f"{test_name} ZC").eq(n_rep))
        .then(pl.lit("less"))
        .otherwise(pl.lit("two-sided"))
        .alias("Alternative Hypothesis")
    )

    return df


def _impute_aon_intensities(df: pl.DataFrame,
                            ctrl_name: str,  test_name: str,
                            ctrl_ints: list[str],  test_ints: list[str],
                            n_rep: int,
                            rcParams: dict, **kwargs) -> pl.DataFrame:


    # Add imputation variables
    loc = rcParams.get("ion.aon_impute_loc", 1e4)
    scale = rcParams.get("ion.aon_impute_scale", 1e4)
    imp_df = pl.from_numpy(np.random.normal(loc=loc, scale=scale, size=(df.height, n_rep)), schema=["IMP"])
    df = df.hstack(imp_df)

    df = \
    df.with_columns(# Impute only on AON ions in control conditions
        pl.when(pl.col(f"{ctrl_name} ZC").eq(n_rep) & pl.col(f"{test_name} ZC").eq(0))
        .then(pl.concat_list("IMP"))
        .otherwise(pl.concat_list(ctrl_ints))
        .alias(f"{ctrl_name} Intensity"),
        # Impute only on AON ions in test conditions
        pl.when(pl.col(f"{ctrl_name} ZC").eq(0) & pl.col(f"{test_name} ZC").eq(n_rep))
        .then(pl.concat_list("IMP"))
        .otherwise(pl.concat_list(test_ints))
        .alias(f"{test_name} Intensity")
    ).with_columns(# Seperate out imputation (if perfromed) vars for control conditions
        pl.col(f"{ctrl_name} Intensity").list.get(i).alias(col) for i, col in enumerate(ctrl_ints)
    ).with_columns(# Seperate out imputation (if perfromed) vars for test conditions
        pl.col(f"{test_name} Intensity").list.get(i).alias(col) for i, col in enumerate(test_ints)
    ).with_columns(# Drop NaNs for further computations, this column will be droped down-stream
        pl.col(f"{ctrl_name} Intensity").list.eval(pl.element().drop_nans()).alias(f"{ctrl_name} Intensity"),
        pl.col(f"{test_name} Intensity").list.eval(pl.element().drop_nans()).alias(f"{test_name} Intensity")
    ).drop(["IMP", f"{ctrl_name} ZC", f"{test_name} ZC"])

    return df


def _add_ttest(df: pl.DataFrame,
               ctrl_name: str,  test_name: str,
               n_rep: int,
               rcParams: dict, **kwargs) -> pl.DataFrame:

    df = \
    df.with_columns(# Prepare the descriptive stats for `.ttest_ind_from_stats()`
        pl.col(f"{ctrl_name} Intensity").list.mean().alias(f"{ctrl_name} Mean"),
        pl.col(f"{ctrl_name} Intensity").list.std().alias(f"{ctrl_name} Std"),

        pl.col(f"{test_name} Intensity").list.mean().alias(f"{test_name} Mean"),
        pl.col(f"{test_name} Intensity").list.std().alias(f"{test_name} Std"),
    ).with_columns(
        pl.struct(# Create struct with `.ttest_ind_from_stats()` var names
            mean1=pl.col(f"{ctrl_name} Mean"),
            std1=pl.col(f"{ctrl_name} Std"),
            nobs1=n_rep,

            mean2=pl.col(f"{test_name} Mean"),
            std2=pl.col(f"{test_name} Std"),
            nobs2=n_rep,
            alternative=pl.col("Alternative Hypothesis")
        ).map_elements(# Perform T-test, slow step
            lambda x: sp.stats.ttest_ind_from_stats(**x, equal_var=False),
            return_dtype=pl.List(pl.Float64)
        ).alias("Stats")
    ).with_columns(# Seperate out T-test vars
        pl.col("Stats").list.first().alias("T-test"),
        pl.col("Stats").list.last().alias("P-value"),
    ).drop(["Stats", f"{ctrl_name} Intensity", f"{test_name} Intensity"])

    return df


def _add_fdr(df: pl.DataFrame, rcParams: dict, **kwargs) -> pl.DataFrame:

    # Sort on P-value
    df = df.sort(by=["Protein ID", "P-value"], descending=[False, False])

    df = \
    df.hstack(
        df.group_by(by="Protein ID", maintain_order=True)
        .agg(# Calculate Adj. P-value
            pl.map_groups(
                exprs="P-value",
                function=lambda x: sp.stats.false_discovery_control(x[0]),
                return_dtype=pl.List(pl.Float64)
            ).alias("Adj. P-value"),
        )# Add Adj. P-value to the original `DataFrame`
        .explode(pl.col("Adj. P-value"))
        .select(pl.col("Adj. P-value"))
    )

    return df


def _add_ratio(df: pl.DataFrame,
               ctrl_name: str,  test_name: str,
               rcParams: dict, **kwargs) -> pl.DataFrame:

    df = \
    df.with_columns(# Calculate FC and CV
        (pl.col(f"{test_name} Mean") / pl.col(f"{ctrl_name} Mean")).alias("FC"),
        (pl.col(f"{test_name} Std") / pl.col(f"{test_name} Mean")).alias("CV")
    ).drop([f"{ctrl_name} Mean", f"{ctrl_name} Std", f"{test_name} Mean", f"{test_name} Std"])

    return df


def _normalize_ratios(df: pl.DataFrame, trp_norm: pl.DataFrame, rcParams: dict, **kwargs) -> pl.DataFrame:

    trp_prot_fc = rcParams.get("trp_protein.fc_sig_tresh", 1.0)
    trp_prot_pval = rcParams.get("trp_protein.pval_sig_tresh", 0.01)

    trp_norm = trp_norm.select(pl.col("Protein ID", "P-value", "FC"))
    trp_norm = trp_norm.select(# Select only the significantly different proteins
        pl.col("Protein ID"),
        pl.when(pl.col("P-value").le(trp_prot_pval) & pl.col("FC").log(base=2).abs().ge(trp_prot_fc))
        .then(pl.col("FC"))
        .otherwise(pl.lit(0.0))
        .alias("Normalization Factor")
    )

    df = df.join(trp_norm, on="Protein ID", how="left")

    df = df.with_columns(# Apply normalization
        pl.col("Normalization Factor").fill_null(strategy="zero"),
        pl.when(pl.col("Normalization Factor").gt(0))
        .then(pl.col("FC") / pl.col("Normalization Factor"))
        .otherwise(pl.col("FC"))
        .alias("Normalized FC")
    )

    return df


def _add_start_end_aa(df: pl.DataFrame, rcParams: dict, **kwargs) -> pl.DataFrame:

    df = \
    df.with_columns(# Parse the Starting and Ending AA of the peptide
        pl.col("Peptide Sequence").str.head(1).alias("Start AA"),
        pl.col("Peptide Sequence").str.tail(1).alias("End AA"),
    )

    return df


def _add_half_trpytic(df: pl.DataFrame, rcParams: dict, **kwargs) -> pl.DataFrame:

    df = \
    df.with_columns(
        pl.when(# [K/R/-].~~~~~(K/R).[X]
            pl.col("Prev AA").is_in(["K", "R", "-"])
            & pl.col("End AA").is_in(["K", "R"])
            & ~pl.col("Next AA").is_in(["-"]))
        .then(pl.lit("FULL_TRP"))
        .when(# [K/R/-].~~~~~(X).[-]
            pl.col("Prev AA").is_in(["K", "R"])
            & pl.col("Next AA").is_in(["-"]))
        .then(pl.lit("FULL_TRP"))
        .when(# [M].2~~~~~(K/R).[X]
            pl.col("Prev AA").is_in(["M"])
            & pl.col("Start").eq(2)
            & pl.col("End AA").is_in(["K", "R"]))
        .then(pl.lit("FULL_TRP"))
        .when(# [K/R/-].~~~~~(X).[X]
            pl.col("Prev AA").is_in(["K", "R", "-"])
            & ~pl.col("End AA").is_in(["K", "R"])
            & ~pl.col("Next AA").is_in(["-"]))
        .then(pl.lit("C_SEMI"))
        .when(# [M].2~~~~~(X).[X]
            pl.col("Prev AA").is_in(["M"])
            & pl.col("Start").eq(2)
            & ~pl.col("End AA").is_in(["K", "R"]))
        .then(pl.lit("C_SEMI"))
        .when(# [X].~~~~~(K/R).[X]
            ~pl.col("Prev AA").is_in(["K", "R", "-"])
            & pl.col("End AA").is_in(["K", "R"])
            & ~pl.col("Next AA").is_in(["-"]))
        .then(pl.lit("N_SEMI"))
        .when(# [X].~~~~~(X).[-]
            ~pl.col("Prev AA").is_in(["K", "R", "-"])
            & pl.col("Next AA").is_in(["-"]))
        .then(pl.lit("N_SEMI"))
        .otherwise(pl.lit(None))
        .alias("Cleavage Type")
    ).filter(# Filter out over-digested peptides
        ~pl.col("Cleavage Type").is_null()
    ).with_columns(# Fill in `Half Tryptic`
        pl.col("Cleavage Type").ne("FULL_TRP").alias("Half Tryptic")
    )

    return df


def _add_cut_sites(df: pl.DataFrame, rcParams: dict, **kwargs) -> pl.DataFrame:

    df = \
    df.with_columns(# Adding unique identifier for all `Cleavage Type`
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
    ).with_columns(
        pl.format("{}_{}", pl.col("Protein ID"), pl.col("Cut Site")).alias("Cut Site ID")
    )

    return df


def _log2(df: pl.DataFrame, col: str) -> pl.DataFrame:

    df = df.with_columns(
        pl.when(pl.col(col).eq(0))
        .then(0.0)
        .otherwise(pl.col(col).log(base=2))
        .alias(f"Log2 {col}")
    )

    return df


def _log10(df: pl.DataFrame, col: str) -> pl.DataFrame:
    return df.with_columns(-pl.col(col).log10().alias(f"-Log10 {col}"))
