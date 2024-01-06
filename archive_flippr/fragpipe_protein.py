# ---------------------------------------------------------------------
# FragPipe Limited-Proteolysis Processor by Edgar Manriquez-Sandoval
# Stephen D. Fried Lab @ The Johns Hopkins University
# ---------------------------------------------------------------------

import numpy as np
from scipy.stats import ttest_ind, variation, norm

from _flippr_io._params import Params
from _flippr_io._fragpipe_parser import FragPipeObject


def validate_max_lfq_abundances(protein: FragPipeObject, cond: str, params: Params):
    cont_reps = params.norm_control_replicates
    cond_reps = params.norm_condition_replicates[cond]

    control_abundances = np.array(
        [
            protein[f"{r} MaxLFQ Intensity"]
            for r in cont_reps
            if protein[f"{r} MaxLFQ Intensity"]
        ],
        dtype=np.float32,
    )
    condition_abundances = np.array(
        [
            protein[f"{r} MaxLFQ Intensity"]
            for r in cond_reps
            if protein[f"{r} MaxLFQ Intensity"]
        ],
        dtype=np.float32,
    )

    control_cnt = np.count_nonzero(control_abundances)
    condition_cnt = np.count_nonzero(condition_abundances)

    replicates = params.replicates

    # set t-test H1
    protein.ttest_alt = "two-sided"

    if (control_cnt + condition_cnt) == 0:
        return

    # If all values are non-zero, do nothing
    if (control_cnt + condition_cnt) == 2 * replicates:
        return protein

    # If the control has all zero values, return control as 1000's
    if control_cnt == 0:
        protein.ttest_alt = "greater"
        protein.update({f"{r} MaxLFQ Intensity": norm.rvs(1e4, 1e3) for r in cont_reps})

    # If the experimental condition has all zero values, return the experimental condition as 1000's
    if condition_cnt == 0:
        protein.ttest_alt = "less"
        protein.update({f"{r} MaxLFQ Intensity": norm.rvs(1e4, 1e3) for r in cond_reps})

    if 1 <= control_cnt < replicates:
        protein.update(
            {
                f"{r} MaxLFQ Intensity": None
                for r in cont_reps
                if not protein[f"{r} MaxLFQ Intensity"]
            }
        )

    if 1 <= condition_cnt < replicates:
        protein.update(
            {
                f"{r} MaxLFQ Intensity": None
                for r in cond_reps
                if not protein[f"{r} MaxLFQ Intensity"]
            }
        )

    return protein


def calculate_max_lfq_ratio(protein: FragPipeObject, cond: str, params: Params):
    cont_reps = params.norm_control_replicates
    cond_reps = params.norm_condition_replicates[cond]

    control_abundances = np.array(
        [
            protein[f"{r} MaxLFQ Intensity"]
            for r in cont_reps
            if protein[f"{r} MaxLFQ Intensity"]
        ],
        dtype=np.float32,
    )
    condition_abundances = np.array(
        [
            protein[f"{r} MaxLFQ Intensity"]
            for r in cond_reps
            if protein[f"{r} MaxLFQ Intensity"]
        ],
        dtype=np.float32,
    )

    ttest, pval = ttest_ind(
        control_abundances,
        condition_abundances,
        alternative=protein.ttest_alt,
        equal_var=False,
    )

    ratio = np.log2(np.mean(condition_abundances) / np.mean(control_abundances))

    var = variation(condition_abundances)

    log10pval = np.abs(np.log10(pval)) if not np.isnan(pval) else 0

    if np.abs(ratio) >= 2 and log10pval >= 1:
        norm_ratio = ratio
    else:
        # this the normalization is not significant, then the normalization factor should be zero
        norm_ratio = 0

    protein.update(
        {
            "Log2 FC": norm_ratio,
            "P-Value": pval,
            "Log10 P-Value": log10pval,
            "T-test": ttest,
            "CV": var,
        }
    )

    return protein


def generate_normalization_ratios(proteins):
    norm_proteins = {}
    for i, protein in enumerate(proteins):
        norm_proteins.setdefault(protein["Protein ID"], protein)

    return norm_proteins
