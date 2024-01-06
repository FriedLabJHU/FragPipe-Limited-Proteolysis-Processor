# ---------------------------------------------------------------------
# FragPipe Limited-Proteolysis Processor by Edgar Manriquez-Sandoval
# Stephen D. Fried Lab @ The Johns Hopkins University
# ---------------------------------------------------------------------

import numpy as np
from scipy.stats import ttest_ind, variation, norm

from _flippr_io._params import Params
from _flippr_io._fragpipe_parser import FragPipeObject


def validate_ion_abundances(ion: FragPipeObject, cond: str, params: Params):
    cont_reps = params.control_replicates
    cond_reps = params.condition_replicates[cond]

    # Separating out the ion abundances
    control_abundances = np.array(
        [ion[f"{r} Intensity"] for r in cont_reps], dtype=np.float32
    )
    condition_abundances = np.array(
        [ion[f"{r} Intensity"] for r in cond_reps], dtype=np.float32
    )

    control_cnt = np.count_nonzero(control_abundances)
    condition_cnt = np.count_nonzero(condition_abundances)

    replicates = params.replicates

    # set t-test H1
    ion.ttest_alt = "two-sided"

    # If all values are non-zero, do nothing
    if (control_cnt + condition_cnt) == 2 * replicates:
        return ion

    # If the control has all zero values and the experimental condition has all non-zero values, return control as 1000's
    if control_cnt == 0 and condition_cnt == replicates:
        ion.ttest_alt = "greater"
        ion.update({f"{r} Intensity": norm.rvs(1e4, 1e3) for r in cont_reps})
        return ion

    # If the control has all non-zero values and the experimental condition has all zero values, return the experimental condition as 1000's
    if control_cnt == replicates and condition_cnt == 0:
        ion.ttest_alt = "less"
        ion.update({f"{r} Intensity": norm.rvs(1e4, 1e3) for r in cond_reps})
        return ion

    # only apply the following imputation if there are more than 3 values
    if replicates > 2:
        # If the control has only one zero value and the the experimental condition has none, return control with only non-zero values
        if control_cnt == replicates - 1 and condition_cnt == replicates:
            ion.update(
                {f"{r} Intensity": None for r in cont_reps if not ion[f"{r} Intensity"]}
            )
            return ion

        # If the the experimental condition has only one zero value and the control has none, return the experimental condition with only non-zero values
        if control_cnt == replicates and condition_cnt == replicates - 1:
            ion.update(
                {f"{r} Intensity": None for r in cond_reps if not ion[f"{r} Intensity"]}
            )
            return ion

    # If none of the above conditions are met, then don't use this ion
    return None


def calculate_ion_ratio(ion: FragPipeObject, cond: str, params: Params):
    cont_reps = params.control_replicates
    cond_reps = params.condition_replicates[cond]

    control_abundances = np.array(
        [ion[f"{r} Intensity"] for r in cont_reps if ion[f"{r} Intensity"]],
        dtype=np.float32,
    )
    condition_abundances = np.array(
        [ion[f"{r} Intensity"] for r in cond_reps if ion[f"{r} Intensity"]],
        dtype=np.float32,
    )

    ratio = np.log2(np.mean(condition_abundances) / np.mean(control_abundances))

    var = variation(condition_abundances)

    ttest, pval = ttest_ind(
        condition_abundances,
        control_abundances,
        alternative=ion.ttest_alt,
        equal_var=False,
    )

    ion.update(
        {
            "Log2 FC": ratio,
            "P-Value": pval,
            "Log10 P-Value": np.abs(np.log10(pval)),
            "T-test": ttest,
            "CV": var,
        }
    )

    return ion
