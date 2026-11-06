rcParams and configuration
==========================

FLiPPR exposes a small configuration dictionary, `rcParams`, defined in `src/flippr/parameters.py`.

Common configuration keys include:

- `ion.missing_intensity_thresh` : threshold for missing ion intensities
- `ion.aon_impute_loc`, `ion.aon_impute_scale` : parameters for AON imputation
- `ion.aon_impute_seed` : RNG seed for AON imputation; an integer (default `0`) makes
  imputation reproducible across runs, or set to `None` for non-deterministic draws
- `trp_protein.intensity_value` : which TrP protein intensity column to use (e.g. "MaxLFQ Intensity")
- significance thresholds for proteins and TrP-derived normalization:
  - `trp_protein.fc_sig_thresh`, `trp_protein.pval_sig_thresh`
  - `protein.fc_sig_thresh`, `protein.pval_sig_thresh`, `protein.adj_pval_sig_thresh`

Modify `rcParams` before running a study, for example:

.. code-block:: python

    from flippr.parameters import rcParams
    rcParams["ion.missing_intensity_thresh"] = 2

Refer to `src/flippr/parameters.py` for the full list of keys and defaults.
