API Reference
=============

This page summarizes the public API used by typical FLiPPR workflows.

Study
-----

Class: `Study`

- Constructor: `Study(lip: str | Path, trp: Optional[str | Path] = None, method: str = "dda")`
- Properties: `samples` (dict)
- Methods: `add_process(pid, lip_ctrl, lip_test, n_rep, trp_ctrl=None, trp_test=None, trp_n_rep=None)`, `run()`

Process & Result
------------------

`Process` and `Result` are types used internally and returned by `Study.run()`.

- `Result` exposes dataframes and summary tables:
  - `ion` : ion-level `polars.DataFrame`
  - `peptide` : peptide-level `polars.DataFrame`
  - `modified_peptide` : modified peptide-level `polars.DataFrame`
  - `cut_site` : cut-site-level `polars.DataFrame`
  - `protein_summary` : protein summary `polars.DataFrame`
  - `name` : human-readable process name

Combine helpers
---------------

The module-level helpers in `combine.py` are convenient for downstream aggregation:

- `combine_by(df: polars.DataFrame, by: str, fc: str) -> polars.DataFrame`
- `summary_by(df: polars.DataFrame, by: str, fc: str, rcParams: dict) -> polars.DataFrame`

Notes
-----

Many functions under `src/flippr` use a leading underscore to indicate internal utilities. The public surface is primarily the `Study` class and the `combine` helpers.
