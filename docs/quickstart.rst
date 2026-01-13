Quickstart
==========

This quickstart shows the minimal steps to run a FLiPPR analysis programmatically.

Installation
------------

Install from source in a modern Python (see `pyproject.toml`):

.. code-block:: bash

    python -m pip install .

Basic example
-------------

The typical workflow uses `Study` to register FragPipe output directories, add one or more processes, and run them.

.. code-block:: python

    from flippr import Study

    # Create a Study for a LiP dataset (optional TrP normalization can be provided)
    study = Study(lip="/path/to/LiP_LFQ", trp=None, method="dda")

    # Add a process: pid, control name, test name, replicates
    study.add_process(pid="exp1", lip_ctrl="CTRL", lip_test="TREAT", n_rep=3)

    # Run all processes
    results = study.run()

    # Access a result and its tables
    res = results["exp1"]
    ion_df = res.ion  # FLiPPR ion-level dataframe (polars.DataFrame)
    peptides = res.peptide
    proteins = res.protein_summary

Replicate argument
------------------

The `n_rep` (and `trp_n_rep`) argument supports several shapes:

- An integer for equal numbers of control/test replicates: `3`
- A 2-tuple for different counts: `(2, 3)`
- A pair of index tuples to reference specific replicate numbers: `((1,2), (1,2,3))`
