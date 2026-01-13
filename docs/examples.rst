Examples & Use Cases
=====================

The examples below show common FLiPPR workflows found in public use-cases.

Simple fold-change (LiP only)
-----------------------------

.. code-block:: python

    from flippr import Study

    study = Study(lip="/data/PXD/LiP_LFQ")
    study.add_process("ctrl_v_treat", "CTRL", "TREAT", 3)
    results = study.run()
    res = results["ctrl_v_treat"]
    # Inspect top peptides
    print(res.peptide.head(10))

Fold-change with TrP normalization
---------------------------------

If you have matched TrP data (no broad protease), provide it to `Study` and pass TrP args to `add_process`:

.. code-block:: python

    study = Study(lip="/data/PXD/LiP_LFQ", trp="/data/PXD/TrP_LFQ")
    study.add_process("normed", "CTRL", "TREAT", 3, trp_ctrl="CTRL_TrP", trp_test="TREAT_TrP", trp_n_rep=3)
    results = study.run()

Multiple-dose experiment
------------------------

Add multiple processes to the same `Study` and run them together:

.. code-block:: python

    study.add_process("dose_lo", "CTRL", "DRUG_LO", 3, trp_ctrl="CTRL_TrP", trp_test="DRUG_TrP", trp_n_rep=3)
    study.add_process("dose_hi", "CTRL", "DRUG_HI", 3, trp_ctrl="CTRL_TrP", trp_test="DRUG_TrP", trp_n_rep=3)
    results = study.run()
