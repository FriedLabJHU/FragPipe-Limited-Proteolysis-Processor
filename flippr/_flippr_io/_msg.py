# ---------------------------------------------------------------------
# FragPipe Limited-Proteolysis Processor by Edgar Manriquez-Sandoval
# Stephen D. Fried Lab @ The Johns Hopkins University
# ---------------------------------------------------------------------

import pandas as pd
from pathlib import Path


def error(s, *args):
    print("ERROR:")
    print(f"\t{s}")
    if args:
        for sarg in args:
            print(f"\t{sarg}")
    print("Exiting program")
    exit(1)


def notice(s, *args):
    print("NOTICE:")
    print(f"\t{s}")
    if args:
        for sarg in args:
            print(f"\t{sarg}")


def msg(s, *args):
    print(s)
    if args:
        for sarg in args:
            print(sarg)


def print_to_excel(data, headers, write_path):
    df = pd.json_normalize(data)
    df = df[headers]
    df.to_excel(write_path, index=False)


def output_ion_data(ions: list, cond: str, params: object):
    cond_reps = params.condition_replicates[cond]
    out_path = params.cond_out_paths[cond]

    header_keys = [
        "Peptide Sequence",
        "Modified Sequence",
        "Prev AA",
        "Next AA",
        "Start",
        "End",
        "Peptide Length",
        "M/Z",
        "Charge",
        "Compensation Voltage",
        "Assigned Modifications",
        "Protein",
        "Protein ID",
        "Entry Name",
        "Gene",
        "Protein Description",
        "Mapped Genes",
        "Mapped Proteins",
    ]

    rep_keys = [rep for rep in params.control_replicates]
    rep_keys += [rep for rep in cond_reps]

    rep_keys = [
        f"{rep} {hk}"
        for hk in ["Spectral Count", "Apex Retention Time", "Intensity", "Match Type"]
        for rep in rep_keys
    ]

    header_keys += rep_keys
    header_keys += [
        "Log2 FC",
        "Log10 P-Value",
        "Log10 Adj. P-Value",
        "P-Value",
        "Adj. P-Value",
        "CV",
    ]

    write_path = Path(
        out_path, f"{params.prefix}{cond}_v_{params.control_annotation}_ions.xlsx"
    )
    print_to_excel(ions, header_keys, write_path)


def output_norm_protein_data(proteins: list, cond: str, params: object):
    cont_reps = params.norm_control_replicates
    cond_reps = params.norm_condition_replicates[cond]
    out_path = params.cond_out_paths[cond]

    header_keys = [
        "Protein",
        "Protein ID",
        "Entry Name",
        "Gene",
        "Protein Length",
        "Organism",
        "Protein Existence",
        "Description",
        "Protein Probability",
        "Top Peptide Probability",
        "Combined Total Peptides",
        "Indistinguishable Proteins",
    ]

    rep_keys = [rep for rep in cont_reps]
    rep_keys += [rep for rep in cond_reps]

    rep_keys = [
        f"{rep} {hk}"
        for hk in ["Spectral Count", "MaxLFQ Intensity"]
        for rep in rep_keys
    ]

    header_keys += rep_keys
    header_keys += ["Log2 FC", "Log10 P-Value", "CV"]

    write_path = Path(
        out_path,
        f"{params.prefix}{cond}_v_{params.control_annotation}_normalization_factors.xlsx",
    )
    print_to_excel(proteins.values(), header_keys, write_path)


def output_peptide_data(
    peptides: list, fixed_headers: list, filename: str, cond: str, params: object
):
    out_path = params.cond_out_paths[cond]

    write_path = Path(
        out_path, f"{params.prefix}{cond}_v_{params.control_annotation}_{filename}.xlsx"
    )
    header_keys = (
        [
            "Item ID",
        ]
        + fixed_headers
        + [
            "No. of Ions",
            "Log2 FC",
            "Normalized Log2 FC",
            "Log10 P-Value",
            "Log10 Adj. P-Value",
            "P-Value",
            "Adj. P-Value",
            "CV",
        ]
    )

    print_to_excel(peptides, header_keys, write_path)


def output_protein_data(proteins: list, cond: str, params: object):
    out_path = params.cond_out_paths[cond]

    write_path = Path(
        out_path,
        f"{params.prefix}{cond}_v_{params.control_annotation}_protein_summary.xlsx",
    )

    header_keys = [
        "Protein ID",
        "Gene",
        "No. of Ions",
        "No. of Modified Peptides",
        "No. of Modified Peptides (Valid P-Value)",
        "No. of Modified Peptides (Significant P-Value)",
        "No. of Modified Peptides (Valid Adj. P-Value)",
        "No. of Modified Peptides (Significant Adj. P-Value)",
        "No. of Peptides",
        "No. of Peptides (Valid P-Value)",
        "No. of Peptides (Significant P-Value)",
        "No. of Peptides (Valid Adj. P-Value)",
        "No. of Peptides (Significant Adj. P-Value)",
        "No. of Cut Sites",
        "No. of Cut Sites (Valid P-Value)",
        "No. of Cut Sites (Significant P-Value)",
        "No. of Cut Sites (Valid Adj. P-Value)",
        "No. of Cut Sites (Significant Adj. P-Value)",
    ]

    print_to_excel(proteins, header_keys, write_path)
