from typing import Any

rcParams: dict[str, Any] = {
    "ion.missing_intensity_thresh": 1,
    "ion.aon_impute_type": "gaussian", # unused for now
    "ion.aon_impute_loc": 1e4,
    "ion.aon_impute_scale": 1e3,
    "trp_protein.intensity_value": "MaxLFQ Intensity", # unused for dia methods
    "trp_protein.fc_sig_tresh": 1.0,
    "trp_protein.pval_sig_tresh": 0.01,
    "protein.fc_sig_thresh": 1.0,
    "protein.pval_sig_thresh": 0.01,
    "protein.adj_pval_sig_thresh": 0.05,
}

_DDA_FP_FILES: list[str] = [
    "combined_ion.tsv",
    "combined_protein.tsv",
    "experiment_annotation.tsv",
]

_DDA_FP_CONSTANT_ION_COLUMNS: list[str] = [
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

_DDA_FP_VARIABLE_ION_COLUMNS: list[str] = [
    "Spectral Count",
    "Apex Retention Time",
    "Intensity",
    "Match Type",
]

_DDA_FP_CONSTANT_PROTEIN_COLUMNS: list[str] = [
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
    "Combined Spectral Count",
    "Combined Unique Spectral Count",
    "Combined Total Spectral Count",
]

_DDA_FP_VARIABLE_PROTEIN_COLUMNS: list[str] = [
    "Spectral Count",
    "Intensity",
    "MaxLFQ Intensity",
]

_DIA_FP_FILES: list[str] = [
    "ion.tsv",
    "dia-quant-output/report.pr_matrix.tsv",
    "dia-quant-output/report.pg_matrix.tsv",
    "experiment_annotation.tsv",
]

_DIA_FP_CONSTANT_ION_COLUMNS: list[str] = [
    "Protein ID",
    "Peptide Sequence",
    "Prev AA",
    "Next AA",
    "Protein Start",
    "Protein End",
]

_DIA_DIANN_CONSTANT_ION_COLUMNS: list[str] = [
    "Protein.Group",
    "Stripped.Sequence",
    "Modified.Sequence",
    "Precursor.Charge",
    "Protein.Ids",
    "Protein.Names",
    "Genes",
    "First.Protein.Description",
]

_DIA_DIANN_CONSTANT_PROTEIN_COLUMNS: list[str] = [
    "Protein.Group",
    "Protein.Names",
    "Genes",
    "First.Protein.Description",
]

_DIA_RENAME_FP_ION: dict[str, str] = {
    "Protein Start": "Start",
    "Protein End": "End",
}

_DIA_RENAME_DIANN_ION: dict[str, str] = {
    "Protein.Group": "Protein ID",
    "Stripped.Sequence": "Peptide Sequence",
    "Modified.Sequence": "Modified Sequence",
    "Precursor.Charge": "Charge",
    "Protein.Ids": "Protein",
    "Protein.Names": "Entry Name",
    "Genes": "Gene",
    "First.Protein.Description": "Protein Description",
}

_DIA_RENAME_DIANN_PROTEIN: dict[str, str] = {
    "Protein.Group": "Protein ID",
    "Protein.Names": "Entry Name",
    "Genes": "Gene",
    "First.Protein.Description": "Protein Description",
}

_FLIPPR_ION_COLUMNS: list[str] = [
    "Protein ID",
    "Gene",
    "Entry Name",
    "Protein Description",
    "Peptide Sequence",
    "Modified Sequence",
    "Prev AA",
    "Next AA",
    "Start",
    "End",
]

_FLIPPR_PROTEIN_COLUMNS: list[str] = [
    "Protein ID",
    "Gene",
    "Entry Name",
    "Protein Probability",
]

_FLIPPR_PROTEIN_SUMMARY_COLUMNS: list[str] = [
    "Gene",
    "Entry Name",
]

_FLIPPR_CUT_SITE_COLUMNS: list[str] = [
    "Cut Site",
    "Gene",
    "Entry Name",
    "Protein Description",
    "Half Tryptic",
]

_FLIPPR_PEPTIDE_COLUMNS: list[str] = [
    "Prev AA",
    "Start AA",
    "End AA",
    "Next AA",
    "Start",
    "End",
    "Gene",
    "Entry Name",
    "Protein Description",
    "Half Tryptic",
    "Cleavage Type",
]

_FLIPPR_COMBINE_KEY: dict[str, list[str]] = {
    "CUT SITE": _FLIPPR_CUT_SITE_COLUMNS,
    "PEPTIDE": _FLIPPR_PEPTIDE_COLUMNS,
    "MODIFIED PEPTIDE": _FLIPPR_PEPTIDE_COLUMNS,
}

# Thank you Holehouse lab!
_STANDARD_AA_CONVERSION: dict[str, str] = {
    "B": "N",
    "U": "C",
    "X": "G",
    "Z": "Q",
    "*": "",
    "-": "",
    " ": "",
}
