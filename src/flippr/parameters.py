rcParams: dict = {
    "ion.missing_intensity_thresh": 1,
    "ion.impute_type": "gaussian",
    "ion.aon_impute_loc": 1e4,
    "ion.aon_impute_scale": 1e3,
    "trp_protein.intensity_value": "MaxLFQ Intensity",
    "trp_protein.fc_sig_tresh": 1.0,
    "trp_protein.pval_sig_tresh": 0.01,
    "protein.fc_sig_thresh": 1.0,
    "protein.pval_sig_thresh": 0.01,
    "protein.adj_pval_sig_thresh": 0.05,
    "protein.high_fc_sig_thresh": 6.0,
    "protein.high_pval_sig_thresh": 0.016,
}

_LFQ_FP_FILES: list[str] = [
    "combined_ion.tsv",
    "combined_protein.tsv",
    "experiment_annotation.tsv",
]

_LFQ_FP_CONSTANT_ION_COLUMNS: list[str] = [
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

_LFQ_FP_VARIABLE_ION_COLUMNS: list[str] = [
    "Spectral Count",
    "Apex Retention Time",
    "Intensity",
    "Match Type",
]

_LFQ_FP_CONSTANT_PROTEIN_COLUMNS: list[str] = [
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

_LFQ_FP_VARIABLE_PROTEIN_COLUMNS: list[str] = [
    "Spectral Count",
    "Intensity",
    "MaxLFQ Intensity",
]

_FLIPPR_PROTEIN_COLUMNS: list[str] = [
    "Protein",
    "Protein ID",
    "Entry Name",
    "Gene",
    "Protein Length",
    "Organism",
    "Protein Existence",
    "Description",
    "Protein Probability",
]


_FLIPPR_PROTEIN_SUMMARY_COLUMNS: list[str] = [
    "Protein",
    "Entry Name",
    "Gene",
    "Protein Description",
    "Mapped Genes",
    "Mapped Proteins",
]

_FLIPPR_CUT_SITE_COLUMNS: list[str] = [
    "Cut Site",
    "Gene",
    "Entry Name",
    "Protein Description",
    "Mapped Genes",
    "Mapped Proteins",
    "Half Tryptic",
]

_FLIPPR_PEPTIDE_COLUMNS: list[str] = [
    "Prev AA",
    "Start AA",
    "End AA",
    "Next AA",
    "Start",
    "End",
    "Peptide Length",
    "Gene",
    "Entry Name",
    "Protein Description",
    "Mapped Genes",
    "Mapped Proteins",
    "Half Tryptic",
    "Cleavage Type",
]

_FLIPPR_MODIFIED_PEPTIDE_COLUMNS: list[str] = [
    "Prev AA",
    "Start AA",
    "End AA",
    "Next AA",
    "Start",
    "End",
    "Peptide Length",
    "Gene",
    "Entry Name",
    "Protein Description",
    "Mapped Genes",
    "Mapped Proteins",
    "Half Tryptic",
    "Cleavage Type",
]

_FLIPPR_COMBINE_KEY: dict[str, list] = {
    "CUT SITE": _FLIPPR_CUT_SITE_COLUMNS,
    "PEPTIDE": _FLIPPR_PEPTIDE_COLUMNS,
    "MODIFIED PEPTIDE": _FLIPPR_MODIFIED_PEPTIDE_COLUMNS,
}

# Thank you Holehouse lab!
_STANDARD_AA_CONVERSION = {
    "B": "N",
    "U": "C",
    "X": "G",
    "Z": "Q",
    "*": "",
    "-": "",
    " ": "",
}
