rcParams: dict = {
    "ion.missing_intensity_thresh": 1,
    "ion.aon_impute_loc": 1e4,
    "ion.aon_impute_scale": 1e3,
    "protein.fc_sig_thresh": 1.0,
    "protein.pval_sig_thresh": 0.01,
    "protein.adj_pval_sig_thresh": 0.05,
    "protein.fc_high_sig_thresh": 6.0,
    "protein.pval_high_sig_thresh": 0.016,
}

LFQ_FP_FILES: list[str] = [
    "combined_ion.tsv",
    "combined_protein.tsv",
    "experiment_annotation.tsv",
]

LFQ_FP_CONSTANT_ION_COLUMNS: list[str] = [
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

LFQ_FP_VARIABLE_ION_COLUMNS: list[str] = [
    "Spectral Count",
    "Apex Retention Time",
    "Intensity",
    "Match Type",
]

LFQ_FP_CONSTANT_PROTEIN_COLUMNS: list[str] = [
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

LFQ_FP_VARIABLE_PROTEIN_COLUMNS: list[str] = [
    "Spectral Count",
    "Intensity",
    "MaxLFQ Intensity",
]

FLIPPR_PROTEIN_COLUMNS: list[str] = [
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


FLIPPR_PROTEIN_SUMMARY_COLUMNS: list[str] = [
    "Protein ID",
    "Protein",
    "Entry Name",
    "Gene",
    "Protein Description",
    "Mapped Genes",
    "Mapped Proteins",
]

FLIPPR_CUT_SITE_COLUMNS: list[str] = [
    "Cut Site",
    "Gene",
    "Entry Name",
    "Protein Description",
    "Mapped Genes",
    "Mapped Proteins",
    "Half Tryptic",
]

FLIPPR_PEPTIDE_COLUMNS: list[str] = [
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

FLIPPR_MODIFIED_PEPTIDE_COLUMNS: list[str] = [
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

FLIPPR_COMBINE_KEY: dict[str, list] = {
    "CUT SITE": FLIPPR_CUT_SITE_COLUMNS,
    "PEPTIDE": FLIPPR_PEPTIDE_COLUMNS,
    "MODIFIED PEPTIDE": FLIPPR_MODIFIED_PEPTIDE_COLUMNS,
}

# Thank you Holehouse lab!
STANDARD_AA_CONVERSION = {
    "B": "N",
    "U": "C",
    "X": "G",
    "Z": "Q",
    "*": "",
    "-": "",
    " ": "",
}

# TODO: Add SILAC compatibility
# SILAC_FP_FILES = None
# SILAC_FP_COLUMNS = None
