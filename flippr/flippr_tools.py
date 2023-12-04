# ---------------------------------------------------------------------
# FragPipe Limited-Proteolysis Processor by Edgar Manriquez-Sandoval
# Stephen D. Fried Lab @ The Johns Hopkins University
# ---------------------------------------------------------------------

import numpy as np
from statistics import mode
from scipy.stats import combine_pvalues
from scipy.stats import false_discovery_control
from _flippr_io._params import CUTSITE_IDS, PEPTIDE_IDS


class Peptide(dict):
    def __init__(self, id: str, headers: list, ions: list) -> None:
        self.ions = ions
        for key in headers:
            self[key] = self.ions[0].get(key, None)

        self["No. of Ions"] = len(self.ions)
        self["Item ID"] = id

        self.update(
            {
                "Log2 FC": None,
                "Normalized Log2 FC": None,
                "P-Value": None,
                "Log10 P-Value": None,
                "Adj. P-Value": None,
                "Log10 Adj. P-Value": None,
                "CV": None,
            }
        )

        self.calculate_ratio()

    @property
    def ratio(self):
        return self.get("Log2 FC", None)

    @property
    def norm_ratio(self):
        return self.get("Normalized Log2 FC", None)

    @property
    def pval(self):
        return self.get("P-Value", None)

    @property
    def log10pval(self):
        return self.get("Log10 P-Value", None)

    @property
    def adjpval(self):
        return self.get("Adj. P-Value", None)

    @property
    def log10adjpval(self):
        return self.get("Log10 Adj. P-Value", None)

    def normalize_ratio(self, norm_proteins: dict):
        key = self["Protein ID"]
        ratio = self.ratio

        if key in norm_proteins.keys():
            norm_protein = norm_proteins.get(key)
            ratio = self.ratio - norm_protein.ratio

        self.update({"Normalized Log2 FC": ratio})

    def calculate_ratio(self) -> None:
        ratio = np.array([ion.ratio for ion in self.ions])
        pval = np.array([ion.pval for ion in self.ions])
        adjpval = np.array([ion.adjpval for ion in self.ions])
        var = np.array([ion.var for ion in self.ions])
        ttest = np.array([ion.ttest for ion in self.ions])
        ttest_sign = np.sign(ttest)
        ttest_sign_mode = 0

        # when there are 3 or more tests to compare...
        if len(ttest) > 2:
            # ... ensure that the mode will not throw a StatisticsError
            # this happens for even len(ttest) with equal number of signs + and -
            if sum(ttest_sign) != 0:
                ttest_sign_mode = mode(ttest_sign)

        # if there is a single t-test, accept the p-value
        if len(ttest) == 1:
            median_ratio, max_var, pval, adjpval = ratio[0], var[0], pval[0], adjpval[0]
            log10pval, log10adjpval = np.abs(np.log10(pval)), np.abs(np.log10(adjpval))
        else:
            median_ratio, max_var = np.median(ratio), np.max(var)
            if len(ttest) == 2:
                if all(t > 0 for t in ttest) or all(t < 0 for t in ttest):
                    chisq, pval = combine_pvalues(pval, method="fisher")
                    chisq, adjpval = combine_pvalues(adjpval, method="fisher")
                    log10pval, log10adjpval = (
                        np.abs(np.log10(pval)),
                        np.abs(np.log10(adjpval)),
                    )
                else:
                    pval, adjpval = 1, 1
                    log10pval, log10adjpval = 0, 0
            elif len(ttest) > 2 and ttest_sign_mode != 0:
                filtered_vals = [
                    (r, p, q, v)
                    for r, p, q, v, t in zip(ratio, pval, adjpval, var, ttest)
                    if np.sign(t) == ttest_sign_mode
                ]
                ratio, pval, adjpval, var = zip(*filtered_vals)
                median_ratio, max_var = np.median(ratio), np.max(var)
                chisq, pval = combine_pvalues(pval, method="fisher")
                chisq, adjpval = combine_pvalues(adjpval, method="fisher")
                log10pval, log10adjpval = (
                    np.abs(np.log10(pval)),
                    np.abs(np.log10(adjpval)),
                )
            else:
                pval, adjpval = 1.0, 1.0
                log10pval, log10adjpval = 0.0, 0.0

        self.update(
            {
                "Log2 FC": median_ratio,
                "P-Value": pval,
                "Log10 P-Value": log10pval,
                "Adj. P-Value": adjpval,
                "Log10 Adj. P-Value": log10adjpval,
                "CV": max_var,
            }
        )


class Protein(dict):
    def __init__(
        self, id: str, modified_peptides: list, peptides: list, cutsites: list
    ) -> None:
        self.modified_peptides = modified_peptides
        self.peptides = peptides
        self.cutsites = cutsites

        self["Protein ID"] = id
        self["Gene"] = self.peptides[0]["Gene"]
        self["No. of Ions"] = sum(pep["No. of Ions"] for pep in self.peptides)
        self["No. of Modified Peptides"] = len(self.modified_peptides)
        self["No. of Peptides"] = len(self.peptides)
        self["No. of Cut Sites"] = len(self.cutsites)

    def normalize_ratios(self, norm_proteins: dict):
        for mod_pep in self.modified_peptides:
            mod_pep.normalize_ratio(norm_proteins)

        for pep in self.peptides:
            pep.normalize_ratio(norm_proteins)

        for cutsite in self.cutsites:
            cutsite.normalize_ratio(norm_proteins)

    def count_significants(self):
        sig_count = {
            "No. of Modified Peptides (Significant P-Value)": 0,
            "No. of Modified Peptides (Valid P-Value)": 0,
            "No. of Peptides (Significant P-Value)": 0,
            "No. of Peptides (Valid P-Value)": 0,
            "No. of Cut Sites (Significant P-Value)": 0,
            "No. of Cut Sites (Valid P-Value)": 0,
            "No. of Modified Peptides (Significant Adj. P-Value)": 0,
            "No. of Modified Peptides (Valid Adj. P-Value)": 0,
            "No. of Peptides (Significant Adj. P-Value)": 0,
            "No. of Peptides (Valid Adj. P-Value)": 0,
            "No. of Cut Sites (Significant Adj. P-Value)": 0,
            "No. of Cut Sites (Valid Adj. P-Value)": 0,
        }

        def p_value_counter(peptides, sig, nonsig):
            for peptide in peptides:
                if peptide.log10pval != 0:
                    sig_count[nonsig] += 1
                    # confidence ratio cut-off R >= 1 P >= 2; high confidence ratio cut-off R >= 6 P >= 1.8
                    if (np.abs(peptide.ratio) >= 1 and peptide.log10pval >= 2) or (
                        np.abs(peptide.ratio) >= 6 and peptide.log10pval >= 1.8
                    ):
                        sig_count[sig] += 1

        def adj_p_value_counter(peptides, sig, nonsig):
            for peptide in peptides:
                if peptide.adjpval <= 0.10:
                    sig_count[nonsig] += 1
                    # confidence ratio cut-off R >= 1 Adj. P <= 5% FDR
                    if np.abs(peptide.ratio) >= 1 and peptide.adjpval <= 0.05:
                        sig_count[sig] += 1

        p_value_counter(
            self.modified_peptides,
            "No. of Modified Peptides (Significant P-Value)",
            "No. of Modified Peptides (Valid P-Value)",
        )

        p_value_counter(
            self.peptides,
            "No. of Peptides (Significant P-Value)",
            "No. of Peptides (Valid P-Value)",
        )

        p_value_counter(
            self.cutsites,
            "No. of Cut Sites (Significant P-Value)",
            "No. of Cut Sites (Valid P-Value)",
        )

        adj_p_value_counter(
            self.modified_peptides,
            "No. of Modified Peptides (Significant Adj. P-Value)",
            "No. of Modified Peptides (Valid Adj. P-Value)",
        )

        adj_p_value_counter(
            self.peptides,
            "No. of Peptides (Significant Adj. P-Value)",
            "No. of Peptides (Valid Adj. P-Value)",
        )

        adj_p_value_counter(
            self.cutsites,
            "No. of Cut Sites (Significant Adj. P-Value)",
            "No. of Cut Sites (Valid Adj. P-Value)",
        )

        self.update(sig_count)


def __generate_key_sorted_ions(ions: list, key: str, headers: list):
    key_sorted_ions = {}
    for ion in ions:
        t = key_sorted_ions.setdefault(ion[key], [])
        t.append(ion)

    key_sorted_ions = [
        Peptide(key, headers, ions) for key, ions in key_sorted_ions.items()
    ]

    return key_sorted_ions


def __sort_by_protein_id(dict_objects: list):
    sorted_dict_object = {}
    for obj in dict_objects:
        t = sorted_dict_object.setdefault(obj["Protein ID"], [])
        t.append(obj)

    return sorted_dict_object


def generate_proteins(ions: list):
    ions_per_protein = __sort_by_protein_id(ions)
    for protein in ions_per_protein.keys():
        ions = ions_per_protein[protein]
        ion_p_vals = [ion.pval for ion in ions]
        ion_adj_p_vals = false_discovery_control(ion_p_vals)
        for ion, pval in zip(ions, ion_adj_p_vals):
            ion.update(
                {"Adj. P-Value": pval, "Log10 Adj. P-Value": np.abs(np.log10(pval))}
            )

    ions = [ion for ions in ions_per_protein.values() for ion in ions]

    modified_peptides = __generate_key_sorted_ions(
        ions, "Modified Sequence", PEPTIDE_IDS
    )
    peptides = __generate_key_sorted_ions(ions, "Peptide Sequence", PEPTIDE_IDS)
    cutsites = __generate_key_sorted_ions(ions, "Cut Site ID", CUTSITE_IDS)

    prot_modified_peptides = __sort_by_protein_id(modified_peptides)
    prot_peptides = __sort_by_protein_id(peptides)
    prot_cutsites = __sort_by_protein_id(cutsites)

    proteins = []
    for protein_id in prot_peptides.keys():
        modified_peptides = prot_modified_peptides[protein_id]
        peptides = prot_peptides[protein_id]
        cutsites = prot_cutsites[protein_id]
        proteins.append(Protein(protein_id, modified_peptides, peptides, cutsites))

    return proteins
