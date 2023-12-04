# ---------------------------------------------------------------------
# FragPipe Limited-Proteolysis Processor by Edgar Manriquez-Sandoval
# Stephen D. Fried Lab @ The Johns Hopkins University
# ---------------------------------------------------------------------

import argparse
import multiprocessing

from _flippr_io._fragpipe_parser import FragPipeFileParser
from _flippr_io._params import Params
import _flippr_io._msg as flippr_io
import flippr_tools as flippr
import fragpipe_ion as fp_ion
import fragpipe_protein as fp_prot
from _flippr_io._params import CUTSITE_IDS, PEPTIDE_IDS


parser = argparse.ArgumentParser()

parser.add_argument("-n", "--replicates", nargs=1, type=int, required=True)
parser.add_argument("-f", "--filename", type=str, required=True)
parser.add_argument("-c", "--control", nargs=1, type=str, required=True)
parser.add_argument(
    "-d", "--conditions", nargs="+", action="append", type=str, required=True
)
parser.add_argument("-z", "--norm_filename", type=str, required=False)
parser.add_argument("-zc", "--norm_control", nargs=1, type=str, required=False)
parser.add_argument(
    "-zd", "--norm_conditions", nargs="+", action="append", type=str, required=False
)
parser.add_argument("-p", "--prefix", type=str, required=False)


# for cond in params.condition_annotations:
def main(cond, params):
    # supppppppppper slow step for very large files :(
    # hope you have a lot of memory, will improve later
    ion_data = FragPipeFileParser(params.ion_file)

    ion_data = list(
        filter(
            lambda x: x is not None,
            [
                fp_ion.validate_ion_abundances(ion, cond=cond, params=params)
                for ion in ion_data
            ],
        )
    )

    ion_data = [
        fp_ion.calculate_ion_ratio(ion, cond=cond, params=params) for ion in ion_data
    ]

    proteins = flippr.generate_proteins(ion_data)

    # normalize proteins
    if params.norm_file is not None:
        norm_protein_data = list(FragPipeFileParser(params.norm_file))

        norm_protein_data = list(
            filter(
                lambda x: x is not None,
                [
                    fp_prot.validate_max_lfq_abundances(
                        norm_protein, cond=cond, params=params
                    )
                    for norm_protein in norm_protein_data
                ],
            )
        )

        norm_protein_data = [
            fp_prot.calculate_max_lfq_ratio(norm_protein, cond=cond, params=params)
            for norm_protein in norm_protein_data
        ]

        norm_proteins = fp_prot.generate_normalization_ratios(norm_protein_data)

        [protein.normalize_ratios(norm_proteins) for protein in proteins]

        flippr_io.output_norm_protein_data(norm_proteins, cond, params)

        del norm_proteins

    [protein.count_significants() for protein in proteins]

    flippr_io.output_ion_data(ion_data, cond, params)

    flippr_io.output_protein_data(proteins, cond, params)

    collect_peptides = [pep for prot in proteins for pep in prot.modified_peptides]
    flippr_io.output_peptide_data(
        collect_peptides, PEPTIDE_IDS, "modified_peptides", cond, params
    )

    collect_peptides = [pep for prot in proteins for pep in prot.peptides]
    flippr_io.output_peptide_data(
        collect_peptides, PEPTIDE_IDS, "peptides", cond, params
    )

    collect_peptides = [pep for prot in proteins for pep in prot.cutsites]
    flippr_io.output_peptide_data(
        collect_peptides, CUTSITE_IDS, "cutsites", cond, params
    )

    del ion_data

    del proteins


if __name__ == "__main__":
    args = parser.parse_args()

    params = Params(args)
    params.validate_args()

    p = multiprocessing.Pool(4)

    for cond in params.condition_annotations:
        p.apply_async(main, args=(cond, params))

    p.close()
    p.join()
