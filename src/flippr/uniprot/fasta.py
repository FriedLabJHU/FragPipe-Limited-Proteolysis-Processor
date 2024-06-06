import re
from pathlib import Path

import polars as pl
from Bio.SeqIO.FastaIO import SimpleFastaParser

from flippr.parameters import _STANDARD_AA_CONVERSION

# UniProtKB accession number validation regex (https://www.uniprot.org/help/accession_numbers)
ACCESSION_RE = re.compile(
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
)


def _read_fasta(fasta_path: str | Path) -> pl.DataFrame:
    """
    This function reads and parses Uniprot FASTA files and returns a Polars dataframe.

    Args:
        fasta_path: Path to a Uniprot FASTA file.

    Returns:
        fasta_df: Polars dataframe containing Uniprot FASTA entry data.

    Raises:
        FileNotFoundError: If the `fasta_path` does not exsist.

        ValueError: If the `fasta_path` does not point to a file with `.fasta` or `.faa` extensions.
    """

    fasta_dict = dict()
    with open(fasta_path) as open_fasta:
        fasta_gen = SimpleFastaParser(open_fasta)

        for faa in fasta_gen:
            acc = __extract_accession(faa[0])
            seq = faa[1].upper()

            for old_aa, new_aa in _STANDARD_AA_CONVERSION.items():
                seq = seq.replace(old_aa, new_aa)

            fasta_dict.update({acc: seq})

    fasta_df = pl.from_dict(
        {"Protein ID": list(fasta_dict.keys()), "Sequence": list(fasta_dict.values())}
    )

    fasta_df = fasta_df.with_columns(
        pl.col("Protein ID").map_elements(__validate_accession, return_dtype=pl.Boolean).alias("Valid ID")
    )

    return fasta_df


def __extract_accession(fasta_header: str) -> str:
    """
    This function extracts the accession number from a Uniprot FASTA header line.

    e.g.
    The Uniprot FASTA header...
        >sp|P08200|IDH_ECOLI Isocitrate dehydrogenase
    produced the accession number...
        P08200

    Most other FASTA headers will return the whole header string.

    Args:
        fasta_header: First line from a two-line Uniprot FASTA entry.

    Returns:
        Accession number string | FASTA header string.
    """

    if "|" not in fasta_header:
        return fasta_header

    try:
        return fasta_header.split("|")[1]

    except IndexError:
        return fasta_header


def __validate_accession(accession: str) -> bool:
    """
    This function will return `True` if an accession number fully matches the Uniprot accession number regex.

    For more information, visit https://www.uniprot.org/help/accession_numbers

    Args:
        accession: Accession number string | FASTA header string.

    Returns:
        True | False
    """

    if re.fullmatch(ACCESSION_RE, accession) is not None:
        return True

    return False
