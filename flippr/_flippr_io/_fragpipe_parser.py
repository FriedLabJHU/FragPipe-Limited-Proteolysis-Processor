# ---------------------------------------------------------------------
# FragPipe Limited-Proteolysis Processor by Edgar Manriquez-Sandoval
# Stephen D. Fried Lab @ The Johns Hopkins University
# ---------------------------------------------------------------------

from csv import reader
from pathlib import Path
from ast import literal_eval


class FragPipeObject(dict):
    """
    `FragPipeObject` is an inherited `dict` object with easily accessible data
    from FragPipe output files.
    """

    def __init__(self, header, fp_tsv_data) -> None:
        """
        The dictionary requires as input the header annotations & data rows from a
        FragPipe tab separated variable file.

        Type casting will be automatically performed for valid data types; otherwise
        `None` will be returned.

        Parameters
        ----------
        header : <str>
            FragPipe tab separated variable file header annotations.
        fp_tsv_data : <list>
            FragPipe tab separated variables parsed into a list.
        """

        # dictionary is being created with header annotation & tsv data
        # the tsv data list is being type cast with the internal method `__type_cast`
        for key, val in zip(header, map(self.__type_cast, fp_tsv_data)):
            self[key] = val

        # making a convenient set of variables that are easily accessible
        self.seq = self.get("Peptide Sequence", None)
        self.mod_seq = self.get("Modified Sequence", None)
        self.mods = self.get("Assigned Modifications", None)
        self.prev_aa = self.get("Prev AA", None)
        self.next_aa = self.get("Next AA", None)
        self.start = self.get("Start", None)
        self.end = self.get("End", None)
        self.pep_len = self.get("Peptide Length", None)
        self.protein = self.get("Protein", None)
        self.protein_id = self.get("Protein ID", None)
        self.gene = self.get("Gene", None)
        self.ttest_alt = None

        # Only perform this for peptide-like FragPipeObjects
        if self.seq is not None:
            self.__set_trypticity()
            self.isht = self.get("Half Tryptic", None)
            self.cutsite = self.get("Cut Site", None)
            self.cutsite_id = self.get("Cut Site ID", None)

        # Adding for general use
        self.update(
            {
                "Log2 FC": None,
                "P-Value": None,
                "Log10 P-Value": None,
                "Adj. P-Value": None,
                "Log10 Adj. P-Value": None,
                "T-test": None,
                "CV": None,
            }
        )

    @property
    def ratio(self):
        return self.get("Log2 FC")

    @property
    def pval(self):
        return self.get("P-Value")

    @property
    def log10pval(self):
        return self.get("Log10 P-Value")

    @property
    def adjpval(self):
        return self.get("Adj. P-Value")

    @property
    def log10adjpval(self):
        return self.get("Log10 Adj. P-Value")

    @property
    def ttest(self):
        return self.get("T-test")

    @property
    def var(self):
        return self.get("CV")

    def __set_trypticity(self) -> None:
        """
        This function will determine if the precursor peptide is half tryptic or not.
        """

        def __is_K_or_R(aminoacid):
            return aminoacid == "K" or aminoacid == "R"

        def __return_results(isht, cutsite_pos, cutsite_res, cutsite):
            self.update(
                {
                    "Half Tryptic": isht,
                    "Cut Site Position": cutsite_pos,
                    "Cut Site Residue": cutsite_res,
                    "Cut Site": cutsite,
                    "Cut Site ID": f"{self.protein_id}_{cutsite}",
                }
            )
            return

        isht, cutsite_pos, cutsite_res, cutsite = None, None, None, None

        main_pep = self.seq
        n_aa = main_pep[0]
        c_aa = main_pep[-1]

        nn_aa = self.prev_aa
        cc_aa = self.next_aa

        if nn_aa == "-" and __is_K_or_R(c_aa):
            isht = False
            cutsite = f"[{self.start},{self.end}]"

            __return_results(isht, cutsite_pos, cutsite_res, cutsite)

        elif nn_aa == "-" and not __is_K_or_R(c_aa):
            isht = True
            cutsite_pos = self.end + 1
            cutsite_res = cc_aa
            cutsite = f"{cc_aa}{self.end+1}"

            __return_results(isht, cutsite_pos, cutsite_res, cutsite)

        elif cc_aa == "-" and __is_K_or_R(nn_aa):
            isht = False
            cutsite = f"[{self.start},{self.end}]"

            __return_results(isht, cutsite_pos, cutsite_res, cutsite)

        elif cc_aa == "-" and not __is_K_or_R(nn_aa):
            isht = True
            cutsite_pos = self.start
            cutsite_res = n_aa
            cutsite = f"{n_aa}{self.start}"

            __return_results(isht, cutsite_pos, cutsite_res, cutsite)

        elif not __is_K_or_R(c_aa) and __is_K_or_R(nn_aa):
            isht = True
            cutsite_pos = self.end + 1
            cutsite_res = cc_aa
            cutsite = f"{cc_aa}{self.end+1}"

            __return_results(isht, cutsite_pos, cutsite_res, cutsite)

        elif __is_K_or_R(c_aa) and not __is_K_or_R(nn_aa):
            isht = True
            cutsite_pos = self.start
            cutsite_res = n_aa
            cutsite = f"{n_aa}{self.start}"

            __return_results(isht, cutsite_pos, cutsite_res, cutsite)

        else:
            isht = False
            cutsite = f"[{self.start},{self.end}]"

            __return_results(isht, cutsite_pos, cutsite_res, cutsite)

    def __type_cast(self, obj: str) -> object:
        """
        Perform type casting on the given object.

        This method attempts to convert the object to its appropriate type.
        If the object can be evaluated as a literal the converted value is returned.
        If the object cannot be evaluated or is empty, the original object is returned.
        If the object is empty, None is returned.

        Parameters
        ----------
        obj <object> :
            The object to be type casted.

        Returns
        -------
        <object> :
            The type-casted object, or the original object if type casting
            is not possible or if the object is empty.
        """

        if obj:
            try:
                return literal_eval(obj)
            except:
                return obj
        return None


class FragPipeFileParser(object):
    """
    `FragPipeFileParser` is an iterator for reading FragPipe files.
    """

    def __init__(self, filename: str | Path) -> None:
        """
        Initialize the `FragPipeFileParser`.

        Parameters
        ----------
        filename <str | Path> :
            The filename or path to the FragPipe file.
        """

        self.filename = Path(filename)
        self.open_file = reader(open(filename, "r"), delimiter="\t")
        self.header = next(self.open_file)

    def __iter__(self) -> object:
        """
        Return the iterator object.

        Returns
        -------
        <FragPipeFileParser> :
            The iterator object.
        """
        return self

    def __next__(self) -> FragPipeObject:
        """
        Get the next `FragPipeObject`.

        Returns
        -------
        <FragPipeObject> :
            The next `FragPipeObject`.
        """

        return FragPipeObject(self.header, next(self.open_file))
