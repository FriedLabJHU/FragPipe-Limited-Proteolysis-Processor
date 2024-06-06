from pathlib import Path

from .parameters import _LFQ_FP_FILES

# from flippr.parameters import SILAC_FP_FILES


def _validate_study(
    lip: str | Path, trp: str | Path | None, fasta: str | Path | None, method: str
) -> tuple[Path, Path | None, Path | None, str]:
    method = __validate_method(method)

    lip = __validate_fragpipe_path(lip, "lip", method)

    if trp is not None:
        trp = __validate_fragpipe_path(trp, "trp", method)

    if fasta is not None:
        fasta = __validate_fasta_file(fasta)

    return lip, trp, fasta, method


def __validate_method(method: str) -> str:
    """
    Validate the quantification method.
    """

    if not isinstance(method, str):
        raise TypeError(
            f'`method` was provided: "{method}" with type `{type(method)}`. The type `{type(method)}`is not recognized. Set `method` to "lfq" or "silac".'
        )

    if method not in ["lfq", "silac"]:
        raise ValueError(
            f'`method` was provided: "{method}". "{method}" is not recognized. Set `method` to "lfq" or "silac".'
        )

    return method


def __validate_fragpipe_path(path: str | Path, liptrp: str, method: str) -> Path:
    """
    Validate FragPipe output directory paths.
    """

    if not isinstance(path, (str, Path)):
        raise TypeError(
            f'`{liptrp}` was provided: "{path}" with type `{type(path)}`. The type `{type(path)}` is not recognized. Set `{liptrp}` to a FragPipe output path.'
        )

    # type cast to `Path`
    # this should have no effect on `Path`
    path = Path(path)

    # invokes `.exists` and checks for directory-ness
    if not path.is_dir():
        raise ValueError(
            f'`{liptrp}` was provided: "{path}". "{path}" is not a directory path. Set `{liptrp}` to a FragPipe output directory path.'
        )

    __validate_fragpipe_files(path, method)

    return Path(path)


def __validate_fragpipe_files(path: Path, method: str) -> None:
    """
    Validate FragPipe output directory path contains all the necessay files.
    """

    if method == "lfq":
        if not all([path.joinpath(f).exists() for f in _LFQ_FP_FILES]):
            raise FileNotFoundError(
                f'Files not found in "{path}". The FragPipe output directory path should minimally contain: '
                + "\t\n".join([f"`{f}`" for f in _LFQ_FP_FILES])
            )

    if method == "silac":
        raise ValueError("This feature does not work yet")
        # if not all([path.joinpath(f).exists() for f in SILAC_FP_FILES]):
        #     raise FileNotFoundError(f"Files not found in \"{path}\". The FragPipe output directory path should minimally contain: " + "\t\n".join([f"`{f}`" for f in SILAC_FP_FILES]))


def __validate_fasta_file(fasta: str | Path) -> Path:
    """
    Validate Uniprot FASTA file.
    """

    if not isinstance(fasta, (str, Path)):
        raise TypeError(
            f'`fasta` was provided: "{fasta}" with type `{type(fasta)}`. The type `{type(fasta)}` is not recognized. Set `fasta` to a Uniprot FASTA file.'
        )

    # type cast to `Path`
    # this should have no effect on `Path`
    fasta = Path(fasta)

    if not fasta.exists():
        raise FileNotFoundError(
            f'"{fasta.name}" was not found in {fasta.parent.resolve()}. Check that your file extension is `.fasta` or `.faa`. Set `fasta` to a Uniprot FASTA file.'
        )

    if fasta.suffix.lower() not in [".fasta", ".faa"]:
        raise ValueError(
            f'"{fasta.name}" has extension `{fasta.suffix}`. Check that your file extension is `.fasta` or `.faa`. Set `fasta` to a Uniprot FASTA file.'
        )

    return Path(fasta)
