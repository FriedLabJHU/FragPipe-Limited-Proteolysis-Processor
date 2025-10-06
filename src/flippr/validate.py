from pathlib import Path
from warnings import warn
from typing import Optional, Literal, cast

from .parameters import _DDA_FP_FILES, _DIA_FP_FILES

def _validate_study(
    lip: str | Path, 
    trp: Optional[str | Path], 
    method: str
) -> tuple[Path, Optional[Path], str]:
    method = __validate_method(method)

    lip = __validate_fragpipe_path(lip, "lip", method)

    if trp is not None:
        trp = __validate_fragpipe_path(trp, "trp", method)

    return lip, trp, method


def __validate_method(method: str) -> str:
    """
    Validate the quantification method.

    """

    if not isinstance(method, str):
        raise TypeError(
            f'`method` was provided: "{method}" with type `{type(method)}`. The type `{type(method)}`is not recognized. Set `method` to "dda" or "dia".'
        )

    if method not in ["dda", "dia"]:
        raise ValueError(
            f'`method` was provided: "{method}". "{method}" is not recognized. Set `method` to "dda" or "dia".'
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

    if method == "dda":
        if not all([path.joinpath(f).exists() for f in _DDA_FP_FILES]):
            raise FileNotFoundError(
                f'Files not found in "{path}". The FragPipe output directory path should minimally contain: \n'
                + "\t\n".join([f"`{f}`" for f in _DDA_FP_FILES])
            )

    if method == "dia":
        if not all([path.joinpath(f).exists() for f in _DIA_FP_FILES]):
            raise FileNotFoundError(
                f'Files not found in "{path}". The FragPipe output directory path should minimally contain: \n'
                + "\t\n".join([f"`{f}`" for f in _DIA_FP_FILES])
            )


def _validate_replicate(replicate: int | tuple[int, int] | tuple[tuple[int, ...], tuple[int, ...]]) -> Literal["int", "tuple", "tuple_tuple"] | None:
    """
    Validate the replicate inputs.
    
    """

    def __int_validation(i: int) -> None:
        if i <= 1:
            raise ValueError(f'Number of replicates was set to `{i}`. Replicate values must be positive and greater than or equal to `2`.')
        
        if i == 2:
            warn("WARNING: Using two replicates will result in an under-powered study.", UserWarning)
            warn("WARNING: Set flippr.rcParams `ion.missing_intensity_thresh` to `0` for the best results.", UserWarning)

        
    if isinstance(replicate, int):
        __int_validation(replicate)
        return "int"

    if isinstance(replicate, tuple):
        replicate = cast(tuple, replicate)

        if len(replicate) != 2:
            raise ValueError(f'Replicate input contains `{len(replicate)}` `{type(replicate)}`. Only `int`, `tuple[int, int]`, or `tuple[tuple[int, ...], tuple[int, ...]]` are allowed.')

        if all(isinstance(i, int) for i in replicate):
            replicate = cast(tuple[int, int], replicate)
            
            (__int_validation(i) for i in replicate)
            return "tuple"

            
        if all(isinstance(i, tuple) for i in replicate):
            replicate = cast(tuple[tuple[int, ...], tuple[int, ...]], replicate)

            if all(all(isinstance(i, int) for i in tup) for tup in replicate):
                (__int_validation(len(tup)) for tup in replicate)
                return "tuple_tuple"
        
    raise TypeError(f'Replicate input is of type `{type(replicate)}`. Only `int`, `tuple[int, int]`, or `tuple[tuple[int, ...], tuple[int, ...]]` are allowed.')
