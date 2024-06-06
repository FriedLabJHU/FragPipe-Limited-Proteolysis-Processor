import asyncio

import httpx
import nest_asyncio
import polars as pl
from Bio.SeqUtils.ProtParam import ProteinAnalysis

nest_asyncio.apply()

def n_batch(x: pl.Series, n: int):
    for i in range(0, len(x), n):
        yield x[i : i + n]


async def uniprot_api(acc: str) -> str:
    async with httpx.AsyncClient(
        base_url="https://rest.uniprot.org/uniprotkb", timeout=None
    ) as client:
        response = await client.get(f"{acc}.json")
        return response.text


async def uniprotkb(accessions: pl.Series) -> list:
    tasks = [uniprot_api(acc) for acc in accessions]
    return await asyncio.gather(*tasks)


# Parse JSON inside of this function
def uniprot_json(accessions: pl.Series):
    for acc in n_batch(accessions, 200):
        accessions_json = asyncio.run(uniprotkb(acc))


def bio_metadata(fasta: pl.DataFrame) -> pl.DataFrame:
    return fasta.with_columns(
        pl.col("Sequence").map_elements(lambda x: ProteinAnalysis(x), return_dtype=pl.Object)
        .alias("Analysis"),
    ).with_columns(
        pl.col("Analysis").map_elements(lambda x: x.isoelectric_point(), return_dtype=pl.Float64)
        .alias("pI"),
        pl.col("Analysis").map_elements(lambda x: x.molecular_weight(), return_dtype=pl.Float64)
        .alias("MW (Da)"),
        pl.col("Sequence").str.len_chars()
        .alias("Len (aa)"),
    ).drop(["Analysis", "Valid ID"])
