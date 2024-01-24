import polars as pl
import asyncio
import httpx
import json
import Bio.SeqUtils.ProtParam
import pprint


async def uniprot_api(acc):
    async with httpx.AsyncClient(
        base_url="https://rest.uniprot.org/uniprotkb", timeout=None
    ) as client:
        response = await client.get(f"{acc}.json")
        return response.text


async def uniprotkb(accessions):
    tasks = [uniprot_api(acc) for acc in accessions]
    return await asyncio.gather(*tasks)


accessions = ["P04949", "P04994", "P05523", "P06611", "U3PVA8"]

aa = asyncio.run(uniprotkb(accessions))

pprint.pprint(json.loads(aa[1]))
