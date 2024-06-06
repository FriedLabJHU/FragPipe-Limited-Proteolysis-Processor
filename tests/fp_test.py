from pathlib import Path

import flippr as fp

lip = Path("LFQ/LiP")
trp = Path("LFQ/TrP")
fasta = Path("UP000000625_83333.fasta")

fp.rcParams.update({"ion.aon_impute_loc": 1e55})
study = fp.Study(lip=lip, trp=trp, fasta=fasta)

print(study.samples)

study.add_process(1  , "Native", "Refolded_001_min", 3, "Native", "Refolded", 3)
study.add_process(5  , "Native", "Refolded_005_min", 3, "Native", "Refolded", 3)
study.add_process(120, "Native", "Refolded_120_min", 3, "Native", "Refolded", 3)

study.run()

for result in study.results.values():
    name = result.name

    result.protein_summary.write_parquet(f"{name}_protein_summary.parquet")
    result.ion.write_parquet(f"{name}_ion.parquet")
    result.peptide.write_parquet(f"{name}_peptide.parquet")
    result.modified_peptide.write_parquet(f"{name}_modified_peptide.parquet")
    result.cut_site.write_parquet(f"{name}_cut_site.parquet")
