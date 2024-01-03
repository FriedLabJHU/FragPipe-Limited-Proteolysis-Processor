from pathlib import Path
import flippr as fp

lip = Path("LFQ/LiP")
trp = Path("LFQ/Control")

study = fp.Study(lip=lip, trp=trp)

print(study.samples)

study.add_process(1  , "Native", "Refolded_001_min", 3, "Native", "Refolded", 3)
study.add_process(5  , "Native", "Refolded_005_min", 3, "Native", "Refolded", 3)
study.add_process(120, "Native", "Refolded_120_min", 3, "Native", "Refolded", 3)

study.run()

for result in study.results.values():
    res = \
    [result.ion,
    result.peptide,
    result.modified_peptide,
    result.cut_site,
    result.protein_summary]

    print(res)