from pathlib import Path
import flippr as fp

lip = Path("LFQ/LiP")
trp = Path("LFQ/TrP")

study = fp.Study(lip=lip, trp=trp)

print(study.samples)

study.add_process(1  , "Native", "Refolded_001_min", 3, "Native", "Refolded", 3)
study.add_process(5  , "Native", "Refolded_005_min", 3, "Native", "Refolded", 3)
study.add_process(120, "Native", "Refolded_120_min", 3, "Native", "Refolded", 3)

study.run()

for result in study.results.values():
    name = result.name
 
    result.ion.write_excel(f"{name}_ion.xlsx", freeze_panes = (1,0))
    result.peptide.write_excel(f"{name}_peptide.xlsx", freeze_panes = (1,0))
    result.cut_site.write_excel(f"{name}_cut_site.xlsx", freeze_panes = (1,0))