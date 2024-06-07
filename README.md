![FLiPPR Logo](/assets/images/flippr_long_logo.png)

## Table of Contents

- [Overview](#overview)
- [Support](#support)
- [Installation](#installation)
- [Usage](#usage)
- [License](#license)
- [Citation](#citation)

## Overview

FLiPPR (**F**ragPipe **Li**mited-**P**roteolysis **Pr**ocessor) is modular, fast, and easy-to-use data processing tool for LiP-MS & LFQ-based experiments analyzed in FragPipe.

FLiPPR is:
* **Super Fast**: FLiPPR utilizes [Polars](https://pola.rs/) in the back-end to ensure all CPU cores contribute to your data processing.
* **Convinent to Use**: [FragPipe](https://fragpipe.nesvilab.org/) produces standarized outputs, FLiPPR takes full advantage of this feat and integrates seemlessly with any FragPipe LFQ-MBR anaylsis.
* **Flexible and Expandable**: Experiments introduce unique and novel variables, FLiPPR ensures capatibility with all experimental setups and gives you full control of the data processing pipeline.

## Support
FLiPPR has been tested with Python 3.10 – 3.12
To learn more about how FLiPPR can help you analyze your LiP-MS or LFQ data, head over to the [FLiPPR docs](https://flippr.readthedocs.io)

## Installation

```bash
    # Initial install
    python -m pip install flippr

    # Update to latest release
    python -m pip install -U flippr
```


## Usage

1. Start a `Study`

```python
import flippr as fp
from flippr import Study

# Pass the FragPipe output directory path from your LiP-MS experiment
study = Study(lip = "path/to/lip")

# Include protein normalization factors from a Trypsin-only experiment
study = Study(lip = "path/to/lip", trp = "path/to/trp")

# Including the protein FASTA file will add metadata and remove contaminants
study = Study(lip = "path/to/lip", trp = "path/to/trp", fasta = "UP000000625_83333.fasta")
```

2. View the experimental sample annotations

```python
print(study.samples)
# >>> { 'LiP': {'Native', 'Refolded_001_min', 'Refolded_005_min', 'Refolded_120_min'}, 'TrP': {'Refolded', 'Native'}}
```

3. Add an experimental process

```python
# Process without normalization
study.add_process(1, "Native", "Refolded_001_min", 3)

# Process with normalizations (only when `trp` is included in the study)
study.add_process(5, "Native", "Refolded_005_min", 3, "Native", "Refolded", 3)

# Multiple processes within one study
study.add_process(1, "Native", "Refolded_001_min", 3, "Native", "Refolded", 3)
study.add_process(5, "Native", "Refolded_005_min", 3, "Native", "Refolded", 3)
study.add_process(120, "Native", "Refolded_120_min", 3, "Native", "Refolded", 3)
```

4. Run your study

```python
# Running a study returns results and populates the `Study().results` dictionary
results = study.run()
# >>> {1: <flippr.Results>, 5: <flippr.Results>, 120: <flippr.Results>}

# Change analysis parameters within `fp.rcParams` before running for more control
fp.rcParams["protein.fc_sig_thresh"] = 2.0
results = study.run()
```

5. View or save your results as Polars DataFrames

```python
# 1 min time point results
results_1_min = results[1]

# View all ions included in the 1 min experimental process
results_1_min.ion

# View higher-order results
results_1_min.modified_peptide
results_1_min.peptide
results_1_min.cut_site

# View protein-level summary of all data
results_1_min.protein_summary
```

| Protein ID | Protein                     | Entry Name  | Gene   | ... | No. of Valid Cut Sites | No. of Significant Cut Sites | No. of Significant Cut Sites |
|------------|-----------------------------|-------------|--------|-----|------------------------|------------------------------|------------------------------|
| P00350     | sp\|P00350\|6PGD_ECOLI      | 6PGD_ECOLI  | gnd    | ... | 111                    | 40                           | 50                           |
| P00363     | sp\|P00363\|FRDA_ECOLI      | FRDA_ECOLI  | frdA   | ... | 15                     | 2                            | 2                            |
| Q57261     | sp\|Q57261\|TRUD_ECOLI      | TRUD_ECOLI  | truD   | ... | 24                     | 7                            | 8                            |
| Q59385     | sp\|Q59385\|COPA_ECOLI      | COPA_ECOLI  | copA   | ... |13                      | 2                            | 1                            |
| Q7DFV3     | sp\|Q7DFV3\|YMGG_ECOLI      | YMGG_ECOLI  | ymgG   | ... | 4                      | 2                            | 4                            |
| Q93K97     | sp\|Q93K97\|ADPP_ECOLI      | ADPP_ECOLI  | nudF   | ... | 9                      | 4                            | 5                            |

## License

This project is licensed under the CC BY-NC-ND 4.0 License. Feel free to use the code according to the terms specified in the license.
This project is licensed under the CC BY-NC-ND 4.0 License. Feel free to use the code according to the terms specified in the license.

Thank you for your interest in FLiPPR!
If you encounter any issues or have suggestions, please open an issue.
We appreciate your feedback!

## Citation
If you found this work helpful in your research, please cite us!  
 
**FLiPPR: A Processor for Limited Proteolysis (LiP) Mass Spectrometry Data Sets Built on FragPipe**  
Edgar Manriquez-Sandoval, Joy Brewer, Gabriela Lule, Samanta Lopez, and Stephen D. Fried  
*Journal of Proteome Research*  
DOI: ***10.1021/acs.jproteome.3c00887***
