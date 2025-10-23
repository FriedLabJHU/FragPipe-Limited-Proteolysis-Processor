![FLiPPR Logo](docs/images/flippr_long_logo.png)

## Table of Contents

- [Overview](#overview)
- [Support](#support)
- [Installation](#installation)
- [Usage](#usage)
- [License](#license)
- [Citation](#citation)

## Overview

FLiPPR (**F**ragPipe **Li**mited-**P**roteolysis **Pr**ocessor) is modular, fast, and easy-to-use data processing tool for LiP-MS experiments analyzed in FragPipe.

FLiPPR is:
* **Super Fast**: FLiPPR utilizes [Polars](https://pola.rs/) in the back-end to ensure all CPU cores contribute to your data processing.
* **Convinent to Use**: [FragPipe](https://fragpipe.nesvilab.org/) produces standarized outputs, FLiPPR takes full advantage of this feat and integrates seemlessly with any FragPipe LFQ DDA or DIA anaylsis.
* **Flexible and Expandable**: Experiments introduce unique and novel variables, FLiPPR ensures compatibility with all experimental setups and gives you full control of the data processing pipeline.

## Support
FLiPPR requires python ≥3.10
To learn more about how FLiPPR can help you analyze your LiP-MS or LFQ data, head over to the [FLiPPR docs](https://fragpipe-limited-proteolysis-processor.readthedocs.io)

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
import flippr

study = flippr.Study(lip = "path/to/fragpipe_output", method="dda")

```

2. View the experimental sample annotations

```python
print(study.samples)
# >>> { 'LiP': {'WT', '1x_Drug', '5x_Drug'}

```

3. Add an experimental process

```python
study.add_process(pid="1x", lip_ctrl="WT", lip_test="1x_Drug", n_rep=3)

study.add_process(pid="5x", lip_ctrl="WT", lip_test="5x_Drug", n_rep=3)

```

4. Run your study

```python
results = study.run()

results
# >>> {'1x': <flippr.Results>, '5x': <flippr.Results>}
```

5. View or save your results as Polars DataFrames

```python
results["1x"].ion

# View higher-order results
results["1x"].modified_peptide
results["1x"].peptide
results["1x"].cut_site

# View protein-level summary of all data
results["1x"].protein_summary

```

Learn about all the ways FLiPPR can analyze your data by heading over to the [FLiPPR docs](https://flippr.readthedocs.io)

## License

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
