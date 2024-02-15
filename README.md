# FLiPPR v0.0.5

## Overview

> [!IMPORTANT]  
> FLiPPR currently does not incorporate metadata. This feature is coming in the next release.

Welcome to FLiPPR! Please note that this version is an early release, and the code base is incomplete. We are actively working on enhancing and expanding the features.

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Installation

To get started with FLiPPR, follow these steps:

1. Clone the repository to your local machine:

    ```
    git clone https://https://github.com/FriedLabJHU/FragPipe-Limited-Proteolysis-Processor.git
    ```

2. Navigate to the project directory:

    ```
    cd FragPipe-Limited-Proteolysis-Processor
    ```

3. Install FLiPPR:

    ```
    pip install .
    ```

## Usage

1. Start a Study

```python
from flippr import Study

# Pass the FragPipe output directory path from your LiP-MS study
study = Study(lip = "path/to/lip")

# Include protein normalization factors from a Trypsin-only study
study = Study(lip = "path/to/lip", trp = "path/to/trp")
```

2. View the experimental sample annotations

```python
print(study.samples)
# > {
#    'LiP': 
#       {'Refolded_005_min', 'Native', 'Refolded_120_min', 'Refolded_001_min'},
#    'TrP': 
#       {'Refolded', 'Native'}
#   }
```

3. Add a Process

```python
# Process without normalization
study.add_process(1  , "Native", "Refolded_001_min")

# Process with normalizations (only when `trp` is included in the study)
study.add_process(5  , "Native", "Refolded_005_min", 3, "Native", "Refolded", 3)
study.add_process(120, "Native", "Refolded_120_min", 3, "Native", "Refolded", 3)

# Process with different normalizations (only when `trp` is included in the study & contains all normalization conditions)
study.add_process(5  , "Native", "Refolded_005_min", 3, "Native", "Refolded_005_min", 3)
study.add_process(120, "Native", "Refolded_120_min", 3, "Native", "Refolded_120_min", 3)
```

4. Run your study

```python
results = study.run()
# > {1: Results<Refolded_001_min_v_Native>,
# >  5: Results<Refolded_005_min_v_Native>,
# >  120: Results<Refolded_120_min_v_Native>}


# Change the numer of missing values tolerated by the study (only recommended for studies with +4 replicates in all conditions)
results = study.run(max_missing_values = 2)

# Change the All-or-Nothing Gaussian imputation parameters
# Defaults: aon_mean = 1e4, aon_std = 1e3
results = study.run(aon_mean = 1e6, aon_std = 1e2)

```

5. View your results

```python
# Viewing polars dataframes
# 1 min time point ions
results[1].ion

# 5 min time point cut-sites
results[5].cut_site

# 120 min time point summary with metadata
results[120].protein_summary
```

## Contributing

:safety_vest: :hammer_and_wrench: Under construction :hammer_and_wrench: :safety_vest:	

## License

This project is licensed under the MIT License. Feel free to use, modify, and distribute the code according to the terms specified in the license.

Thank you for your interest in FLiPPR! If you encounter any issues or have suggestions, please open an issue. We appreciate your feedback!

## Citation

:safety_vest: :hammer_and_wrench: Under construction :hammer_and_wrench: :safety_vest:	
