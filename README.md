# FLiPPR v0.0.3

## Overview

> [!IMPORTANT]  
> FLiPPR currently does not incorporate metadata. This feature is coming in the next release.

Welcome to the v0.0.2 release of FLiPPR! Please note that this version is an early release, and the code base is incomplete. We are actively working on enhancing and expanding the features.

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

> [!NOTE]  
> Documentation is coming, I promise!

1. Start a Study

```python
import flippr as fp

# Pass the FragPipe output directory path from your LiP-MS study
study = fp.Study(lip = "path/to/lip")

# Include protein normalization factors from a Trypsin-only study
study = fp.Study(lip = "path/to/lip", trp = "path/to/trp")
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
```

4. Run your study & view your results

```python
results = study.run()
# > {1: Results<Refolded_001_min_v_Native>,
# >  5: Results<Refolded_005_min_v_Native>,
# >  120: Results<Refolded_120_min_v_Native>}

# view polars dataframes
results[1].ion

results[5].cut_site

results[120].protein_summary
```

5. Save to Excel

```python
results[120].protein_summary.write_excel(f"{results[120].name}_flippr_summary.xlsx")
```

## Contributing

:safety_vest: :hammer_and_wrench: Under construction :hammer_and_wrench: :safety_vest:	

## License

This project is licensed under the MIT License. Feel free to use, modify, and distribute the code according to the terms specified in the license.

Thank you for your interest in FLiPPR! If you encounter any issues or have suggestions, please open an issue. We appreciate your feedback!

## Citation

:safety_vest: :hammer_and_wrench: Under construction :hammer_and_wrench: :safety_vest:	
