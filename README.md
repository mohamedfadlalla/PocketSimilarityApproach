# Pocket Similarity Approach

## Introduction
The Pocket Similarity Approach is a set of Python scripts for automating pocket matching in protein structures. This repository supports the reproduction of results from a manuscript that explores drug repurposing through pocket similarity approaches.

## Technologies
- Python 3.7.1
- Biopython
- Pandas
- Numpy

## Installation
To set up this project locally:
1. Clone the repository:
   ```bash
   git clone https://github.com/mohamedfadlalla/PocketSimilarityApproach.git
   ```
2. Install the required dependencies:
   ```bash
   pip install biopython pandas numpy
   ```

## Usage
The repository includes scripts for preparing drug and COVID-19 protein pockets for analysis:
- `Drug_pocket_prepration.py`: Extracts the pocket of a drug.
- `COVID19_pocket_prepration.py`: Handles the extraction of the COVID-19 protein pocket.

Run the scripts as follows:
```bash
python Drug_pocket_prepration.py
python COVID19_pocket_prepration.py
```

Example output files for Fpocket results on COVID19 and drug proteins structures are provided in the respective `.rar` files.

## Features
- Automated extraction and processing of protein pockets.
- Support for analyzing approved drug data (`HET_of_approved_drugs.txt`).

## Contributing
Contributions to the Pocket Similarity Approach are welcome. Please fork the repository and submit a pull request with your proposed changes.

## License
This project is licensed under the terms of the [MIT License](LICENSE.md).

## Authors
- Mohamed Fadlalla - *Initial work* - [mohamedfadlalla](https://github.com/mohamedfadlalla)

## Acknowledgments
- The manuscript associated with this work can be found [here](https://chemrxiv.org/articles/preprint/COVID19_Approved_Drug_Repurposing_Pocket_Similarity_Approach/12722483).
