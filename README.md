# DESeq2 Example in Python

This repository contains an example of using DESeq2 in Python for differential gene expression analysis.

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Files](#files)
- [References](#references)
- [Contributing](#contributing)
- [Contact](#contact)
- [Acknowledgements](#acknowledgements)
- [License](#license)

## Introduction

[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) is a widely used tool for differential gene expression (DGE) analysis in RNA-Seq data, and it is mainly implemented in [R](https://cran.r-project.org/).
However, you can use [rpy2](https://pypi.org/project/rpy2/), a Python library that allows you to run R code. If you're interested in learning more about RNA-seq analysis, I highly recommend visiting [RNA-seqlopedia](https://rnaseq.uoregon.edu/). This comprehensive resource provides valuable insights into RNA-seq, covering everything from basic concepts to advanced techniques. It's an excellent starting point for anyone looking to deepen their understanding of RNA-seq analysis.
## Installation 

To run the DESeq2 example in Python, Ensure you have installed `R` on your system, and install the required packages. Follow these steps:


1. **Install required packages:**
    ```bash
    pip install -r requirements.txt
    ```

2. **Install DESeq2:**
    ```r
    # Open R and run the following command to install DESeq2
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("DESeq2")
    ```


## Usage

To run the DESeq2 analysis, follow these steps:

1. **Prepare your count data:**
    Ensure your count data is in a CSV file format. An example file is provided in the repository.

2. **Run the DESeq2 script:**
    ```bash
    python run_deseq2.py 
    ```

3. **Analyze the results:**
    The results will be saved in your provided path or directory. You can use the provided `analyze_results.py` script to generate plots and interpret results.
    ```bash
    python analyze_results.py 
    ```


## Files

- `Data/filtered_gene_expression_counts.csv`: Example count data for DESeq2 analysis.
- `Data/deseq2_results.csv`: Example of the result.
- `Code/run_deseq2.py`: Script to run DESeq2 analysis.
- `Code/analyze_results.py`: Script to analyze DESeq2 results and generate plots.
- `requirements.txt`: List of Python dependencies.
- `README.md`: This file.

## References

- [DESeq2 Documentation](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- [Bioconductor](https://bioconductor.org/)
- [rpy2](https://pypi.org/project/rpy2/)

## Contributing
If you would like to contribute to this project, please fork the repository and submit a pull request. For major changes, please open an issue to discuss what you would like to change.

## Contact
Thanks for visiting this repository! If you have any questions or feedback, feel free to contact me at [ptlsugam@gmail.com].

## Acknowledgements
Thank you for visiting this repository. Your feedback and contributions are greatly appreciated. Feel free to reach out if you have any suggestions for improvements or new features.

## License

This repository is licensed under the MIT License. See the [LICENSE](LICENSE.md) file for details.
