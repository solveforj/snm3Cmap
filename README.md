
## Overview

This pipeline contains tools for mapping and processing Hi-C/3C data to generate chromatin contacts. It is compatible with both non-bisulfite-converted and bisulfite-converted DNA. Additionally, it has a specialized module for mapping and processing snm3C-seq data, which has the additional functionality of calling genome-wide methylation in addition to chromatin contacts.

## Installation

```
conda env create -n snm3Cmap --file snm3Cmap.yml
conda activate snm3Cmap
pip install git+https://github.com/solveforj/snm3Cmap.git
```

## Usage

Please use the snm3Cmap conda environment.
```
conda activate snm3Cmap
```

Task-specific documentation can be found at the following locations:

* [snm3C-seq](docs/snm3Cseq/README.md) - For demultiplexing and mapping snm3C-seq data
