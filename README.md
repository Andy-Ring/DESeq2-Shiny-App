# DESeq2 Analysis and Visualization Shiny App

## Overview
This Shiny app provides an interactive interface for performing differential gene expression analysis using DESeq2. Users can upload count matrices, select experimental conditions, perform pairwise comparisons, and visualize results through PCA plots, heatmaps, and volcano plots.

## Features
- **Data Input:** Upload count data in CSV format or Kallisto abundance files in TSV format.
- **Metadata Handling:** Import metadata files and select design variables for DESeq2 analysis.
- **Exploratory Analysis:** Generate PCA plots and heatmaps to assess sample clustering.
- **Pairwise Comparisons:** Perform differential expression analysis and visualize results with volcano plots.
- **Gene-Level Analysis:** Inspect individual gene expression patterns across conditions.
- **Downloadable Results:** Save high-resolution plots and significant gene lists for further analysis.

## Web Access
This app can be run on the web at https://0195156d-49f7-f338-c8bf-7b1165e8e62d.share.connect.posit.cloud/

## Running App Locally
### Installation
To run this Shiny app, ensure you have R and the required packages installed. You can install the necessary dependencies using:
```r
install.packages(c("shiny", "shinydashboard", "bslib", "DESeq2", "tximport", "ggplot2", "pheatmap", "ggrepel", "matrixStats", "EnhancedVolcano", "org.Hs.eg.db", "dplyr", "rlang"))
```

### Running the App
1. Clone this repository or download the script files.
2. Open R or RStudio and navigate to the app directory.
3. Run the following command:
   ```r
   shiny::runApp()
   ```

## User Guide
### 1. Data Input
- Choose between **Kallisto** (tximport) or **CSV Count Matrix**.
- Upload the corresponding files.
- Download a metadata template if needed.
- Upload the metadata file
- Select the design variables.
- Click **Run DESeq2 Analysis** to process the data.

### 2. Exploratory Analysis
- **PCA Plot:** Visualize sample clustering based on selected variables.
- **Heatmap:** Display the top 50 most variable genes with optional clustering.

### 3. Pairwise Comparison
- Select a contrast variable and groups for differential analysis.
- Adjust p-value and log2 fold-change cutoffs.
- Generate a **volcano plot** to visualize significant genes.
- Download significant gene lists and plots.

### 4. Gene-Level Analysis
- Select a specific gene to view expression trends.
- Generate a bar plot showing normalized counts across conditions.

## Downloads
- **PCA Plot** (PNG)
- **Heatmap** (PNG)
- **Volcano Plot** (PNG)
- **Gene Expression Plot** (PNG)
- **Significant Genes List** (CSV)

## Credits
- **Author:** Andy Ring
- **Version:** 1.0.1
- **Date:** February 17th, 2025

## License
This project is licensed under the MIT License. Feel free to use and modify it for your research and analysis needs.

