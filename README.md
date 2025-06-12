[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15642732.svg)](https://doi.org/10.5281/zenodo.15642732)

# Spatial transcriptomics of Adenomyosis lesions using Nanostring GeoMX

This repo contains the software and data required to repeat the analysis in the paper titled: 'Spatially-resolved transcriptomics uncovers the hybrid molecular identify, ciliated phenotype, and immune signature of adenomyosis lesions'. Analysis is contained in reproducible quarto markdown (.qmd) documents. The documents are numbered so they can be ran in the correct order. Details are below. 

![image](https://github.com/CBFLivUni/AdenomyosisGeoMX/blob/main/supporting_images/experimental_design.png)



## Running the analysis
### Using RScript

1.	Install the needed R packages
    ```console
     RScript install/install.R
    ```
2.	Run the analysis and render the html notebooks
    ```console
     RScript run_analysis.R
    ```
