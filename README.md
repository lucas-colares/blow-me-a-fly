![alt text](https://i.imgur.com/AGtNp4d.png)
# BlowMeAFly repository 🪰🌳
## Overview

This repository contains the data and R code used in the study **"Changes in blowfly (Diptera: Calliphoridae) wing morphology during succession in rat carcasses across forest and grassland habitats in South Brazil"** by Colares et al. (2024). The study investigates the relationship between dispersal capacity and the timing of species arrival during ecological succession, focusing on carrion blowflies.

# Data outline
In this repository, you will find:
- Landmark data that represent the shape of (A) wing (wingFixed2.TPS), (B) dorsal view of the thorax (toraxdorsalFixed2.TPS), and (C) lateral view of the thorax (toraxlateralFixed2.TPS) of 355 individuals of carrion Blowflies, distributed across 11 species.
![alt text](https://i.imgur.com/f6fsR8U.png)
- Abundance data from 44 samples (rat carcasses) across space and time (in 4 spatial samples collected throughout the whole decomposition process of rats; Blow_Occ.csv). Also, you will find data of nine functional traits: (i) wing shape; (ii) thorax shape in dorsal view; (iii) thorax shape in lateral view; (iv) wing size; (v) thorax size in dorsal view; (vi) thorax area in lateral view; (vii) thorax volume; (viii) thorax volume/wing size ratio; and (ix) thorax height. Finally, you will find the environmental data corresponding to each sample as measures of temperature, precipitation and seasons (specimenID.csv).
![alt text](https://i.imgur.com/nY4unsp.jpg)
- An R code (BlowflyAnalysis_Jan2023.R) to reproduce our morphological and statistical analysis of traits and the graphs presented in the main paper and supporting information.

## Repository Structure

- **Data Files**:
  - `wingFixed2.TPS`: Landmark data representing the wing shapes of 355 individual blowflies across 11 species.
  - `toraxdorsalFixed2.TPS`: Landmark data representing the dorsal view of the thorax.
  - `toraxlateralFixed2.TPS`: Landmark data representing the lateral view of the thorax.
  - `Blow_Occ.csv`: Abundance data from 44 samples (rat carcasses) collected throughout the decomposition process.
  - `specimenID.csv`: Metadata linking specimens to their respective data points.

- **R Script**:
  - `BlowflyAnalysis_Jan2023.R`: Contains the R code to reproduce the analyses and generate the graphs presented in the study.

## Usage

To reproduce the analyses:

1. **Prerequisites**: Ensure you have R installed on your system. Required R packages are specified within the script.

2. **Running the Analysis**:
   - Clone this repository:
   ```bash  
   git clone https://github.com/lucas-colares/blow-me-a-fly.git  
   ``` 
   - Open the downloaded folder, then open the `BlowflyAnalysis_Jan2023.R` script in R or RStudio.
   - Execute the script to perform the analyses and generate the graphs.

## License

This project is licensed under the **MIT License**. You are free to use, modify, and distribute the code with proper attribution. 

## Citation

If you use this repository, please cite:

> Colares, L., Herdina, A., Bender M. & Dambros, C. (2024). *Changes in blowfly (Diptera: Calliphoridae) wing morphology during succession in rat carcasses across forest and grassland habitats in South Brazil*. *Insect Science*. [https://doi.org/10.1111/1744-7917.13485](https://doi.org/10.1111/1744-7917.13485)
