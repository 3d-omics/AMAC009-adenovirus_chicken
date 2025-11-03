# Adenovirus Microbiome Analysis - R Bookdown Project

This repository contains the complete R analysis code for the adenovirus trial microbiome study, compiled as an R Bookdown document for easy sharing and reproducibility.

**📖 View the compiled bookdown**: [https://3d-omics.github.io/AMAC009-adenovirus_chicken/](https://3d-omics.github.io/AMAC009-adenovirus_chicken/)

## Project Overview

This project analyzes adenovirus presence and its effects on the chicken microbiome using metagenomic sequencing data. The analysis compares two chicken breeds (layers and broilers) under different experimental treatments across multiple time points following infection.

## Repository Structure

This repository contains the complete R analysis code for the adenovirus trial microbiome study. The analysis is organized into chapters covering:

- **Data Preparation** (Chapters 01A, 01B, 02): Data import, processing, and preparation for analysis
- **MAG Catalogue** (Chapter 03): Overview of the metagenome-assembled genome catalogue
- **Adenovirus Analysis** (Chapters 04A-04G): Comprehensive analyses including:
  - Adenovirus presence over time
  - Sequencing statistics
  - Antimicrobial effects
  - Alpha and beta diversity analyses
  - Diversity partitioning
- **Sample Selection** (Chapters 05A-05C): Identification of samples for downstream sequencing (PacBio, LMD, Metatranscriptomics)

The compiled HTML bookdown document provides a complete view of all analyses with code, plots, and interpretations.

## Data Files

The repository includes processed data files (`.Rdata`) used in the analyses:

- **`data/macro/plot_data.Rdata`**: Processed plotting data including sample metadata, sequencing statistics, and tidy-format genome count tables
- **`data/MAG_catalogue/data.Rdata`**: MAG catalogue information including phylogenetic trees, taxonomic data, and functional annotations
- **`data/macro/data.Rdata`**: Core sample metadata and sequencing statistics
- **`data/data_colors.Rdata`**: Color schemes and plotting themes used throughout the analyses
- **`data/macro/extract_concentration_data.Rdata`**: DNA extraction concentration data for sample selection
- **`data/MAG_catalogue/prep_*.Rdata`**: Complete preparation data files for additional analyses (taxonomy, quality metrics, annotations)


## Contact

For questions about the analysis or data access, please contact:
- Amalia Bogri: amalia.bogri@sund.ku.dk
- Antton Alberdi: antton.alberdi@sund.ku.dk
- Jorge Langa: jorge.langa@sund.ku.dk

## Citation

If you use this code in your research, please cite the repository:

Bogri, A., Alberdi, A., & Langa, J. (2024). Adenovirus Microbiome Analysis. GitHub repository. https://github.com/3d-omics/AMAC009-adenovirus_chicken

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
