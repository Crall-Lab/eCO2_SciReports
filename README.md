# eCO2_SciReports
Code and data files used for analyses and figure making for Bernauer et al. in prep. - <i>"Elevated atmospheric CO<sub>2</sub> has inconsistent effects on pollen chemistry and plant growth across flowering plant species"</i>.  
Data and code for this manuscript can also be found on Dryad: Effects of eCO<sub>2</sub> on plant growth and pollen chemistry in 14 angiosperms at <https://doi.org/10.5061/dryad.70rxwdc5n>

# Datasets
Here is a list of the included data sets included here and used for analyses in the manuscript. 

1. <b>Metadata.xlsx</b> is the the metadata for all subsequent data files. Each file name is listed in bold with all column headings in column A and column descriptions, including units in column B. 
2. <b>CO2_log_exp1.csv</b> is the data file for the CO<sub>2</sub> logger data from experiment 1 containing time-stamped CO<sub>2</sub> and temperature readings for ambient and elevated CO<sub>2</sub> conditions. 
3. <b>CO2_log_exp2.csv</b> is the data file for the CO<sub>2</sub> logger data from experiment 2 containing time-stamped CO<sub>2</sub>, humidity, and temperature readings for ambient and elevated CO<sub>2</sub> conditions in rounds 1 and 2. Outliers (high and low) are noted in a column so they can be removed for analyses.  
4. <b>Pollen_metabolomics_exp1.csv</b> is the data file containing the targeted metabolomics data from experiment 1. For each pollen sample, this data file records, the unique sample identifier, plant species, CO<sub>2</sub> treatment, sample type (pollen vs. anther), unique file name, the quanitities of 18 targeted compounds in each sample, and the experimental round each sample was from. 
5. <b>Pollen_chemistry_exp1.csv</b> contains pollen chemistry (%N, %C, C:N) data for experiment 1. For each pollen sample, this data file records the unique sample identifier, the plant species the sample came from, the CO<sub>2</sub> treatment, experimental round, %N and %C present in each sample, notes about the sample and a omit column to remove samples that were too low to be detected or duplicate samples. 
6. <b>Pollen_chemistry_exp2.csv</b> contains pollen chemistry (%N, %C, C:N) data for experiment 2. For each pollen sample, this data file records the unique sample identifier, the plant species the sample came from, the CO<sub>2</sub> treatment, experimental round, %N and %C present in each sample, notes about the sample and a omit column to remove samples that were too low to be detected or duplicate samples. 
7. [<b>Biomass_exp2.csv</b>](https://github.com/Crall-Lab/eCO2_SciReports/blob/979f5964d17bb705b5a8e535265b44bec3a39363/Biomass_exp2.csv) is the data file recording biomass data for experiment 2. Here, above-ground, dry biomass for each plant sample is reported along with plant species, unique plant identifier, greenhouse chamber, and CO<sub>2</sub> treatment. Biomass data was only collected for round 1 of experiment 2.  
8. <b>Flower_size_exp2.csv</b> contains the flower size (measured as diameter) for experiment 2. Along with flower size (as measured by flower diameter as a mean of three diameter measuremenets on the same flower), this data set contains the associated greenhouse chamber, CO<sub>2</sub> treatment, experimental round, unique plant identifier, and plant species. 
9. <b>Flowering_initiation_exp2.csv</b> is the data set containing the date each plant that flowered in experiment 2, initiated flowering. In addition to flowering initiation (flowering start date) this data set reportsthe associated green house chamber, CO<sub>2</sub> treatment, experimental round, plant species, and unique plant identifier. 
10. <b>Plant_growth_exp2.csv</b> contains the plant growth data for experiment 2. This data set reports the date the data was recorded, the associated greenhouse chamber, CO<sub>2</sub> treatment, and experimental round along with the plant species and unique plant identifier and the number of leaves (Leaf_no), number of flowers (Flower_no), and plant's height (in cm). There is also a notes column. 

# Code
## Analysis and visualization
Here is a list of the R (v 4.2.2) code files used to analyse the above data and generate the figures used in the manuscript. All necessary packages are listed at the top of each code file. 
1. <b>CO2_plotting.R</b> contains the script file used to create Figure S1 which depicts the CO<sub>2</sub> levels across experimental rounds for both experiments 1 and 2. This code file uses data files: <b>CO2_log_exp1.csv</b> and <b>CO2_log_exp2.csv</b>.
2. <b>Pollen_chemistry.R</b> contains the script file used to create Figures 1 and 3 and to analyze the pollen chemistry data (%N, %C, and C:N) for both experiments 1 and 2. This code file uses data files: <b>Pollen_chemistry_exp1.csv</b> and <b>Pollen_chemistry_exp2.csv</b>.
3. <b>Targeted_metabolomics_exp1.R</b> contains the script file used to create Figure 2 and to analyze the targeted metabolomics data from experiment 2 using the <b>Pollen_metabolomics_exp1.csv</b> data file. 
4. <b>Plant_growth.R</b> contains the script file used to create Figures 4-8 and to analyze the plant growth data (biomass, flower size, flowering initiation, height, flower number and leaf number) for experiment 2. This code file uses the data files: [<b>Biomass_exp2.csv</b>](https://github.com/Crall-Lab/eCO2_SciReports/blob/979f5964d17bb705b5a8e535265b44bec3a39363/Biomass_exp2.csv), <b>Flower_size_exp2.csv</b>, <b>Flowering_initiation_exp2.csv</b>, and <b>Plant_growth_exp2.csv</b>.

## Greenhouse CO<sub>2</sub> control
To control the CO2 loggers in the greenhouses 

# Access information
This data and code can be found on github (<https://github.com/Crall-Lab/eCO2_SciReports.git>) and dryad at <https://doi.org/10.5061/dryad.70rxwdc5n>.
