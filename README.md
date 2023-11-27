# CORPSE-LIDET
CORPSE model with LIDET-derived litter decomposition parameters 

[![DOI](https://zenodo.org/badge/523086274.svg)](https://zenodo.org/badge/latestdoi/523086274)

**Description:**<br> This is a version of the CORPSE model (Carbon, Organisms, Rhizosphere and Protection in the Soil Environment, Sulman et al. 2014, Sulman et al. 2017) that uses litter decomposition parameters derived from a modified Monte Carlo simulation using the LIDET litter decomposition dataset (Long-term Intersite Decomposition Experiment Team, Harmon 2013). The code also includes the Baseline parameters, and the 8 other best parameter sets identified in the modified Monte Carlo simulation.

**Creators:**<br> Stephanie Juice, Joanna Ridgeway, Melannie Hartman, William Parton, Danielle Berardi, Benjamin Sulman, Kara Allen, and Edward Brzostek

**Contact information:**<br> stephanie.juice@mail.wvu.edu, erbrzostek@mix.wvu.edu 

**Software:** R

**Related publication:**<br>
Juice, S.M., Ridgeway, J.R., Hartman, M.D., Parton, W.J., Berardi, D.M., Sulman, B.N., Allen, K.E., & Brzostek, E.R. Reparameterizing litter decomposition using a simplified Monte Carlo method improves litter decay simulated by a microbial model and alters bioenergy soil carbon estimates.

**Accompanying files**:<br> The folder "Input Files" contains one folder for each LIDET site with data necessary to run the model. *Note that "(site)" in the filenames below indicates where the LIDET site code appears (see Table 1 for site codes).* Data streams include: 

1. **CORPSE_full_spinup_litter.csv, CORPSE_full_spinup_rhizo.csv, CORPSE_full_spinup_bulk.csv, litterbag_init_100g_6spp.csv:** initial C and N (kg C or N/m<sup>2</sup>) pool values for each soil layer, the litterbag_init_100_6spp.csv file is for the litterbag layer and is the same file for all sites.

	All initial C and N files have the same columns:

| **Column**    | **Description**                       | **Units**                  |
| -------------- | ------------------------------------- | -------------------------- |
| uFastC         | Unprotected fast decomposing carbon   | kg  carbon/m<sup>2</sup>   |
| uSlowC         | Unprotected slow decomposing carbon   | kg  carbon/m<sup>2</sup>   |
| uNecroC        | Unprotected necromass carbon          | kg  carbon/m<sup>2</sup>   |
| pFastC         | Protected fast decomposing carbon     | kg  carbon/m<sup>2</sup>   |
| pSlowC         | Protected slow decomposing carbon     | kg  carbon/m<sup>2</sup>   |
| pNecroC        | Protected necromass carbon            | kg  carbon/m<sup>2</sup>   |
| livingMicrobeC | Carbon in living microbial biomass    | kg  carbon/m<sup>2</sup>   |
| uFastN         | Unprotected fast decomposing nitrogen | kg  nitrogen/m<sup>2</sup> |
| uSlowN         | Unprotected slow decomposing nitrogen | kg  nitrogen/m<sup>2</sup> |
| uNecroN        | Unprotected necromass nitrogen        | kg  nitrogen/m<sup>2</sup> |
| pFastN         | Protected fast decomposing nitrogen   | kg nitrogen/m<sup>2</sup>  |
| pSlowN         | Protected slow decomposing nitrogen   | kg  nitrogen/m<sup>2</sup> |
| pNecroN        | Protected necromass nitrogen          | kg  nitrogen/m<sup>2</sup> |
| inorganicN     | Inorganic nitrogen                    | kg  nitrogen/m<sup>2</sup> |
| CO2            | Carbon in carbon dioxide              | kg  carbon/m<sup>2</sup>   |
| livingMicrobeN | Nitrogen in living microbial biomass  | kg  nitrogen/m<sup>2</sup> |

 

2. **soilT (site) DOY274start.csv**: Average daily soil temperature (<sup>o</sup>C) interpolated from previously calculated monthly values used in DayCent LIDET simulations (Bonan et al., 2013).

3. **soilT (site) DOY274start.csv**: Average daily soil volumetric water content (VWC) scalar interpolated from previously calculated monthly values used in DayCent LIDET simulations (Bonan et al., 2013).

4. **litter production.csv**: Average daily litter production values for each site, data sources listed in Table S3 of related publication. 
	
5. **litter (site) CN.csv**: C:N ratio for each species from LIDET dataset (Table 2, Harmon 2013). 

6. **(site).csv**: Table indicating number of observations for each species decomposed at each site. 

\
\
**Instructions:** 

1) Save the model code (“CORPSE_LIDET.R”) and "Input Files" folder in the same folder. Also make a folder for the model output (e.g., "results_Baseline") in the same folder. 

2) Set the working directory (setwd) in the model code to the folder with the files saved in step #1.
   
4) Select the parameter set to use for the litter and litterbag compartments, comment out all other parameter sets. 

5) Run code. Output will be saved in the folder made in step 1. Output destination can be changed as necessary in code section called "Running the model."

\
\
**Table 1** LIDET sites and site codes used in model files. 

| **Site Code**    | **Site** | 
| -------------- | -----------|
|AND| H.J. Andrews Experimental Forest	|
|BNZ| Bonanza Creek Experimental Forest	
|BSF| Blodgett Research Forest	
|CDR| Cedar Creek Natural History Area	
|CPR| Central Plains Experimental Range	
|HBR| Hubbard Brook Experimental Forest	
|HFR| Harvard Forest	
|JUN| Juneau	
|KBS| Kellogg Biological Station	
|KNZ| Konza Prairie Research Natural Area	
|NWT| Niwot Ridge/Green Lakes Valley	
|OLY| Olympic National Park	OLY	Conifer forest	
|SEV| Sevilleta National Wildlife Refuge	
|SMR| Santa Margarita Ecological Reserve	
|UFL| University of Florida	
|VCR| Virginia Coast Reserve

\
\
**Table 2** LIDET species and species codes used in model files (6 common species). 

| **Species**    | **Species Code** | 
| -------------- | -----------|
Sugar maple (Acer saccharum)	|	ACSA
Drypetes (Drypetes glauca)	|	DRGL
Red pine (Pinus resinosa)	|	PIRE
Chestnut oak (Quercus prinus)	|	QUPR
Western redcedar (Thuja plicata)	|	THPL
Wheat (Triticum aestivum)	|	TRAE

\
\
**References:**<br>
Bonan, G. B., Hartman, M. D., Parton, W. J., & Wieder, W. R. (2013). Evaluating litter decomposition in earth system models with long-term litterbag experiments: an example using the Community Land Model version 4 (CLM4). Global Change Biology, 19(3), 957-974. https://doi.org/https://doi.org/10.1111/gcb.12031 

Harmon, M. (2013). LTER Intersite Fine Litter Decomposition Experiment (LIDET), 1990 to 2002. Long-Term Ecological Research. Forest Science Data Bank, Corvallis, OR. [Data set]. Accessed http://andlter.forestry.oregonstate.edu/data/abstract.aspx?dbcode=TD023. https://doi.org/10.6073/pasta/f35f56bea52d78b6a1ecf1952b4889c5.

Sulman, B. N., Phillips, R. P., Oishi, A. C., Shevliakova, E., & Pacala, S. W. (2014). Microbe-driven turnover offsets mineral-mediated storage of soil carbon under elevated CO2. Nature Climate Change, 4, 1099 - 1102. https://doi.org/10.1038/nclimate2436

Sulman, B. N., Brzostek, E. R., Medici, C., Shevliakova, E., Menge, D. N., & Phillips, R. P. (2017). Feedbacks between plant N demand and rhizosphere priming depend on type of mycorrhizal association. Ecology letters, 20(8), 1043-1053.
