# Data from: Habitat suitability ensembles of genotype and disease resistance of the endangered Juglans cinerea tree to assist restoration efforts

These methods and data were developed to build regression and machine-learning based predictive models that were integerated into ensembles. The goal of the project was to understand what factors influence the distrubtion of sub-populations of a native tree (butternut). The butternut tree is currently imperilled by a non-native disease (Butternut canker disease) and hybridizes with a non-native congener (Japanese heartnut). In addition to identifying key drivers of population distrubtion, we sought to highlight areas where conservation efforts should be focus; both in an effor to find additional families within each sub-population, but also to identify areas climatically amenable to restoration effots. 


##Description of the data and file structure

**Data acquisition and pre-processing**
	Family BCD resistance scores and species/hybrid categorizations were obtained as part of ongoing research at the Hardwood Tree Improvement and Regeneration Center (Conrad et al., 2025). Accuracy and precision of family location varied across the dataset. Therefore, we elected to aggregate these data, and all data thereafter, at the county level. Here, “butternut” refers to genotypes that are unambiguously assigned to the species Juglans cinerea. Meanwhile “hybrid” refers to genotypes with evidence of Japanese walnut introgression. Families were already binned to “butternut” or “hybrid” based on single nucleotide polymorphism-based genotyping and morphology. We binned butternut into Oc-j “more resistant” or “less resistant” based on the quantitative BCD resistance phenotype provided by Conrad et al., 2025. For the purposes of this study, “more resistant” families had a positive breeding value, while “less resistant” families were negative. Importantly, “more resistant” trees could still develop symptoms of BCD, although “more resistant” families generally contained trees with less severe disease symptoms. Counties were then split into “butternut”, “hybrid”, “more resistant” and “less resistant” datasets to build an ensemble of habitat suitability models for each category (Figure 1). For each, counties were subjected to a random 70-30 training-test split (Figure 1A). For the training data, additional counties were assigned as pseudo-absences at a ratio of three absences for every one presence. For the test data, the ratio was 1:1 as has been suggested (Allouche et al. 2006). 
We sourced the climate data from the Parameter-elevation Regressions on Independent Slopes Model (PRISM) (Daly et al. 1994; PRISM 2024). We used the thirty year (1980-2020) normal for each variable at a four-kilometer resolution. We included both monthly and annual averages. We calculated the climate variable average for every county or county-equivalent in all states east of Minnesota, Iowa, Missouri, Arkansas, and Louisiana, inclusive, using custom R (v4.5.0 (2025-04-11) commands (Team 2025)), using R Studio (v2025.5.0.496, (Team 2025)) with the tigris (Walker 2015), raster (Hijmans 2010), sf (Pebesma 2018; Pebesma and Bivand 2023), and dplyr (Wickham et al. 2023) packages (Figure 1B).
Organic soil carbon and butternut basal area estimates were derived from the US Department of Agriculture Forest Service Forest Inventory and Analysis program data base (FIADB v2.0.6, (Bechtold and Patterson 2005)) using custom SQL queries. We queried the FIADB in R using the RSQLite package (Müller et al. 2024). SQL queries that required lengthy computational time utilized the beepr package (Bååth 2024) to inform the user when the query had finished. Custom SQL queries were modified from those provided in the SQL_QUERY_SE column of the REF_POP_ATTRIBUTE table within the FIADB. The queries for organic soil carbon and basal area were based off SQL_QUERY_SE numbers 52 (“Carbon in organic soil, in short tons, on forested land”), and 1004 (“Basal area of live trees (at least 1 inch d.b.h/d.r.c), in square feet, on forested land”), respectively. For soil organic carbon, we restricted our query to only include plots designated by the US Forest Service as occurring on accessible forested land (n = 77,360 plots). For butternut basal area, we restricted our query to include plots from accessible forested land that had live butternut stems (n = 220 plots). Estimates were validated with the US Forest Service’s EVALIDator tool (https://apps.fs.usda.gov/fiadb-api/evalidator) before converting to metric units. Both estimates used the data from the 2020 evaluation cycle. In 2023, Connecticut overhauled their county system in favor of county-equivalent planning regions with boundaries that differed from previous counties. Averages made from the climate data used the planning region boundaries. Therefore, we assigned the state-wide average of county organic soil carbon value to each planning region to integrate both datasets (Fig S1). Additionally, due to FIA plot sampling intensity and variations in forest cover, several counties lacked an estimate for organic soil carbon. In these instances, we assigned a county an organic soil carbon value calculated from the average of adjacent counties (Fig S1). 
To facilitate the production of parsimonious models, we manually removed variables that had a Pearson’s correlation coefficient -0.7 < r > 0.7 within the training data. When possible, we preferentially retained annual over monthly means (Figure 1C). Consistently, temperature variables (Tmin, Tmean, and Tmax) had strong correlation values and were removed. All predictors used in the final composite models are include in Supplemental Dataset 1.


**Citizen-science resources for ensemble filters and ensemble validation**
	We utilized two citizen-science databases to filter the merged ensemble models by counties where there have been confirmed sightings of butternut. TreeSnap (https://treesnap.org/) was built specifically to engage the public in efforts to find and record imperiled tree species such as ash (Fraxinus spp.), American elm (Ulmus americana), and American chestnut to assist conservation and breeding efforts (Crocker et al. 2020). TreeSnap data were accessed on June 11th, 2025, and include 289 observations of butternut from October 11th, 2020 through June 1st, 2025. Observations were binned into counties (n = 98) which were subsequently used as a filter for the ensembles. The Global Biodiversity Information Facility (GBIF, https://www.gbif.org) is a digital world-wide repository of information with the goal of enhancing interoperability between sources of biodiversity data (Edwards 2001). GBIF collates data from herbaria and fossil collections through genome sequences and citizen science platforms. GBIF was accessed on June 13th, 2025 and the downloaded data included all mentions of Juglans cinerea in North America. The GBIF data were then curated to include only live butternut trees in the United States. The resulting dataset includes 2,587 observations from May 2000 through January 2025. 98.15% of observations downloaded from GBIF were categorized as originating from iNaturalist (https://www.inaturalist.org), a nature sharing citizen science social network. Observations from GBIF were then binned into counties (n = 478) and those counties were used as a filter for the ensembles.
	Butternut basal area estimates from the FIADB include 220 plots across 141 counties inventoried from 2015 through 2020. Counties with both an ensemble score and empirical basal area estimate were then used as an independent assessment of the ensemble. We performed a single-order linear regression of basal area by ensemble score using base R.


###Files and variables

####File: AllCounties_climate_OSC.csv
Description: Climate and organic soil carbon (OSC) values aggregated at the county level.

Variables

GEOID: 5-digit unique identifier for each county or county-equivilent in the study area. 

mean_norm_ppt_01-12: Mean 30-year normal precipitation for each month (01 - Jan, through 12 - December). Units are in millimeters.

mean_norm_ppt_annual: Mean 30-year annual normal precipitation. Units are in millimeters.

mean_norm_solclear_01-12: Mean 30-year normal total solar radiation incident on a horizontal surface at the surface of the earth when the cloud cover is less than 10%. Units are MJ m^-2 day^-1. (01 - Jan, through 12 - December). 

mean_norm_solclear_annual: Mean 30-year normal annual total solar radiation incident on a horizontal surface at the surface of the earth when the cloud cover is less than 10%. Units are MJ m^-2 day^-1.

mean_norm_soltotal_01-12: Mean 30-year monthly normal amount of the total solar radiation incident on a horizontal surface at the surface of the earth.  Units are MJ m^-2 day^-1 (01 - Jan, through 12 - December).

mean_norm_soltotal_annual: Mean 30-year normal annual amount of the total solar radiation incident on a horizontal surface at the surface of the earth.  Units are MJ m^-2 day^-1

mean_norm_soltrans_01-12: Mean 30-year monthly normal Cloud Transmittance, (ratio of cloud cover, and the units are as such) (01 - Jan, through 12 - December).

mean_norm_soltrans_annual:  Mean 30-year annual normal Cloud Transmittance, (ratio of cloud cover, and the units are as such).

mean_norm_tdmean_01-12: Mean 30-year monthly dew point temperature (01 - Jan, through 12 - December). Units are in degrees Celcius.

mean_norm_tdmean_annual: Mean 30-year annual dew point temperature. Units are in degrees Celcius.

mean_norm_tmax_01-12: Mean 30-year monthly maximum temperature (01 - Jan, through 12 - December). Units are in degrees Celcius.

mean_norm_tmax_annual: Mean 30-year annual maximum temperature. Units are in degrees Celcius.

mean_norm_tmax_01-12: Mean 30-year monthly mean temperature (01 - Jan, through 12 - December). Units are in degrees Celcius.

mean_norm_tmax_annual: Mean 30-year annual mean temperature. Units are in degrees Celcius.

mean_norm_vpdmax_01-12: Mean 30-year monthly maximum vapor pressure deficient (01 - Jan, through 12 - December). Units are kiloPascals (kPa).

mean_norm_vpdmax_annual: Mean 30-year annual maximum vapor pressure deficient. Units are kiloPascals (kPa).

mean_norm_vpdmin_01-12: Mean 30-year monthly minimum vapor pressure deficient (01 - Jan, through 12 - December). Units are kiloPascals (kPa).

mean_norm_vpdmin_annual: Mean 30-year annual minimum vapor pressure deficient. Units are kiloPascals (kPa).

osc_mtu_per_ha: Organic soil carbon, in metric tonnes per hectare, in forested land.

####File: BN_BA_CountyLvlError.R
Description: SQLite code used to query the FIADB to generate basal area estimates of butternut, with sample error, at the county level.

####File: BN_Ensemble_Validation.R
Description: Code used to comapre ensembles with butternut basal area estimates and comparing genotype and phenotype ensembles.

####File: BN_HSM_Butternut_LessResistant.R
Description: Code to generate habitat suitability ensembles for less BCD resistant butternut.

####File: BN_HSM_Butternut_MoreResistant.R
Description: Code to generate habitat suitability ensembles for more BCD resistant butternut.

####File: BN_HSM_Butternut.R
Description: Code to generate habitat suitability ensembles for the butternut genotype.

####File: BN_HSM_Genotype_Ensemble.R
Description: Code to merge habitat suitability ensembles of butternut and hybrid genotypes.

####File: BN_HSM_Hybrid.R
Description: Code to generate habitat suitability ensembles for the hybrid genotype.

####File: BN_HSM_Phenotype_Ensemble.R
Description: Code to merge habitat suitability ensembles of more resistant and less resistant phenotypes.

####File: Butternut_Families.csv
Description: Genotype, phenotype, and location data used to train the HSMs.

Variables

Breeding_Value: BCD resistance phenotype value assigned to a family. Negative is less resistant, positive is more resistant.

species: Genotype, JC is Juglans cinerea (butternut). JXC indicates a hybrid.

latitude: Latitude of family origin location

longitude: Longitude of famiy origin location

state_of_origin: Two-letter abbreviation indicating the state from which the family originiated.

fips:  5-digit unique code for the originating county for each family.

####File:FIA_CountyLevel_OSC.R
Description: SQLite code used to query the FIA database to generate organic soil carbon estimates for forested land in each county.


# Butternut_HSM
