library(tmap)
library(usmap)
library(sf)
library(tigris)
library(ggplot2)
library(gridExtra)
library(ggVennDiagram)


setwd("/Volumes/Onni_Drive/Fearer_Projects/Butternut_Landscape_Genomics")
ba_data <- read.csv("./Results/20250207_PlotLvlError/County_Level_Error_BasalArea_2020_v2.csv", header=TRUE)
ba_data$county <- sprintf("%03d", ba_data$county)
ba_data$GEOID <- as.numeric(paste(ba_data$state, ba_data$county, sep=""))
ba_data$GEOID <- sprintf("%05d", ba_data$GEOID)
ba_data <- ba_data[,c("GEOID","ba_msqrd_per_ha","se_ba_msqrd_per_ha")]
nrow(ba_data)
expvol_FIA <- merge(biomass_data, ba_data, by = "GEOID")
nrow(expvol_FIA)


#Getting numbers from the ensembles####
RnS_ensemble <- read.csv("./Results/20250523_HabitSuitMod_2/GenoJC_RnS_ensemble.csv", header=TRUE)
nrow(RnS_ensemble[!(is.na(RnS_ensemble$ensemble_w_tss_adj_merge)),])
"1040 counties have a score of some kind"
nrow(RnS_ensemble[!(is.na(RnS_ensemble$ensemble_w_tss_adj_merge)) 
                      & RnS_ensemble$ensemble_w_tss_adj_merge > 0,])
"457 counties in the final ensemble were predicted to be more suitable to resilient than 
susceptible pure butternut"
nrow(RnS_ensemble[RnS_ensemble$ensemble_w_tss_adj > 0,])
"718 counties are at least moderately suitable to resilient butternut"
nrow(RnS_ensemble[RnS_ensemble$ensemble_w_tss_adj == 1,])
"65 counties have the highest degree of confidence for resilient butternut suitability"
nrow(RnS_ensemble[RnS_ensemble$ensemble_w_tss_adj_inverted < 0,])
"786 counties are at least moderately suitable to susceptible butternut"
nrow(RnS_ensemble[RnS_ensemble$ensemble_w_tss_adj_inverted == -1,])



Hybrids_ensemble <- read.csv("./Results/20250523_HabitSuitMod_2/JC_JXC_ensemble.csv", header=TRUE)
nrow(Hybrids_ensemble[!(is.na(Hybrids_ensemble$ensemble_w_tss_adj_merge)),])
"1105 counties with a score of some kind"
nrow(Hybrids_ensemble[!(is.na(Hybrids_ensemble$ensemble_w_tss_adj_merge)) 
                      & Hybrids_ensemble$ensemble_w_tss_adj_merge > 0,])
"365 counties have a posiive (pure butternut) score"
nrow(Hybrids_ensemble[Hybrids_ensemble$ensemble_w_tss_adj > 0,])
"829 counties are at least moderately suitable to pure butternut"
nrow(Hybrids_ensemble[Hybrids_ensemble$ensemble_w_tss_adj == 1,])
"151 counties have the highest degree of confidence of butternut suitability"

nrow(Hybrids_ensemble[Hybrids_ensemble$ensemble_w_tss_adj_inverted < 0,])
"934 counties are at least moderately suitable to hybrid butternut"
nrow(Hybrids_ensemble[Hybrids_ensemble$ensemble_w_tss_adj_inverted == -1,])
"220 counties have the highest degree of confidience of hybrid suitability"



##Empirical FIA data and the Phenotype Ensemble####
RnS_ensemble <- read.csv("./Results/20250523_HabitSuitMod_2/GenoJC_RnS_ensemble.csv", header=TRUE)
RnS_ensemble$GEOID <- sprintf("%05d", RnS_ensemble$GEOID)
nrow(RnS_ensemble)

RnS_ensemble <- merge(RnS_ensemble, expvol_FIA, by = "GEOID", all = TRUE)

data_fia_filtered <- RnS_ensemble[RnS_ensemble$GEOID %in% c(expvol_FIA$GEOID),]
nrow(data_fia_filtered[!is.na(data_fia_filtered$ensemble_w_tss_adj_merge),])

lm_2 <- lm(ba_msqrd_per_ha ~ ensemble_w_tss_adj_merge, data = RnS_ensemble )
summary(lm_2)
summary(lm_2)$r.squared


BA_plot <- ggplot(data = RnS_ensemble) +
  geom_errorbar(data = RnS_ensemble, aes(x = ensemble_w_tss_adj_merge, 
                                         ymin = ba_msqrd_per_ha - se_ba_msqrd_per_ha, 
                                         ymax = ba_msqrd_per_ha + se_ba_msqrd_per_ha),
                width = 0.01, color = "black") +
 labs(title = " ",
       y = "JC Basal Area \n (m^2/ha)",
       x = "Merged Ensemble, \nTSS Weighted") +
  geom_point(aes(x = ensemble_w_tss_adj_merge, y = ba_msqrd_per_ha)) +
  geom_smooth(method = "lm", data = RnS_ensemble, 
              aes(x = ensemble_w_tss_adj_merge, y = ba_msqrd_per_ha),
              color = "red",
              se = ifelse(summary(lm_2)$coefficients[2, 4] > 0.05, FALSE, TRUE), 
              linetype = ifelse(summary(lm_2)$coefficients[2, 4] > 0.05, "dashed", "solid")) +
  theme_minimal()

file_name3 <- paste("./Results/20250523_HabitSuitMod_2/pngs/",
                    "HSM_RnSvFIA_scatterplots",".png", sep = "")
file_name4 <- paste("./Results/20250523_HabitSuitMod_2/svgs/",
                    "HSM_RnSvFIA_scatterplots",".svg", sep = "")

ggsave(file_name3, BA_plot, width = 7, height = 3.5, units = "in", dpi = 300)
ggsave(file_name4, BA_plot, width = 7, height = 3.5, units = "in", dpi = 300)


##Empirical FIA and Genotype Models#####
Hybrids_ensemble <- read.csv("./Results/20250523_HabitSuitMod_2/JC_JXC_ensemble.csv", header=TRUE)
Hybrids_ensemble$GEOID <- sprintf("%05d", Hybrids_ensemble$GEOID)
data <- Hybrids_ensemble
Hybrids_ensemble <- merge(Hybrids_ensemble, expvol_FIA, by = "GEOID", all = TRUE)
Hybrids_ensemble[!(is.na(Hybrids_ensemble$ensemble_w_tss_adj_merge)),]

lm_4 <- lm(ba_msqrd_per_ha ~ ensemble_w_tss_adj_merge, data = Hybrids_ensemble )
summary(lm_4)


BA_plot <- ggplot(data = Hybrids_ensemble) +
  geom_errorbar(data = Hybrids_ensemble, aes(x = ensemble_w_tss_adj_merge, 
                                         ymin = ba_msqrd_per_ha - se_ba_msqrd_per_ha, 
                                         ymax = ba_msqrd_per_ha + se_ba_msqrd_per_ha),
                width = 0.01, color = "black") +
  labs(title = " ",
       y = "JC Basal Area \n (m^2/ha)",
       x = "Merged Ensemble, \nTSS Weighted") +
  geom_point(aes(x = ensemble_w_tss_adj_merge, y = ba_msqrd_per_ha)) +
  geom_smooth(method = "lm", data = Hybrids_ensemble, 
              aes(x = ensemble_w_tss_adj_merge, y = ba_msqrd_per_ha),
              color = "red",
              se = ifelse(summary(lm_4)$coefficients[2, 4] > 0.05, FALSE, TRUE), 
              linetype = ifelse(summary(lm_4)$coefficients[2, 4] > 0.05, "dashed", "solid")) +
  theme_minimal()

arranged_plots <- grid.arrange(Biomass_plot,BA_plot, ncol = 2)

file_name3 <- paste("./Results/20250523_HabitSuitMod_2/pngs/",
                    "HSM_HybvFIA_scatterplots",".png", sep = "")
file_name4 <- paste("./Results/20250523_HabitSuitMod_2/svgs/",
                    "HSM_HybvFIA_scatterplots",".svg", sep = "")
ggsave(file_name3, BA_plot, width = 7, height = 3.5, units = "in", dpi = 300)
ggsave(file_name4, BA_plot, width = 7, height = 3.5, units = "in", dpi = 300)





##Align after TreeSnap filter
us_counties_shapefile <- counties(cb = TRUE)
us_counties_shapefile <- us_counties_shapefile[us_counties_shapefile$STUSPS %in% c("ME","NH","VT","RI","CT","MA","NY","PA",
                                                                                   "DE","WV","NJ","MD","VA","OH","SC","NC",
                                                                                   "TN","KY","GA","FL","AL","MS","IN","IL",
                                                                                   "MI","MN","IA","AR","MO","LA","WI"),]
us_counties_shapefile_simpl <- st_simplify(us_counties_shapefile, dTolerance = 1000)
us_states_shapefile <- states(cb = TRUE)
us_states_shapefile <- us_states_shapefile[us_states_shapefile$STUSPS %in% c("ME","NH","VT","RI","CT","MA","NY","PA",
                                                                             "DE","WV","NJ","MD","VA","OH","SC","NC",
                                                                             "TN","KY","GA","FL","AL","MS","IN","IL",
                                                                             "MI","MN","IA","AR","MO","LA","WI"),]
us_states_shapefile_simpl <- st_simplify(us_states_shapefile, dTolerance = 1000)

#ensemble data
RnS_ensemble <- read.csv("./Results/20250523_HabitSuitMod_2/GenoJC_RnS_ensemble.csv", header=TRUE)
RnS_ensemble$GEOID <- sprintf("%05d", RnS_ensemble$GEOID)
RnS_ensemble <- RnS_ensemble[,c("GEOID","ensemble_w_tss_adj_merge")]
RnS_ensemble$RnS_merge <- RnS_ensemble$ensemble_w_tss_adj_merge
RnS_ensemble$ensemble_w_tss_adj_merge <- NULL

Hybrids_ensemble <- read.csv("./Results/20250523_HabitSuitMod_2/JC_JXC_ensemble.csv", header=TRUE)
Hybrids_ensemble$GEOID <- sprintf("%05d", Hybrids_ensemble$GEOID)
Hybrids_ensemble <- Hybrids_ensemble[,c("GEOID","ensemble_w_tss_adj_merge")]
Hybrids_ensemble$SvH_merge <- Hybrids_ensemble$ensemble_w_tss_adj_merge
Hybrids_ensemble$ensemble_w_tss_adj_merge <- NULL

data <- merge(RnS_ensemble, Hybrids_ensemble, by = "GEOID")
data$color_guide <- data$RnS_merge + data$SvH_merge

#TreeSnap Filter####
TreeSnap_Data <- read.csv("./CitSci_JCobvs/TreeSnap_Butternut_06_11_2025.csv", header=TRUE)
TreeSnap_Data <- TreeSnap_Data[,c("Unique.ID","Latitude","Longitude")]
TreeSnap_Data_sf <- st_as_sf(TreeSnap_Data, coords = c("Longitude", "Latitude"), crs = st_crs(us_counties_shapefile))
TreeSnap_Data_object <- st_join(TreeSnap_Data_sf, us_counties_shapefile)
TreeSnap_Data_Counties <- as.data.frame(TreeSnap_Data_object[,c("GEOID")])
TreeSnap_Data_Counties$geometry <- NULL
TreeSnap_Data_Counties_Mappable <- merge(TreeSnap_Data_Counties,us_counties_shapefile_simpl, by = "GEOID" )


data_treesnap_filtered <- data[data$GEOID %in% c(unique(TreeSnap_Data_Counties$GEOID)),]
data_treesnap_filtered <- data_treesnap_filtered[complete.cases(data_treesnap_filtered$color_guide),]
nrow(data_treesnap_filtered)
data_3 <- data_treesnap_filtered[order(data_treesnap_filtered$color_guide,decreasing = TRUE),]
head(data_3)

#GBIF Filter####
GBIF_Data <- read.csv("./CitSci_JCobvs/GBIF_Jcinerea.csv", header=TRUE)
GBIF_Data <- GBIF_Data[complete.cases(GBIF_Data$decimalLatitude),c("gbifID","decimalLatitude","decimalLongitude","institutionCode")]
GBIF_Data$longitude <- GBIF_Data$decimalLongitude
GBIF_Data$latitude <- GBIF_Data$decimalLatitude
GBIF_Data$decimalLatitude <- NULL
GBIF_Data$decimalLongitude <- NULL

nrow(GBIF_Data[GBIF_Data$institutionCode == "iNaturalist",])/nrow(GBIF_Data)*100

GBIF_Data_sf <- st_as_sf(GBIF_Data, coords = c("longitude", "latitude"), crs = st_crs(us_counties_shapefile_simpl))
GBIF_Data_object <- st_join(GBIF_Data_sf, us_counties_shapefile)
#as.data.frame(GBIF_Data_object[is.na(GBIF_Data_object$STATEFP),])
GBIF_Data_Counties <- GBIF_Data_object[!(GBIF_Data_object$GEOID %in% as.character(c(51013,51510,51520,51530,51540,51570,51580,51590,51595,51620,51630,51640,51650,
                                                                                    51660,51670,51678,51680,51683,51690,51700,51710,51720,51730,51735,51740,51750,51760,
                                                                                    51770,51775,51790,51800,51810,51820,51830,51840)))
                                       , c("GEOID")]
GBIF_Data_Counties$geometry <- NULL
GBIF_Data_Counties_Mappable <- merge(GBIF_Data_Counties,us_counties_shapefile, by = "GEOID" )
length(unique(GBIF_Data_Counties$GEOID))

data_gbif_filtered <- data[data$GEOID %in% c(GBIF_Data_Counties$GEOID),]
data_gbif_filtered <- data_gbif_filtered[complete.cases(data_gbif_filtered$color_guide),]
data_4 <- data_gbif_filtered[order(data_gbif_filtered$color_guide,decreasing = TRUE),]
head(data_4)

#Do the ensembles align?####
RnS_ensemble <- read.csv("./Results/20250523_HabitSuitMod_2/GenoJC_RnS_ensemble.csv", header=TRUE)
RnS_ensemble$GEOID <- sprintf("%05d", RnS_ensemble$GEOID)
RnS_ensemble <- RnS_ensemble[,c("GEOID","ensemble_w_tss_adj_merge")]
RnS_ensemble$RnS_merge <- RnS_ensemble$ensemble_w_tss_adj_merge
RnS_ensemble$ensemble_w_tss_adj_merge <- NULL

Hybrids_ensemble <- read.csv("./Results/20250523_HabitSuitMod_2/JC_JXC_ensemble.csv", header=TRUE)
Hybrids_ensemble$GEOID <- sprintf("%05d", Hybrids_ensemble$GEOID)
Hybrids_ensemble <- Hybrids_ensemble[,c("GEOID","ensemble_w_tss_adj_merge")]
Hybrids_ensemble$SvH_merge <- Hybrids_ensemble$ensemble_w_tss_adj_merge
Hybrids_ensemble$ensemble_w_tss_adj_merge <- NULL

data <- merge(RnS_ensemble, Hybrids_ensemble, by = "GEOID")

data$color_guide <- data$RnS_merge + data$SvH_merge

TreeSnap_Data_Counties$geometry <- NULL
GBIF_Data_Counties

data$TreeSnap_Dataset <- ""
data$TreeSnap_Dataset <- ifelse(data$GEOID %in% c(TreeSnap_Data_Counties$GEOID), TRUE, FALSE)
data$GBIF_Dataset <- ifelse(data$GEOID %in% c(GBIF_Data_Counties$GEOID), TRUE, FALSE)

county_names <- us_counties_shapefile[,c("GEOID","NAME","STATE_NAME")]
county_names$geometry <- NULL

data <- merge(data, county_names, by = "GEOID")

write.csv(data, "./Results/20250523_HabitSuitMod_2/Combined_Enembles_Filters.csv", row.names=FALSE)

lm <- lm(SvH_merge ~ RnS_merge, data = data)
summary(lm)
summary(lm)$fstatistic
summary(lm)$df[2]
summary(lm)$coefficients[2,4]
summary(lm)$r.squared
summary(lm)$coefficients[2,1]
summary(lm)$coefficients[2,2]

head(data[order(data$color_guide,decreasing = TRUE),])
data[order(data$RnS_merge,decreasing = TRUE),]


ensemble_alignment <- ggplot() +
  geom_point(data = data, aes(x = RnS_merge, y = SvH_merge, color = color_guide, fill = color_guide), alpha = 0.7) +
  scale_fill_gradient2(low = "navy",
                       mid = "lightgrey",
                       high = "gold",
                       midpoint = 0,
                       na.value = "grey") +
  scale_color_gradient2(low = "navy",
                       mid = "lightgrey",
                       high = "gold",
                       midpoint = 0,
                       na.value = "grey") +
  geom_smooth(method = "lm", data = data, aes(x = RnS_merge, y = SvH_merge),
              se = ifelse(summary(lm)$coefficients[2, 4] > 0.05, FALSE, TRUE), 
              linetype = ifelse(summary(lm)$coefficients[2, 4] > 0.05, "dashed", "solid")) +
  labs( y = "Hybrid <---------> Pure Butternut",
        x = "Suscept <---------> Resili",
        color = "Ensemble Sum",) +
  theme_light()

file_name1 <- paste("./Results/20250523_HabitSuitMod_2/pngs/",
                    "HSM_ensemble_alignment",".png", sep = "")
file_name2 <- paste("./Results/20250523_HabitSuitMod_2/svgs/",
                    "HSM_ensemble_alignment",".svg", sep = "")
ggsave(file_name1, ensemble_alignment, width = 5, height = 3, units = "in", dpi = 300)
ggsave(file_name2, ensemble_alignment, width = 5, height = 3, units = "in", dpi = 300)

RnS_ensemble <- read.csv("./Results/20250523_HabitSuitMod_2/GenoJC_RnS_ensemble.csv", header=TRUE)
RnS_ensemble$GEOID <- sprintf("%05d", RnS_ensemble$GEOID)
RnS_ensemble <- RnS_ensemble[,c("GEOID","ensemble_w_tss_adj_merge")]
RnS_ensemble$RnS_merge <- RnS_ensemble$ensemble_w_tss_adj_merge
RnS_ensemble$ensemble_w_tss_adj_merge <- NULL

Hybrids_ensemble <- read.csv("./Results/20250523_HabitSuitMod_2/JC_JXC_ensemble.csv", header=TRUE)
Hybrids_ensemble$GEOID <- sprintf("%05d", Hybrids_ensemble$GEOID)
Hybrids_ensemble <- Hybrids_ensemble[,c("GEOID","ensemble_w_tss_adj_merge")]
Hybrids_ensemble$SvH_merge <- Hybrids_ensemble$ensemble_w_tss_adj_merge
Hybrids_ensemble$ensemble_w_tss_adj_merge <- NULL

data <- merge(RnS_ensemble, Hybrids_ensemble, by = "GEOID")
data$SvH_merge <- data$SvH_merge*-1

data$color_guide <- data$RnS_merge + data$SvH_merge
us_counties_shapefile <- counties(cb = TRUE)
us_counties_shapefile <- us_counties_shapefile[us_counties_shapefile$STUSPS %in% c("ME","NH","VT","RI","CT","MA","NY","PA",
                                                                                   "DE","WV","NJ","MD","VA","OH","SC","NC",
                                                                                   "TN","KY","GA","FL","AL","MS","IN","IL",
                                                                                   "MI","MN","IA","AR","MO","LA","WI"),]
us_counties_shapefile_simpl <- st_simplify(us_counties_shapefile, dTolerance = 1000)
us_states_shapefile <- states(cb = TRUE)
us_states_shapefile <- us_states_shapefile[us_states_shapefile$STUSPS %in% c("ME","NH","VT","RI","CT","MA","NY","PA",
                                                                             "DE","WV","NJ","MD","VA","OH","SC","NC",
                                                                             "TN","KY","GA","FL","AL","MS","IN","IL",
                                                                             "MI","MN","IA","AR","MO","LA","WI"),]
us_states_shapefile_simpl <- st_simplify(us_states_shapefile, dTolerance = 1000)
Ensmble_Merge_Mappable <- merge(data,us_counties_shapefile_simpl, by = "GEOID" )

ensemble_alignment_map <- ggplot() +
  geom_sf(data = us_counties_shapefile_simpl, aes(geometry = geometry), fill = "grey", color = "NA", show.legend = TRUE) +
  geom_sf(data = Ensmble_Merge_Mappable, aes(geometry = geometry,fill = color_guide), color = "NA", show.legend = TRUE) +
  geom_sf(data = us_states_shapefile_simpl, aes(geometry = geometry), fill = "NA", color = "dark grey", show.legend = TRUE) +
  theme_light() +
  labs(title = " ",
       subtitle = " ",
       fill = "Ensemble Sum") +
  scale_fill_gradient2(low = "navy",
                       mid = "whitesmoke",
                       high = "gold",
                       midpoint = 0,
                       na.value = "grey") +
  theme(
    plot.title = element_text(size = 12),  
    plot.subtitle = element_text(size = 12),
    axis.title.x = element_text(size = 10),  
    axis.text.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.y = element_text(size = 10)) 

file_name1 <- paste("./Results/20250523_HabitSuitMod_2/pngs/",
                    "HSM_ensemble_alignment_map",".png", sep = "")
file_name2 <- paste("./Results/20250523_HabitSuitMod_2/svgs/",
                    "HSM_ensemble_alignment_map",".svg", sep = "")
ggsave(file_name1, ensemble_alignment_map, width = 5, height = 3, units = "in", dpi = 300)
ggsave(file_name2, ensemble_alignment_map, width = 5, height = 3, units = "in", dpi = 300)
