#install.packages("gbm")
library(sf)
library(dplyr)
library(raster)
library(tmap)
library(usmap)
library(tigris)
library(ggplot2)
library(maps)
library(corrplot)
setwd("/Volumes/Onni_Drive/Fearer_Projects/Butternut_Landscape_Genomics")

#model performance data
model_scores <- read.csv("./Results/20250523_HabitSuitMod_2/HSM_Stats.csv", header = TRUE)

#JC Geno data
JC_Models <- read.csv("./Results/20250501_HabitSuitMod/GenoJCAll_ModelsOutputContinuous.csv", header=TRUE)
JC_Models$GEOID <- sprintf("%05d", JC_Models$GEOID)

JC_Models <- JC_Models[,c("GEOID","glm_binary","gam_binary","rf_binary","brt_binary","me_binary")]
JC_Models_scores <- model_scores[model_scores$DataSet == "JCGenoAll",]
w_auc_glm <- JC_Models_scores[1,4]/sum(JC_Models_scores$AUC)
w_auc_gam <- JC_Models_scores[2,4]/sum(JC_Models_scores$AUC)
w_auc_rf <- JC_Models_scores[3,4]/sum(JC_Models_scores$AUC)
w_auc_brt <- JC_Models_scores[4,4]/sum(JC_Models_scores$AUC)
w_auc_me <- JC_Models_scores[5,4]/sum(JC_Models_scores$AUC)

w_tss_glm <- JC_Models_scores[1,7]/sum(JC_Models_scores$TSS)
w_tss_gam <- JC_Models_scores[2,7]/sum(JC_Models_scores$TSS)
w_tss_rf <- JC_Models_scores[3,7]/sum(JC_Models_scores$TSS)
w_tss_brt <- JC_Models_scores[4,7]/sum(JC_Models_scores$TSS)
w_tss_me <- JC_Models_scores[5,7]/sum(JC_Models_scores$TSS)

JC_Models$glm_weighted_auc_adj <- JC_Models$glm_binary*w_auc_glm
JC_Models$glm_weighted_tss_adj <- JC_Models$glm_binary*w_tss_glm

JC_Models$gam_weighted_auc_adj <- JC_Models$gam_binary*w_auc_gam
JC_Models$gam_weighted_tss_adj <- JC_Models$gam_binary*w_tss_gam

JC_Models$rf_weighted_auc_adj <- JC_Models$rf_binary*w_auc_rf
JC_Models$rf_weighted_tss_adj <- JC_Models$rf_binary*w_tss_rf

JC_Models$brt_weighted_auc_adj <- JC_Models$brt_binary*w_auc_brt
JC_Models$brt_weighted_tss_adj <- JC_Models$brt_binary*w_tss_brt

JC_Models$me_weighted_auc_adj <- JC_Models$me_binary*w_auc_me
JC_Models$me_weighted_tss_adj <- JC_Models$me_binary*w_tss_me

JC_Models$ensemble_w_auc_adj <- rowSums(JC_Models[, c("glm_weighted_auc_adj", "gam_weighted_auc_adj", "rf_weighted_auc_adj", "brt_weighted_auc_adj", "me_weighted_auc_adj")], na.rm = TRUE)
JC_Models$ensemble_w_tss_adj <- rowSums(JC_Models[, c("glm_weighted_tss_adj", "gam_weighted_tss_adj", "rf_weighted_tss_adj", "brt_weighted_tss_adj", "me_weighted_tss_adj")], na.rm = TRUE)

#JXC Geno data
JXC_Models <- read.csv("./Results/20250501_HabitSuitMod/GenoJXCAll_ModelsOutputContinuous.csv", header=TRUE)
JXC_Models <- JXC_Models[,c("GEOID","glm_binary","gam_binary","rf_binary","brt_binary","me_binary")]
JXC_Models$GEOID <- sprintf("%05d", JXC_Models$GEOID)
JXC_Models_scores <- model_scores[model_scores$DataSet == "JXCGenoAll",]
w_auc_glm <- JXC_Models_scores[1,4]/sum(JXC_Models_scores$AUC)
w_auc_gam <- JXC_Models_scores[2,4]/sum(JXC_Models_scores$AUC)
w_auc_rf <- JXC_Models_scores[3,4]/sum(JXC_Models_scores$AUC)
w_auc_brt <- JXC_Models_scores[4,4]/sum(JXC_Models_scores$AUC)
w_auc_me <- JXC_Models_scores[5,4]/sum(JXC_Models_scores$AUC)
w_tss_glm <- JXC_Models_scores[1,7]/sum(JXC_Models_scores$TSS)
w_tss_gam <- JXC_Models_scores[2,7]/sum(JXC_Models_scores$TSS)
w_tss_rf <- JXC_Models_scores[3,7]/sum(JXC_Models_scores$TSS)
w_tss_brt <- JXC_Models_scores[4,7]/sum(JXC_Models_scores$TSS)
w_tss_me <- JXC_Models_scores[5,7]/sum(JXC_Models_scores$TSS)

JXC_Models$glm_weighted_auc_adj <- JXC_Models$glm_binary*w_auc_glm
JXC_Models$glm_weighted_tss_adj <- JXC_Models$glm_binary*w_tss_glm

JXC_Models$gam_weighted_auc_adj <- JXC_Models$gam_binary*w_auc_gam
JXC_Models$gam_weighted_tss_adj <- JXC_Models$gam_binary*w_tss_gam

JXC_Models$rf_weighted_auc_adj <- JXC_Models$rf_binary*w_auc_rf
JXC_Models$rf_weighted_tss_adj <- JXC_Models$rf_binary*w_tss_rf

JXC_Models$brt_weighted_auc_adj <- JXC_Models$brt_binary*w_auc_brt
JXC_Models$brt_weighted_tss_adj <- JXC_Models$brt_binary*w_tss_brt

JXC_Models$me_weighted_auc_adj <- JXC_Models$me_binary*w_auc_me
JXC_Models$me_weighted_tss_adj <- JXC_Models$me_binary*w_tss_me

JXC_Models$ensemble_w_auc_adj <- rowSums(JXC_Models[, c("glm_weighted_auc_adj", "gam_weighted_auc_adj", "rf_weighted_auc_adj", "brt_weighted_auc_adj", "me_weighted_auc_adj")], na.rm = TRUE)
JXC_Models$ensemble_w_tss_adj <- rowSums(JXC_Models[, c("glm_weighted_tss_adj", "gam_weighted_tss_adj", "rf_weighted_tss_adj", "brt_weighted_tss_adj", "me_weighted_tss_adj")], na.rm = TRUE)


JXC_Models$ensemble_w_auc_adj_inverted <- -1*JXC_Models$ensemble_w_auc_adj
JXC_Models$ensemble_w_tss_adj_inverted <- -1*JXC_Models$ensemble_w_tss_adj

sub.data1 <- JXC_Models[,c("GEOID","ensemble_w_auc_adj_inverted","ensemble_w_tss_adj_inverted")]
sub.data1[sub.data1$GEOID == "23005",]
JXC_Model_Zeros <- sub.data1[sub.data1$ensemble_w_auc_adj_inverted == 0, c("GEOID")]
sub.data2 <- JC_Models[,c("GEOID","ensemble_w_auc_adj","ensemble_w_tss_adj")]
sub.data2[sub.data2$GEOID == "23005",]
JC_Model_Zeros <- sub.data2[sub.data2$ensemble_w_auc_adj == 0, c("GEOID")]
common_zeros <- intersect(JXC_Model_Zeros,JC_Model_Zeros)
data <- merge(sub.data1, sub.data2, by = "GEOID")
data$ensemble_w_auc_adj_merge <- (data$ensemble_w_auc_adj + data$ensemble_w_auc_adj_inverted)
data$ensemble_w_auc_adj_merge <- ifelse(data$GEOID %in% c(common_zeros), NA , data$ensemble_w_auc_adj_merge)
data$ensemble_w_tss_adj_merge <- (data$ensemble_w_tss_adj + data$ensemble_w_tss_adj_inverted)
data$ensemble_w_tss_adj_merge <- ifelse(data$GEOID %in% c(common_zeros), NA , data$ensemble_w_tss_adj_merge)

write.csv(data, "./Results/20250523_HabitSuitMod_2/JC_JXC_ensemble.csv", row.names=FALSE)

#preprocssed data#### 
data <- read.csv("./Results/20250523_HabitSuitMod_2/JC_JXC_ensemble.csv", header=TRUE)
data$GEOID <- sprintf("%05d", data$GEOID)
nrow(data[!is.na(data$ensemble_w_tss_adj_merge) & data$ensemble_w_tss_adj_merge >= 0,])

JC_LittleRange <- read_sf("/Volumes/Onni_Drive/Fearer_Projects/Butternut_Landscape_Genomics/LittleShapefile/litt601av/litt601av.shp")
littles_proj4 <- "+proj=aea +lat_1=38 +lat_2=42 +lat_0=40 +lon_0=-82 +x_0=0 +y_0=0 +ellps=clrk66 +units=m +no_defs"
st_crs(JC_LittleRange) <- littles_proj4

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

Uamericana_range_transformed <- st_transform(JC_LittleRange, st_crs(us_states_shapefile_simpl))

data_to_map <- merge(data,us_counties_shapefile_simpl, by = "GEOID" )

plot_1 <- ggplot() +
  geom_sf(data = us_counties_shapefile_simpl, aes(geometry = geometry), fill = "grey", color = "NA", show.legend = TRUE) +
  geom_sf(data = data_to_map, aes(geometry = geometry, fill = ensemble_w_tss_adj_merge), color = "NA", show.legend = TRUE) +
  geom_sf(data = us_states_shapefile_simpl, aes(geometry = geometry), fill = "NA", color = "darkgrey", show.legend = TRUE) +
  geom_sf(data = JC_LittleRange, aes(geometry = geometry), color = "blue", fill = "NA", show.legend = TRUE) +
  theme_light() +
  labs(title = "Pure Species & Hybrid Suitability Ensemble Merge",
       subtitle = "Weighted TSS adjustment") +
  scale_fill_gradient2(name = "",
                       low = "darkorchid",
                       mid = "white",
                       high = "limegreen",
                       na.value = "grey",
                       labels = c("Hybrid", "", "", "", "Pure"),
                       breaks = c(min(data_to_map$ensemble_w_tss_adj_merge, na.rm = TRUE),
                                  quantile(data_to_map$ensemble_w_tss_adj_merge, 0.25, na.rm = TRUE),
                                  0,
                                  quantile(data_to_map$ensemble_w_tss_adj_merge, 0.75, na.rm = TRUE),
                                  max(data_to_map$ensemble_w_tss_adj_merge, na.rm = TRUE))) +
  theme(
    plot.title = element_text(size = 12),  
    plot.subtitle = element_text(size = 12),
    axis.title.x = element_text(size = 10),  
    axis.text.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.y = element_text(size = 10))

plot_3 <- ggplot() +
  geom_sf(data = us_counties_shapefile_simpl, aes(geometry = geometry), fill = "grey", color = "NA", show.legend = TRUE) +
  geom_sf(data = data_to_map, aes(geometry = geometry, fill = (-1*ensemble_w_tss_adj_inverted)), color = "NA", show.legend = TRUE) +
  geom_sf(data = us_states_shapefile_simpl, aes(geometry = geometry), fill = "NA", color = "dark grey", show.legend = TRUE) +
  geom_sf(data = JC_LittleRange, aes(geometry = geometry), color = "blue", fill = "NA", show.legend = TRUE) +
  theme_light() +
  labs(title = "Hybrid Suitability Ensemble",
      subtitle = "Weighted TSS adjustment") +
  scale_fill_gradient(name = "",
                       low = "white",
                       high = "darkorchid",
                       #midpoint = 0,
                       na.value = "grey",
                       labels = c("Unsuitable", "", "", "", "Suitable")) +#,
                       #breaks = c(min(data_to_map$ensemble_w_tss_adj_inverted, na.rm = TRUE),
                                  #quantile(data_to_map$ensemble_w_tss_adj_inverted, 0.25, na.rm = TRUE),
                                  #0,
                                  #quantile(data_to_map$ensemble_w_tss_adj_inverted, 0.75, na.rm = TRUE),
                                  #max(data_to_map$ensemble_w_tss_adj_inverted, na.rm = TRUE))) +
  theme(
    plot.title = element_text(size = 12),  
    plot.subtitle = element_text(size = 12),
    axis.title.x = element_text(size = 10),  
    axis.text.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.y = element_text(size = 10))



file_name1 <- paste("./Results/20250523_HabitSuitMod_2/pngs/",
                    "HSM_GenoJC_JXC_TSS_W8d_EnsembleModel",".png", sep = "")
file_name2 <- paste("./Results/20250523_HabitSuitMod_2/pngs/",
                    "HSM_Hybrid_Geno_TSS_W8d_EnsembleModel",".png", sep = "")

file_name3 <- paste("./Results/20250523_HabitSuitMod_2/svgs/",
                    "HSM_GenoJC_JXC_TSS_W8d_EnsembleModel",".svg", sep = "")
file_name4 <- paste("./Results/20250523_HabitSuitMod_2/svgs/",
                    "HSM_Hybrid_Geno_TSS_W8d_EnsembleModel",".svg", sep = "")

ggsave(file_name1, plot_1, width = 7, height = 7, units = "in", dpi = 300)
ggsave(file_name2, plot_3, width = 7, height = 7, units = "in", dpi = 300)
ggsave(file_name3, plot_1, width = 7, height = 7, units = "in", dpi = 300)
ggsave(file_name4, plot_3, width = 7, height = 7, units = "in", dpi = 300)


#Filter by iNaturalist, TreeSnap####
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

#TreeSnap Filter
TreeSnap_Data <- read.csv("./CitSci_JCobvs/TreeSnap_Butternut_06_11_2025.csv", header=TRUE)
TreeSnap_Data <- TreeSnap_Data[,c("Unique.ID","Latitude","Longitude")]

TreeSnap_Data_sf <- st_as_sf(TreeSnap_Data, coords = c("Longitude", "Latitude"), crs = st_crs(us_counties_shapefile))
TreeSnap_Data_object <- st_join(TreeSnap_Data_sf, us_counties_shapefile)
TreeSnap_Data_Counties <- as.data.frame(TreeSnap_Data_object[,c("GEOID")])
length(unique(TreeSnap_Data_Counties$GEOID))
TreeSnap_Data_Counties$geometry <- NULL
TreeSnap_Data_Counties_Mappable <- merge(TreeSnap_Data_Counties,us_counties_shapefile_simpl, by = "GEOID" )

data <- read.csv("./Results/20250523_HabitSuitMod_2/JC_JXC_ensemble.csv", header=TRUE)
data$GEOID <- sprintf("%05d", data$GEOID)
data_treesnap_filtered <- data[data$GEOID %in% c(TreeSnap_Data_Counties$GEOID),]

nrow(data_treesnap_filtered[!is.na(data_treesnap_filtered$ensemble_w_tss_adj_merge),])
data_treesnap_filtered[is.na(data_treesnap_filtered$ensemble_w_tss_adj_merge),]

treesnap_mappable <- merge(data_treesnap_filtered, us_counties_shapefile, by = "GEOID")
  
#map filtered data
tree_snap <- ggplot() +
  geom_sf(data = us_counties_shapefile_simpl, aes(geometry = geometry), fill = "grey", color = "NA", show.legend = TRUE) +
  geom_sf(data = treesnap_mappable, aes(geometry = geometry,fill = ensemble_w_tss_adj_merge), color = "NA", show.legend = TRUE) +
  geom_sf(data = us_states_shapefile_simpl, aes(geometry = geometry), fill = "NA", color = "dark grey", show.legend = TRUE) +
  theme_light() +
  labs(title = "Species vs Hybrid Merged Ensemble \nTreesnap Filter",
       subtitle = "Weighted TSS adjustment") +
  scale_fill_gradient2(name = "",
                       low = "darkorchid",
                       high = "limegreen",
                       midpoint = 0,
                       na.value = "black",
                       labels = c("Hybrid", "", "", "", "Species"),
                       breaks = c(min(treesnap_mappable$ensemble_w_tss_adj_merge, na.rm = TRUE),
                                  quantile(treesnap_mappable$ensemble_w_tss_adj_merge, 0.25, na.rm = TRUE),
                                  0,
                                  quantile(treesnap_mappable$ensemble_w_tss_adj_merge, 0.75, na.rm = TRUE),
                                  max(treesnap_mappable$ensemble_w_tss_adj_merge, na.rm = TRUE))) +
  theme(
    plot.title = element_text(size = 12),  
    plot.subtitle = element_text(size = 12),
    axis.title.x = element_text(size = 10),  
    axis.text.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.y = element_text(size = 10)) 

file_name1 <- paste("./Results/20250523_HabitSuitMod_2/pngs/",
                    "HSM_GenoJC_JXC_TSS_W8d_EnsembleModel_Merge_TreeSnap",".png", sep = "")
file_name2 <- paste("./Results/20250523_HabitSuitMod_2/svgs/",
                    "HSM_GenoJC_JXC_TSS_W8d_EnsembleModel_Merge_TreeSnap",".svg", sep = "")
ggsave(file_name1, tree_snap, width = 7, height = 7, units = "in", dpi = 300)
ggsave(file_name2, tree_snap, width = 7, height = 7, units = "in", dpi = 300)


data_treesnap_filtered[order(data_treesnap_filtered$ensemble_w_tss_adj_merge,decreasing = TRUE),c("GEOID","ensemble_w_tss_adj_merge")]



#iNaturalist Filter
iNaturalist_Data <- read.csv("./CitSci_JCobvs/iNaturalist/observations-585230.csv", header=TRUE)
iNaturalist_Data <- iNaturalist_Data[,c("id","latitude","longitude")]
iNaturalist_Data_sf <- st_as_sf(iNaturalist_Data, coords = c("longitude", "latitude"), crs = st_crs(us_counties_shapefile_simpl))
iNaturalist_Data_object <- st_join(iNaturalist_Data_sf, us_counties_shapefile_simpl)
iNaturalist_Data_Counties <- as.data.frame(iNaturalist_Data_object[,c("GEOID")])
iNaturalist_Data_Counties$geometry <- NULL
iNaturalist_Data_Counties_Mappable <- merge(iNaturalist_Data_Counties,us_counties_shapefile_simpl, by = "GEOID" )

data_inaturalist_filtered <- data[data$GEOID %in% c(iNaturalist_Data_Counties$GEOID),]

inaturalist_mappable <- merge(data_inaturalist_filtered, us_counties_shapefile, by = "GEOID")


iNaturalist <- ggplot() +
  geom_sf(data = us_counties_shapefile_simpl, aes(geometry = geometry), fill = "grey", color = "NA", show.legend = TRUE) +
  geom_sf(data = inaturalist_mappable, aes(geometry = geometry,fill = ensemble_w_tss_adj_merge), color = "NA", show.legend = TRUE) +
  geom_sf(data = us_states_shapefile_simpl, aes(geometry = geometry), fill = "NA", color = "dark grey", show.legend = TRUE) +
  theme_light() +
  labs(title = "Species vs Hybrid Merged Ensemble \niNaturalist Filter",
       subtitle = "Weighted TSS adjustment") +
  scale_fill_gradient2(name = "",
                       low = "darkorchid",
                       high = "limegreen",
                       midpoint = 0,
                       na.value = "black",
                       labels = c("Hybrid", "", "", "", "Species"),
                       breaks = c(min(inaturalist_mappable$ensemble_w_tss_adj_merge, na.rm = TRUE),
                                  quantile(inaturalist_mappable$ensemble_w_tss_adj_merge, 0.25, na.rm = TRUE),
                                  0,
                                  quantile(inaturalist_mappable$ensemble_w_tss_adj_merge, 0.75, na.rm = TRUE),
                                  max(inaturalist_mappable$ensemble_w_tss_adj_merge, na.rm = TRUE))) +
  theme(
    plot.title = element_text(size = 12),  
    plot.subtitle = element_text(size = 12),
    axis.title.x = element_text(size = 10),  
    axis.text.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.y = element_text(size = 10)) 

file_name1 <- paste("./Results/20250523_HabitSuitMod_2/pngs/",
                    "HSM_GenoJC_JXC_TSS_W8d_EnsembleModel_Merge_iNaturalist",".png", sep = "")
file_name2 <- paste("./Results/20250523_HabitSuitMod_2/svgs/",
                    "HSM_GenoJC_JXC_TSS_W8d_EnsembleModel_Merge_iNaturalist",".svg", sep = "")
ggsave(file_name1, iNaturalist, width = 7, height = 7, units = "in", dpi = 300)
ggsave(file_name2, iNaturalist, width = 7, height = 7, units = "in", dpi = 300)

#GBIF filter####
GBIF_Data <- read.csv("./CitSci_JCobvs/GBIF_Jcinerea.csv", header=TRUE)
GBIF_Data <- GBIF_Data[complete.cases(GBIF_Data$decimalLatitude),c("gbifID","decimalLatitude","decimalLongitude")]
GBIF_Data$longitude <- GBIF_Data$decimalLongitude
GBIF_Data$latitude <- GBIF_Data$decimalLatitude
GBIF_Data$decimalLatitude <- NULL
GBIF_Data$decimalLongitude <- NULL

GBIF_Data_sf <- st_as_sf(GBIF_Data, coords = c("longitude", "latitude"), crs = st_crs(us_counties_shapefile))
GBIF_Data_object <- st_join(GBIF_Data_sf, us_counties_shapefile_simpl)
GBIF_Data_Counties <- GBIF_Data_object[!(GBIF_Data_object$GEOID %in% as.character(c(51013,51510,51520,51530,51540,51570,51580,51590,51595,51620,51630,51640,51650,
                                                                                    51660,51670,51678,51680,51683,51690,51700,51710,51720,51730,51735,51740,51750,51760,
                                                                                    51770,51775,51790,51800,51810,51820,51830,51840)))
                                       , c("GEOID")]
length(unique(GBIF_Data_Counties$GEOID))

GBIF_Data_Counties$geometry <- NULL
GBIF_Data_Counties_Mappable <- merge(GBIF_Data_Counties,us_counties_shapefile_simpl, by = "GEOID" )


data <- read.csv("./Results/20250523_HabitSuitMod_2/JC_JXC_ensemble.csv", header=TRUE)
data$GEOID <- sprintf("%05d", data$GEOID)
nrow(data[!is.na(data$ensemble_w_tss_adj_merge),])

data_gbif_filtered <- data[data$GEOID %in% c(GBIF_Data_Counties$GEOID),]

data_3 <- data_gbif_filtered[order(data_gbif_filtered$ensemble_w_tss_adj_merge,decreasing = TRUE),
                                 c("GEOID","ensemble_w_tss_adj_merge")]
head(data_3)


nrow(data_gbif_filtered[!is.na(data_gbif_filtered$ensemble_w_tss_adj_merge),])

gbif_mappable <- merge(data_gbif_filtered, us_counties_shapefile, by = "GEOID")


gbif <- ggplot() +
  geom_sf(data = us_counties_shapefile_simpl, aes(geometry = geometry), fill = "grey", color = "NA", show.legend = TRUE) +
  geom_sf(data = gbif_mappable, aes(geometry = geometry,fill = ensemble_w_tss_adj_merge), color = "NA", show.legend = TRUE) +
  geom_sf(data = us_states_shapefile_simpl, aes(geometry = geometry), fill = "NA", color = "dark grey", show.legend = TRUE) +
  theme_light() +
  labs(title = "Species vs Hybrid Merged Ensemble \nGBIF Filter",
       subtitle = "Weighted TSS adjustment") +
  scale_fill_gradient2(name = "",
                       low = "darkorchid",
                       high = "limegreen",
                       midpoint = 0,
                       na.value = "black",
                       labels = c("Hybrid", "", "", "", "Species"),
                       breaks = c(min(gbif_mappable$ensemble_w_tss_adj_merge, na.rm = TRUE),
                                  quantile(gbif_mappable$ensemble_w_tss_adj_merge, 0.25, na.rm = TRUE),
                                  0,
                                  quantile(gbif_mappable$ensemble_w_tss_adj_merge, 0.75, na.rm = TRUE),
                                  max(gbif_mappable$ensemble_w_tss_adj_merge, na.rm = TRUE))) +
  theme(
    plot.title = element_text(size = 12),  
    plot.subtitle = element_text(size = 12),
    axis.title.x = element_text(size = 10),  
    axis.text.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.y = element_text(size = 10)) 

file_name1 <- paste("./Results/20250523_HabitSuitMod_2/pngs/",
                    "HSM_GenoJC_JXC_TSS_W8d_EnsembleModel_Merge_GBIF",".png", sep = "")
file_name2 <- paste("./Results/20250523_HabitSuitMod_2/svgs/",
                    "HSM_GenoJC_JXC_TSS_W8d_EnsembleModel_Merge_GBIF",".svg", sep = "")
ggsave(file_name1, gbif, width = 7, height = 7, units = "in", dpi = 300)
ggsave(file_name2, gbif, width = 7, height = 7, units = "in", dpi = 300)

data_gbif_filtered[order(data_gbif_filtered$ensemble_w_tss_adj_merge, decreasing = TRUE), c("GEOID","ensemble_w_tss_adj_merge")]

#FIA filter####
FIA_Data <- read.csv("./CitSci_JCobvs/GBIF_Jcinerea.csv", header=TRUE)
FIA_Data <- read.csv("./Results/20250207_PlotLvlError/County_Level_Error_Biomass_2020_v2.csv", header=TRUE)
FIA_Data$county <- sprintf("%03d", FIA_Data$county)
FIA_Data$GEOID <- as.numeric(paste(FIA_Data$state, FIA_Data$county, sep=""))
FIA_Data$GEOID <- sprintf("%05d", FIA_Data$GEOID)
FIA_Data <- FIA_Data[,c("GEOID","mtu_per_ha","se_mtu_per_ha")]

length(unique(FIA_Data$GEOID))

FIA_Data_Counties_Mappable <- merge(FIA_Data,us_counties_shapefile_simpl, by = "GEOID" )

data <- read.csv("./Results/20250523_HabitSuitMod_2/JC_JXC_ensemble.csv", header=TRUE)
data$GEOID <- sprintf("%05d", data$GEOID)
nrow(data[!is.na(data$ensemble_w_tss_adj_merge),])

data_fia_filtered <- data[data$GEOID %in% c(FIA_Data$GEOID),]
nrow(data_fia_filtered[!is.na(data_fia_filtered$ensemble_w_tss_adj_merge),])

fia_mappable <- merge(data_fia_filtered, us_counties_shapefile, by = "GEOID")


fia <- ggplot() +
  geom_sf(data = us_counties_shapefile_simpl, aes(geometry = geometry), fill = "grey", color = "NA", show.legend = TRUE) +
  geom_sf(data = fia_mappable, aes(geometry = geometry,fill = ensemble_w_tss_adj_merge), color = "NA", show.legend = TRUE) +
  geom_sf(data = us_states_shapefile_simpl, aes(geometry = geometry), fill = "NA", color = "dark grey", show.legend = TRUE) +
  theme_light() +
  labs(title = "Final Genotype Ensemble \nFIA Filter",
       subtitle = "Weighted TSS adjustment") +
  scale_fill_gradient2(name = "",
                       low = "darkorchid",
                       high = "limegreen",
                       midpoint = 0,
                       na.value = "black",
                       labels = c("Hybrid", "", "", "", "Species"),
                       breaks = c(min(fia_mappable$ensemble_w_tss_adj_merge, na.rm = TRUE),
                                  quantile(fia_mappable$ensemble_w_tss_adj_merge, 0.25, na.rm = TRUE),
                                  0,
                                  quantile(fia_mappable$ensemble_w_tss_adj_merge, 0.75, na.rm = TRUE),
                                  max(fia_mappable$ensemble_w_tss_adj_merge, na.rm = TRUE))) +
  theme(
    plot.title = element_text(size = 12),  
    plot.subtitle = element_text(size = 12),
    axis.title.x = element_text(size = 10),  
    axis.text.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.y = element_text(size = 10)) 

file_name1 <- paste("./Results/20250523_HabitSuitMod_2/pngs/",
                    "HSM_GenoJC_JXC_TSS_W8d_EnsembleModel_Merge_FIA",".png", sep = "")
file_name2 <- paste("./Results/20250523_HabitSuitMod_2/svgs/",
                    "HSM_GenoJC_JXC_TSS_W8d_EnsembleModel_Merge_FIA",".svg", sep = "")
ggsave(file_name1, fia, width = 7, height = 7, units = "in", dpi = 300)
ggsave(file_name2, fia, width = 7, height = 7, units = "in", dpi = 300)
