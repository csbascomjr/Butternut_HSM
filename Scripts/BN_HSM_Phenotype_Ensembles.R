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
model_scores <- read.csv("./Results/20250523_HabitSuitMod_2/HSM_Performance_Stats.csv", header = TRUE)

#Resistant JC Geno data
Resist_Models <- read.csv("./Results/20250501_HabitSuitMod/GenoJCResist_ModelsOutputContinuous.csv", header=TRUE)
Resist_Models <- Resist_Models[,c("GEOID","glm_binary","gam_binary","rf_binary","brt_binary","me_binary")]
Resist_Model_scores <- model_scores[model_scores$DataSet == "JCGenoResist",]
w_auc_glm <- Resist_Model_scores[1,4]/sum(Resist_Model_scores$AUC)
w_auc_gam <- Resist_Model_scores[2,4]/sum(Resist_Model_scores$AUC)
w_auc_rf <- Resist_Model_scores[3,4]/sum(Resist_Model_scores$AUC)
w_auc_brt <- Resist_Model_scores[4,4]/sum(Resist_Model_scores$AUC)
w_auc_me <- Resist_Model_scores[5,4]/sum(Resist_Model_scores$AUC)

w_tss_glm <- Resist_Model_scores[1,7]/sum(Resist_Model_scores$TSS)
w_tss_gam <- Resist_Model_scores[2,7]/sum(Resist_Model_scores$TSS)
w_tss_rf <- Resist_Model_scores[3,7]/sum(Resist_Model_scores$TSS)
w_tss_brt <- Resist_Model_scores[4,7]/sum(Resist_Model_scores$TSS)
w_tss_me <- Resist_Model_scores[5,7]/sum(Resist_Model_scores$TSS)

Resist_Models$glm_weighted_auc_adj <- Resist_Models$glm_binary*w_auc_glm
Resist_Models$glm_auc_adj <- Resist_Models$glm_binary*Resist_Model_scores[1,4]
Resist_Models$glm_weighted_tss_adj <- Resist_Models$glm_binary*w_tss_glm
Resist_Models$glm_tss_adj <- Resist_Models$glm_binary*Resist_Model_scores[1,7]

Resist_Models$gam_weighted_auc_adj <- Resist_Models$gam_binary*w_auc_gam
Resist_Models$gam_auc_adj <- Resist_Models$gam_binary*Resist_Model_scores[1,4]
Resist_Models$gam_weighted_tss_adj <- Resist_Models$gam_binary*w_tss_gam
Resist_Models$gam_tss_adj <- Resist_Models$gam_binary*Resist_Model_scores[1,7]

Resist_Models$rf_weighted_auc_adj <- Resist_Models$rf_binary*w_auc_rf
Resist_Models$rf_auc_adj <- Resist_Models$rf_binary*Resist_Model_scores[1,4]
Resist_Models$rf_weighted_tss_adj <- Resist_Models$rf_binary*w_tss_rf
Resist_Models$rf_tss_adj <- Resist_Models$rf_binary*Resist_Model_scores[1,7]

Resist_Models$brt_weighted_auc_adj <- Resist_Models$brt_binary*w_auc_brt
Resist_Models$brt_auc_adj <- Resist_Models$brt_binary*Resist_Model_scores[1,4]
Resist_Models$brt_weighted_tss_adj <- Resist_Models$brt_binary*w_tss_brt
Resist_Models$brt_tss_adj <- Resist_Models$brt_binary*Resist_Model_scores[1,7]

Resist_Models$me_weighted_auc_adj <- Resist_Models$me_binary*w_auc_me
Resist_Models$me_auc_adj <- Resist_Models$me_binary*Resist_Model_scores[1,4]
Resist_Models$me_weighted_tss_adj <- Resist_Models$me_binary*w_tss_me
Resist_Models$me_tss_adj <- Resist_Models$me_binary*Resist_Model_scores[1,7]

Resist_Models$ensemble_w_auc_adj <- rowSums(Resist_Models[, c("glm_weighted_auc_adj", "gam_weighted_auc_adj", "rf_weighted_auc_adj", "brt_weighted_auc_adj", "me_weighted_auc_adj")], na.rm = TRUE)
Resist_Models$ensemble_auc_adj <- rowSums(Resist_Models[, c("glm_auc_adj", "gam_auc_adj", "rf_auc_adj", "brt_auc_adj", "me_auc_adj")], na.rm = TRUE)
Resist_Models$ensemble_tss_adj <- rowSums(Resist_Models[, c("glm_tss_adj", "gam_tss_adj", "rf_tss_adj", "brt_tss_adj", "me_tss_adj")], na.rm = TRUE)
Resist_Models$ensemble_w_tss_adj <- rowSums(Resist_Models[, c("glm_weighted_tss_adj", "gam_weighted_tss_adj", "rf_weighted_tss_adj", "brt_weighted_tss_adj", "me_weighted_tss_adj")], na.rm = TRUE)

#Susceptible JC Geno data
Suscept_Models <- read.csv("./Results/20250501_HabitSuitMod/GenoJCSuscept_ModelsOutputContinuous.csv", header=TRUE)
Suscept_Models <- Suscept_Models[,c("GEOID","glm_binary","gam_binary","rf_binary","brt_binary","me_binary")]
Suscept_Models_scores <- model_scores[model_scores$DataSet == "JCGenoSuscept",]
w_auc_glm <- Suscept_Models_scores[1,4]/sum(Suscept_Models_scores$AUC)
w_auc_gam <- Suscept_Models_scores[2,4]/sum(Suscept_Models_scores$AUC)
w_auc_rf <- Suscept_Models_scores[3,4]/sum(Suscept_Models_scores$AUC)
w_auc_brt <- Suscept_Models_scores[4,4]/sum(Suscept_Models_scores$AUC)
w_auc_me <- Suscept_Models_scores[5,4]/sum(Suscept_Models_scores$AUC)
w_tss_glm <- Suscept_Models_scores[1,7]/sum(Suscept_Models_scores$TSS)
w_tss_gam <- Suscept_Models_scores[2,7]/sum(Suscept_Models_scores$TSS)
w_tss_rf <- Suscept_Models_scores[3,7]/sum(Suscept_Models_scores$TSS)
w_tss_brt <- Suscept_Models_scores[4,7]/sum(Suscept_Models_scores$TSS)
w_tss_me <- Suscept_Models_scores[5,7]/sum(Suscept_Models_scores$TSS)

Suscept_Models$glm_weighted_auc_adj <- Suscept_Models$glm_binary*w_auc_glm
Suscept_Models$glm_auc_adj <- Suscept_Models$glm_binary*Suscept_Models_scores[1,4]
Suscept_Models$glm_weighted_tss_adj <- Suscept_Models$glm_binary*w_tss_glm
Suscept_Models$glm_tss_adj <- Suscept_Models$glm_binary*Suscept_Models_scores[1,7]

Suscept_Models$gam_weighted_auc_adj <- Suscept_Models$gam_binary*w_auc_gam
Suscept_Models$gam_auc_adj <- Suscept_Models$gam_binary*Suscept_Models_scores[1,4]
Suscept_Models$gam_weighted_tss_adj <- Suscept_Models$gam_binary*w_tss_gam
Suscept_Models$gam_tss_adj <- Suscept_Models$gam_binary*Suscept_Models_scores[1,7]

Suscept_Models$rf_weighted_auc_adj <- Suscept_Models$rf_binary*w_auc_rf
Suscept_Models$rf_auc_adj <- Suscept_Models$rf_binary*Suscept_Models_scores[1,4]
Suscept_Models$rf_weighted_tss_adj <- Suscept_Models$rf_binary*w_tss_rf
Suscept_Models$rf_tss_adj <- Suscept_Models$rf_binary*Suscept_Models_scores[1,7]

Suscept_Models$brt_weighted_auc_adj <- Suscept_Models$brt_binary*w_auc_brt
Suscept_Models$brt_auc_adj <- Suscept_Models$brt_binary*Suscept_Models_scores[1,4]
Suscept_Models$brt_weighted_tss_adj <- Suscept_Models$brt_binary*w_tss_brt
Suscept_Models$brt_tss_adj <- Suscept_Models$brt_binary*Suscept_Models_scores[1,7]

Suscept_Models$me_weighted_auc_adj <- Suscept_Models$me_binary*w_auc_me
Suscept_Models$me_auc_adj <- Suscept_Models$me_binary*Suscept_Models_scores[1,4]
Suscept_Models$me_weighted_tss_adj <- Suscept_Models$me_binary*w_tss_me
Suscept_Models$me_tss_adj <- Suscept_Models$me_binary*Suscept_Models_scores[1,7]

Suscept_Models$ensemble_w_auc_adj <- rowSums(Suscept_Models[, c("glm_weighted_auc_adj", "gam_weighted_auc_adj", "rf_weighted_auc_adj", "brt_weighted_auc_adj", "me_weighted_auc_adj")], na.rm = TRUE)
Suscept_Models$ensemble_auc_adj <- rowSums(Suscept_Models[, c("glm_auc_adj", "gam_auc_adj", "rf_auc_adj", "brt_auc_adj", "me_auc_adj")], na.rm = TRUE)
Suscept_Models$ensemble_w_tss_adj <- rowSums(Suscept_Models[, c("glm_weighted_tss_adj", "gam_weighted_tss_adj", "rf_weighted_tss_adj", "brt_weighted_tss_adj", "me_weighted_tss_adj")], na.rm = TRUE)
Suscept_Models$ensemble_tss_adj <- rowSums(Suscept_Models[, c("glm_tss_adj", "gam_tss_adj", "rf_tss_adj", "brt_tss_adj", "me_tss_adj")], na.rm = TRUE)


Suscept_Models$ensemble_w_auc_adj_inverted <- -1*Suscept_Models$ensemble_w_auc_adj
Suscept_Models$ensemble_auc_adj_inverted <- -1*Suscept_Models$ensemble_auc_adj
Suscept_Models$ensemble_w_tss_adj_inverted <- -1*Suscept_Models$ensemble_w_tss_adj
Suscept_Models$ensemble_tss_adj_inverted <- -1*Suscept_Models$ensemble_tss_adj
write.csv(Suscept_Models, "./output1.csv",row.names = FALSE)

summary(Suscept_Models$ensemble_w_auc_adj_inverted)
#ensamble_w_auc_adj has a healthy variety of numbers
"Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.0000  0.0000  0.0000  0.1716  0.1998  1.0000 
"


sub.data1 <- Suscept_Models[,c("GEOID","ensemble_w_auc_adj_inverted","ensemble_auc_adj_inverted","ensemble_w_tss_adj_inverted","ensemble_tss_adj_inverted")]
sub.data1[sub.data1$GEOID == "23005",]
Suscept_Model_Zeros <- sub.data1[sub.data1$ensemble_w_auc_adj_inverted == 0, c("GEOID")]
sub.data2 <- Resist_Models[,c("GEOID","ensemble_w_auc_adj","ensemble_auc_adj","ensemble_w_tss_adj","ensemble_tss_adj")]
sub.data2[sub.data2$GEOID == "23005",]
Resist_Model_Zeros <- sub.data2[sub.data2$ensemble_w_auc_adj == 0, c("GEOID")]
common_zeros <- intersect(Resist_Model_Zeros,Suscept_Model_Zeros)
data <- merge(sub.data1, sub.data2, by = "GEOID")
data$ensemble_w_auc_adj_merge <- (data$ensemble_w_auc_adj + data$ensemble_w_auc_adj_inverted)
data$ensemble_w_auc_adj_merge <- ifelse(data$GEOID %in% c(common_zeros), NA , data$ensemble_w_auc_adj_merge)
data$ensemble_auc_adj_merge <- (data$ensemble_auc_adj + data$ensemble_auc_adj_inverted)
data$ensemble_auc_adj_merge <- ifelse(data$GEOID %in% c(common_zeros), NA , data$ensemble_auc_adj_merge)
data$ensemble_tss_adj_merge <- (data$ensemble_tss_adj + data$ensemble_tss_adj_inverted)
data$ensemble_tss_adj_merge <- ifelse(data$GEOID %in% c(common_zeros), NA , data$ensemble_tss_adj_merge)
data$ensemble_w_tss_adj_merge <- (data$ensemble_w_tss_adj + data$ensemble_w_tss_adj_inverted)
data$ensemble_w_tss_adj_merge <- ifelse(data$GEOID %in% c(common_zeros), NA , data$ensemble_w_tss_adj_merge)

write.csv(data, "./Results/20250523_HabitSuitMod_2/GenoJC_RnS_ensemble.csv", row.names=FALSE)

#preprocssed data 
data <- read.csv("./Results/20250523_HabitSuitMod_2/GenoJC_RnS_ensemble.csv", header=TRUE)
data$GEOID <- sprintf("%05d", data$GEOID)
nrow(data[!is.na(data$ensemble_w_tss_adj_merge),])
nrow(data[!is.na(data$ensemble_w_tss_adj_merge) & data$ensemble_w_tss_adj_merge >= 0,])

data[data$GEOID == "09110",]


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


data_to_map <- merge(data,us_counties_shapefile_simpl, by = "GEOID" )

plot_1 <- ggplot() +
  geom_sf(data = data_to_map, aes(geometry = geometry, fill = ensemble_w_tss_adj), color = "NA", show.legend = TRUE) +
  geom_sf(data = us_states_shapefile_simpl, aes(geometry = geometry), fill = "NA", color = "dark grey", show.legend = TRUE) +
  theme_light() +
  labs(title = "Resistant Butternut Ensemble",
       subtitle = "Weighted TSS adjustment") +
  scale_fill_gradient(name = "",
                       #low = "#d8b365",
                       low = "white",
                       high = "#5ab4ac",
                       #midpoint = 0,
                       na.value = "grey",
                       labels = c("Unsuitable", "", "", "", "Suitable"))+#,
                       #breaks = c(min(data_to_map$ensemble_tss_adj, na.rm = TRUE),
                      #            quantile(data_to_map$ensemble_tss_adj, 0.25, na.rm = TRUE),
                      #            0,
                      #            quantile(data_to_map$ensemble_tss_adj, 0.75, na.rm = TRUE),
                      #            max(data_to_map$ensemble_tss_adj, na.rm = TRUE))) +
  theme(
    plot.title = element_text(size = 12),  
    plot.subtitle = element_text(size = 12),
    axis.title.x = element_text(size = 10),  
    axis.text.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.y = element_text(size = 10))

plot_3 <- ggplot() +
  geom_sf(data = data_to_map, aes(geometry = geometry, fill = (-1*ensemble_w_tss_adj_inverted)), color = "NA", show.legend = TRUE) +
  geom_sf(data = us_states_shapefile_simpl, aes(geometry = geometry), fill = "NA", color = "dark grey", show.legend = TRUE) +
  theme_light() +
  labs(title = "Susceptible Butternut Ensemble",
       subtitle = "Weighted TSS adjustment") +
  scale_fill_gradient(name = "",
                       low = "white",
                       #mid = "#f5f5f5",
                       high = "#d8b365",
                       #midpoint = 0,
                       na.value = "grey",
                       labels = c("Unsuitable", "", "", "", "Suitable"))+#,
                       #breaks = c(min(data_to_map$ensemble_w_tss_adj_merge, na.rm = TRUE),
                      #            quantile(data_to_map$ensemble_w_tss_adj_merge, 0.25, na.rm = TRUE),
                      #            0,
                      #            quantile(data_to_map$ensemble_w_tss_adj_merge, 0.75, na.rm = TRUE),
                      #            max(data_to_map$ensemble_w_tss_adj_merge, na.rm = TRUE))) +
  theme(
    plot.title = element_text(size = 12),  
    plot.subtitle = element_text(size = 12),
    axis.title.x = element_text(size = 10),  
    axis.text.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.y = element_text(size = 10))



file_name1 <- paste("./Results/20250523_HabitSuitMod_2/pngs/",
                    "HSM_Res_Tolerance_TSS_W8d_EnsembleModel",".png", sep = "")
file_name2 <- paste("./Results/20250523_HabitSuitMod_2/pngs/",
                    "HSM_Sucs_Tolerance_TSS_W8d_EnsembleModel",".png", sep = "")

file_name3 <- paste("./Results/20250523_HabitSuitMod_2/svgs/",
                    "HSM_Res_Tolerance_TSS_W8d_EnsembleModel",".svg", sep = "")
file_name4 <- paste("./Results/20250523_HabitSuitMod_2/svgs/",
                    "HSM_Sucs_Tolerance_TSS_W8d_EnsembleModel",".svg", sep = "")

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
as.data.frame(TreeSnap_Data_object[is.na(TreeSnap_Data_object$STATEFP),])

TreeSnap_Data_Counties <- as.data.frame(TreeSnap_Data_object[,c("GEOID")])
length(unique(TreeSnap_Data_Counties$GEOID))
TreeSnap_Data_Counties$geometry <- NULL
TreeSnap_Data_Counties_Mappable <- merge(TreeSnap_Data_Counties,us_counties_shapefile_simpl, by = "GEOID" )

data <- read.csv("./Results/20250523_HabitSuitMod_2/GenoJC_RnS_ensemble.csv", header=TRUE)
data$GEOID <- sprintf("%05d", data$GEOID)

data_treesnap_filtered <- data[data$GEOID %in% c(unique(TreeSnap_Data_Counties$GEOID)),]
nrow(data_treesnap_filtered)
nrow(data_treesnap_filtered[!is.na(data_treesnap_filtered$ensemble_w_tss_adj_merge),])

treesnap_mappable <- merge(data_treesnap_filtered, us_counties_shapefile, by = "GEOID")


data_3 <- data_treesnap_filtered[order(data_treesnap_filtered$ensemble_w_tss_adj_merge,decreasing = TRUE),
                                 c("GEOID","ensemble_w_tss_adj_merge")]
head(data_3)




#map filtered data
tree_snap <- ggplot() +
  geom_sf(data = us_counties_shapefile_simpl, aes(geometry = geometry), fill = "grey", color = "NA", show.legend = TRUE) +
  geom_sf(data = treesnap_mappable, aes(geometry = geometry,fill = ensemble_w_tss_adj_merge), color = "NA", show.legend = TRUE) +
  geom_sf(data = us_states_shapefile_simpl, aes(geometry = geometry), fill = "NA", color = "dark grey", show.legend = TRUE) +
  theme_light() +
  labs(title = "Resistant vs Susceptible Merged Ensemble \nTreesnap Filter",
       subtitle = "Weighted TSS adjustment") +
  scale_fill_gradient2(name = "",
                       low = "#d8b365",
                       mid = "#f5f5f5",
                       high = "#5ab4ac",
                       midpoint = 0,
                       na.value = "black",
                       labels = c("Susceptible", "", "", "", "Resistant"),
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
                    "HSM_GenoJC_RvS_TSS_W8d_EnsembleModel_Merge_TreeSnap",".png", sep = "")
file_name2 <- paste("./Results/20250523_HabitSuitMod_2/svgs/",
                    "HSM_GenoJC_RvS_TSS_W8d_EnsembleModel_Merge_TreeSnap",".svg", sep = "")
ggsave(file_name1, tree_snap, width = 7, height = 7, units = "in", dpi = 300)
ggsave(file_name2, tree_snap, width = 7, height = 7, units = "in", dpi = 300)

data_treesnap_filtered[order(data_treesnap_filtered$ensemble_w_tss_adj_merge,decreasing = TRUE), c("GEOID","ensemble_w_tss_adj_merge")]


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
  labs(title = "Resistant vs Susceptible Merged Ensemble \niNaturalist Filter",
       subtitle = "Weighted TSS adjustment") +
  scale_fill_gradient2(name = "",
                       low = "#d8b365",
                       mid = "#f5f5f5",
                       high = "#5ab4ac",
                       midpoint = 0,
                       na.value = "black",
                       labels = c("Susceptible", "", "", "", "Resistant"),
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
                    "HSM_GenoJC_RvS_TSS_W8d_EnsembleModel_Merge_iNaturalist",".png", sep = "")
file_name2 <- paste("./Results/20250523_HabitSuitMod_2/svgs/",
                    "HSM_GenoJC_RvS_TSS_W8d_EnsembleModel_Merge_iNaturalist",".svg", sep = "")
ggsave(file_name1, iNaturalist, width = 7, height = 7, units = "in", dpi = 300)
ggsave(file_name2, iNaturalist, width = 7, height = 7, units = "in", dpi = 300)

#GBIF filter####
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



data <- read.csv("./Results/20250523_HabitSuitMod_2/GenoJC_RnS_ensemble.csv", header=TRUE)
data$GEOID <- sprintf("%05d", data$GEOID)

setdiff(GBIF_Data_Counties$GEOID, data$GEOID)

nrow(data[!is.na(data$ensemble_w_tss_adj_merge),])
nrow(data[data$GEOID %in% c(unique(GBIF_Data_Counties$GEOID)),])


data_gbif_filtered <- data[data$GEOID %in% c(GBIF_Data_Counties$GEOID),]
nrow(data_gbif_filtered)

nrow(data_gbif_filtered[!is.na(data_gbif_filtered$ensemble_w_tss_adj_merge),])
nrow(data_gbif_filtered[is.na(data_gbif_filtered$ensemble_w_tss_adj_merge),])


data_3 <- data_gbif_filtered[order(data_gbif_filtered$ensemble_w_tss_adj_merge,decreasing = TRUE),
                             c("GEOID","ensemble_w_tss_adj_merge")]
head(data_3)


gbif_mappable <- merge(data_gbif_filtered, us_counties_shapefile, by = "GEOID")


gbif <- ggplot() +
  geom_sf(data = us_counties_shapefile_simpl, aes(geometry = geometry), fill = "grey", color = "NA", show.legend = TRUE) +
  geom_sf(data = gbif_mappable, aes(geometry = geometry,fill = ensemble_w_tss_adj_merge), color = "NA", show.legend = TRUE) +
  geom_sf(data = us_states_shapefile_simpl, aes(geometry = geometry), fill = "NA", color = "dark grey", show.legend = TRUE) +
  theme_light() +
  labs(title = "Resistant vs Susceptible Merged Ensemble \nGBIF Filter",
       subtitle = "Weighted TSS adjustment") +
  scale_fill_gradient2(name = "",
                       low = "#d8b365",
                       mid = "#f5f5f5",
                       high = "#5ab4ac",
                       midpoint = 0,
                       na.value = "black",
                       labels = c("Susceptible", "", "", "", "Resistant"),
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
                    "HSM_GenoJC_RvS_TSS_W8d_EnsembleModel_Merge_GBIF",".png", sep = "")
file_name2 <- paste("./Results/20250523_HabitSuitMod_2/svgs/",
                    "HSM_GenoJC_RvS_TSS_W8d_EnsembleModel_Merge_GBIF",".svg", sep = "")
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

data <- read.csv("./Results/20250523_HabitSuitMod_2/GenoJC_RnS_ensemble.csv", header=TRUE)
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
  labs(title = "Final Tolerance Ensemble \nFIA Filter",
       subtitle = "Weighted TSS adjustment") +
  scale_fill_gradient2(name = "",
                       low = "#d8b365",
                       mid = "#f5f5f5",
                       high = "#5ab4ac",
                       na.value = "black",
                       labels = c("Susceptible", "", "", "", "Resistant"),
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
                    "HSM_GenoJC_RvS_TSS_W8d_EnsembleModel_Merge_FIA",".png", sep = "")
file_name2 <- paste("./Results/20250523_HabitSuitMod_2/svgs/",
                    "HSM_GenoJC_RvS_TSS_W8d_EnsembleModel_Merge_FIA",".svg", sep = "")
ggsave(file_name1, fia, width = 7, height = 7, units = "in", dpi = 300)
ggsave(file_name2, fia, width = 7, height = 7, units = "in", dpi = 300)
