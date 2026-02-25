##Adapted from https://rpubs.com/JoshCarrell/HSM
#install.packages("flexsdm")
library(sf)
library(dplyr)
library(tmap)
library(usmap)
library(ggplot2)
library(pracma)
library(randomForest)
library(caret)
library(mgcv)
library(maps)
library(PresenceAbsence)
library(DAAG)
library(corrplot)
library(terra)
library(pROC)
library(raster)
library(tigris)
library(pander)
library(gridExtra)
library(gbm)
library(dismo)
library(maxnet)
library(matrixStats)
setwd("/Volumes/Onni_Drive/******_Projects/Butternut_Landscape_Genomics")

Location_Data <- read.csv("./Results/20240819_AccessionLocations/bnut_families_wBV_sp.csv", header=TRUE)
#Subset family data to just butternut families that are more resistant (resistant) to BCD
Location_Data_Resistant_JC <- Location_Data[Location_Data$species == "JC" 
                                            & Location_Data$bays_bcd_bv >= 0, c("bays_bcd_bv","fips")]
length(unique(Location_Data_Resistant_JC[!is.na(Location_Data_Resistant_JC$bays_bcd_bv), c("fips")]))
Location_Data_Resistant_JC$GEOID <- as.character(Location_Data_Resistant_JC$fips)
Location_Data_Resistant_JC$fips <- NULL
Location_Data_Resistant_JC[Location_Data$GEOID == "17085",]
Location_Data_Resistant_JC <- na.omit(Location_Data_Resistant_JC)
length(unique(Location_Data_Resistant_JC$GEOID))


VA_City_GEOIDS <- as.character(c(51510,51520,51530,51540,51550,51570,51580,51590,51595,51620,51630,51640,51650,
                                 51660,51670,51680,51683,51690,51700,51710,51720,51730,51735,51740,51750,51760,
                                 51770,51775,51790,51800,51810,51820,51830,51840))


PRISM_Data <- read.csv("./Results/20250415_Accessions_Climate/AllCounties_climate_OSC_elevation_v3.csv", header=TRUE)
PRISM_Data$GEOID <-sprintf("%05d", PRISM_Data$GEOID)
PRISM_Data <- PRISM_Data[ , !grepl("^sd_", names(PRISM_Data))]

us_counties_shapefile <- counties(cb = TRUE)
us_counties_shapefile <- us_counties_shapefile[us_counties_shapefile$STUSPS %in% c("ME","NH","VT","RI","CT","MA","NY","PA",
                                                                                   "DE","WV","NJ","MD","VA","OH","SC","NC",
                                                                                   "TN","KY","GA","FL","AL","MS","IN","IL",
                                                                                   "MI","MN","IA","AR","MO","LA","WI"),]
us_counties_shapefile_simpl <- st_simplify(us_counties_shapefile, dTolerance = 1000)

####move to actually making the HSM
##Generate pseduo-absense data###

set.seed(6644)

df <- data.frame(
  GEOID = unique(Location_Data_Resistant_JC$GEOID),
  ind = sample(2, length(unique(Location_Data_Resistant_JC$GEOID)), replace = TRUE, prob = c(0.7, 0.3))
)
nrow(df[df$ind == 2,])


occurance_geoids <- df[df$ind == 1, c("GEOID")]
unavailable_GEOIDS <- c(VA_City_GEOIDS, as.character(Location_Data_Resistant_JC$GEOID))
available_GEOIDS <- PRISM_Data[!PRISM_Data$GEOID %in% unavailable_GEOIDS, "GEOID"]
train.absence_GEOIDS <- sample(available_GEOIDS, 
                              size = nrow(df[df$ind == 1,])*3)


unavailable_GEOIDS <- c(train.absence_GEOIDS, VA_City_GEOIDS, as.character(Location_Data_Resistant_JC$GEOID))
available_GEOIDS <- PRISM_Data[!PRISM_Data$GEOID %in% unavailable_GEOIDS, "GEOID"]
test.absence_GEOIDS <- sample(available_GEOIDS, 
                              size = nrow(df[df$ind == 2,]))

test.occurance <- data.frame(
  GEOID = as.character(df[df$ind == 2, c("GEOID")]),
  PA = 1)
test.absence <- data.frame(
  GEOID = test.absence_GEOIDS,
  Resistance = 0.002,
  PA = 0
)
nrow(test.occurance)


test.occurance <- test.occurance %>%
  left_join(Location_Data_Resistant_JC %>%
              group_by(GEOID) %>%
              summarise(bays_bcd_bv = mean(bays_bcd_bv, na.rm = TRUE)),
            by = "GEOID")
colnames(test.occurance) <- c("GEOID","PA","Resistance")

test.location <- bind_rows(test.occurance,test.absence)
nrow(test.location[test.location$PA == "0",])
nrow(test.location)


train.occurance <- data.frame(
  GEOID = occurance_geoids,
  PA = 1)

train.occurance <- train.occurance %>%
  left_join(Location_Data_Resistant_JC %>%
              group_by(GEOID) %>%
              summarise(bays_bcd_bv = mean(bays_bcd_bv, na.rm = TRUE)),
            by = "GEOID")
colnames(train.occurance) <- c("GEOID","PA","Resistance")


train.absence <- data.frame(
  GEOID = train.absence_GEOIDS,
  PA = 0,
  Resistance = 0.002)

train.occurance$GEOID <- as.character(train.occurance$GEOID)
train.absence$GEOID <- as.character(train.absence$GEOID)
train.location <- bind_rows(train.occurance, train.absence)
nrow(train.location[train.location$PA == "0",])
nrow(train.location)

#Check to make sure the location makes sense
train.data_to_map <- merge(train.location,us_counties_shapefile_simpl, by = "GEOID")
test.data_to_map <- merge(test.location,us_counties_shapefile_simpl, by = "GEOID")

ggplot() +
  geom_sf(data = us_counties_shapefile_simpl, aes(geometry = geometry), fill = "grey", color = "NA", show.legend = TRUE) +
  geom_sf(data = train.data_to_map, aes(geometry = geometry, fill = as.factor(PA)), color = "NA", show.legend = TRUE) +
  theme_light() +
  scale_fill_manual(name = "Presense/\nAbsense", values = c("0" = "deepskyblue", "1" = "darkgreen", "NA"="grey"))+
  theme(
    plot.title = element_text(size = 12),  
    plot.subtitle = element_text(size = 12),
    axis.title.x = element_text(size = 10),  
    axis.text.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.y = element_text(size = 10)) 

ggplot() +
  geom_sf(data = us_counties_shapefile_simpl, aes(geometry = geometry), fill = "grey", color = "NA", show.legend = TRUE) +
  geom_sf(data = test.data_to_map, aes(geometry = geometry), fill = "firebrick", color = "NA", show.legend = TRUE) +
  theme_light() +
  theme(
    plot.title = element_text(size = 12),  
    plot.subtitle = element_text(size = 12),
    axis.title.x = element_text(size = 10),  
    axis.text.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10),  
    axis.text.y = element_text(size = 10)) 

predictor_table <- merge(train.location, PRISM_Data, by = "GEOID", all.x=TRUE)

str(predictor_table)
#Check for correlations among predictors
spec.cor.topo <- predictor_table[3:134] # subset dataframe to show just extracted predictor info
topo.cor <- cor(spec.cor.topo) # correlation
high.cor.pairs <- which(abs(topo.cor) > 0.7 & abs(topo.cor) < 1, arr.ind = TRUE)
vars.to.remove <- unique(rownames(high.cor.pairs)[high.cor.pairs[,1] < high.cor.pairs[,2]])
topo.reduced <- spec.cor.topo[ , !(names(spec.cor.topo) %in% vars.to.remove)]
topo.reduced.cor <- cor(topo.reduced)
corrplot.mixed(topo.reduced.cor, lower = 'number', upper = 'square')
corrplot.mixed(topo.reduced.cor)
rownames(topo.reduced.cor)
write.csv(topo.cor, "./Results/20250501_HabitSuitMod/PRISM_Data_GenoJCResist_Correlations_v1.csv", row.names=TRUE)

#Reiterate with potentially non-correlary predictors
spec.cor.topo_v2 <- predictor_table[,c("Resistance",
                                       "mean_norm_ppt_05",
                                       'mean_norm_ppt_06',
                                       'mean_norm_ppt_09',
                                       "mean_norm_ppt_10",
                                       'mean_norm_ppt_annual',
                                       'mean_norm_solclear_05',
                                       'mean_norm_solclear_07',
                                       'mean_norm_soltotal_03',
                                       "mean_norm_soltotal_06",
                                       'mean_norm_soltotal_07',
                                       'mean_norm_soltrans_02',
                                       'mean_norm_soltrans_09',
                                       'mean_norm_vpdmax_08',
                                       'mean_norm_vpdmin_annual',
                                       "osc_mtu_per_ha")]
topo.cor_v2 <- cor(spec.cor.topo_v2) # correlation
corrplot.mixed(topo.cor_v2, lower = 'number', upper = 'square')
write.csv(topo.cor_v2, "./Results/20250501_HabitSuitMod/PRISM_Data_GenoJCResist_Correlations_v2.csv", row.names=TRUE)

predictor_table_final <- predictor_table[,c("GEOID","PA",
                                            "Resistance",
                                            "mean_norm_ppt_05",
                                            'mean_norm_ppt_06',
                                            'mean_norm_ppt_09',
                                            "mean_norm_ppt_10",
                                            'mean_norm_ppt_annual',
                                            'mean_norm_solclear_05',
                                            'mean_norm_solclear_07',
                                            'mean_norm_soltotal_03',
                                            "mean_norm_soltotal_06",
                                            'mean_norm_soltotal_07',
                                            'mean_norm_soltrans_02',
                                            'mean_norm_soltrans_09',
                                            'mean_norm_vpdmax_08',
                                            'mean_norm_vpdmin_annual',
                                            "osc_mtu_per_ha")]
write.csv(predictor_table_final, "./output.csv", row.names = FALSE)
topo.cor_v2 <- cor(predictor_table_final[3:14]) # correlation
corrplot.mixed(topo.cor_v2)

write.csv(predictor_table_final, "./Results/20250501_HabitSuitMod/PRISM_Data_GenoJCResist_PredictorTableFinal_v1.csv", row.names=FALSE)
#Bring up final predictor table, if wanted
predictor_table_final <- read.csv("./Results/20250501_HabitSuitMod/PRISM_Data_GenoJCResist_PredictorTableFinal_v1.csv", header=TRUE)
nrow(predictor_table_final)
names(predictor_table_final)
#Generalized Linear Model####
mod.glm.1 <- glm(as.factor(PA) ~ 
                 mean_norm_ppt_05 +
                 mean_norm_ppt_06 +
                 mean_norm_ppt_09 +
                 mean_norm_ppt_10 +
                 mean_norm_ppt_annual +
                 mean_norm_solclear_05 +
                 mean_norm_solclear_07 +
                 mean_norm_soltotal_03 +
                 mean_norm_soltotal_06 +
                 mean_norm_soltotal_07 +
                 mean_norm_soltrans_02 +
                 mean_norm_soltrans_09 +
                 mean_norm_vpdmax_08 +
                 mean_norm_vpdmin_annual +
                 osc_mtu_per_ha,
                 family = binomial, data = predictor_table_final)
summary(mod.glm.1)
mod.glm.2 <- glm(as.factor(PA) ~
                   mean_norm_ppt_05 +
                   mean_norm_soltotal_07 +
                   mean_norm_soltrans_02 +
                   mean_norm_vpdmax_08,
                 family = binomial, data = predictor_table_final)
summary(mod.glm.2)

mod.glm.3 <- glm(as.factor(PA) ~ 
                   mean_norm_ppt_10 +
                   mean_norm_ppt_annual +
                   mean_norm_solclear_05 +
                   mean_norm_soltotal_07 +
                   mean_norm_vpdmax_08 +
                   mean_norm_vpdmin_annual,
                 family = binomial, data = predictor_table_final)
summary(mod.glm.3)
AIC(mod.glm.1,mod.glm.2,mod.glm.3)
mod1.fit <- 100*(1-mod.glm.1$deviance/mod.glm.1$null.deviance)
mod2.fit <- 100*(1-mod.glm.2$deviance/mod.glm.2$null.deviance)
mod3.fit <- 100*(1-mod.glm.3$deviance/mod.glm.3$null.deviance)
print(mod1.fit)
print(mod2.fit)
print(mod3.fit)

mod.glm.final <- mod.glm.1

pred_probs <- predict(mod.glm.final, type = "response")
roc_curve <- roc(predictor_table_final$PA, pred_probs)
auc(roc_curve) 

##Cross validation of GLM####
perform_glm_cv_caret <- function(predictor_table_final, k = 5) {
  train_control <- trainControl(
    method = "cv",
    number = k,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final"
  )
  
  data_caret <- predictor_table_final
  data_caret$PA <- factor(data_caret$PA, levels = c(0, 1), labels = c("Absent", "Present"))
  
  cv_model <- train(
    PA ~ mean_norm_ppt_10 + mean_norm_ppt_annual + mean_norm_solclear_05 + 
      mean_norm_soltotal_07 + mean_norm_soltrans_08 + mean_norm_soltrans_annual,
    data = data_caret,
    method = "glm",
    family = "binomial",
    trControl = train_control,
    metric = "ROC"
  )
  
  return(cv_model)
}

cv_results_caret <- perform_glm_cv_caret(train.location, k = 5)
print(cv_results_caret)

#Generalized Additive Model####
mod.gam.1 <- gam(as.factor(PA) ~ 
                   s(mean_norm_ppt_05) +
                   s(mean_norm_ppt_06) +
                   s(mean_norm_ppt_09) +
                   s(mean_norm_ppt_10) +
                   s(mean_norm_ppt_annual) +
                   s(mean_norm_solclear_05) +
                   s(mean_norm_solclear_07) +
                   s(mean_norm_soltotal_03) +
                   s(mean_norm_soltotal_06) +
                   s(mean_norm_soltotal_07) +
                   s(mean_norm_soltrans_02) +
                   s(mean_norm_soltrans_09) +
                   s(mean_norm_vpdmax_08) +
                   s(mean_norm_vpdmin_annual) +
                   s(osc_mtu_per_ha),
                  family = binomial, 
                  select = TRUE,
                  data = predictor_table_final)
summary(mod.gam.1)
mod.gam.2 <- gam(as.factor(PA) ~ 
                   s(mean_norm_ppt_05, k = 3) +
                   s(mean_norm_ppt_06, k = 3) +
                   s(mean_norm_ppt_09, k = 3) +
                   s(mean_norm_ppt_10, k = 3) +
                   s(mean_norm_ppt_annual, k = 3) +
                   s(mean_norm_solclear_05, k = 3) +
                   s(mean_norm_solclear_07, k = 3) +
                   s(mean_norm_soltotal_03, k = 3) +
                   s(mean_norm_soltotal_06, k = 3) +
                   s(mean_norm_soltotal_07, k = 3) +
                   s(mean_norm_soltrans_02, k = 3) +
                   s(mean_norm_soltrans_09, k = 3) +
                   s(mean_norm_vpdmax_08, k = 3) +
                   s(mean_norm_vpdmin_annual, k = 3) +
                   s(osc_mtu_per_ha, k = 3),
                 family = binomial, 
                 select = TRUE,
                 data = predictor_table_final)
summary(mod.gam.2)
mod.gam.3 <- gam(as.factor(PA) ~ 
                   #s(mean_norm_ppt_05, k = 3) +
                   #s(mean_norm_ppt_06, k = 3) +
                   s(mean_norm_ppt_09, k = 3) +
                   #s(mean_norm_ppt_10, k = 3) +
                   #s(mean_norm_ppt_annual, k = 3) +
                   #s(mean_norm_solclear_05, k = 3) +
                   s(mean_norm_solclear_07, k = 3) +
                   s(mean_norm_soltotal_03, k = 3) +
                   #s(mean_norm_soltotal_06, k = 3) +
                   #s(mean_norm_soltotal_07, k = 3) +
                   s(mean_norm_soltrans_02, k = 3) +
                   #s(mean_norm_soltrans_09, k = 3) +
                   s(mean_norm_vpdmax_08, k = 3),# +
                   #s(mean_norm_vpdmin_annual, k = 3) +
                   #s(osc_mtu_per_ha, k = 3),
                 family = binomial, 
                 select = TRUE,
                 data = predictor_table_final)
gam.check(mod.gam.3)
summary(mod.gam.3)

AIC(mod.gam.1, mod.gam.2, mod.gam.3)
mod1.fit <- 100*(1-mod.gam.1$deviance/mod.gam.1$null.deviance)
mod2.fit <- 100*(1-mod.gam.2$deviance/mod.gam.2$null.deviance)
mod3.fit <- 100*(1-mod.gam.3$deviance/mod.gam.3$null.deviance)
print(mod1.fit)
print(mod2.fit)
print(mod3.fit)
#model 3 is the winner.

mod.gam.final <- mod.gam.3

#Random Forest####
set.seed(9242)
mod.rf.1 <- randomForest(as.factor(PA) ~ 
                           mean_norm_ppt_05 +
                           mean_norm_ppt_06 +
                           mean_norm_ppt_09 +
                           mean_norm_ppt_10 +
                           mean_norm_ppt_annual +
                           mean_norm_solclear_05 +
                           mean_norm_solclear_07 +
                           mean_norm_soltotal_03 +
                           mean_norm_soltotal_06 +
                           mean_norm_soltotal_07 +
                           mean_norm_soltrans_02 +
                           mean_norm_soltrans_09 +
                           mean_norm_vpdmax_08 +
                           mean_norm_vpdmin_annual +
                           osc_mtu_per_ha,
                           data = predictor_table_final, 
                           sampsize = c('0'=20,'1'=20), #RF benefits from even sampling (allegedly)
                           ntree = 1000,
                           mtry = 4,
                           importance = TRUE,
                           proximity=TRUE) 
mod.rf.1
randomForest::importance(mod.rf.1)
varImpPlot(mod.rf.1)


mod.rf.2 <- randomForest(as.factor(PA) ~ 
                           mean_norm_soltrans_02 +
                           mean_norm_soltotal_03 +
                           mean_norm_ppt_09,
                         data = predictor_table_final, 
                         sampsize = c('0'=20,'1'=20), #RF benefits from even sampling (allegedly)
                         ntree = 1000,
                         mtry = 3,
                         importance = TRUE,
                         proximity=TRUE) 
mod.rf.2
importance(mod.rf.2)
varImpPlot(mod.rf.2)

mod.rf.3 <-  randomForest(as.factor(PA) ~ 
                            mean_norm_ppt_05 +
                            #mean_norm_ppt_06 +
                            #mean_norm_ppt_08 +
                            #mean_norm_ppt_10 +
                            #mean_norm_ppt_annual +
                            #mean_norm_solclear_05 +
                            #mean_norm_solclear_07 +
                            #mean_norm_soltotal_06 +
                            #mean_norm_soltotal_07 +
                            #mean_norm_soltrans_08 +
                            mean_norm_soltrans_annual,# +
                          #mean_norm_vpdmin_annual +
                          #osc_mtu_per_ha,
                          data = predictor_table_final, 
                          sampsize = c('0'=30,'1'=30), #RF benefits from even sampling (allegedly)
                          ntree = 1000, 
                          mtry = 2,
                          importance = TRUE,
                          proximity=TRUE) 
mod.rf.3
importance(mod.rf.3)
varImpPlot(mod.rf.3)

###Code to systematically reduce OOB
var_imp <- importance(mod.rf.1)[, "MeanDecreaseAccuracy"]
varImpPlot(mod.rf.1, sort = TRUE, main = "Variable Importance")
threshold <- quantile(var_imp, 0.7)  # You can adjust this (0.7 = top 30%)
top_vars <- names(var_imp[var_imp > threshold])
# Print selected variables
cat("Top predictors:\n")
print(top_vars)

# Create formula dynamically
rf_formula <- as.formula(paste("as.factor(PA) ~", paste(top_vars, collapse = " + ")))
set.seed(9242)
mod.rf.4 <- randomForest(
  formula = rf_formula,
  data = predictor_table_final,
  ntree = 1000,
  sampsize = c('0'=30,'1'=30), #RF benefits from even sampling (allegedly)
  importance = TRUE,
  proximity = TRUE
)
mod.rf.4
importance(mod.rf.4)
varImpPlot(mod.rf.4)


pred_probs <- predict(mod.rf.4, type = "prob")[,2]
roc_curve <- roc(predictor_table_final$PA, pred_probs)
auc(roc_curve) 


mod.rf.final <- mod.rf.1

#Boosted Regression Tree (BRT) Model ####
mod.BRT.1_weighted <- gbm(formula = PA ~
                            mean_norm_ppt_05 +
                            mean_norm_ppt_06 +
                            mean_norm_ppt_09 +
                            mean_norm_ppt_10 +
                            mean_norm_ppt_annual +
                            mean_norm_solclear_05 +
                            mean_norm_solclear_07 +
                            mean_norm_soltotal_03 +
                            mean_norm_soltotal_06 +
                            mean_norm_soltotal_07 +
                            mean_norm_soltrans_02 +
                            mean_norm_soltrans_09 +
                            mean_norm_vpdmax_08 +
                            mean_norm_vpdmin_annual +
                            osc_mtu_per_ha,
                        weights = predictor_table_final$Resistance,
                        data=predictor_table_final,
                        distribution = "bernoulli",
                        n.trees = 3000,
                        cv.folds = 10,
                        interaction.depth = 4,
                        shrinkage = 0.001,
                        n.minobsinnode = 10)
summary(mod.BRT.1_weighted)
best_iter <- gbm.perf(mod.BRT.1_weighted, method = "cv")
print(best_iter)

mod.BRT.2_weighted <- gbm(formula = PA ~
                            mean_norm_ppt_05 +
                            mean_norm_ppt_06 +
                            mean_norm_ppt_09 +
                            mean_norm_ppt_10 +
                            mean_norm_ppt_annual +
                            mean_norm_solclear_05 +
                            mean_norm_solclear_07 +
                            mean_norm_soltotal_03 +
                            mean_norm_soltotal_06 +
                            mean_norm_soltotal_07 +
                            mean_norm_soltrans_02 +
                            mean_norm_soltrans_09 +
                            mean_norm_vpdmax_08 +
                            mean_norm_vpdmin_annual +
                            osc_mtu_per_ha,
                          weights = predictor_table_final$Resistance,
                          data=predictor_table_final,
                          distribution = "bernoulli",
                          n.trees = best_iter,
                          cv.folds = 10,
                          interaction.depth = 4,
                          shrinkage = 0.001,
                          n.minobsinnode = 10)
summary(mod.BRT.2_weighted)
best_iter <- gbm.perf(mod.BRT.2_weighted, method = "cv")
print(best_iter)

gbm.perf(mod.BRT.2_weighted, method = "cv")

mod.BRT.3_weighted <- gbm(formula = PA ~
                            mean_norm_ppt_05 +
                            mean_norm_ppt_06 +
                            #mean_norm_ppt_08 +
                            mean_norm_ppt_10 +
                            mean_norm_ppt_annual +
                            #mean_norm_solclear_05 +
                            mean_norm_solclear_07 +
                            mean_norm_soltotal_06 +
                            #mean_norm_soltotal_07 +
                            mean_norm_soltrans_annual +
                            mean_norm_vpdmin_annual +
                            osc_mtu_per_ha,
                          weights = predictor_table_final$Resistance,
                          data=predictor_table_final,
                          distribution = "bernoulli",
                          n.trees = 5000,
                          cv.folds = 10,
                          interaction.depth = 4,
                          shrinkage = 0.001,
                          n.minobsinnode = 10)
summary(mod.BRT.3_weighted)
gbm.perf(mod.BRT.3_weighted, method = "cv")

pred_probs <- predict(mod.BRT.2_weighted, type = "response")
roc_curve <- roc(predictor_table_final$PA, pred_probs)
auc(roc_curve) 
coords_best <- coords(roc_curve, "best", ret = "threshold", transpose = FALSE)
best_thresh <- coords_best$threshold
pred_class <- ifelse(pred_probs >= best_thresh, 1, 0)
table(Predicted = pred_class, Actual = predictor_table_final$PA)

# Number 2 has the best AUC score. Over fitted? They're all 'great'

mod.brt.final <- mod.BRT.2_weighted


#Maximum Entropy (MaxEnt) Model ####
maxent_predictors_1 <- predictor_table_final[, c(
  "mean_norm_ppt_05",
  'mean_norm_ppt_06',
  'mean_norm_ppt_09',
  "mean_norm_ppt_10",
  'mean_norm_ppt_annual',
  'mean_norm_solclear_05',
  'mean_norm_solclear_07',
  'mean_norm_soltotal_03',
  "mean_norm_soltotal_06",
  'mean_norm_soltotal_07',
  'mean_norm_soltrans_02',
  'mean_norm_soltrans_09',
  'mean_norm_vpdmax_08',
  'mean_norm_vpdmin_annual',
  "osc_mtu_per_ha"
)]
#Run MaxEnt
mod.me.1 <- maxnet(p = predictor_table_final$PA, data = maxent_predictors_1)
coefficients <- mod.me.1$betas
var_contributions <- abs(coefficients)
var_contrib_percent <- (var_contributions / sum(var_contributions)) * 100

contrib_summary <- data.frame(
  Variable = names(var_contrib_percent),
  Coefficient = coefficients,
  Contribution_Percent = as.numeric(var_contrib_percent)
) %>%
  arrange(desc(Contribution_Percent))

print(contrib_summary)

extract_base_varname <- function(feature_name) {
  return(feature_name)
}
contrib_by_var <- contrib_summary %>%
  mutate(Base_Variable = sapply(Variable, extract_base_varname)) %>%
  group_by(Base_Variable) %>%
  summarise(
    Total_Contribution = sum(Contribution_Percent),
    N_Features = n(),
    .groups = 'drop'
  ) %>%
  arrange(desc(Total_Contribution))

print(as.data.frame(contrib_by_var))


roc_obj <- roc(
  predictor_table_final$PA,
  predict(mod.me.1, newdata = maxent_predictors_1, type = "cloglog")
)
auc(roc_obj)
plot(roc_obj)
pred_class <- ifelse(pred_probs >= 0.4, 1, 0)
coords_best <- coords(roc_obj, "best", ret = "threshold", transpose = FALSE)
best_thresh <- coords_best$threshold
pred_class <- ifelse(pred_probs >= best_thresh, 1, 0)
table(Predicted = pred_class, Actual = predictor_table_final$PA)


maxent_predictors_2 <- predictor_table_final[, c(
  "mean_norm_ppt_05",
  #"mean_norm_ppt_06",
  "mean_norm_ppt_08",
  #"mean_norm_ppt_10",
  #"mean_norm_ppt_annual",
  "mean_norm_solclear_05",
  "mean_norm_solclear_07",
  #"mean_norm_soltotal_06",
  "mean_norm_soltotal_07",
  "mean_norm_soltrans_08",
  "mean_norm_soltrans_annual",
  "mean_norm_vpdmin_annual"
  #"osc_mtu_per_ha"
)]
mod.ME.2 <- maxnet(p = predictor_table_final$PA, data = maxent_predictors_2)


roc_obj <- roc(
  predictor_table_final$PA,
  predict(mod.ME.2, newdata = maxent_predictors_2, type = "cloglog")
)
auc(roc_obj)
plot(roc_obj)

mod.me.final <- mod.me.1

#Evaluate
test_predictor_table <- merge(test.location, PRISM_Data, by = "GEOID", all.x=TRUE)
test_predictor_table <- test_predictor_table[,c("GEOID","PA", "mean_norm_ppt_05" , 
                                                "mean_norm_ppt_05",
                                                "mean_norm_ppt_06",
                                                "mean_norm_ppt_08",
                                                "mean_norm_ppt_10",
                                                "mean_norm_ppt_annual",
                                                "mean_norm_solclear_05",
                                                "mean_norm_solclear_07",
                                                "mean_norm_soltotal_06",
                                                "mean_norm_soltotal_07",
                                                "mean_norm_soltrans_08",
                                                "mean_norm_soltrans_annual",
                                                "mean_norm_vpdmin_annual",
                                                "osc_mtu_per_ha")]

test_predictors <- test_predictor_table[,c(         "mean_norm_ppt_05",
                                                    "mean_norm_ppt_06",
                                                    "mean_norm_ppt_08",
                                                    "mean_norm_ppt_10",
                                                    "mean_norm_ppt_annual",
                                                    "mean_norm_solclear_05",
                                                    "mean_norm_solclear_07",
                                                    "mean_norm_soltotal_06",
                                                    "mean_norm_soltotal_07",
                                                    "mean_norm_soltrans_08",
                                                    "mean_norm_soltrans_annual",
                                                    "mean_norm_vpdmin_annual",
                                                    "osc_mtu_per_ha")]
test_PA <- test_predictor_table$PA
predicted <- predict(mod.ME.1, test_predictors, type = "cloglog")
eval_result <- dismo::evaluate(p = predicted[test_PA == 1],
                        a = predicted[test_PA == 0])

print(eval_result)


"-----------------------------------------------------------------------------------------------------"
#Use models to predict across whole study area####
"-----------------------------------------------------------------------------------------------------"
test_predictor_table <- merge(test.location, PRISM_Data, by = "GEOID", all.x=TRUE)
write.csv(test_predictor_table, "./Results/20250501_HabitSuitMod/PRISM_Data_GenoJCResist_TestDataPredictorTable.csv", row.names=FALSE)

test_predictor_table <- test_predictor_table[,c("GEOID","PA","Resistance",
                                                "mean_norm_ppt_05",
                                                'mean_norm_ppt_06',
                                                'mean_norm_ppt_09',
                                                "mean_norm_ppt_10",
                                                'mean_norm_ppt_annual',
                                                'mean_norm_solclear_05',
                                                'mean_norm_solclear_07',
                                                'mean_norm_soltotal_03',
                                                "mean_norm_soltotal_06",
                                                'mean_norm_soltotal_07',
                                                'mean_norm_soltrans_02',
                                                'mean_norm_soltrans_09',
                                                'mean_norm_vpdmax_08',
                                                'mean_norm_vpdmin_annual',
                                                "osc_mtu_per_ha")]



test_predictor_table$mod.glm.predict <- 100*predict(mod.glm.final, newdata = test_predictor_table, type = "response") 
test_predictor_table$mod.gam.predict <- 100*predict(mod.gam.final, newdata = test_predictor_table, type = "response") 
test_predictor_table$mod.rf.predict  <- 100*predict(mod.rf.final, newdata = test_predictor_table, type = "prob")[,2] 
test_predictor_table$mod.brt.predict <- 100*predict(mod.brt.final, newdata = test_predictor_table, type = "response",  n.trees = best_iter)

predictor_vars <- colnames(maxent_predictors_1)
test_data_subset <- test_predictor_table[, predictor_vars, drop = FALSE]
test_predictor_table$mod.me.predict <- 100*predict(mod.me.final, test_data_subset, type = "cloglog")



#Determine threshold using training data
mod.glm.predict <- predict(mod.glm.final, type = "response") 
mod.gam.predict <- predict(mod.gam.final, type = "response") 
mod.rf.predict  <- predict(mod.rf.final, type = "prob")[,2] 
mod.brt.predict <- predict(mod.brt.final, type = "response",  n.trees = best_iter)
mod.me.predict <-  predict(mod.me.final, maxent_predictors_1, type = "cloglog")


modl.glm <- "mod.glm.final"
modl.gam <- "mod.gam.final"
modl.rf <- "mod.rf.final"
modl.brt <- "mod.brt.final"
modl.me <- "mod.me.final"

data.glm <- cbind(modl.glm, predictor_table_final[2], mod.glm.predict)
data.gam <- cbind(modl.gam, predictor_table_final[2], mod.gam.predict)
data.rf <- cbind(modl.rf, predictor_table_final[2], mod.rf.predict)
data.brt <- cbind(modl.brt, predictor_table_final[2], mod.brt.predict)
data.me <- cbind(modl.me, predictor_table_final[2], mod.me.predict)
#Threshold determination
glm.mod.cut <- optimal.thresholds(data.glm, opt.methods = c("MaxKappa"))
gam.mod.cut <- optimal.thresholds(data.gam, opt.methods = c("MaxKappa"))
rf.mod.cut <- optimal.thresholds(data.rf, opt.methods = c("MaxKappa"))
brt.mod.cut <- optimal.thresholds(data.brt, opt.methods = c("MaxKappa"))
me.mod.cut <- optimal.thresholds(data.me, opt.methods = c("MaxKappa"))
glm.mod.cut
gam.mod.cut
rf.mod.cut[,2]
brt.mod.cut
me.mod.cut

#Assess strength of models on test data
# Binary predictions
test_predictor_table$me.test.pred.binary <- ifelse(test_predictor_table$mod.me.predict >= 100*me.mod.cut[,2], 1, 0)
# Confusion matrix and metrics
conf.mat <- table(Predicted = test_predictor_table$me.test.pred.binary, Observed = test_predictor_table$PA)
print(conf.mat)

# Accuracy
accuracy <- sum(diag(conf.mat)) / sum(conf.mat)
print(accuracy*100)

"
GLM accuracy: %
GAM accuracy: %
 RF accuracy: %
BRT accuracy: %
 ME accuracy: %
"

#Training Data
mod.gam.final.predict <- predict(mod.gam.final, type = "response") 
#mod.brt.final.predict <- predict(mod.brt.final, type = "response", n.trees = best_iter) 
#mod.me.final.predict <-  predict(mod.me.final, maxent_predictors_1, type = "cloglog")
mod.gam.final.acc <- presence.absence.accuracy(data.gam, threshold = gam.mod.cut$mod.gam.predict,  st.dev = F)
tss <- mod.gam.final.acc$sensitivity + mod.gam.final.acc$specificity - 1
mod.gam.final.acc <- cbind(mod.gam.final.acc[1:7], tss)
pander::pander(mod.gam.final.acc[c(1,4:5,7:8)])

#Test Data
mod.predict.test <- predict(mod.gam.final, newdata = test_predictor_table, type = "response") 
mod.test.acc <- presence.absence.accuracy(data.gam, threshold = gam.mod.cut$mod.gam.predict,  st.dev = F)
tss <- mod.test.acc$sensitivity + mod.test.acc$specificity - 1
mod.test.acc <- cbind(mod.test.acc[1:7], tss)
pander::pander(mod.test.acc[c(1,4:5,7:8)])





###Pull in all climate data, predict####
PRISM_Data <- read.csv("./Results/20250415_Accessions_Climate/AllCounties_climate_OSC_elevation_v3.csv", header=TRUE)
PRISM_Data_4_glm <- PRISM_Data[,c("GEOID",
                                  "mean_norm_ppt_05",
                                  'mean_norm_ppt_06',
                                  'mean_norm_ppt_09',
                                  "mean_norm_ppt_10",
                                  'mean_norm_ppt_annual',
                                  'mean_norm_solclear_05',
                                  'mean_norm_solclear_07',
                                  'mean_norm_soltotal_03',
                                  "mean_norm_soltotal_06",
                                  'mean_norm_soltotal_07',
                                  'mean_norm_soltrans_02',
                                  'mean_norm_soltrans_09',
                                  'mean_norm_vpdmax_08',
                                  'mean_norm_vpdmin_annual',
                                  "osc_mtu_per_ha")]
PRISM_Data_4_glm$predicted_prob_glm <- predict(mod.glm.final, newdata = PRISM_Data_4_glm, type = "response")

PRISM_Data_4_gam <- PRISM_Data[,c("GEOID","mean_norm_ppt_05",
                                  'mean_norm_ppt_06',
                                  'mean_norm_ppt_09',
                                  "mean_norm_ppt_10",
                                  'mean_norm_ppt_annual',
                                  'mean_norm_solclear_05',
                                  'mean_norm_solclear_07',
                                  'mean_norm_soltotal_03',
                                  "mean_norm_soltotal_06",
                                  'mean_norm_soltotal_07',
                                  'mean_norm_soltrans_02',
                                  'mean_norm_soltrans_09',
                                  'mean_norm_vpdmax_08',
                                  'mean_norm_vpdmin_annual',
                                  "osc_mtu_per_ha")]
PRISM_Data_4_gam$predicted_prob_gam <- predict(mod.gam.final, newdata = PRISM_Data_4_gam, type = "response")

PRISM_Data_4_rf <- PRISM_Data[,c("GEOID", "mean_norm_ppt_05",
                                 'mean_norm_ppt_06',
                                 'mean_norm_ppt_09',
                                 "mean_norm_ppt_10",
                                 'mean_norm_ppt_annual',
                                 'mean_norm_solclear_05',
                                 'mean_norm_solclear_07',
                                 'mean_norm_soltotal_03',
                                 "mean_norm_soltotal_06",
                                 'mean_norm_soltotal_07',
                                 'mean_norm_soltrans_02',
                                 'mean_norm_soltrans_09',
                                 'mean_norm_vpdmax_08',
                                 'mean_norm_vpdmin_annual',
                                 "osc_mtu_per_ha")]
PRISM_Data_4_rf$predicted_prob_rf <- predict(mod.rf.final, newdata = PRISM_Data_4_rf, type = "prob")[,2]


PRISM_Data_4_brt <- PRISM_Data[,c("GEOID","mean_norm_ppt_05",
                                  'mean_norm_ppt_06',
                                  'mean_norm_ppt_09',
                                  "mean_norm_ppt_10",
                                  'mean_norm_ppt_annual',
                                  'mean_norm_solclear_05',
                                  'mean_norm_solclear_07',
                                  'mean_norm_soltotal_03',
                                  "mean_norm_soltotal_06",
                                  'mean_norm_soltotal_07',
                                  'mean_norm_soltrans_02',
                                  'mean_norm_soltrans_09',
                                  'mean_norm_vpdmax_08',
                                  'mean_norm_vpdmin_annual',
                                  "osc_mtu_per_ha")]
PRISM_Data_4_brt <- merge(PRISM_Data_4_brt, train.location, by = "GEOID", all.x = TRUE)
PRISM_Data_4_brt$Resistance <- ifelse(is.na(PRISM_Data_4_brt$Resistance), 0.002, PRISM_Data_4_brt$Resistance)
PRISM_Data_4_brt$predicted_prob_brt <- predict(mod.brt.final, newdata = PRISM_Data_4_brt, type = "response", n.trees = best_iter)

PRISM_Data_4_me <- PRISM_Data[complete.cases(PRISM_Data$osc_mtu_per_ha),c("GEOID","mean_norm_ppt_05",
                                                                          'mean_norm_ppt_06',
                                                                          'mean_norm_ppt_09',
                                                                          "mean_norm_ppt_10",
                                                                          'mean_norm_ppt_annual',
                                                                          'mean_norm_solclear_05',
                                                                          'mean_norm_solclear_07',
                                                                          'mean_norm_soltotal_03',
                                                                          "mean_norm_soltotal_06",
                                                                          'mean_norm_soltotal_07',
                                                                          'mean_norm_soltrans_02',
                                                                          'mean_norm_soltrans_09',
                                                                          'mean_norm_vpdmax_08',
                                                                          'mean_norm_vpdmin_annual',
                                                                          "osc_mtu_per_ha")]
temp_geoids <- PRISM_Data_4_me$GEOID
me_predictors <- PRISM_Data_4_me[,c("mean_norm_ppt_05",
                                    'mean_norm_ppt_06',
                                    'mean_norm_ppt_09',
                                    "mean_norm_ppt_10",
                                    'mean_norm_ppt_annual',
                                    'mean_norm_solclear_05',
                                    'mean_norm_solclear_07',
                                    'mean_norm_soltotal_03',
                                    "mean_norm_soltotal_06",
                                    'mean_norm_soltotal_07',
                                    'mean_norm_soltrans_02',
                                    'mean_norm_soltrans_09',
                                    'mean_norm_vpdmax_08',
                                    'mean_norm_vpdmin_annual',
                                    "osc_mtu_per_ha")]
predicted_prob_me <-  predict(mod.me.final, x = me_predictors, newdata = me_predictors,  type = "cloglog")
PRISM_Data_4_me$predicted_prob_me <- predicted_prob_me



Models_DF <- merge(PRISM_Data_4_glm,PRISM_Data_4_gam, by = "GEOID")
Models_DF <- merge(Models_DF,PRISM_Data_4_rf, by = "GEOID" )
Models_DF <- merge(Models_DF,PRISM_Data_4_brt, by = "GEOID" )
Models_DF <- merge(Models_DF,PRISM_Data_4_me, by = "GEOID" )
Models_DF <- Models_DF[,c("GEOID","predicted_prob_glm","predicted_prob_gam","predicted_prob_rf","predicted_prob_brt", "predicted_prob_me")]

#Models_DF <- Models_DF %>%
#  mutate(
#    mod_avg = rowMeans(dplyr::select(., predicted_prob_glm, predicted_prob_gam, predicted_prob_rf, predicted_prob_brt), na.rm = TRUE),
#    mod_stdev = rowSds(as.matrix(dplyr::select(., predicted_prob_glm, predicted_prob_gam, predicted_prob_rf, predicted_prob_brt), na.rm = TRUE)),
#    mod_variance = (rowSds(as.matrix(dplyr::select(., predicted_prob_glm, predicted_prob_gam, predicted_prob_rf, predicted_prob_brt), na.rm = TRUE)))^2
#  )

write.csv(Models_DF, "./Results/20250501_HabitSuitMod/GenoJCResist_ModelsOutputContinuous.csv", row.names=FALSE)
data <- read.csv("./Results/20250501_HabitSuitMod/GenoJCResist_ModelsOutputContinuous.csv", header=TRUE)
glm.mod.cut <- 0.44
'0.44'
gam.mod.cut <- 0.335
'0.335'
rf.mod.cut <- 0.58
'0.58'
brt.mod.cut <- 0.98
'0.98'
me.mod.cut <- 0.61
'0.61'

#map it####

data <- Models_DF
data$glm_binary <- ifelse(data$predicted_prob_glm >= glm.mod.cut,1,0)
data$gam_binary <- ifelse(data$predicted_prob_gam >= gam.mod.cut,1,0)
data$rf_binary <- ifelse(data$predicted_prob_rf >= rf.mod.cut,1,0)
data$brt_binary <- ifelse(data$predicted_prob_brt >= brt.mod.cut,1,0)
data$me_binary <- ifelse(data$predicted_prob_me >= me.mod.cut,1,0)
data$model_concensus <- rowSums(data[, c("glm_binary", "gam_binary", "rf_binary", "brt_binary", "me_binary")], na.rm = TRUE)

write.csv(data, "./Results/20250501_HabitSuitMod/GenoJCResist_ModelsOutputContinuous.csv", row.names=FALSE)

