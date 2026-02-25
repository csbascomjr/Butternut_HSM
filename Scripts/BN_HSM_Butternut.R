##Adapted from https://rpubs.com/JoshCarrell/HSM
library(ENMeval)
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
library(ranger)
library(pander)
library(gridExtra)
library(gbm)
library(dismo)
library(maxnet)
library(matrixStats)
setwd("/Volumes/Onni_Drive/Fearer_Projects/Butternut_Landscape_Genomics")

###Uselocations for all JC Genotyped Accession ####
Location_Data <- read.csv("./Results/20240819_AccessionLocations/bnut_families_wBV_sp.csv", header=TRUE)
#Subset location data to include just butternut genotyped families (both more and less resistant)
Location_Data_JC <- Location_Data[Location_Data$species == "JC", c("family","GEOID","bays_bcd_bv")]
Location_Data$GEOID <- sprintf("%05d", Location_Data$fips)

VA_City_GEOIDS <- as.character(c(51510,51520,51530,51540,51550,51570,51580,51590,51595,51620,51630,51640,51650,
                                 51660,51670,51680,51683,51690,51700,51710,51720,51730,51735,51740,51750,51760,
                                 51770,51775,51790,51800,51810,51820,51830,51840))

PRISM_Data <- read.csv("./Results/20250415_Accessions_Climate/AllCounties_climate_OSC_elevation_v3.csv", header=TRUE)
PRISM_Data$GEOID <- sprintf("%05d", PRISM_Data$GEOID)
PRISM_Data <- PRISM_Data[ , !grepl("^sd_", names(PRISM_Data))]
PRISM_Data[PRISM_Data$GEOID == "09003",]
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
  GEOID = unique(Location_Data_JC$GEOID),
  ind = sample(2, length(unique(Location_Data_JC$GEOID)), replace = TRUE, prob = c(0.7, 0.3))
)

nrow(df[df$ind == 1,])
nrow(df[df$ind == 2,])

occurance_geoids <- df[df$ind == 1, c("GEOID")]

unavailable_GEOIDS <- c(VA_City_GEOIDS, as.character(Location_Data_JC$GEOID))
available_GEOIDS <- PRISM_Data[!PRISM_Data$GEOID %in% unavailable_GEOIDS, "GEOID"]
train.absence_GEOIDS <- sample(available_GEOIDS, 
                               size = nrow(df[df$ind == 1,])*3)


unavailable_GEOIDS <- c(train.absence_GEOIDS, VA_City_GEOIDS, as.character(Location_Data_JC$GEOID))
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
  left_join(Location_Data %>%
              group_by(GEOID) %>%
              summarise(bays_bcd_bv = mean(bays_bcd_bv, na.rm = TRUE)),
            by = "GEOID")
colnames(test.occurance) <- c("GEOID","PA","Resistance")

test.location <- bind_rows(test.occurance,test.absence)
nrow(test.location[test.location$PA == "0",])


train.occurance <- data.frame(
  GEOID = occurance_geoids,
  PA = 1)

train.absence <- data.frame(
  GEOID = train.absence_GEOIDS,
  PA = 0)

train.occurance$GEOID <- as.character(train.occurance$GEOID)
train.absence$GEOID <- as.character(train.absence$GEOID)
train.location <- bind_rows(train.occurance, train.absence)
nrow(train.location[train.location$PA == "0",])
nrow(train.location)

#Check chosen counties to see if they make sense
train.data_to_map <- merge(train.location,us_counties_shapefile_simpl, by = "GEOID")
test.data_to_map <- merge(test.location,us_counties_shapefile_simpl, by = "GEOID")


ggplot() +
  geom_sf(data = us_counties_shapefile_simpl, aes(geometry = geometry), fill = "grey", color = "NA", show.legend = FALSE) +
  geom_sf(data = train.data_to_map, aes(geometry = geometry, fill = as.factor(PA)), alpha = 0.7, color = "NA", show.legend = FALSE) +
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
  geom_sf(data = test.data_to_map, aes(geometry = geometry), fill = "firebrick", alpha = 0.7, color = "NA", show.legend = TRUE) +
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

#Automatic way that might actually be better but I'm not convinced
high.cor.pairs <- which(abs(topo.cor) > 0.7 & abs(topo.cor) < 1, arr.ind = TRUE)
vars.to.remove <- unique(rownames(high.cor.pairs)[high.cor.pairs[,1] < high.cor.pairs[,2]])
topo.reduced <- spec.cor.topo[ , !(names(spec.cor.topo) %in% vars.to.remove)]
topo.reduced.cor <- cor(topo.reduced)
corrplot.mixed(topo.reduced.cor, lower = 'number', upper = 'square')
corrplot.mixed(topo.reduced.cor)
rownames(topo.reduced.cor)
write.csv(topo.reduced.cor, "./PRISM_Data_Correaltions_v4.csv", row.names=TRUE)

#Manual way that is pretty tedious but has, thus far, provided better results.
spec.cor.topo_v2 <- predictor_table[,c(3:134)]
topo.cor_v2 <- cor(spec.cor.topo_v2) # correlation
corrplot.mixed(topo.cor_v2, lower = 'number', upper = 'square')
write.csv(topo.cor_v2, "./Results/20250501_HabitSuitMod/PRISM_Data_GenoJCAll_Correlations_v1.csv", row.names=TRUE)

predictor_table_final <- predictor_table[,c("GEOID","PA",
                                            "mean_norm_ppt_04",
                                            "mean_norm_ppt_05",
                                            "mean_norm_ppt_10",
                                            'mean_norm_ppt_annual',
                                            'mean_norm_solclear_05',
                                            "mean_norm_solclear_07",
                                            "mean_norm_soltotal_06",
                                            "mean_norm_soltotal_07",
                                            "mean_norm_soltrans_09",
                                            "mean_norm_soltrans_annual",
                                            "mean_norm_vpdmax_07",
                                            'mean_norm_vpdmin_annual',
                                            "osc_mtu_per_ha")]
spec.cor.topo_v2 <- predictor_table_final[,c(3:15)]
topo.cor_v2 <- cor(spec.cor.topo_v2) # correlation
corrplot.mixed(topo.cor_v2, lower = 'number', upper = 'square')

write.csv(predictor_table_final, "./Results/20250501_HabitSuitMod/PRISM_Data_GenoJCAll_PredictorTableFinal_v1.csv", row.names=FALSE)
#Bring up final predictor table, if wanted
predictor_table_final <- read.csv("./Results/20250501_HabitSuitMod/PRISM_Data_GenoJCAll_PredictorTableFinal_v1.csv", header=TRUE)

#Generalized Linear Model####
mod.glm.1 <- glm(as.factor(PA) ~ 
                 mean_norm_ppt_04 +
                 mean_norm_ppt_05 +
                 mean_norm_ppt_10 +
                 mean_norm_ppt_annual +
                 mean_norm_solclear_05 +
                 mean_norm_solclear_07 +
                 mean_norm_soltotal_06 +
                 mean_norm_soltotal_07 +
                 mean_norm_soltrans_09 +
                 mean_norm_soltrans_annual +
                 mean_norm_vpdmax_07 +
                 mean_norm_vpdmin_annual +
                 osc_mtu_per_ha,
                 family = binomial, data = predictor_table_final)
summary(mod.glm.1)
mod.glm.2 <- glm(as.factor(PA) ~ 
                   mean_norm_ppt_04 +
                   mean_norm_ppt_05 +
                   mean_norm_ppt_10 +
                   mean_norm_ppt_annual +
                   #mean_norm_solclear_05 +
                   #mean_norm_solclear_07 +
                   #mean_norm_soltotal_06 +
                   #mean_norm_soltotal_07 +
                   #mean_norm_soltrans_09 +
                   mean_norm_soltrans_annual, # +
                   #mean_norm_vpdmax_07 +
                   #mean_norm_vpdmin_annual +
                   #osc_mtu_per_ha,
                 family = binomial, data = predictor_table_final)
summary(mod.glm.2)
mod.glm.3 <- glm(as.factor(PA) ~ 
                   #mean_norm_ppt_04 +
                   mean_norm_ppt_05 +
                   #mean_norm_ppt_10 +
                   #mean_norm_ppt_annual +
                   #mean_norm_solclear_05 +
                   #mean_norm_solclear_07 +
                   #mean_norm_soltotal_06 +
                   #mean_norm_soltotal_07 +
                   #mean_norm_soltrans_09 +
                   mean_norm_soltrans_annual, # +
                 #mean_norm_vpdmax_07 +
                 #mean_norm_vpdmin_annual +
                 #osc_mtu_per_ha,
                 family = binomial, data = predictor_table_final)
summary(mod.glm.3)
AIC(mod.glm.1,mod.glm.2,mod.glm.3)
mod1.fit <- 100*(1-mod.glm.1$deviance/mod.glm.1$null.deviance)
mod2.fit <- 100*(1-mod.glm.2$deviance/mod.glm.2$null.deviance)
mod3.fit <- 100*(1-mod.glm.3$deviance/mod.glm.3$null.deviance)
print(mod1.fit)
print(mod2.fit)
print(mod3.fit)

roc_obj <- roc(predictor_table_final$PA, predict(mod.glm.2, type = "response"))
auc(roc_obj)

mod.glm.final <- mod.glm.2

#Generalized Additive Model####
mod.gam.1 <- gam(as.factor(PA) ~ 
                   s(mean_norm_ppt_04) +
                   s(mean_norm_ppt_05) +
                   s(mean_norm_ppt_10) +
                   s(mean_norm_ppt_annual) +
                   s(mean_norm_solclear_05) +
                   s(mean_norm_solclear_07) +
                   s(mean_norm_soltotal_06) +
                   s(mean_norm_soltotal_07) +
                   s(mean_norm_soltrans_09) +
                   s(mean_norm_soltrans_annual) +
                   s(mean_norm_vpdmax_07) +
                   s(mean_norm_vpdmin_annual) +
                   s(osc_mtu_per_ha),
                  family = binomial, 
                  select = TRUE,
                  data = predictor_table_final)
summary(mod.gam.1)
mod.gam.2 <- gam(as.factor(PA) ~ 
                   s(mean_norm_ppt_04) +
                   s(mean_norm_ppt_05) +
                   #s(mean_norm_ppt_10) +
                   #s(mean_norm_ppt_annual) +
                   #s(mean_norm_solclear_05) +
                   s(mean_norm_solclear_07) +
                   #s(mean_norm_soltotal_06) +
                   s(mean_norm_soltotal_07) +
                   #s(mean_norm_soltrans_09) +
                   s(mean_norm_soltrans_annual) +
                   #s(mean_norm_vpdmax_07) +
                   #s(mean_norm_vpdmin_annual) +
                   s(osc_mtu_per_ha),
                  family = binomial, 
                  select = TRUE,
                  data = predictor_table_final)
summary(mod.gam.2)
dfs <- sum(influence(mod.gam.2))
mod.gam.3 <- gam(as.factor(PA) ~
                   s(mean_norm_ppt_04, k = 3) +
                   s(mean_norm_ppt_05, k = 3) +
                   #s(mean_norm_ppt_10) +
                   #s(mean_norm_ppt_annual) +
                   #s(mean_norm_solclear_05) +
                   s(mean_norm_solclear_07, k = 3) +
                   #s(mean_norm_soltotal_06) +
                   s(mean_norm_soltotal_07, k = 3) +
                   #s(mean_norm_soltrans_09) +
                   s(mean_norm_soltrans_annual, k = 3) +
                   #s(mean_norm_vpdmax_07) +
                   #s(mean_norm_vpdmin_annual) +
                   s(osc_mtu_per_ha, k = 3),
                 family = binomial, 
                 select = TRUE,
                 data = predictor_table_final)
summary(mod.gam.3)
gam.check(mod.gam.3)

mod.gam.4 <- gam(as.factor(PA) ~ 
                   s(mean_norm_ppt_04, k = 3) +
                   s(mean_norm_ppt_05, k = 3) +
                   #s(mean_norm_ppt_10) +
                   #s(mean_norm_ppt_annual) +
                   #s(mean_norm_solclear_05) +
                   #s(mean_norm_solclear_07, k = 3) +
                   #s(mean_norm_soltotal_06) +
                   #s(mean_norm_soltotal_07, k = 3) +
                   #s(mean_norm_soltrans_09) +
                   s(mean_norm_soltrans_annual, k = 3) +
                   #s(mean_norm_vpdmax_07) +
                   #s(mean_norm_vpdmin_annual) +
                   s(osc_mtu_per_ha, k = 3),
                 family = binomial, 
                 select = TRUE,
                 data = predictor_table_final)
summary(mod.gam.4)
gam.check(mod.gam.4)
AIC(mod.gam.1, mod.gam.2, mod.gam.3,mod.gam.4)


roc_obj <- roc(predictor_table_final$PA, predict(mod.gam.2, type = "response"))
auc(roc_obj)

mod1.fit <- 100*(1-mod.gam.1$deviance/mod.gam.1$null.deviance)
mod2.fit <- 100*(1-mod.gam.2$deviance/mod.gam.2$null.deviance)
mod3.fit <- 100*(1-mod.gam.3$deviance/mod.gam.3$null.deviance)
print(mod1.fit)
print(mod2.fit)
print(mod3.fit)

#Stick with model 3

mod.gam.final <- mod.gam.2


#Random Forest####
set.seed(9242)
mod.rf.1 <- randomForest(as.factor(PA) ~ 
                           mean_norm_ppt_04 +
                           mean_norm_ppt_05 +
                           mean_norm_ppt_10 +
                           mean_norm_ppt_annual +
                           mean_norm_solclear_05 +
                           mean_norm_solclear_07 +
                           mean_norm_soltotal_06 +
                           mean_norm_soltotal_07 +
                           mean_norm_soltrans_09 +
                           mean_norm_soltrans_annual +
                           mean_norm_vpdmax_07 +
                           mean_norm_vpdmin_annual +
                           osc_mtu_per_ha,
                           data = predictor_table_final, 
                           ntree = 3000, 
                           mtry = 5,
                           sampsize = c('0'=68,'1'=68),
                           #classwt = c("0" = 1, "1" = 5),
                           importance = TRUE,
                           proximity=TRUE) 
print(mod.rf.1)
pred_probs <- predict(mod.rf.1, type = "prob")[,2]
roc_curve <- roc(predictor_table_final$PA, pred_probs)
auc(roc_curve) 

importance(mod.rf.1)
varImpPlot(mod.rf.1)

nrow(predictor_table_final[predictor_table_final$PA == "1",])

mod.rf.2 <- randomForest(as.factor(PA) ~ 
                           mean_norm_ppt_04 +
                           mean_norm_ppt_05 +
                           mean_norm_ppt_10 +
                           mean_norm_ppt_annual +
                           #mean_norm_solclear_05 +
                           mean_norm_solclear_07 +
                           #mean_norm_soltotal_06 +
                           #mean_norm_soltotal_07 +
                           mean_norm_soltrans_09 +
                           mean_norm_soltrans_annual +
                           mean_norm_vpdmax_07 +
                           mean_norm_vpdmin_annual +
                           osc_mtu_per_ha,
                         data = predictor_table_final, 
                         ntree = 3000, 
                         mtry = 5,
                         sampsize = c('0'=68,'1'=68),
                           #classwt = c("0" = 1, "1" = 10),
                           importance = TRUE,
                           proximity=TRUE)
print(mod.rf.2)
mod.rf.2$call
summary(mod.rf.2)
randomForest::importance(mod.rf.2)
varImpPlot(mod.rf.2)

pred_probs <- predict(mod.rf.2, type = "prob")[,2]
roc_curve <- roc(predictor_table_final$PA, pred_probs)
auc(roc_curve) 

mod.rf.3 <- randomForest(as.factor(PA) ~ 
                           #mean_norm_ppt_05 +
                           #mean_norm_ppt_09 +
                           #mean_norm_ppt_10 +
                           mean_norm_ppt_annual +
                           #mean_norm_solclear_05 +
                           #mean_norm_solclear_07 +
                           #mean_norm_soltotal_06 +
                           #mean_norm_soltrans_07 +
                           mean_norm_soltrans_annual +
                           #mean_norm_vpdmin_02 +
                           #mean_norm_vpdmin_annual +
                           osc_mtu_per_ha,
                         data = predictor_table_final, 
                         ntree = 4000, 
                         mtry = 1,
                         sampsize = c('0'=50,'1'=50),
                         nodesize = 1,
                        #classwt = c("0" = 1, "1" = 5),
                         importance = TRUE,
                         proximity=TRUE) 
print(mod.rf.3)
pred_probs <- predict(mod.rf.3, type = "prob")[,2]
roc_curve <- roc(predictor_table_final$PA, pred_probs)
auc(roc_curve) 

predictor_table_final$mean_norm_ppt_soltrans_interaction <- predictor_table_final$mean_norm_ppt_annual*predictor_table_final$mean_norm_soltrans_annual
predictor_table_final$osc_mtu_per_ha_squared <- predictor_table_final$osc_mtu_per_ha^2

mod.rf.3_engineered <- randomForest(as.factor(PA) ~ 
                           #mean_norm_ppt_05 +
                           #mean_norm_ppt_09 +
                           #mean_norm_ppt_10 +
                           mean_norm_ppt_annual +
                           mean_norm_ppt_soltrans_interaction +
                           osc_mtu_per_ha_squared + 
                           #mean_norm_solclear_05 +
                           #mean_norm_solclear_07 +
                           #mean_norm_soltotal_06 +
                           #mean_norm_soltrans_07 +
                           mean_norm_soltrans_annual +
                           #mean_norm_vpdmin_02 +
                           #mean_norm_vpdmin_annual +
                           osc_mtu_per_ha,
                         data = predictor_table_final, 
                         ntree = 4000, 
                         mtry = 4,
                         sampsize = c('0'=68,'1'=68),
                         nodesize = 1,
                         #classwt = c("0" = 1, "1" = 5),
                         importance = TRUE,
                         proximity=TRUE) 
print(mod.rf.3_engineered)

#

mod.rf.final <- mod.rf.2

#Boosted Regression Tree (BRT) Model ####
mod.BRT.1 <- gbm(formula = PA ~ 
                   mean_norm_ppt_04 +
                   mean_norm_ppt_05 +
                   mean_norm_ppt_10 +
                   mean_norm_ppt_annual +
                   mean_norm_solclear_05 +
                   mean_norm_solclear_07 +
                   mean_norm_soltotal_06 +
                   mean_norm_soltotal_07 +
                   mean_norm_soltrans_09 +
                   mean_norm_soltrans_annual +
                   mean_norm_vpdmax_07 +
                   mean_norm_vpdmin_annual +
                   osc_mtu_per_ha,
                        data = predictor_table_final,
                        distribution = "bernoulli",
                        n.trees = 4000,
                        cv.folds = 10,
                        interaction.depth = 4,
                        shrinkage = 0.001,
                        n.minobsinnode = 10)
summary(mod.BRT.1)
best_iter <- gbm.perf(mod.BRT.1, method = "cv")
print(best_iter)
pred_probs <- predict(mod.BRT.1, type = "response", n.trees = best_iter)
roc_curve <- roc(predictor_table_final$PA, pred_probs)
plot(roc_curve, main = "ROC Curve", col = "blue")
pROC::auc(roc_curve)

mod.BRT.2 <- gbm(formula = PA ~ 
                   mean_norm_ppt_04 +
                   mean_norm_ppt_05 +
                   mean_norm_ppt_10 +
                   mean_norm_ppt_annual +
                   mean_norm_solclear_05 +
                   mean_norm_solclear_07 +
                   #mean_norm_soltotal_06 +
                   mean_norm_soltotal_07 +
                   mean_norm_soltrans_09 +
                   mean_norm_soltrans_annual +
                   mean_norm_vpdmax_07 +
                   mean_norm_vpdmin_annual +
                   osc_mtu_per_ha,
                 data = predictor_table_final,
                 distribution = "bernoulli",
                 n.trees = best_iter,
                 cv.folds = 10,
                 interaction.depth = 4,
                 shrinkage = 0.001,
                 n.minobsinnode = 10)
summary(mod.BRT.2)
best_iter <- gbm.perf(mod.BRT.2, method = "cv")
print(best_iter)
pred_probs <- predict(mod.BRT.2, type = "response", n.trees = best_iter)
roc_curve <- roc(predictor_table_final$PA, pred_probs)
plot(roc_curve, main = "ROC Curve", col = "blue")
pROC::auc(roc_curve)
coords_best <- coords(roc_curve, "best", ret = "threshold", transpose = FALSE)
best_thresh <- coords_best$threshold
pred_class <- ifelse(pred_probs >= best_thresh, 1, 0)
table(Predicted = pred_class, Actual = predictor_table_final$PA)


mod.brt.final <- mod.BRT.1

#Maximum Entropy (MaxEnt) Model ####
maxent_predictors_1 <- predictor_table_final[, c(
  "mean_norm_ppt_04",
  "mean_norm_ppt_05",
  "mean_norm_ppt_10",
  'mean_norm_ppt_annual',
  'mean_norm_solclear_05',
  "mean_norm_solclear_07",
  "mean_norm_soltotal_06",
  "mean_norm_soltotal_07",
  "mean_norm_soltrans_09",
  "mean_norm_soltrans_annual",
  "mean_norm_vpdmax_07",
  'mean_norm_vpdmin_annual',
  "osc_mtu_per_ha"
)]

# Run MaxEnt
mod.me.1 <- maxnet(p = predictor_table_final$PA, data = maxent_predictors_1)
str(mod.me.1)
print(mod.me.1)
plot(mod.me.1)
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
str(contrib_summary)

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
pred_class <- ifelse(pred_probs >= 0.3, 1, 0)
table(Predicted = pred_class, Actual = predictor_table_final$PA)
plot(mod.ME.1)
response(mod.ME.1)

mod.me.final <- mod.me.1


"-----------------------------------------------------------------------------------------------------"
#Use models to predict across whole study area####
"-----------------------------------------------------------------------------------------------------"
test_predictor_table <- merge(test.location, PRISM_Data, by = "GEOID", all.x=TRUE)
test_predictor_table <- test_predictor_table[,c("GEOID","PA",
                                                "mean_norm_ppt_04",
                                                "mean_norm_ppt_05",
                                                "mean_norm_ppt_10",
                                                'mean_norm_ppt_annual',
                                                'mean_norm_solclear_05',
                                                "mean_norm_solclear_07",
                                                "mean_norm_soltotal_06",
                                                "mean_norm_soltotal_07",
                                                "mean_norm_soltrans_09",
                                                "mean_norm_soltrans_annual",
                                                "mean_norm_vpdmax_07",
                                                'mean_norm_vpdmin_annual',
                                                "osc_mtu_per_ha")]
test_predictor_table$mean_norm_ppt_soltrans_interaction <- test_predictor_table$mean_norm_ppt_annual*test_predictor_table$mean_norm_soltrans_annual
test_predictor_table$mean_norm_ppt_soltrans_squared <- test_predictor_table$mean_norm_soltrans_annual^2

write.csv(test_predictor_table, "./Results/20250501_HabitSuitMod/PRISM_Data_GenoJCAll_TestDataPredictorTable.csv", row.names=FALSE)


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

#Training Data model assessment
mod.glm.final.predict <- predict(mod.glm.final, type = "response") 
#mod.brt.final.predict <- predict(mod.brt.final, type = "response", n.trees = best_iter) 
#mod.me.final.predict <-  predict(mod.me.final, maxent_predictors_1, type = "cloglog")
mod.final.acc <- presence.absence.accuracy(data.glm, threshold = glm.mod.cut$mod.glm.predict,  st.dev = F)
tss <- mod.final.acc$sensitivity + mod.final.acc$specificity - 1
mod.final.acc <- cbind(mod.final.acc[1:7], tss)
pander::pander(mod.final.acc[c(1,4:5,7:8)])


###Pull in all climate data, predict####
PRISM_Data <- read.csv("./Results/20250415_Accessions_Climate/AllCounties_climate_OSC_elevation_v3.csv", header=TRUE)
PRISM_Data_4_glm <- PRISM_Data[,c("GEOID","mean_norm_ppt_04",
                                  "mean_norm_ppt_05",
                                  "mean_norm_ppt_10",
                                  'mean_norm_ppt_annual',
                                  'mean_norm_solclear_05',
                                  "mean_norm_solclear_07",
                                  "mean_norm_soltotal_06",
                                  "mean_norm_soltotal_07",
                                  "mean_norm_soltrans_09",
                                  "mean_norm_soltrans_annual",
                                  "mean_norm_vpdmax_07",
                                  'mean_norm_vpdmin_annual',
                                  "osc_mtu_per_ha")]
PRISM_Data_4_glm$predicted_prob_glm <- predict(mod.glm.final, newdata = PRISM_Data_4_glm, type = "response")

PRISM_Data_4_gam <- PRISM_Data[,c("GEOID","mean_norm_ppt_04",
                                  "mean_norm_ppt_05",
                                  "mean_norm_ppt_10",
                                  'mean_norm_ppt_annual',
                                  'mean_norm_solclear_05',
                                  "mean_norm_solclear_07",
                                  "mean_norm_soltotal_06",
                                  "mean_norm_soltotal_07",
                                  "mean_norm_soltrans_09",
                                  "mean_norm_soltrans_annual",
                                  "mean_norm_vpdmax_07",
                                  'mean_norm_vpdmin_annual',
                                  "osc_mtu_per_ha")]
PRISM_Data_4_gam$predicted_prob_gam <- predict(mod.gam.final, newdata = PRISM_Data_4_gam, type = "response")

PRISM_Data_4_rf <- PRISM_Data[,c("GEOID","mean_norm_ppt_04",
                                 "mean_norm_ppt_05",
                                 "mean_norm_ppt_10",
                                 'mean_norm_ppt_annual',
                                 'mean_norm_solclear_05',
                                 "mean_norm_solclear_07",
                                 "mean_norm_soltotal_06",
                                 "mean_norm_soltotal_07",
                                 "mean_norm_soltrans_09",
                                 "mean_norm_soltrans_annual",
                                 "mean_norm_vpdmax_07",
                                 'mean_norm_vpdmin_annual',
                                 "osc_mtu_per_ha")]

PRISM_Data_4_rf$mean_norm_ppt_soltrans_interaction <- PRISM_Data_4_rf$mean_norm_ppt_annual*PRISM_Data_4_rf$mean_norm_soltrans_annual
PRISM_Data_4_rf$mean_norm_ppt_soltrans_squared <- PRISM_Data_4_rf$mean_norm_soltrans_annual^2

PRISM_Data_4_rf$predicted_prob_rf <- predict(mod.rf.final, newdata = PRISM_Data_4_rf, type = "prob")[,2]


PRISM_Data_4_brt <- PRISM_Data[,c("GEOID","mean_norm_ppt_04",
                                  "mean_norm_ppt_05",
                                  "mean_norm_ppt_10",
                                  'mean_norm_ppt_annual',
                                  'mean_norm_solclear_05',
                                  "mean_norm_solclear_07",
                                  "mean_norm_soltotal_06",
                                  "mean_norm_soltotal_07",
                                  "mean_norm_soltrans_09",
                                  "mean_norm_soltrans_annual",
                                  "mean_norm_vpdmax_07",
                                  'mean_norm_vpdmin_annual',
                                  "osc_mtu_per_ha")]
PRISM_Data_4_brt$predicted_prob_brt <- predict(mod.brt.final, newdata = PRISM_Data_4_brt, type = "response", n.trees = best_iter)

PRISM_Data_4_me <- PRISM_Data[complete.cases(PRISM_Data$osc_mtu_per_ha),c("GEOID","mean_norm_ppt_04",
                                                                          "mean_norm_ppt_05",
                                                                          "mean_norm_ppt_10",
                                                                          'mean_norm_ppt_annual',
                                                                          'mean_norm_solclear_05',
                                                                          "mean_norm_solclear_07",
                                                                          "mean_norm_soltotal_06",
                                                                          "mean_norm_soltotal_07",
                                                                          "mean_norm_soltrans_09",
                                                                          "mean_norm_soltrans_annual",
                                                                          "mean_norm_vpdmax_07",
                                                                          'mean_norm_vpdmin_annual',
                                                                          "osc_mtu_per_ha")]
temp_geoids <- PRISM_Data_4_me$GEOID
me_predictors <- PRISM_Data_4_me[,c("mean_norm_ppt_04",
                                    "mean_norm_ppt_05",
                                    "mean_norm_ppt_10",
                                    'mean_norm_ppt_annual',
                                    'mean_norm_solclear_05',
                                    "mean_norm_solclear_07",
                                    "mean_norm_soltotal_06",
                                    "mean_norm_soltotal_07",
                                    "mean_norm_soltrans_09",
                                    "mean_norm_soltrans_annual",
                                    "mean_norm_vpdmax_07",
                                    'mean_norm_vpdmin_annual',
                                    "osc_mtu_per_ha")]
predicted_prob_me <-  predict(mod.me.final, x = me_predictors, newdata = me_predictors,  type = "cloglog")
PRISM_Data_4_me$predicted_prob_me <- predicted_prob_me


Models_DF <- merge(PRISM_Data_4_glm,PRISM_Data_4_gam, by = "GEOID")
Models_DF <- merge(Models_DF,PRISM_Data_4_rf, by = "GEOID" )
Models_DF <- merge(Models_DF,PRISM_Data_4_brt, by = "GEOID" )
Models_DF <- merge(Models_DF,PRISM_Data_4_me, by = "GEOID" )
Models_DF <- Models_DF[,c("GEOID","predicted_prob_glm","predicted_prob_gam","predicted_prob_rf","predicted_prob_brt", "predicted_prob_me")]

##Models_DF <- Models_DF %>%
#  mutate(
##    mod_avg = rowMeans(dplyr::select(., predicted_prob_glm, predicted_prob_gam, predicted_prob_rf, predicted_prob_brt), na.rm = TRUE),
#    mod_stdev = rowSds(as.matrix(dplyr::select(., predicted_prob_glm, predicted_prob_gam, predicted_prob_rf, predicted_prob_brt), na.rm = TRUE)),
#    mod_variance = (rowSds(as.matrix(dplyr::select(., predicted_prob_glm, predicted_prob_gam, predicted_prob_rf, predicted_prob_brt), na.rm = TRUE)))^2
#  )

data <- Models_DF
data$glm_binary <- ifelse(data$predicted_prob_glm >= glm.mod.cut[,2],1,0)
data$gam_binary <- ifelse(data$predicted_prob_gam >= gam.mod.cut[,2],1,0)
data$rf_binary <- ifelse(data$predicted_prob_rf >= rf.mod.cut[,2],1,0)
data$brt_binary <- ifelse(data$predicted_prob_brt >= brt.mod.cut[,2],1,0)
data$me_binary <- ifelse(data$predicted_prob_me >= me.mod.cut[,2],1,0)
data$model_concensus <- rowSums(data[, c("glm_binary", "gam_binary", "rf_binary", "brt_binary", "me_binary")], na.rm = TRUE)

write.csv(data, "./Results/20250501_HabitSuitMod/GenoJCAll_ModelsOutputcontinuous.csv", row.names=FALSE)






