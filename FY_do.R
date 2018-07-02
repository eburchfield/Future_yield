source("FY_load.R")
source("FY_func.R")

# scale features, divide by max or mean normalization - 
# wonky cotton results with negative yield for one county
# wwheat gdd/sdd
# better visualizations
# finish writing up methods
# read more on predict() function

###################################################################################################################
# Create historical panels
###################################################################################################################

# construct_historical_data()
# construct_future_data()
# future_interpolation()
# future_panels()

###################################################################################################################
# Load historical panels
###################################################################################################################

corn <- readRDS("./out/panels/corn_PRISM.RDS")
corn$YEAR <- as.numeric(corn$YEAR)
cotton <- readRDS("./out/panels/cotton_PRISM.RDS")
cotton$Yield <- cotton$Yield/32 # from pounds/acre to bushels/acre
cotton$YEAR <- as.numeric(cotton$YEAR)
soy <- readRDS("./out/panels/soy_PRISM.RDS")
soy$YEAR <- as.numeric(soy$YEAR)
wheat <- readRDS("./out/panels/wheat_PRISM.RDS")
wheat$YEAR <- as.numeric(wheat$YEAR)
wwheat <- readRDS("./out/panels/wwheat_PRISM.RDS")
wwheat$YEAR <- as.numeric(wwheat$YEAR)

###################################################################################################################
# Calibrate MFP
###################################################################################################################

voi <- c("GDD", "SDD", "TP", "YEAR", "pH_sur", "AWC_sur", "Slope") 
set.seed(100)

corn <- corn %>% filter(!is.na(Yield))
random_rn <- sample(nrow(corn), ceiling(nrow(corn)*.25))
HO <- corn[random_rn,]
sample <- corn[-random_rn,]
fm <- "Yield ~ fp(GDD) + fp(SDD) + fp(TP) + fp(pH_sur) +AWC_sur + fp(Slope, df=2) + YEAR"
corn_mod <- model_results(df = sample, fm = fm)
corn_data <- predict_data_PRISM(corn_mod, voi, res = 200)
ho_predict <- predict(corn_mod[["model"]], newdata = HO, type = "response")
rmse <- RMSE(HO$Yield, ho_predict)
rmse_null <- RMSE(HO$Yield, mean(HO$Yield, na.rm = T))
corn_mod[["RMSE"]] <- rmse
corn_mod[["RMSENull"]] <- rmse_null

soy <- soy %>% filter(!is.na(Yield))
random_rn <- sample(nrow(soy), ceiling(nrow(soy)*.25))
HO <- soy[random_rn,]
sample <- soy[-random_rn,]
fm <- "Yield ~ fp(GDD) + fp(SDD) + fp(TP) + fp(pH_sur) +AWC_sur + fp(Slope, df = 2) + YEAR"
soy_mod <- model_results(df = sample, fm = fm)
soy_data <- predict_data_PRISM(soy_mod, voi, res = 200)
ho_predict <- predict(soy_mod[["model"]], newdata = HO, type = "response")
rmse <- RMSE(HO$Yield, ho_predict)
rmse_null <- RMSE(HO$Yield, mean(HO$Yield, na.rm = T))
soy_mod[["RMSE"]] <- rmse
soy_mod[["RMSENull"]] <- rmse_null

cotton <- cotton %>% filter(!is.na(Yield))
random_rn <- sample(nrow(cotton), ceiling(nrow(cotton)*.25))
HO <- cotton[random_rn,]
sample <- cotton[-random_rn,]
fm <- "Yield ~ fp(GDD, df = 4) + fp(SDD) + fp(TP) + fp(pH_sur) + AWC_sur + fp(Slope, df=2) + YEAR"
cotton_mod <- model_results(df = sample, fm = fm)
cotton_data <- predict_data_PRISM(cotton_mod, voi, res = 200)
ho_predict <- predict(cotton_mod[["model"]], newdata = HO, type = "response")
rmse <- RMSE(HO$Yield, ho_predict)
rmse_null <- RMSE(HO$Yield, mean(HO$Yield, na.rm = T))
cotton_mod[["RMSE"]] <- rmse
cotton_mod[["RMSENull"]] <- rmse_null

wheat <- wheat %>% filter(!is.na(Yield))
random_rn <- sample(nrow(wheat), ceiling(nrow(wheat)*.25))
HO <- wheat[random_rn,]
sample <- wheat[-random_rn,]
fm <- "Yield ~ fp(GDD) + fp(SDD) + fp(TP) + fp(pH_sur, df=2) +AWC_sur + fp(Slope, df=2) + YEAR"
wheat_mod <- model_results(df = sample, fm = fm)
wheat_data <- predict_data_PRISM(wheat_mod, voi, res = 200)
ho_predict <- predict(wheat_mod[["model"]], newdata = HO, type = "response")
rmse <- RMSE(HO$Yield, ho_predict)
rmse_null <- RMSE(HO$Yield, mean(HO$Yield, na.rm = T))
wheat_mod[["RMSE"]] <- rmse
wheat_mod[["RMSENull"]] <- rmse_null

wwheat <- wwheat %>% filter(!is.na(Yield))
random_rn <- sample(nrow(wwheat), ceiling(nrow(wwheat)*.25))
HO <- wwheat[random_rn,]
sample <- wwheat[-random_rn,]
fm <- "Yield ~ fp(GDD) + fp(SDD) + fp(TP) + fp(pH_sur) +AWC_sur + fp(Slope, df=2) + YEAR"
wwheat_mod <- model_results(df = sample, fm = fm)
wwheat_data <- predict_data_PRISM(wwheat_mod, voi, res = 200)
ho_predict <- predict(wwheat_mod[["model"]], newdata = HO, type = "response")
rmse <- RMSE(HO$Yield, ho_predict)
rmse_null <- RMSE(HO$Yield, mean(HO$Yield, na.rm = T))
wwheat_mod[["RMSE"]] <- rmse
wwheat_mod[["RMSENull"]] <- rmse_null

saveRDS(corn_data,"./out/model_results/corn_PRISM.RDS")
saveRDS(cotton_data,"./out/model_results/cotton_PRISM.RDS")
saveRDS(soy_data,"./out/model_results/soy_PRISM.RDS")
saveRDS(wheat_data,"./out/model_results/wheat_PRISM.RDS")
saveRDS(wwheat_data,"./out/model_results/wwheat_PRISM.RDS")

saveRDS(corn_mod,"./out/model_results/corn_mod_PRISM.RDS")
saveRDS(cotton_mod,"./out/model_results/cotton_mod_PRISM.RDS")
saveRDS(soy_mod,"./out/model_results/soy_mod_PRISM.RDS")
saveRDS(wheat_mod,"./out/model_results/wheat_mod_PRISM.RDS")
saveRDS(wwheat_mod,"./out/model_results/wwheat_mod_PRISM.RDS")

###################################################################################################################
# Load MFP results
###################################################################################################################

corn_data <- readRDS("./out/model_results/corn_PRISM.RDS")
cotton_data <- readRDS("./out/model_results/cotton_PRISM.RDS")
soy_data <- readRDS("./out/model_results/soy_PRISM.RDS")
wheat_data <- readRDS("./out/model_results/wheat_PRISM.RDS")
wwheat_data <- readRDS("./out/model_results/wwheat_PRISM.RDS")

corn_mod <- readRDS("./out/model_results/corn_mod_PRISM.RDS")
cotton_mod <- readRDS("./out/model_results/cotton_mod_PRISM.RDS")
soy_mod <- readRDS("./out/model_results/soy_mod_PRISM.RDS")
wheat_mod <- readRDS("./out/model_results/wheat_mod_PRISM.RDS")
wwheat_mod <- readRDS("./out/model_results/wwheat_mod_PRISM.RDS")

corn_future <- predict_yield_int(corn_mod, "corn", yearc = F)
cotton_future <- predict_yield_int(cotton_mod, "cotton", yearc=F)
soy_future <- predict_yield_int(soy_mod, "soy", yearc=F)
wwheat_future <- predict_yield_int(wwheat_mod, "wwheat", yearc=F)
wheat_future <- predict_yield_int(wheat_mod, "wheat", yearc=F)

###################################################################################################################
# Visualize historical response curves
###################################################################################################################

plot_data <- soy_data
soy_gdd <- plot_crop_response(plot_data, "GDD")
soy_sdd <- plot_crop_response(plot_data, "SDD")
soy_tp <- plot_crop_response(plot_data, "TP")
soy_pHsur <- plot_crop_response(plot_data, "pH_sur")
soy_AWC <- plot_crop_response(plot_data, "AWC_sur")
soy_Slope <- plot_crop_response(plot_data, "Slope")
soy_Y <- plot_crop_response(plot_data, "YEAR")

plot_data <- corn_data
corn_gdd <- plot_crop_response(plot_data, "GDD")
corn_sdd <- plot_crop_response(plot_data, "SDD")
corn_tp <- plot_crop_response(plot_data, "TP")
#corn_efp <- plot_crop_response(plot_data, "EfP")
#corn_exp <- plot_crop_response(plot_data, "ExP")
corn_pHsur <- plot_crop_response(plot_data, "pH_sur")
corn_AWC <- plot_crop_response(plot_data, "AWC_sur")
corn_Slope <- plot_crop_response(plot_data, "Slope")
corn_Y <- plot_crop_response(plot_data, "YEAR")

plot_data <- cotton_data
cotton_gdd <- plot_crop_response(plot_data, "GDD")
cotton_sdd <- plot_crop_response(plot_data, "SDD")
#cotton_efp <- plot_crop_response(plot_data, "EfP")
#cotton_exp <- plot_crop_response(plot_data, "ExP")
cotton_tp <- plot_crop_response(plot_data, "TP")
cotton_pHsur <- plot_crop_response(plot_data, "pH_sur")
cotton_AWC <- plot_crop_response(plot_data, "AWC_sur")
cotton_Slope <- plot_crop_response(plot_data, "Slope")
cotton_Y <- plot_crop_response(plot_data, "YEAR")

plot_data <- wheat_data
wheat_gdd <- plot_crop_response(plot_data, "GDD")
wheat_sdd <- plot_crop_response(plot_data, "SDD")
#wheat_efp <- plot_crop_response(plot_data, "EfP")
#wheat_exp <- plot_crop_response(plot_data, "ExP")
wheat_tp <- plot_crop_response(plot_data, "TP")
wheat_pHsur <- plot_crop_response(plot_data, "pH_sur")
wheat_AWC <- plot_crop_response(plot_data, "AWC_sur")
wheat_Slope <- plot_crop_response(plot_data, "Slope")
wheat_Y <- plot_crop_response(plot_data, "YEAR")

plot_data <- wwheat_data
wwheat_gdd <- plot_crop_response(plot_data, "GDD")
wwheat_sdd <- plot_crop_response(plot_data, "SDD")
#wwheat_efp <- plot_crop_response(plot_data, "EfP")
#wwheat_exp <- plot_crop_response(plot_data, "ExP")
wwheat_tp <- plot_crop_response(plot_data, "TP")
wwheat_pHsur <- plot_crop_response(plot_data, "pH_sur")
wwheat_AWC <- plot_crop_response(plot_data, "AWC_sur")
wwheat_Slope <- plot_crop_response(plot_data, "Slope")
wwheat_Y <- plot_crop_response(plot_data, "YEAR")

gdd_all <- plot_all_response(v = "GDD")
sdd_all <- plot_all_response(v = "SDD")
#efp_all <- plot_all_response(v = "EfP")
#exp_all <- plot_all_response(v = "ExP")
tp_all <- plot_all_response(v = "TP")
year_all <- plot_all_response(v = "YEAR")

