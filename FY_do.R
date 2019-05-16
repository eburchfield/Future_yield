source("FY_load.R")
source("FY_func.R")

# re-run historical models with scaling vector, linear year, perhaps center on scale paraneter with lower start value for scale vector - play with this
# negative slopes in the future?  this may mean 85 is better b/c of CO2 fertilization 
# 0.92, what's direction of scaling (tech trend?) - think through what tscale, tyear plots mean
# finalize historical models and clean up scripts/data outputs to this point
# update predictions with scaled data vector - may be applied post models, think on this
# generate future prediction maps and plots


###################################################################################################################
# Create data
###################################################################################################################

# construct_historical_data()
# add_hist_co2()
# build_ensembles()
# construct_future_data()
# add_future_co2()

###################################################################################################################
# Load historical data
###################################################################################################################

# FY_sensitivity testing documents model construction
# Response_curves_C02 documents model runs including CO2 predictor

corn <- readRDS("./out/panels/corn_PRISM_30.RDS")
cotton <- readRDS("./out/panels/cotton_PRISM_30.RDS")
cotton$YIELD <- cotton$YIELD/32 # from pounds/acre to bushels/acre
cotton <- gam_data_prep(cotton)
cotton <- cotton %>% filter(SDD < 100)
soy <- readRDS("./out/panels/soy_PRISM_30.RDS")
wwheat <- readRDS("./out/panels/wwheat_PRISM_30.RDS")
wheat <- readRDS("./out/panels/wheat_PRISM_30.RDS")


###################################################################################################################
# Run GAM models
###################################################################################################################

set.seed(0)
library(mgcv)

model_name = "gamly"
bs <- "ps" 
fm <- as.formula("YIELD ~ s(GDD, bs=bs) + s(SDD, bs=bs)+ s(ExP, bs=bs) + s(EfP, bs=bs) + YEAR + s(GEOID, bs='re')")

corn_mod <- run_gam(formula = fm, model_name = model_name, crop = "corn")
soy_mod <- run_gam(formula = fm, model_name = model_name, crop = "soy")
cotton_mod <- run_gam(formula = fm, model_name = model_name, crop = "cotton")
gam <- bam(fm, data=cotton, family = "gaussian", method="REML") # deal with extreme SDD values
saveRDS(gam, "./out/gam/gamly_cotton.RDS")
wheat_mod <- run_gam(formula = fm, model_name = model_name, crop = "wheat")
wwheat_mod <- run_gam(formula = fm, model_name = model_name, crop = "wwheat")

c <- run_gam_RMSE(fm, "corn") # fm, model_name, crop
c
# 17.45, 36.80
co <- run_gam_RMSE(fm, "cotton") # fm, model_name, crop
co
# 145.18, 234.51
s <- run_gam_RMSE(fm, "soy") # fm, model_name, crop
s
# 5.28, 10.37
w <- run_gam_RMSE(fm, "wheat") # fm, model_name, crop
w
#7.37, 10.31
ww <- run_gam_RMSE(fm, "wwheat") # fm, model_name, crop
ww
# 8.46, 13.92

###################################################################################################################
# Run GAM projections
###################################################################################################################

bs <- "ps" 
fm <- as.formula("YIELD ~ s(GDD, bs=bs) + s(SDD, bs=bs)+ s(ExP, bs=bs) + s(EfP, bs=bs) + offset(I(scale*YEAR)) + s(GEOID, bs='re')")

# High
scenario <- "High"
model_name <- "gamly_High"
scale <- scenario_dvt("corn")[,scenario]
corn_mod <- run_gam(formula = fm, model_name = model_name, crop = "corn")
scale <- scenario_dvt("soy")[,scenario]
soy_mod <- run_gam(formula = fm, model_name = model_name, crop = "soy")
scale <- scenario_dvt("cotton")[,scenario]
cotton_mod <- run_gam(formula = fm, model_name = model_name, crop = "cotton")
scale <- scenario_dvt("wheat")[,scenario]
wheat_mod <- run_gam(formula = fm, model_name = model_name, crop = "wheat")
scale <- scenario_dvt("wwheat")[,scenario]
wwheat_mod <- run_gam(formula = fm, model_name = model_name, crop = "wwheat")

# Average
scenario <- "Average"
model_name <- "gamly_Average"
scale <- scenario_dvt("corn")[,scenario]
corn_mod <- run_gam(formula = fm, model_name = model_name, crop = "corn")
scale <- scenario_dvt("soy")[,scenario]
soy_mod <- run_gam(formula = fm, model_name = model_name, crop = "soy")
scale <- scenario_dvt("cotton")[,scenario]
cotton_mod <- run_gam(formula = fm, model_name = model_name, crop = "cotton")
scale <- scenario_dvt("wheat")[,scenario]
wheat_mod <- run_gam(formula = fm, model_name = model_name, crop = "wheat")
scale <- scenario_dvt("wwheat")[,scenario]
wwheat_mod <- run_gam(formula = fm, model_name = model_name, crop = "wwheat")

# Low
scenario <- "Low"
model_name <- "gamly_Low"
scale <- scenario_dvt("corn")[,scenario]
corn_mod <- run_gam(formula = fm, model_name = model_name, crop = "corn")
scale <- scenario_dvt("soy")[,scenario]
soy_mod <- run_gam(formula = fm, model_name = model_name, crop = "soy")
scale <- scenario_dvt("cotton")[,scenario]
cotton_mod <- run_gam(formula = fm, model_name = model_name, crop = "cotton")
scale <- scenario_dvt("wheat")[,scenario]
wheat_mod <- run_gam(formula = fm, model_n5212 Old Main Hillame = model_name, crop = "wheat")
scale <- scenario_dvt("wwheat")[,scenario]
wwheat_mod <- run_gam(formula = fm, model_name = model_name, crop = "wwheat")

# None
scenario <- "None"
model_name <- "gamly_None"
scale <- scenario_dvt("corn")[,scenario]
corn_mod <- run_gam(formula = fm, model_name = model_name, crop = "corn")
scale <- scenario_dvt("soy")[,scenario]
soy_mod <- run_gam(formula = fm, model_name = model_name, crop = "soy")
scale <- scenario_dvt("cotton")[,scenario]
cotton_mod <- run_gam(formula = fm, model_name = model_name, crop = "cotton")
scale <- scenario_dvt("wheat")[,scenario]
wheat_mod <- run_gam(formula = fm, model_name = model_name, crop = "wheat")
scale <- scenario_dvt("wwheat")[,scenario]
wwheat_mod <- run_gam(formula = fm, model_name = model_name, crop = "wwheat")

###################################################################################################################
# GAM scenarios
###################################################################################################################

corn <- all_maps("corn")
#cotton <- all_maps("cotton")
soy <- all_maps("soy")
#wheat <- all_maps("wheat")
wwheat <- all_maps("wwheat")

###################################################################################################################
# GAM prediction
###################################################################################################################


# create visuals











