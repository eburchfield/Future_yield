library(tidyverse)
library(sp)
library(raster)
library(spdplyr)
library(visreg)
library(mgcv)
library(lme4)
library(gstat)
library(rasterVis)
select <- dplyr::select

# PRISM data constructed in stack_data.R on D drive and run using the extract_panels_GDD script in SESYNC cluster
# Shapefiles, SA_counties_analyzed original dataset, updated with season stop start as coi_startstop.RDS

coi <- readRDS("./data/coi_stopstart.RDS")
county <- readRDS("./data/county.RDS")
coist <- readRDS("./data/states.RDS")

################################################################################################################################
# Station coordinates to shapefile
################################################################################################################################

coords <- read.table("./data/latlon.txt")
coords <- rownames_to_column(coords, "ID")
colnames(coords) <- c("ID", "LAT", "LON")
 
# coords_sp <- SpatialPointsDataFrame(coords = cbind(coords[,"LON"], coords[,"LAT"]),
#                                     data = coords, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# coords_sp <- coords_sp[coi,]
# coords_cty <- over(coords_sp, coi[,"GEOID"], returnList = F) # link coordinates to counties
# coords_cty$ID <- row.names(coords_cty)
# coords_sp <- merge(coords_sp, coords_cty, by.x = "ID", by.y = "ID")
# saveRDS(coords_sp, "./data/coords_sp.RDS")

coords_sp <- readRDS("./data/coords_sp.RDS")

################################################################################################################################
# Historical station data 
################################################################################################################################

# rain <- read.table("./data/master_prcp_214.txt")
# comment(rain) <- "rain"
# tmax <- read.table("./data/master_tmax_214.txt")
# comment(tmax) <- "tmax"

################################################################################################################################
# Irrigation subset
################################################################################################################################

irr <- readRDS("./data/irr.RDS") %>%
  filter(GEOID %in% unique(coi$GEOID)) %>%
  mutate(PERC = IRRIGATED_ACRES/TOTAL_ACRES) %>%
  filter(PERC < 0.5)
irr_list <- unique(irr$GEOID)


################################################################################################################################
# GDD parameters
################################################################################################################################

# gdd range, from p. 8 in Neil's thesis
gdd_range <- list()
gdd_range[["corn"]] <- c(10,30)
gdd_range[["soy"]] <- c(10, 30)
gdd_range[["wheat"]] <- c(0,35)
gdd_range[["wwheat"]] <- c(0, 30)
gdd_range[["cotton"]] <- c(15.6, 37.8)

################################################################################################################################
# MFP visualization parameters
################################################################################################################################

var_range <- list()
var_range[["GDD"]] <- c(0, 5000)
var_range[["SDD"]] <- c(0, 600)
var_range[["EfP"]] <- c(0, 1200)
var_range[["ExP"]] <- c(0, 500)
var_range[["Year"]] <- c(1981, 2010)
var_range[["AWC_sur"]] <- c(0,0.4)
var_range[["pH_sur"]] <- c(1,8)
var_range[["Slope"]] <- c(0,30)
var_range[["SDI"]] <- c(0, 14)

################################################################################################################################
# Historical PRISM data
################################################################################################################################

# ppt <- readRDS("./data/prism_ppt_roi.RDS")
# tmax <- readRDS("./data/prism_tmax_roi.RDS")

################################################################################################################################
# CO2 scaling from Attavanich
################################################################################################################################

# from 1990 - 2010
# amount of tech trend explained by tech
co2 <- list()
co2[["corn"]] <- c(.08)
co2[["soy"]] <- c(.13)
co2[["wheat"]] <- c(.15)
co2[["wwheat"]] <- c(.15)
co2[["cotton"]] <- c(.34)

# 7-22, 4-47, 5-26, 65-96, and 3-35 % for yields of corn, sorghum, soybeans, cotton, and wheat, respectively.

################################################################################################################################
# Model results
################################################################################################################################

# corn <- readRDS("./out/panels/corn_PRISM.RDS")
# cotton <- readRDS("./out/panels/cotton_PRISM.RDS")
# cotton$Yield <- cotton$Yield/32 # from pounds/acre to bushels/acre
# soy <- readRDS("./out/panels/soy_PRISM.RDS")
# wwheat <- readRDS("./out/panels/wwheat_PRISM.RDS")
# corn$Year <- corn$YEAR
# cotton$Year <- cotton$YEAR
# soy$Year <- soy$YEAR
# wwheat$Year <- wwheat$YEAR
# 
# corn_data <- readRDS("./out/model_results/corn_PRISM.RDS")
# cotton_data <- readRDS("./out/model_results/cotton_PRISM.RDS")
# soy_data <- readRDS("./out/model_results/soy_PRISM.RDS")
# wwheat_data <- readRDS("./out/model_results/wwheat_PRISM.RDS")
# corn_mod <- readRDS("./out/model_results/corn_mod_PRISM.RDS")
# cotton_mod <- readRDS("./out/model_results/cotton_mod_PRISM.RDS")
# soy_mod <- readRDS("./out/model_results/soy_mod_PRISM.RDS")
# wwheat_mod <- readRDS("./out/model_results/wwheat_mod_PRISM.RDS")

