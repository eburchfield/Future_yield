source("FY_load.R")
source("FY_func.R")

# This scripts documents any modifications made to input datasets used in the func and do scripts

## COI/season stop/start processing
# coi <- readRDS("./data/SA_counties_analyzed.RDS")  
# coi <- spTransform(coi, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
# coi <- coi[,c("study_ar_5", "State", "County", "study_ar_3", "study_ar_4")]
# colnames(coi@data) <- c("GEOID", "STATE", "COUNTY", "STATEFP", "COUNTYFP")
# corn_plant <- raster("./data/corn.plant.start.asc")
# corn_harvest <- raster("./data/corn.harvest.end.asc")
# cotton_plant <- raster("./data/cotton.plant.start.asc")
# cotton_harvest <- raster("./data/cotton.harvest.end.asc")
# soy_plant <- raster("./data/soy.plant.start.asc")
# soy_harvest <- raster("./data/soy.harvest.end.asc")
# wheat_plant <- raster("./data/wheat.plant.start.asc")
# wheat_harvest <- raster("./data/wheat.harvest.end.asc")
# wwheat_plant <- raster("./data/wwheat.plant.start.asc")
# wwheat_harvest <- raster("./data/wwheat.harvest.end.asc")
# stop_start <- brick(corn_plant, corn_harvest, 
#                     cotton_plant, cotton_harvest, 
#                     soy_plant, soy_harvest, 
#                     wheat_plant, wheat_harvest,
#                     wwheat_plant, wwheat_harvest)
# coi <- extract(stop_start, coi, fun=mean, sp=T)
# saveRDS(coi, "./data/coi_stopstart.RDS")

# add FIPS to wheat control data
# wheat <- read.csv("./data/controls/wheat_aggregated.csv")
# wwheat <- read.csv("./data/controls/wwheat_aggregated.csv")
# wwheat <- wwheat %>% filter(Year == 1980) %>% select(COUNTY_State, FIPS)
# wheat <- merge(wwheat, wheat, by = "COUNTY_State", all.x=T)

################################################################################################################################
# PRISM functions for historical data - making tmean and ppt from original list of data from raster extract
################################################################################################################################

roi <- unique(coi$GEOID)

tmax <- readRDS("./data/tmax_list.RDS")
tmax <- do.call(rbind, tmax)
tmax <- tmax %>% filter(GEOID %in% roi)  %>% filter(YEAR < 2011) # subset to only include counties in ROI
tmax <- merge(tmax, coi, by.x = "GEOID", by.y = "GEOID")  # CI
tmax <- tmax %>% mutate(TMAX = MEAN) %>% dplyr::select(-MEAN)
comment(tmax) <- "TMAX"

ppt <-readRDS("./data/ppt_list.RDS")
ppt <- ppt %>% filter(GEOID %in% roi) %>% filter(YEAR < 2011)
ppt <- merge(ppt, coi, by.x = "GEOID", by.y = "GEOID") # CI

saveRDS(tmax, "./data/prism_tmax_roi.RDS")
saveRDS(ppt, "./data/prism_ppt_roi.RDS")

#################################################################################################################################
# Extract PRISM to county
#################################################################################################################################

# see extract_PRISM_FINAL.R - run on SESYNC cluster
