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

library(lubridate)

ss <- readRDS("./data/coi_stopstart.RDS")
roi <- unique(coi$GEOID)

tmax <- readRDS("C:/Users/A02256433/Desktop/Data/PRISM/daily_county_TMAX.RDS")  # PRISM data folder
colnames(tmax) <- c("GEOID", "YEAR", "MONTH", "DAY", "TMAX")
tmax <- tmax %>% filter(GEOID %in% roi)
tmax_ss <- merge(tmax, ss, by = "GEOID", all = T)
comment(tmax_ss) <- "TMAX"
tmax_ss$DOY <- yday(as.Date(paste(tmax_ss$YEAR,tmax_ss$MONTH,tmax_ss$DAY,sep="-"))) 

ppt <-readRDS("C:/Users/A02256433/Desktop/Data/PRISM/daily_county_PPT.RDS")
colnames(ppt) <- c("GEOID", "YEAR", "MONTH", "DAY", "PPT")
ppt <- ppt %>% filter(GEOID %in% roi)
ppt_ss <- merge(ppt, ss, by = "GEOID", all = T)
comment(ppt_ss) <- "PPT"
ppt_ss$DOY <- yday(as.Date(paste(ppt_ss$YEAR,ppt_ss$MONTH,ppt_ss$DAY,sep="-"))) 

saveRDS(tmax_ss, "./data/prism_tmax_roi.RDS")
saveRDS(ppt_ss, "./data/prism_ppt_roi.RDS")

#################################################################################################################################
# Extract PRISM to county
#################################################################################################################################

# see extract_PRISM_FINAL.R - run on SESYNC cluster

#################################################################################################################################
# Yield
#################################################################################################################################

clist <- Sys.glob("./data/yield/*.csv")

allcrops <- list()
for (i in 1:length(clist)) {
  c <- read.csv(clist[i], stringsAsFactors = F)
  c <- c %>% 
    filter(!is.na(County.ANSI)) %>% 
    dplyr::select(Year, State, State.ANSI, Ag.District, Ag.District.Code, County, County.ANSI, Commodity, Value, Data.Item) %>%
    mutate(State.ANSI = str_pad(State.ANSI, width = 2, side = "left", pad = "0"),
           County.ANSI = str_pad(County.ANSI, width = 3, side = "left", pad = "0"),
           GEOID = paste0(State.ANSI, County.ANSI))
  allcrops[[i]] <- c
}
allcrops <- do.call(rbind, allcrops)
allcrops$Value[allcrops$Value == "999"] <- "NA"
allcrops$Value[allcrops$Value == "0"] <- "NA"  
allcrops$Value[allcrops$Value == "(D)"] <- "NA"  
allcrops$Value <- as.numeric(gsub(",", "", allcrops$Value))  
allcrops$GEOID <- as.factor(allcrops$GEOID)
allcrops$State <- as.factor(allcrops$State)
allcrops$State.ANSI <- as.factor(allcrops$State.ANSI)
allcrops$Ag.District <- as.factor(allcrops$Ag.District)
allcrops$Ag.District.Code <- as.factor(allcrops$Ag.District.Code)
allcrops$County <- as.factor(allcrops$County)
allcrops$County.ANSI <- as.factor(allcrops$County.ANSI)
# remove duplicates from spreadsheet issues/overlap
allcrops <- distinct(allcrops)
allcrops$Commodity[allcrops$Data.Item == "WHEAT, SPRING, (EXCL DURUM) - YIELD, MEASURED IN BU / ACRE"] <- "WHEAT"
allcrops$Commodity[allcrops$Data.Item == "WHEAT, WINTER - YIELD, MEASURED IN BU / ACRE"] <- "WWHEAT"
allcrops$Commodity <- as.factor(allcrops$Commodity)
yield <- allcrops %>% select(Year, GEOID, Commodity, Value)
colnames(yield) <- c("YEAR", "GEOID", "CROP", "YIELD")
saveRDS(yield, "./data/yield/yield.RDS")
