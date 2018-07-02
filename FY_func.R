library(tidyverse)
select <- dplyr::select
library(stringr)
library(lubridate)
library(mfp)
library(spdplyr)
library(gstat)
library(rasterVis)

################################################################################################################################
# Climate metrics
################################################################################################################################

gdd <- function(df, crop, gdd_range) {
  gdd <- ifelse(df$TMAX<= gdd_range[[crop]][1], 0, df$TMAX)
  gdd <- ifelse(df$TMAX >= gdd_range[[crop]][1] & df$TMAX <= gdd_range[[crop]][2], df$TMAX - gdd_range[[crop]][1], gdd)
  gdd <- ifelse(df$TMAX >= gdd_range[[crop]][2], gdd_range[[crop]][2] - gdd_range[[crop]][1], gdd) 
  return(gdd)
}

sdd <- function(df, crop, gdd_range) {
  sdd <- ifelse(df$TMAX >= gdd_range[[crop]][2], df$TMAX - gdd_range[[crop]][2], 0)
  return(sdd)
}

efp <- function(df, th) {
  efp <- ifelse(df$PPT >= th, th, df$PPT)
  return(efp)
}

add_controls <- function(crop) {
  # ch, including NASS QuickStats values
  controls <- read.csv(paste0("./data/controls/", crop, "_aggregated.csv"), stringsAsFactors = F)
  controls$GEOID <- str_pad(controls$FIPS, 5, pad = "0")
  controls <- controls %>%
    select("GEOID", "Year", "Slope", "pH_sur", "Yield", "AWC_sur")
  return(controls)
}

################################################################################################################################
# Historical PRISM
################################################################################################################################

time_subset_PRISM <- function(df, voi, crop)  {
  
  sta_vn <- paste0(crop, ".plant.start")
  sto_vn <- paste0(crop, ".harvest.end")
  if (crop == "wwheat") {
    df[,voi][df$DOY >= df[,sto_vn] & df$DOY <= df[,sta_vn]] <- NA  # planted in fall, harvested in summer/fall
  }
  else {
    df[,voi][df$DOY <= df[,sta_vn]] <- NA
    df[,voi][df$DOY >= df[,sto_vn]] <- NA
  }
  return(df)
}

time_subset <- function(df, crop)  {
  sta_vn <- paste0(crop, ".plant.start")
  sto_vn <- paste0(crop, ".harvest.end")
  
  if (crop == "wwheat") {
    df$measurement[df$DOY > df[,sto_vn] & df$DOY < df[,sta_vn]] <- NA  # planted in fall, harvested in summer/fall
  }
  else {
    df$measurement[df$DOY < df[,sta_vn]] <- NA
    df$measurement[df$DOY > df[,sto_vn]] <- NA
  }
  return(df)
}

prep_climate_PRISM <- function(crop, dfl) {
  
  df <- dfl[[1]]
  df$FID <- paste0(df$GEOID, df$DOY, df$YEAR)
  var_name <- comment(df)
  df <- time_subset_PRISM(df, var_name, crop)
  rf <- df %>% select(PPT, FID)
  
  df <- dfl[[2]]
  df$FID <- paste0(df$GEOID, df$DOY, df$YEAR)
  var_name <- comment(df)
  temp <- time_subset_PRISM(df, var_name, crop)

  final_df <- merge(rf, temp, by = "FID")
  comment(final_df) <- crop
  final_df$gdd <- gdd(final_df, crop, gdd_range[crop])
  final_df$sdd <- sdd(final_df, crop, gdd_range[crop])
  final_df$efp <- efp(final_df, th = 20)
  return(final_df)
}

seasonal_panel_PRISM <- function(df, crop) {
  
    season <- df %>% group_by(YEAR, GEOID) %>% 
      summarize(GDD = sum(gdd, na.rm=T), SDD = sum(sdd, na.rm=T), 
                TP = sum(PPT, na.rm=T), EfP = sum(efp, na.rm = T)) %>%
      mutate(ExP = TP - EfP)
  season <- as.data.frame(season)
  controls <- add_controls(crop)
  controls <- controls %>% filter(Year < 2011 & Year > 1981)
  controls$YEAR <- controls$Year
  season <- merge(season, controls, by = c("YEAR", "GEOID")) 
  comment(season) <- crop
  return(season)
}

construct_hist_data_PRISM <- function() {
  
  corn_daily <- prep_climate_PRISM(crop = "corn", dfl = list(ppt, tmax)) #1 min
  wheat_daily <- prep_climate_PRISM(dfl = list(ppt, tmax), crop = "wheat")
  wwheat_daily <- prep_climate_PRISM(dfl = list(ppt, tmax), crop = "wwheat")
  cotton_daily <- prep_climate_PRISM(dfl = list(ppt, tmax), crop = "cotton")
  soy_daily <- prep_climate_PRISM(dfl = list(ppt, tmax), crop = "soy")
  
  corn_panel <- seasonal_panel_PRISM(corn_daily, "corn")
  wheat_panel <- seasonal_panel_PRISM(wheat_daily, "wheat")
  wwheat_panel <- seasonal_panel_PRISM(wwheat_daily, "wwheat")
  cotton_panel <- seasonal_panel_PRISM(cotton_daily, "cotton")
  soy_panel <- seasonal_panel_PRISM(soy_daily, "soy")
  
  saveRDS(corn_panel, "./out/panels/corn_PRISM.RDS")
  saveRDS(cotton_panel, "./out/panels/cotton_PRISM.RDS")
  saveRDS(soy_panel, "./out/panels/soy_PRISM.RDS")
  saveRDS(wheat_panel, "./out/panels/wheat_PRISM.RDS")
  saveRDS(wwheat_panel, "./out/panels/wwheat_PRISM.RDS")
  
}

################################################################################################################################
# MFP Model and Prediction
################################################################################################################################

model_results <- function(df, fm) {
  model_list <- list()
  mod <- mfp(formula = as.formula(fm),
             data = df, 
             family = gaussian, 
             alpha = 0.05, # p-values for testing between FP models
             select=0.05, # p-values for selection on each parameter
             na.action = na.omit,
             maxits = 10,
             verbose=T,
             y=T, x=T, rescale=T)
  model_list[["data"]]<- df
  model_list[["model"]] <- mod
  model_list[["R2"]] <- round(1 - (mod$deviance/mod$null.deviance), 3)
  return(model_list)
}

predict_data_PRISM <- function(model, voi, res) {
  
  df <- model[["data"]]
  cm <- length(unique(df$Year))
  
  # filter out complete cases of missing data
  id_list <- df %>% group_by(GEOID) %>% 
    summarize(n = sum(is.na(Yield))) %>% 
    filter(n == cm) %>% select(GEOID) 
  df <- df %>% filter(!GEOID %in% id_list$GEOID) 
  
  final_list <- list()
  
  for (i in 1:length(voi)) {
    
    v <- voi[i] 
    
    dfmean <- t(as.data.frame(colMeans(df[,voi], na.rm=T))) 
#   dfmean[1,"YEAR"] <- 2010
    
    rownames(dfmean) <- NULL
    dfm <- as.data.frame(dfmean[rep(seq_len(nrow(dfmean)), res),])
    rownames(dfm) <- NULL
    
    x <- seq(quantile(df[,v], 0.01, na.rm = T), 
            quantile(df[,v], 0.99, na.rm = T), 
            length.out = res)
    dfm[,v] <- x 

    Yield <- predict(model[["model"]], newdata = dfm, type = "response", se.fit = T)
    Ymn <- as.data.frame(Yield[[1]])
    Ysd <- as.data.frame(Yield[[2]])
    Yield <- cbind(Ymn, Ysd)
    colnames(Yield) <- c(paste0("Predict", v), paste0("SD", v))
    dfm$Yield <- Yield[,paste0("Predict",v)]
    dfm$YieldSD <- Yield[,paste0("SD", v)]
    final_list[[v]] <- dfm
    }

  return(final_list)
  
}

RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}

################################################################################################################################
# Visualize results
################################################################################################################################

plot_crop_response <- function(out, v) {
  # ch
  voi_df <- out[[v]]
  voi_df[,v] <- as.factor(voi_df[,v])
  voi_mean <- voi_df %>% group_by(voi_df[,v]) %>% summarize(mYield = mean(Yield))  # summarize values across range of voi values
  colnames(voi_mean) <- c(v, "Yield")
  voi_mean <- as.data.frame(voi_mean)
  
  vec <- as.numeric(as.character(voi_mean[,v]))
  voi_mean[,v] <- vec

  p1 <- ggplot(data = voi_mean, aes(voi_mean[,v], Yield)) +
    geom_point() +
    xlab(v) +
    scale_x_continuous(limits = c(min(voi_mean[,v]), max(voi_mean[,v])))
  
  return(p1)
}

plot_all_response <- function(v) {
  
  corn <- cbind(corn_data[[v]], Crop = rep("Corn", nrow(corn_data[[v]])))
  corn <- corn %>% mutate(YieldZ = (Yield - mean(Yield, na.rm=T))/sd(Yield, na.rm=T))
  cotton <- cbind(cotton_data[[v]], Crop = rep("Cotton", nrow(cotton_data[[v]])))
  cotton <- cotton %>% mutate(YieldZ = (Yield - mean(Yield, na.rm=T))/sd(Yield, na.rm=T))
  soy <- cbind(soy_data[[v]], Crop = rep("Soy", nrow(soy_data[[v]])))
  soy <- soy %>% mutate(YieldZ = (Yield - mean(Yield, na.rm=T))/sd(Yield, na.rm=T))  
  #wheat <- cbind(wheat_data[[v]], Crop = rep("Wheat", nrow(wheat_data[[v]])))
  #wheat <- wheat %>% mutate(YieldZ = (Yield - mean(Yield, na.rm=T))/sd(Yield, na.rm=T))
  wwheat <- cbind(wwheat_data[[v]], Crop = rep("Wwheat", nrow(wwheat_data[[v]])))
  wwheat <- wwheat %>% mutate(YieldZ = (Yield - mean(Yield, na.rm=T))/sd(Yield, na.rm=T))
  
  out <- rbind(corn, cotton, soy, wwheat)
  out[,v] <- as.factor(out[,v])
  voi_mean <- out %>% group_by(Crop, out[,v]) %>% summarize(mYield = mean(YieldZ))
  colnames(voi_mean) <- c("Crop", v, "YieldZ")
  #voi_mean[,"Yield"] <- exp(voi_mean[,"Yield"])
  voi_mean <- as.data.frame(voi_mean)
  
  vec <- as.numeric(as.character(voi_mean[,v]))
  voi_mean[,v] <- vec
  
  p1 <- ggplot(data = voi_mean, aes(voi_mean[,v], YieldZ, color = Crop)) +
    geom_point() +
    xlab(v) +
    ylab("Standardized Yield")
    scale_x_continuous(limits = c(min(voi_mean[,v]), max(voi_mean[,v])))
  
  return(p1)
}

################################################################################################################################
# Create future panels with station data
################################################################################################################################

build_ensembles <-function() {
  # across three models, get average for each RCP
  model_list <- list("./data/models/ipslcm5alr_", "./data/models/mricgcm3_", "./data/models/noresm1_" )
  climate_vars <- c("tmax", "pamt")
  
  for (v in 1:length(climate_vars)) {
    # average across three models for rcp45
    ips <- as.matrix(read.table(paste0(model_list[1], climate_vars[v], "_gen_rcp45_2006_2100.txt")))
    mri <- as.matrix(read.table(paste0(model_list[2], climate_vars[v], "_gen_rcp45_2006_2100.txt"))) 
    nor <- as.matrix(read.table(paste0(model_list[3], climate_vars[v], "_gen_rcp45_2006_2100.txt")))
    ens45 <- (ips[,4:217] + mri[,4:217] + nor[,4:217])/3
    ens45 <- cbind(nor[,1:3], ens45)
    comment(ens45) <- climate_vars[v]
    saveRDS(ens45, paste0("./out/ensembles/rcp45", climate_vars[v], "_ens_2006_2100.RDS"))
    
    # average across three models for rcp85
    ips <- as.matrix(read.table(paste0(model_list[1], climate_vars[v], "_gen_rcp85_2006_2100.txt")))
    mri <- as.matrix(read.table(paste0(model_list[2], climate_vars[v], "_gen_rcp85_2006_2100.txt"))) 
    nor <- as.matrix(read.table(paste0(model_list[3], climate_vars[v], "_gen_rcp85_2006_2100.txt")))
    ens85 <- (ips[,4:217] + mri[,4:217] + nor[,4:217])/3
    ens85 <- cbind(nor[,1:3], ens85)
    comment(ens85) <- climate_vars[v]
    saveRDS(ens85, paste0("./out/ensembles/rcp85", climate_vars[v], "_ens_2006_2100.RDS"))

    }
}

prep_future_climate <- function(crop, voi = c("PPT", "TMAX")) {
  
  rcp_list <- list("./out/ensembles/rcp45", "./out/ensembles/rcp85")
  
  all_models <- list()
  
  for (m in 1:length(rcp_list)) {
    dr <- rcp_list[m]
    final_df <- data.frame()
    
    for (v in 1:length(voi)) {
      var_name <- voi[v]
      df <- readRDS(paste0(dr, voi[v], "_ens_2006_2100.RDS"))
      df <- as.data.frame(df) %>% filter(V1 > 2009)
      colnames(df) <- c("Year", "Month", "Day", coords[,"ID"])  

      # tidy Schoof dataset
      df <- gather(df, station, measurement, 4:217)
      # merge to county shapefile and station shapefile
      df <- merge(df, coords_sp@data, by.x = "station", by.y = "ID")
      coi_crop <- coi[,c("GEOID", "STATE", "COUNTY", paste0(crop, ".plant.start"), paste0(crop,".harvest.end"))]
      df <- merge(df, coi_crop, by.x = "GEOID", by.y = "GEOID")
      
      df$DOY <- df$Day
      df <- time_subset(df, crop)
      names(df)[names(df) == 'measurement'] <- var_name
      
      ifelse(v == 1, final_df <- df, final_df <- cbind(final_df, df[,var_name]))
      names(final_df)[names(final_df) == 'df[, var_name]'] <- var_name
    }
    
    final_df$gdd <- gdd(final_df, crop, gdd_range)
    final_df$sdd <- sdd(final_df, crop, gdd_range)
    final_df$efp <- efp(final_df, th = 20)
    comment(final_df) <- paste0(rcp_list[m], crop)
    all_models[[m]] <-final_df
    
  }
  return(all_models)
}

seasonal_panel_future <- function(mlist, crop) {
  
  dr <-  c("./out/fpanels/rcp45_", "./out/fpanels/rcp85_")
  
  for (l in 1:length(mlist)) {
    df <- mlist[[l]]
    df$YEAR <- df$Year
    
    season <- df %>% group_by(YEAR, GEOID) %>% 
        summarize(GDD = sum(gdd, na.rm=T), SDD = sum(sdd, na.rm=T), 
                  TP = sum(PPT, na.rm=T), EfP = sum(efp, na.rm = T)) %>%
        mutate(ExP = TP - EfP)

    season <- as.data.frame(season)
    season$Year <- season$YEAR
    season <- season %>% select(-YEAR)
    comment(season) <- crop
    saveRDS(season, paste0(dr[l], crop, ".RDS"))
  }}

construct_future_data <- function() {
  
  ens_list <- c("./out/ensembles/rcp45_", "./out/ensembles/rcp85")
  voi <- c("pamt", "tmax")
  
  corn_daily <- prep_future_climate(crop = "corn") 
  wheat_daily <- prep_future_climate(crop = "wheat")
  wwheat_daily <- prep_future_climate(crop = "wwheat")
  cotton_daily <- prep_future_climate(crop = "cotton")
  soy_daily <- prep_future_climate(crop = "soy")
  
  corn_panel <- seasonal_panel_future(corn_daily, "corn")
  wheat_panel <- seasonal_panel_future(wheat_daily, "wheat")
  wwheat_panel <- seasonal_panel_future(wwheat_daily, "wwheat")
  cotton_panel <- seasonal_panel_future(cotton_daily, "cotton")
  soy_panel <- seasonal_panel_future(soy_daily, "soy")
}

###################################################################################################################
# Create future panels with interpolated version of construct_future_data results
###################################################################################################################

future_interpolation <- function(crop)   {
  # interpolate future daily climate data from 102 stations
  ens_list <-  c("./out/fpanels/rcp45_", "./out/fpanels/rcp85_")
  rcp <- c("45", "85")
  
  for (i in 1:2) {
    df <- as.data.frame(readRDS(paste0(ens_list[i], crop, ".RDS")))
    years <- unique(df$Year)
    grd <- SpatialPixels(SpatialPoints(makegrid(coi, n=7000)), proj4string = proj4string(coi))
    grd <- grd[coi,] 
    
    gdd <- stack()
    sdd <- stack()
    efp <- stack()
    exp <- stack()
    tp <- stack()
    
    for (y in 1:length(years)) {
      df_sub <- df %>% filter(Year == years[y])
      df_sp <- merge(coi, df_sub, by = "GEOID")
      
      pcentroids <- as.data.frame(coordinates(df_sp))
      colnames(pcentroids) <- c("lon", "lat") 
      pcentroids <- data.frame("ID" = 1:nrow(pcentroids), pcentroids)
      coordinates(pcentroids) <- c("lon", "lat") 
      proj4string(pcentroids) <- proj4string(df_sp) # assign projection
      pcentroids@data <- sp::over(x = pcentroids, y = df_sp, returnList = FALSE)
      
      c_gdd <- pcentroids %>% filter(!is.na(GDD))
      p.gdd <- idw(formula = GDD ~ 1, c_gdd, grd,  idp = 2, nmax = 12)
      p.gdd <- raster(p.gdd)
      comment(p.gdd) <- paste0("GDD", years[y])
      gdd <- stack(gdd, p.gdd)

      c_sdd <- pcentroids %>% filter(!is.na(SDD))
      p.sdd <- idw(formula = SDD ~ 1, c_sdd, grd,  idp = 2, nmax = 12)
      p.sdd <- raster(p.sdd)
      comment(p.sdd) <- paste0("SDD", years[y])
      sdd <- stack(sdd, p.sdd)

      c_efp <- pcentroids %>% filter(!is.na(EfP))
      p.efp <- idw(formula = EfP ~ 1, c_efp, grd,  idp = 2, nmax = 12)
      p.efp <- raster(p.efp)
      comment(p.efp) <- paste0("EfP", years[y])
      efp <- stack(efp, p.efp)

      c_exp <- pcentroids %>% filter(!is.na(ExP))
      p.exp <- idw(formula = ExP ~ 1, c_exp, grd,  idp = 2, nmax = 12)
      p.exp <- raster(p.exp)
      comment(p.exp) <- paste0("ExP", years[y])
      exp <- stack(exp, p.exp)

      c_tp <- pcentroids %>% filter(!is.na(TP))
      p.tp <- idw(formula = TP ~ 1, c_tp, grd,  idp = 2, nmax = 12)
      p.tp <- raster(p.tp)
      comment(p.tp) <- paste0("TP", years[y])
      tp <- stack(tp, p.tp)
    }
    saveRDS(gdd, paste0("./out/stacks/gdd_int_", rcp[i], "_", crop, ".RDS"))
    saveRDS(sdd, paste0("./out/stacks/sdd_int_", rcp[i], "_", crop, ".RDS"))
    saveRDS(exp, paste0("./out/stacks/exp_int_", rcp[i], "_", crop, ".RDS"))
    saveRDS(efp, paste0("./out/stacks/efp_int_", rcp[i], "_", crop, ".RDS"))
    saveRDS(tp, paste0("./out/stacks/tp_int_", rcp[i], "_", crop, ".RDS"))   
  }
}

# future_interpolation("corn")
# future_interpolation("soy")
# future_interpolation("wheat")
# future_interpolation("wwheat")
# future_interpolation("cotton")

future_panels <- function(crop) {
  # extracts interpolated daily 102 climate station data to county scale
  # takes a sec to run b/c of extract function
  rcp <- c("45", "85")
  gdd <- readRDS(paste0("./out/stacks/gdd_int_", rcp[1], "_", crop, ".RDS"))
  blank_ras <- gdd[[1]]
  blank_ras <- setValues(blank_ras, NA) 
  
  for (i in 1:2) {
    
    gdd <- readRDS(paste0("./out/stacks/gdd_int_", rcp[i], "_", crop, ".RDS"))
    sdd <- readRDS(paste0("./out/stacks/sdd_int_", rcp[i], "_", crop, ".RDS"))
    exp <- readRDS(paste0("./out/stacks/exp_int_", rcp[i], "_", crop, ".RDS"))
    efp <- readRDS(paste0("./out/stacks/efp_int_", rcp[i], "_", crop, ".RDS"))
    tp <- readRDS(paste0("./out/stacks/tp_int_", rcp[i], "_", crop, ".RDS"))
    
    gdd_p <- raster::extract(gdd, coi, fun = mean, df=T)  # extract average to county
    gdd_p$GEOID <- coi$GEOID
    colnames(gdd_p) <- c("ID", 2010:2100, "GEOID")
    gdd_p <- gather(gdd_p, Year, GDD, 2:92)

    sdd_p <- raster::extract(sdd, coi, fun = mean, df=T)
    colnames(sdd_p) <- c("ID", 2010:2100)
    sdd_p <- gather(sdd_p, Year, SDD, 2:92)

    efp_p <- raster::extract(efp, coi, fun = mean, df=T)
    colnames(efp_p) <- c("ID", 2010:2100)
    efp_p <- gather(efp_p, Year, EfP, 2:92)

    exp_p <- raster::extract(exp, coi, fun = mean, df=T)
    colnames(exp_p) <- c("ID", 2010:2100)
    exp_p <- gather(exp_p, Year, ExP, 2:92)
    
    tp_p <- raster::extract(tp, coi, fun = mean, df=T)
    colnames(tp_p) <- c("ID", 2010:2100)
    tp_p <- gather(tp_p, Year, TP, 2:92)
    
    out <- cbind(gdd_p %>% select(Year, GEOID, GDD), sdd_p %>% select(SDD), 
                 efp_p %>% select(EfP), exp_p %>% select(ExP), tp_p %>% select(TP))
    comment(out) <- crop
    saveRDS(out, paste0("./out/intpanels/rcp", rcp[i], "_", crop, ".RDS"))
  }
}

# future_panels("corn")
# future_panels("cotton")
# future_panels("soy")
# # future_panels("wheat")
# future_panels("wwheat")

add_future_controls <- function(crop) {
  
  # fl <- Sys.glob(paste0("./out/intpanels/rcp*", crop, ".RDS"))
  fl <- Sys.glob(paste0("./out/fpanels/rcp*", crop, ".RDS"))
  
  for (i in 1:length(fl)) {
    r <- readRDS(fl[i])
    c <- add_controls(crop)
    c <- c %>% 
      group_by(GEOID) %>% 
      summarize(Slope = mean(Slope, na.rm = T), pH_sur = mean(pH_sur, na.rm = T), AWC_sur = mean(AWC_sur, na.rm=T))
    fp <- merge(r, c, by = "GEOID", all.x = T)
    saveRDS(fp, fl[i])
  }
}

# add_future_controls("corn")
# add_future_controls("soy")
# add_future_controls("cotton")
# # add_future_controls("wheat")
# add_future_controls("wwheat")

###################################################################################################################
# Predict future yield
###################################################################################################################

predict_yield <- function(model, crop, yearc = F) {
  
  mod <- model[["model"]]
  ens_list <-  c("./out/fpanels/rcp45_", "./out/fpanels/rcp85_")
  
  yield_list <- list()
  
  for (i in 1:2) { 
    df <- as.data.frame(readRDS(paste0(ens_list[i], crop, ".RDS"))) #n - 102
    df$YEAR <- df$Year
    
    if (yearc == TRUE) {
      df$Year1 <- as.numeric(df$Year)
      df$YEAR <- 2010
    } else {
      dfYEAR <- as.numeric(df$Year)
    }
    
    Yield <- predict(mod, df, type = "response", se.fit=T)
    Ymn <- as.data.frame(Yield[[1]])
    Ysd <- as.data.frame(Yield[[2]])
    Yield <- cbind(Ymn, Ysd)
    colnames(Yield) <- c("Yield", "YieldSD")
    Yield <- cbind(Yield, df)
    yield_list[[i]] <- Yield
  }
  return(yield_list)
}

predict_yield_int <- function(model, crop, yearc = FALSE) {
  mod <- model[["model"]]
  in_data <- model[["data"]]
  in_data$Year <- in_data$YEAR
  ens_list <-  c("./out/intpanels/rcp45_", "./out/intpanels/rcp85_")
  yield_list <- list()
  
  for (i in 1:2) { 
    future_data <- as.data.frame(readRDS(paste0(ens_list[i], crop, ".RDS")))
    if (yearc == TRUE) {
      future_data$Year1 <- as.numeric(future_data$Year)
      future_data$YEAR <- 2010
    } else {
      future_data$YEAR <- as.numeric(future_data$Year)
    }
    
    Yield <- predict(mod, future_data, type = "response", se.fit=T)
    Ymn <- as.data.frame(Yield[[1]])
    Ysd <- as.data.frame(Yield[[2]])
    Yield <- cbind(Ymn, Ysd)
    colnames(Yield) <- c("Yield", "YieldSD")
    Yield <- cbind(Yield, future_data)
    yield_list[[i]] <- Yield
    
  }
  return(yield_list)
}

spatial_projections <- function(future_data, hist_data, notech = F) {
  fy <- future_data
  py <- hist_data %>% filter(Year %in% 2000:2010) %>% group_by(GEOID) %>% summarize(Yield = mean(Yield, na.rm=T))
  py0 <- hist_data %>% filter(Year %in% 1990:2000) %>% group_by(GEOID) %>% summarize(Yield = mean(Yield, na.rm=T))
  
  if (notech == T) {
    fy[[1]]$Year <- fy[[1]]$Year1
    fy[[2]]$Year <- fy[[2]]$Year1
  }
  fy30_45 <- fy[[1]] %>% filter(Year %in% 2030:2040) %>% group_by(GEOID) %>% 
    summarize(GDD = mean(GDD), SDD = mean(SDD), EfP = mean(EfP), ExP = mean(ExP), Yield = mean(Yield), YieldSD = mean(YieldSD)) %>% 
    mutate(Model = "RCP4.5", Decade = "2030s")
  comment(fy30_45) <- "2030s/RCP4.5"
  
  fy30_85 <- fy[[2]] %>% filter(Year %in% 2030:2040) %>% group_by(GEOID) %>% 
    summarize(GDD = mean(GDD), SDD = mean(SDD), EfP = mean(EfP), ExP = mean(ExP), Yield = mean(Yield), YieldSD = mean(YieldSD)) %>% 
    mutate(Model = "RCP8.5", Decade = "2030s")
  comment(fy30_85) <- "2030s/RCP8.5"
  
  fy60_45 <- fy[[1]] %>% filter(Year %in% 2060:2070) %>% group_by(GEOID) %>% 
    summarize(GDD = mean(GDD), SDD = mean(SDD), EfP = mean(EfP), ExP = mean(ExP), Yield = mean(Yield), YieldSD = mean(YieldSD)) %>% 
    mutate(Model = "RCP4.5", Decade = "2060s")
  comment(fy60_45) <- "2060s/RCP4.5"
  
  fy60_85 <- fy[[2]] %>% filter(Year %in% 2060:2070) %>% group_by(GEOID) %>% 
    summarize(GDD = mean(GDD), SDD = mean(SDD), EfP = mean(EfP), ExP = mean(ExP), Yield = mean(Yield), YieldSD = mean(YieldSD)) %>% 
    mutate(Model = "RCP8.5", Decade = "2060s")
  comment(fy60_85) <- "2060s/RCP8.5"
  
  fy90_45 <- fy[[1]] %>% filter(Year %in% 2090:2100) %>% group_by(GEOID) %>% 
    summarize(GDD = mean(GDD), SDD = mean(SDD), EfP = mean(EfP), ExP = mean(ExP), Yield = mean(Yield), YieldSD = mean(YieldSD)) %>% 
    mutate(Model = "RCP4.5", Decade = "2090s")
  comment(fy90_45) <- "2090s/RCP4.5"
  
  fy90_85 <- fy[[2]] %>% filter(Year %in% 2090:2100) %>% group_by(GEOID) %>% 
    summarize(GDD = mean(GDD), SDD = mean(SDD), EfP = mean(EfP), ExP = mean(ExP), Yield = mean(Yield), YieldSD = mean(YieldSD)) %>% 
    mutate(Model = "RCP8.5", Decade = "2090s")
  comment(fy90_85) <- "2090s/RCP8.5"
  
  rcp <- list(fy30_45, fy30_85, fy60_45, fy60_85, fy90_45, fy90_85)
  rcp_df <- rbind(fy30_45, fy30_85, fy60_45, fy60_85, fy90_45, fy90_85)
  grd <- SpatialPixels(SpatialPoints(makegrid(coi, n=7000)), proj4string = proj4string(coi))
  grd <- grd[coi,] 
  
  pys <- merge(coi, py, by = "GEOID")
  pcentroids <- as.data.frame(coordinates(pys))
  colnames(pcentroids) <- c("lon", "lat") 
  pcentroids <- data.frame("ID" = 1:nrow(pcentroids), pcentroids)
  coordinates(pcentroids) <- c("lon", "lat") 
  proj4string(pcentroids) <- proj4string(pys) # assign projection
  pcentroids@data <- sp::over(x = pcentroids, y = pys, returnList = FALSE)
  pcentroids <- pcentroids %>% filter(!is.na(Yield))
  p.idw <- idw(formula = Yield ~ 1, pcentroids, grd,  idp = 2, nmax = 12)
  p.idw <- raster(p.idw)
  comment(p.idw) <- "2000 - 2010"
  
  pys2 <- merge(coi, py0, by = "GEOID")
  pcentroids2 <- as.data.frame(coordinates(pys2))
  colnames(pcentroids2) <- c("lon", "lat") 
  pcentroids2 <- data.frame("ID" = 1:nrow(pcentroids2), pcentroids2)
  coordinates(pcentroids2) <- c("lon", "lat") 
  proj4string(pcentroids2) <- proj4string(pys2) # assign projection
  pcentroids2@data <- sp::over(x = pcentroids2, y = pys2, returnList = FALSE)
  pcentroids2 <- pcentroids2 %>% filter(!is.na(Yield))
  p.idw2 <- idw(formula = Yield ~ 1, pcentroids2, grd,  idp = 2, nmax = 12)
  p.idw2 <- raster(p.idw2)
  comment(p.idw2) <- "1990 - 2000"
  
  diff_brick <- stack(p.idw2, p.idw)
  
  for (i in 1:length(rcp)) {
    fy_sub <- rcp[[i]]
    fys <- merge(coi, fy_sub, by = "GEOID")
    centroids <- as.data.frame(coordinates(fys))
    colnames(centroids) <- c("lon", "lat") 
    centroids <- data.frame("ID" = 1:nrow(centroids), centroids)
    coordinates(centroids) <- c("lon", "lat") 
    proj4string(centroids) <- proj4string(fys) # assign projection
    centroids@data <- sp::over(x = centroids, y = fys, returnList = FALSE)
    centroids <- centroids %>% filter(!is.na(Yield))
    f.idw  <- idw(formula = Yield ~ 1, centroids, grd,  idp = 2, nmax = 12)
    f.idw <- raster(f.idw)
    comment(f.idw) <- comment(fy_sub)
    names(f.idw) <- comment(fy_sub)
    
    #diff <- raster(f.idw) - p.idw
    #comment(diff) <- comment(fy_sub)
    
    diff_brick <- stack(diff_brick , f.idw)
  }
  
  out <- list(rcp_df, diff_brick)
  return(out)
}

spatial_projections_int <- function(future_data, hist_data, notech = F) {
  # generates rasters for visualization from interpolated 102 station estimates
  
  fy <- future_data
  py <- hist_data %>% filter(GEOID %in% unique(fy[[1]]$GEOID)) %>% filter(Year %in% 2000:2010) %>% group_by(GEOID) %>% summarize(Yield = mean(Yield, na.rm=T))
  py0 <- hist_data %>% filter(GEOID %in% unique(fy[[1]]$GEOID)) %>% filter(Year %in% 1990:2000) %>% group_by(GEOID) %>% summarize(Yield = mean(Yield, na.rm=T))
  
  if (notech == T) {
    fy[[1]]$Year <- fy[[1]]$Year1
    fy[[2]]$Year <- fy[[2]]$Year1
  }
  
  fy30_45 <- fy[[1]] %>% filter(Year %in% 2030:2040) %>% group_by(GEOID) %>% 
    summarize(GDD = mean(GDD), SDD = mean(SDD), EfP = mean(EfP), ExP = mean(ExP), Yield = mean(Yield), YieldSD = mean(YieldSD)) %>% 
    mutate(Model = "RCP4.5", Decade = "2030s")
  comment(fy30_45) <- "2030s/RCP4.5"
  
  fy30_85 <- fy[[2]] %>% filter(Year %in% 2030:2040) %>% group_by(GEOID) %>% 
    summarize(GDD = mean(GDD), SDD = mean(SDD), EfP = mean(EfP), ExP = mean(ExP), Yield = mean(Yield), YieldSD = mean(YieldSD)) %>% 
    mutate(Model = "RCP8.5", Decade = "2030s")
  comment(fy30_85) <- "2030s/RCP8.5"
  
  fy60_45 <- fy[[1]] %>% filter(Year %in% 2060:2070) %>% group_by(GEOID) %>% 
    summarize(GDD = mean(GDD), SDD = mean(SDD), EfP = mean(EfP), ExP = mean(ExP), Yield = mean(Yield), YieldSD = mean(YieldSD)) %>% 
    mutate(Model = "RCP4.5", Decade = "2060s")
  comment(fy60_45) <- "2060s/RCP4.5"
  
  fy60_85 <- fy[[2]] %>% filter(Year %in% 2060:2070) %>% group_by(GEOID) %>% 
    summarize(GDD = mean(GDD), SDD = mean(SDD), EfP = mean(EfP), ExP = mean(ExP), Yield = mean(Yield), YieldSD = mean(YieldSD)) %>% 
    mutate(Model = "RCP8.5", Decade = "2060s")
  comment(fy60_85) <- "2060s/RCP8.5"
  
  fy90_45 <- fy[[1]] %>% filter(Year %in% 2090:2100) %>% group_by(GEOID) %>% 
    summarize(GDD = mean(GDD), SDD = mean(SDD), EfP = mean(EfP), ExP = mean(ExP), Yield = mean(Yield), YieldSD = mean(YieldSD)) %>% 
    mutate(Model = "RCP4.5", Decade = "2090s")
  comment(fy90_45) <- "2090s/RCP4.5"
  
  fy90_85 <- fy[[2]] %>% filter(Year %in% 2090:2100) %>% group_by(GEOID) %>% 
    summarize(GDD = mean(GDD), SDD = mean(SDD), EfP = mean(EfP), ExP = mean(ExP), Yield = mean(Yield), YieldSD = mean(YieldSD)) %>% 
    mutate(Model = "RCP8.5", Decade = "2090s")
  comment(fy90_85) <- "2090s/RCP8.5"
  
  rcp <- list(fy30_45, fy30_85, fy60_45, fy60_85, fy90_45, fy90_85)
  rcp_df <- rbind(fy30_45, fy30_85, fy60_45, fy60_85, fy90_45, fy90_85)
  r <- raster(ncol=150, nrow=150)
  
  pys <- merge(coi, py, by="GEOID", all.x = T)
  py0s <- merge(coi, py0, by = "GEOID", all.x = T)
  extent(r) <- extent(pys)
  pyr <- rasterize(pys, r, 'Yield')
  comment(pyr) <- comment(pys)
  names(pyr) <- comment(pys)
  pyr0 <- rasterize(py0s, r, 'Yield')
  comment(pyr0) <- comment(py0s)
  names(pyr0) <- comment(py0s)
  
  diff_brick <- stack(pyr, pyr0)
  
  for (i in 1:length(rcp)) {
    fy_sub <- rcp[[i]]
    fys <- merge(coi, fy_sub, by = "GEOID", all.x = T)
    extent(r) <- extent(fys)
    rp <- rasterize(fys, r, 'Yield')
    comment(rp) <- comment(fy_sub)
    names(rp) <- comment(fy_sub)
    
    diff_brick <- stack(diff_brick , rp)
  }
  
  out <- list(rcp_df, diff_brick)
  return(out)
}


###################################################################################################################
# Archive
###################################################################################################################

# predict_data <- function(model, voi, res, TFE=F) {
#   
#   df <- model[["data"]]
#   
#   # drop GEOID with all yield missing 
#   id_list <- df %>% group_by(GEOID) %>% summarize(n = sum(is.na(Yield))) %>% filter(n == 31) %>% select(GEOID) 
#   df <- df %>% filter(!GEOID %in% id_list$GEOID) 
#   
#   final_list <- list()
#   for (i in 1:length(voi)) {
#     
#     dfmean <- t(as.data.frame(colMeans(df[,voi], na.rm=T))) 
#     rownames(dfmean) <- NULL
#     v <- voi[i] 
#     
#     if (TFE == TRUE) {
#       geotime <- merge(unique(df$GEOID), unique(df$Year))
#       colnames(geotime) <- c("GEOID", "Year")
#       
#       final <- list()
#       for (r in 1:nrow(geotime)) {
#         GEOID <- rep(df[r,"GEOID"], res)
#         Year <- rep(df[r, "Year"], res)
#         dfm <- as.data.frame(dfmean[rep(seq_len(nrow(dfmean)), res),])
#         
#         x = seq(quantile(df[,v], 0.01, na.rm = T), 
#                 quantile(df[,v], 0.99, na.rm = T), 
#                 length.out = res)
#         dfm[,v] <- x 
#         dfm <- cbind(dfm, GEOID, Year)
#         
#         # predict yield
#         Yield <- predict(model[["model"]], newdata = dfm, type = "response", se.fit = T)
#         Yield <- as.data.frame(Yield)
#         colnames(Yield) <- c(paste0("Predict", v))
#         dfm$Yield <- Yield[,paste0("Predict",v)]
#         final[[r]] <- dfm
#       }
#       
#       final_df <- do.call(rbind, final) 
#       final_list[[v]] <- final_df }
#     
#     if (TFE == F) {
#       geo <- as.data.frame(unique(df$GEOID))
#       colnames(geo) <- c("GEOID")
#       
#       final <- list()
#       for (r in 1:nrow(geo)) {
#         GEOID <- rep(df[r,"GEOID"], res)
#         dfm <- as.data.frame(dfmean[rep(seq_len(nrow(dfmean)), res),])
#         
#         lsd <- mean(df[,v], na.rm=T) - 2*sd(df[,v], na.rm=T)
#         lsd <- ifelse(lsd < 0, 0, lsd)
#         #hsd <- mean(df[,v], na.rm=T) + 2*sd(df[,v], na.rm=T)
#         x <- seq(lsd, var_range[[v]][2], length.out = res)
#         dfm[,v] <- x 
#         dfm <- cbind(dfm, GEOID)
#         
#         # predict yield
#         Yield <- predict(model[["model"]], newdata = dfm, type = "response")
#         Yield <- as.data.frame(Yield)
#         colnames(Yield) <- c(paste0("Predict", v))
#         dfm$Yield <- Yield[,paste0("Predict",v)]
#         final[[r]] <- dfm
#       }
#       
#       final_df <- do.call(rbind, final) 
#       final_list[[v]] <- final_df }
#   }
#   
#   return(final_list)}
# 
# prep_climate <- function(crop, dfl, yield, controls) {
#   final_df <- data.frame()
#   
#   for (i in 1:length(dfl)) {
#     
#     df <- dfl[[i]]
#     # data.frame metadata
#     var_name <- comment(df)
#     colnames(df) <- c("Year", "Month", "Day", coords[,"ID"])
#     
#     # month, day to DOY
#     df$Day <- str_pad(df$Day, 2, pad="0")
#     df$Month <- str_pad(df$Month, 2, pad = "0")
#     date <- ymd(as.character(paste0(df$Year, df$Month, df$Day)))
#     df$DOY <- yday(date)
#     
#     # tidy Schoof dataset
#     df <- gather(df, station, measurement, 4:217)
#     # merge to county shapefile and station shapefile
#     # without all.x = T, this drops all counties but the 102 with stations
#     df <- merge(df, coords_sp@data, by.x = "station", by.y = "ID") 
#     coi_crop <- coi[,c("GEOID", "STATE", "COUNTY", paste0(crop, ".plant.start"), paste0(crop,".harvest.end"))]
#     df <- merge(df, coi_crop, by.x = "GEOID", by.y = "GEOID")
#     
#     df <- time_subset(df, crop)
#     df$measurement[df$measurement == "NaN"] <- NA
#     names(df)[names(df) == 'measurement'] <- var_name
#     
#     ifelse(i == 1, final_df <- df, final_df <- cbind(final_df, df[,var_name]))
#     names(final_df)[names(final_df) == 'df[, var_name]'] <- var_name
#     comment(final_df) <- crop
#   }
#   
#   final_df$gdd <- gdd(final_df, crop, gdd_range)
#   final_df$sdd <- sdd(final_df, crop, gdd_range)
#   final_df$efp <- efp(final_df, th = 20)
#   return(final_df)
# }
# 
# seasonal_panel <- function(df, crop) {
#   season <- df %>% group_by(Year, GEOID) %>% 
#     summarize(GDD = sum(gdd, na.rm=T), SDD = sum(sdd, na.rm=T), 
#               TP = sum(rain, na.rm=T), EfP = sum(efp, na.rm = T)) %>%
#     mutate(ExP = TP - EfP)
#   season <- season %>% filter(Year > 1979)
#   controls <- add_controls(crop)
#   season <- merge(season, controls, by = c("Year", "GEOID"))
#   
#   comment(season) <- crop
#   return(season)
# }
# 
# construct_hist_data <- function(daily = F, save = F) {
#   
#   corn_daily <- prep_climate(dfl = list(rain, tmax), crop = "corn") #1 min
#   wheat_daily <- prep_climate(dfl = list(rain, tmax), crop = "wheat")
#   wwheat_daily <- prep_climate(dfl = list(rain, tmax), crop = "wwheat")
#   cotton_daily <- prep_climate(dfl = list(rain, tmax), crop = "cotton")
#   soy_daily <- prep_climate(dfl = list(rain, tmax), crop = "soy")
#   
#   corn_panel <- seasonal_panel(corn_daily, "corn")
#   wheat_panel <- seasonal_panel(wheat_daily, "wheat")
#   wwheat_panel <- seasonal_panel(wwheat_daily, "wwheat")
#   cotton_panel <- seasonal_panel(cotton_daily, "cotton")
#   soy_panel <- seasonal_panel(soy_daily, "soy")
#   
#   if (daily == T) {
#     return(c(corn_daily, wheat_daily, wwheat_daily, cotton_daily, soy_daily))
#     return(c(corn_panel, wheat_panel, wwheat_panel, cotton_panel, soy_panel))
#   }
#   else {
#     return(c(corn_panel, wheat_panel, wwheat_panel, cotton_panel, soy_panel))
#   }
#   
#   if (save == T) {
#     saveRDS(corn_panel, "./out/panels/corn.RDS")
#     saveRDS(cotton_panel, "./out/panels/cotton.RDS")
#     saveRDS(soy_panel, "./out/panels/soy.RDS")
#     saveRDS(wheat_panel, "./out/panels/wheat.RDS")
#     saveRDS(wwheat_panel, "./out/panels/wwheat.RDS")
#   }
# }

# predict_data_PRISM <- function(model, voi, res) {
#   # ch
#   # predicts for one variable varying from .01 to .99 of dxn while others are held at mean, averaging over all GEOID effects
#   
#   df <- model[["data"]]
#   df$Year <- df$YEAR
#   cm <- length(unique(df$Year))
#   
#   # filter out complete cases of missing data
#   id_list <- df %>% group_by(GEOID) %>% summarize(n = sum(is.na(Yield))) %>% filter(n == cm) %>% select(GEOID) # drop GEOID with all missing Yield
#   id_list_gdd <- df %>% group_by(GEOID) %>% summarize(n = sum(is.na(GDD))) %>% filter(n == cm) %>% select(GEOID)
#   
#   df <- df %>% filter(!GEOID %in% id_list$GEOID)  %>% filter(!GEOID %in% id_list_gdd$GEOID)
#   
#   final_list <- list()
#   
#   for (i in 1:length(voi)) {
#     
#     dfmean <- t(as.data.frame(colMeans(df[,voi], na.rm=T))) # row means
#     dfmean[1,"YEAR"] <- 2010
#     rownames(dfmean) <- NULL
#     v <- voi[i] 
#     
#     geo <- as.data.frame(unique(df$GEOID))
#     colnames(geo) <- c("GEOID")
#     
#     final <- list()
#     
#     #simulate for each GEOID
#     for (r in 1:nrow(geo)) {
#       
#       GEOID <- rep(df[r,"GEOID"], res) #repeat GEOID res times
#       dfm <- as.data.frame(dfmean[rep(seq_len(nrow(dfmean)), res),])  #repeat voi at mean
#       rownames(dfm) <- NULL
#       
#       x = seq(quantile(df[,v], 0.01, na.rm = T), 
#               quantile(df[,v], 0.99, na.rm = T), 
#               length.out = res)
#       dfm[,v] <- x 
#       dfm <- cbind(dfm, GEOID)
#       
#       # predict yield in county GEOID
#       Yield <- predict(model[["model"]], newdata = dfm, type = "response", se.fit = T) 
#       Ymn <- as.data.frame(Yield[[1]])
#       Ysd <- as.data.frame(Yield[[2]])
#       Yield <- cbind(Ymn, Ysd)
#       colnames(Yield) <- c(paste0("Predict", v), paste0("SD", v))
#       dfm$Yield <- Yield[,paste0("Predict",v)]
#       dfm$YieldSD <- Yield[,paste0("SD", v)]
#       final[[r]] <- dfm
#     }
#     final_df <- do.call(rbind, final) 
#     final_list[[v]] <- final_df }
#   
#   return(final_list)
#   
# }
# 



