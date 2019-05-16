# FUNCTIONS

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
  # from yield folder, QuickStats
  yield <- readRDS("./data/yield/yield.RDS")
  
  if (crop == "corn") {
    y <- yield %>% filter(CROP == "CORN")
  }
  if (crop == "wwheat") {
    y <- yield %>% filter(CROP == "WWHEAT")
  }
  if (crop == "soy") {
    y <- yield %>% filter(CROP == "SOYBEANS")
  }
  if (crop == "cotton") {
    y <- yield %>% filter(CROP == "COTTON")
  }
  if (crop == "wheat") {
    y <- yield %>% filter(CROP == "WHEAT")
  }
  y <- y %>% select(-c(CROP))
  return(y)
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
  final_df$efp <- efp(final_df, th = 30)
  return(final_df)
}

seasonal_panel_PRISM <- function(df, crop) {
  
  season <- df %>% group_by(YEAR, GEOID) %>% 
    summarize(GDD = sum(gdd, na.rm=T), SDD = sum(sdd, na.rm=T), 
              TP = sum(PPT, na.rm=T), EfP = sum(efp, na.rm = T)) %>%
    mutate(ExP = TP - EfP)
  season <- as.data.frame(season)
  yield <- add_controls(crop)
  season <- merge(season, yield, by = c("YEAR", "GEOID")) 
  comment(season) <- crop
  return(season)
}

construct_hist_data_PRISM <- function() {
  
  corn_daily <- prep_climate_PRISM(crop = "corn", dfl = list(ppt, tmax)) #1 min
  wwheat_daily <- prep_climate_PRISM(dfl = list(ppt, tmax), crop = "wwheat")
  wheat_daily <- prep_climate_PRISM(dfl = list(ppt, tmax), crop = "wheat")
  cotton_daily <- prep_climate_PRISM(dfl = list(ppt, tmax), crop = "cotton")
  soy_daily <- prep_climate_PRISM(dfl = list(ppt, tmax), crop = "soy")
  
  corn_panel <- seasonal_panel_PRISM(corn_daily, "corn")
  wwheat_panel <- seasonal_panel_PRISM(wwheat_daily, "wwheat")
  wheat_panel <- seasonal_panel_PRISM(wheat_daily, "wheat")
  cotton_panel <- seasonal_panel_PRISM(cotton_daily, "cotton")
  soy_panel <- seasonal_panel_PRISM(soy_daily, "soy")
  
  saveRDS(corn_panel, "./out/panels/corn_PRISM_30.RDS")
  saveRDS(cotton_panel, "./out/panels/cotton_PRISM_30.RDS")
  saveRDS(wheat_panel, "./out/panels/wheat_PRISM_30.RDS")
  saveRDS(soy_panel, "./out/panels/soy_PRISM_30.RDS")
  saveRDS(wwheat_panel, "./out/panels/wwheat_PRISM_30.RDS")
  
  # run add_hist_co2 to add co2 data
  
}

add_hist_co2 <- function(crop) {
  # Meinshausen 2011 et al
  co <- read.csv("./data/co2.csv")
  df <- readRDS(paste0("./out/panels/", crop, "_PRISM_30.RDS"))
  co <- co %>% filter(Year %in% unique(df$YEAR))
  df <- merge(df, co, by.x = "YEAR", by.y = "Year", all=T)
  saveRDS(df, paste0("./out/panels/", crop, "_PRISM_30.RDS"))
}

###################################################################################################################
# GAM functions
###################################################################################################################

gam_data_prep <- function(df) {
  
  df <- df %>% filter(!is.na(YIELD)) %>% filter(GEOID %in% irr_list)
  df <- df[df$YIELD > quantile(df$YIELD, 0.01), ]
  df <- df[df$YIELD < quantile(df$YIELD, 0.99), ]
  df <- df %>% arrange(YEAR, GEOID)
  
  return(df)
}

run_gam <- function(formula, model_name, crop) {
  df <- readRDS(paste0("./out/panels/", crop, "_PRISM_30.RDS"))
  df <- gam_data_prep(df)
  gam <- bam(formula, data=df, family = "gaussian", method="REML")
  saveRDS(gam, paste0("./out/gam/", model_name, "_", crop, ".RDS"))
  return(gam)
}

run_gam_RMSE <- function(formula, crop) {
  df <- readRDS(paste0("./out/panels/", crop, "_PRISM_30.RDS"))
  df <- gam_data_prep(df)
  df$RN <- 1:nrow(df)
  
  t <- df %>% group_by(GEOID) %>%
    summarize(n=n()) %>% filter(n > 3)
  geoid <- unique(t$GEOID)
  
  HO <- list()
  for (i in 1:length(geoid)) {
    sub <- df[df$GEOID == geoid[i],]
    sub <- sub[sample(nrow(sub), 2), ]
    HO[[i]] <- sub
  }
  
  HO <- do.call(rbind.data.frame, HO)
  RMSE <- rmse(HO$YIELD, HO$PRED)
  RMSE_NULL <- rmse(mean(HO$YIELD, na.rm=T), HO$PRED)
  sample <- df %>% filter(!RN %in% unique(HO$RN))
  
  gam <- bam(formula, data=sample, family = "gaussian", method="REML")
  pred <- predict(gam, newdata = HO, type = "response", se=T)
  HO$PRED <- pred[[1]]
  
  rmse <- rmse(HO$YIELD, HO$PRED)
  rmse_n <- rmse(HO$YIELD, mean(HO$YIELD, na.rm = T))
  
  return(list(rmse, rmse_n))
  
}

spatial_effects <- function(gam, crop) {
  
  df <- readRDS(paste0("./out/panels/", crop, "_PRISM_30.RDS"))
  df <- gam_data_prep(df)
  names <- names(coef(gam))
  idx <- which(names == "s(GEOID).1")
  geoids <- levels(df$GEOID) %in% unique(df$GEOID)
  geoids <- levels(df$GEOID)[geoids]
  coef <- cbind.data.frame(coef(gam)[idx:length(names)], geoids)
  colnames(coef) <- c("SE", "GEOID")
  se <- county %>% filter(GEOID %in% unique(df$GEOID)) 
  se@data <- merge(se@data, coef, by="GEOID", all=T)
  spplot(se, "SE")
  
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
  
  # load ensembles built in build_ensembles()
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
      df <- merge(df, coords_sp@data, by.x = "station", by.y = "ID")  # reduces to 102 stations only
      coi_crop <- coi[,c("GEOID", "STATE", "COUNTY", paste0(crop, ".plant.start"), paste0(crop,".harvest.end"))]
      df <- merge(df, coi_crop, by.x = "GEOID", by.y = "GEOID") # n = 102
      
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
  
  ens_list <- c("./out/ensembles/rcp45", "./out/ensembles/rcp85")
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

add_future_co2 <- function(crop) {
  # meinshausen 2011 et al
  co <- read.csv("./data/co2.csv")
  df <- readRDS(paste0("./out/fpanels/rcp45_", crop, ".RDS"))
  co <- co %>% filter(Year %in% unique(df$Year))
  df <- merge(co, df, by = "Year", all = T)
  saveRDS(df, paste0("./out/fpanels/rcp45_", crop, ".RDS"))
  
  co <- read.csv("./data/co2.csv")
  df <- readRDS(paste0("./out/fpanels/rcp85_", crop, ".RDS"))
  co <- co %>% filter(Year %in% unique(df$Year))
  df <- merge(co, df, by = "Year", all = T)
  saveRDS(df, paste0("./out/fpanels/rcp85_", crop, ".RDS"))
}

# add_future_co2("corn")
# add_future_co2("soy")
# add_future_co2("cotton")
# add_future_co2("wheat")
# add_future_co2("wwheat")

###################################################################################################################
# GAM prediction 
###################################################################################################################

gam_predict <- function(gam, crop, scenario) {
  
  df <- readRDS(paste0("./out/panels/", crop, "_PRISM_30.RDS"))
  df <- gam_data_prep(df)
  rcp45 <- co2_scaling(crop, rcp = "45")[[2]] %>% 
    filter(GEOID %in% unique(df$GEOID))
  
  scale <- scenario_dvt(crop)[,scenario]
  rcp45$scale <- rep(scale, nrow(rcp45))
  pred45 <- predict(gam, newdata = rcp45, type = "response", se=T)
  
  rcp85 <- co2_scaling(crop, rcp = "85")[[2]] %>% 
    filter(GEOID %in% unique(df$GEOID))
  rcp85$scale <- rep(scale, nrow(rcp85))
  pred85 <- predict(gam, newdata = rcp85, type = "response", se=T)
  
  rcp45$YIELD <- pred45$fit
  rcp45$YIELD_SE <- pred45$se.fit
  rcp45$RCP <- "RCP4.5"
  rcp45$SCENARIO <- scenario
  rcp85$YIELD <- pred85$fit
  rcp85$YIELD_SE <- pred45$se.fit
  rcp85$RCP <- "RCP8.5"
  rcp85$SCENARIO <- scenario
  
  out <- rbind.data.frame(rcp45, rcp85)
  out$YIELD <- out$YIELD * (1-out$SCALE)
  
  return(out)
}

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

# play with increasing grid size for higher reoslution projections
future_int <- function(pred) {
  
  bl <- pred %>% 
    filter(YEAR %in% 2000:2010) %>% 
    group_by(GEOID) %>% 
    summarize(mnYIELD = mean(YIELD, na.rm=T))
  comment(bl) <- "2000s"
  
  grd <- SpatialPixels(SpatialPoints(makegrid(coi, n=7000)), proj4string = proj4string(coi))
  grd <- grd[coi,] 
  
  bls <- merge(coi, bl, by = "GEOID", all=T)
  pcentroids <- as.data.frame(coordinates(bls))
  colnames(pcentroids) <- c("lon", "lat") 
  pcentroids <- data.frame("ID" = 1:nrow(pcentroids), pcentroids)
  coordinates(pcentroids) <- c("lon", "lat") 
  proj4string(pcentroids) <- proj4string(bls) # assign projection
  pcentroids@data <- sp::over(x = pcentroids, y = bls, returnList = FALSE)
  pcentroids <- pcentroids %>% filter(!is.na(mnYIELD))
  
  p.idw <- idw(formula = mnYIELD ~ 1, pcentroids, grd,  idp = 2, nmax = 12)
  p.idw <- raster(p.idw)
  comment(p.idw) <- "2000s"
  
  rcp_list <- list("RCP4.5", "RCP8.5")
  
  out <- list()
  
  for (i in 1:length(rcp_list)) {
    
    psub <- pred %>% filter(RCP == rcp_list[[i]])
    f30 <- psub %>% filter(YEAR %in% 2030:2040) %>% 
      group_by(GEOID) %>% 
      summarize(mnYIELD = mean(YIELD, na.rm=T))
    comment(f30) <- paste0("2030s/", rcp_list[[i]])
    f60 <- psub %>% filter(YEAR %in% 2060:2070) %>% 
      group_by(GEOID) %>% 
      summarize(mnYIELD = mean(YIELD, na.rm=T))
    comment(f60) <- paste0("2060s/", rcp_list[[i]])
    f90 <- psub %>% filter(YEAR %in% 2090:2100) %>% 
      group_by(GEOID) %>% 
      summarize(mnYIELD = mean(YIELD, na.rm=T))
    comment(f90) <- paste0("2090s/", rcp_list[[i]])
    
    rcp <- list(f30, f60, f90)
    out_brick <- stack(p.idw)
    
    for (s in 1:length(rcp)) {
      fy_sub <- rcp[[s]]
      fys <- merge(coi, fy_sub, by = "GEOID", all.x = T)
      centroids <- as.data.frame(coordinates(fys))
      colnames(centroids) <- c("lon", "lat") 
      centroids <- data.frame("ID" = 1:nrow(centroids), centroids)
      coordinates(centroids) <- c("lon", "lat") 
      proj4string(centroids) <- proj4string(fys) # assign projection
      centroids@data <- sp::over(x = centroids, y = fys, returnList = FALSE)
      centroids <- centroids %>% filter(!is.na(mnYIELD))
      f.idw  <- idw(formula = mnYIELD ~ 1, centroids, grd,  idp = 2, nmax = 12)
      f.idw <- raster(f.idw)
      comment(f.idw) <- comment(fy_sub)
      
      out_brick <- stack(out_brick , f.idw)
      
    }
    out[[i]] <- out_brick
  }
  return(out)
}

future_yield_map <- function(maps, crop, scenario, save = F) {
  
  # maps result of future_int() function
  
  cropn <- ifelse(crop == "wwheat", "Winter wheat", str_to_title(crop))
  cropn <- ifelse(cropn == "Corn", "Maize", str_to_title(crop))
  
  rcp45 <- stack()
  for (i in 2:dim(maps[[1]])[3]) {
    r <- ((maps[[1]][[i]] - maps[[1]][[1]])/maps[[1]][[1]])*100
    rcp45 <- stack(rcp45, r)
  }
  
  rcp85 <- stack()
  for (i in 2:dim(maps[[2]])[3]) {
    r <- ((maps[[2]][[i]] - maps[[2]][[1]])/maps[[2]][[1]])*100
    rcp85 <- stack(rcp85, r)
  }
  
  fy <- stack(rcp45[[1]], rcp85[[1]],
              rcp45[[2]], rcp85[[2]],
              rcp45[[3]], rcp85[[3]])
  
  stack_names <- c("2030s-RCP4.5", "2030s-RCP8.5", 
                   "2060s-RCP4.5", "2060s-RCP8.5",
                   "2090s-RCP4.5", "2090s-RCP8.5")
  
  names(fy) <- stack_names
  
  vmax <- max(maxValue(fy))
  vmin <- min(minValue(fy))
  
  seq_c <- c(-150, -100, -50, 0, 50, 100, 150, 200, 250)
  fy[fy>max(seq_c)] <- max(seq_c) - 1
  fy[fy<min(seq_c)] <- min(seq_c) + 1
  
  coist <- spTransform(coist, fy@crs)
  # col.l <- colorRampPalette(c("#C99C56", "#DABA75", "#F5F5F5",
  #                             "#CAE6E2",
  #                             "#B5DFD8", "#9FD7CF", "#8AD0C5",
  #                             "#74C6B9", "#5DB9AB", "#46AC9C",
  #                             "#2F9F8E"))
  col.l <- colorRampPalette(c("#C99C56", "#DABA75", "#ecdbb7", "#F5F5F5", "#CAE6E2",
                              "#B2DED7", "#90D2C8",
                              "#6DC2B5", 
                              "#49AE9E", "#259987"))
  
  # labels <- c("<-150", "", "",  "0", "", "", "150", "", "", ">300")
  
  # scenario <- ifelse(scenario == "None", "No", scenario)
  rp <- levelplot(fy, 
                  scales=list(draw=FALSE), 
                  col.regions = col.l,
                  at = seq_c,
                  names.attr=stack_names,
                  ylab = "",
                  main = "",
                  layout = c(2, 3),
                  par.strip.text=list(cex=.6, lines=1), xlab=NULL,
                  colorkey = list(space="bottom", at=seq_c))
  
                 # xlab.top = paste(paste(toupper(substr(cropn, 1, 1)), 
                  #                       substr(cropn, 2, nchar(cropn)), sep=""),
                   #                "yields (bu/ac) - ", scenario, "time effect"))
  rpf <- rp + layer(sp.polygons(coist, alpha = 0.5)) 
  
  png(filename=paste0("./fig/", crop, "_", scenario, ".png"), width = 5, height = 5, units = 'in', res = 600)
  #png(filename=paste0("./fig/", crop, "_", scenario, ".png"))
  print(rpf)
  dev.off()
  return(rpf)
}

future_yield_plot <- function(pred, crop) {
  
  ysum <- pred %>% 
    group_by(SCENARIO, RCP, YEAR) %>% 
    dplyr::summarize(mnYIELD = mean(YIELD, na.rm=T), sdYIELD = sd(YIELD, na.rm=T))
  
  ss <- ggplot(ysum, aes(x = YEAR, y = mnYIELD, color = SCENARIO, shape = RCP, group = interaction(SCENARIO, RCP))) +
    #geom_point(size = 2) +
    geom_point(data = pred, aes(x = YEAR, y = YIELD), alpha = 0.03) +
    geom_line(aes(linetype = RCP), size=1, alpha=1) +
    scale_color_manual(values=c("#9999CC", "#CC6666", "#008000", "#FFD700")) +
    xlab("") +
    ylab("") +
    ggtitle(paste0(simpleCap(crop), " yields (bu/ac)")) +
    theme_minimal()
  
  return(ss)
}

scenario_dvt <- function(crop) {
  
  df <- readRDS(paste0("./out/panels/", crop, "_PRISM_30.RDS"))
  if (crop == "cotton") {
    df$YIELD <- df$YIELD/32
  }
  decades <- c(1981, 1991, 2001, 2011)  # change based on Chris' response
  
  fe_list <- list()
  
  for (i in 1:length(decades)) {
    
    dfs <- df %>% filter(YEAR %in% seq(decades[i],(decades[i]+10)))
    fe <- lmer(YIELD ~ GDD + SDD + ExP + EfP + YEAR + (1|GEOID), data=dfs)
    coef <- summary(fe)$coefficient["YEAR", "Estimate"]
    fe_list[[i]] <- c(decades[i], coef)
    
    if (i ==length(decades)) {
      fe <- lmer(YIELD ~ GDD + SDD + ExP + EfP + YEAR + (1|GEOID), data=df)
      coef <- summary(fe)$coefficient["YEAR", "Estimate"]
      fe_list[[i+1]] <- c('ALL', coef)
    }
    
  }
  
  fe <- do.call(rbind.data.frame, fe_list)
  colnames(fe) <- c("Decade", "Coef")
  fe$Coef <- as.numeric(as.character(fe$Coef))
  
  scenarios <- cbind.data.frame(max(fe$Coef[1:length(decades)]), min(fe$Coef[1:length(decades)]),
                                fe$Coef[length(decades)+1], 0.0001)
  colnames(scenarios) <- c("High", "Low", "Average", "None")
  return(scenarios)
  
}

co2_scaling <- function(crop, rcp) {
  
  # rcp = "45"
  df <- readRDS(paste0("./out/panels/", crop, "_PRISM_30.RDS"))
  if (crop == "cotton") {
    df$YIELD <- df$YIELD/32
  }
  f <- readRDS(paste0("./out/fpanels/rcp", rcp, "_", crop, ".RDS")) %>%
    mutate(YEAR = Year, YIELD = NA) %>% select(-c("Slope", "pH_sur", "AWC_sur", "Year")) 
  
  scale <- co2[[crop]]  # use this to drop YEAR effect historically
  
  varname <- paste0("CO2.RCP", rcp)
  df_sub <- df %>% filter(YEAR %in% 1990:2010) 
  hist_co2 <- mean(df_sub[,varname], na.rm=T)
  
  df <- df %>% mutate(SCALE = (scale*df[,varname])/hist_co2) 
  f <- f %>% mutate(SCALE = (scale*f[,varname])/hist_co2) 
  return(list(df, f))
}

all_scenarios <- function(crop) {
  
  #High
  model <- "gamly_High" # High
  scenario <- "High"
  gam <- readRDS(paste0("./out/gam/", model, "_", crop, ".RDS"))
  predH <- gam_predict(gam, crop, scenario = scenario)
  
  #Average
  model <- "gamly_Average" 
  scenario <- "Average"
  gam <- readRDS(paste0("./out/gam/", model, "_", crop, ".RDS"))
  predA <- gam_predict(gam, crop, scenario = scenario)
  
  #Low
  model <- "gamly_Low" 
  scenario <- "Low"
  gam <- readRDS(paste0("./out/gam/", model, "_", crop, ".RDS"))
  predL <- gam_predict(gam, crop, scenario = scenario)
  
  #None
  model <- "gamly_None" 
  scenario <- "None"
  gam <- readRDS(paste0("./out/gam/", model, "_", crop, ".RDS"))
  predN <- gam_predict(gam, crop, scenario = scenario)
  
  pred <- rbind.data.frame(predH, predA, predL, predN)
  plt <- future_yield_plot(pred, crop)
  return(plt)
  
}

all_maps <- function(crop) {
  
  #High
  model <- "gamly_High" # High
  scenario <- "High"
  gam <- readRDS(paste0("./out/gam/", model, "_", crop, ".RDS"))
  predH <- gam_predict(gam, crop, scenario = scenario)
  
  #Average
  model <- "gamly_Average" 
  scenario <- "Average"
  gam <- readRDS(paste0("./out/gam/", model, "_", crop, ".RDS"))
  predA <- gam_predict(gam, crop, scenario = scenario)
  
  #Low
  model <- "gamly_Low" 
  scenario <- "Low"
  gam <- readRDS(paste0("./out/gam/", model, "_", crop, ".RDS"))
  predL <- gam_predict(gam, crop, scenario = scenario)
  
  #None
  model <- "gamly_None" 
  scenario <- "None"
  gam <- readRDS(paste0("./out/gam/", model, "_", crop, ".RDS"))
  predN <- gam_predict(gam, crop, scenario = scenario)
  
  r <- rbind.data.frame(predH, predA, predL, predN)
  
  ### maps ###
  
  scenario <- "High"
  mapsH <- future_int(predH)
  mapH <- future_yield_map(mapsH, crop, scenario, save=T)
  print("Saved high")
  
  scenario <- "Average"
  mapsA <- future_int(predA)
  mapA <- future_yield_map(mapsA, crop, scenario, save=T)
  print("Saved avg")
  
  scenario <- "Low"
  mapsL <- future_int(predL)
  mapL <- future_yield_map(mapsL, crop, scenario, save=T)
  print("Saved low")
  
  scenario <- "None"
  mapsN <- future_int(predN)
  mapN <- future_yield_map(mapsN, crop, scenario, save=T)
  print("Saved none")
  
  all_maps <- list(mapH, mapA, mapL, mapN)
  
  return(all_maps)
}

