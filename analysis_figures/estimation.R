#Calibration for Amazon dynamics
dev.off()
rm(list=ls())
#setwd("")
set.seed(1234)

#load libraries
libraries=c("terra","Rcpp","truncnorm","purrr","readxl","webshot2","writexl","htmlwidgets","deSolve","doParallel","parallel","plotly","lhs","fixest","xtable","ggplot2","lubridate","dplyr","tidyr","knitr","kableExtra","plm","lmtest","sandwich")
lapply(libraries, require, character.only = TRUE) #load libraries
library(patchwork)  

#define direction
mydirection= file.path(getwd(),"analysis_figures/")
mydirection_numericalmodel= file.path(getwd(),"numerical_model/parameters/")

###########################################################################
############### Step 1: load data and prepare panel FE  ###################
###########################################################################
#define parameter
mcwd_to_carbon = 0.05 #estimate from Phillips et al. (2009)

#load Amazon shapefile #from Souza
boundaries <- vect(file.path(mydirection,"data/subregions_souza2016/amazon_subregions.shp"))
boundaries <- project(boundaries, "EPSG:4326")

#load temperature data for each RCP/ESM
#create one stars for each RCP/ESM
#clean to amazon rainforest extent
#compute mean annual temperature from daily mean temperature
#mean temperature is from october n to septembre n+1
#from october 2006 to september 2099
#for reproducibility, here is how we create temp and mcwd anomaly from raw datasets
#for reproducibility of estimates, start from these transformed files directly

MODEL <- c("GFDL", "HadGEM","IPSL","MIROC")
RCP <- c("rcp26", "rcp60", "rcp85")
files_temperature <- list.files(file.path(mydirection,"data/temperature/"), pattern = "\\.nc4$", full.names = TRUE)

get_model_rcp <- function(f) {
  model <- MODEL[sapply(MODEL, function(m) grepl(m, f))]
  rcp <- RCP[sapply(RCP, function(r) grepl(r, f))]
  if (length(model) == 1 && length(rcp) == 1) {return(paste(model, rcp, sep = "_"))} else {return(NA)}}

file_groups_temperature <- split(files_temperature, sapply(files_temperature, get_model_rcp))

raster_groups_temperature <- lapply(file_groups_temperature, function(group_files) {
  r_list <- lapply(group_files, rast)
  names(r_list) <- NULL  
  r <- do.call(c, r_list) 
  return(r)})

rast_temperature <- lapply(raster_groups_temperature, function(r) {
  if (!same.crs(r, boundaries)) {boundaries <- project(boundaries, crs(r))}
  r <- crop(r, boundaries)
  r <- mask(r, boundaries)
  dates <- time(r)
  hydro_year <- ifelse(month(dates) >= 10, year(dates) + 1, year(dates)) #oct to sept
  unique_years <- sort(unique(hydro_year))
  r_annual <- rast()
  for (y in unique_years) {
    idx <- which(hydro_year == y)
    if (length(idx) > 0) {
      r_year <- mean(r[[idx]], na.rm = TRUE)
      names(r_year) <- paste0("year_", y)
      r_annual <- c(r_annual, r_year)}}
  r_annual = r_annual - 273.15 #kelvin to celsius
  return(r_annual)})

names(rast_temperature) <- tolower(names(rast_temperature))
output_dir <- file.path(mydirection, "output/temperature")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

lapply(names(rast_temperature), function(name) {
  r <- rast_temperature[[name]]
  out_path <- file.path(output_dir, paste0(name, ".tif"))
  writeRaster(r, out_path, overwrite = TRUE)})

#load monthly precipitation data
#compute MCWD for historical (1980-2004)
#compute MCWD anomaly wrt baseline for each RCP/ESM
MODEL <- c("gfdl", "hadgem","ipsl","miroc")
RCP <- c("rcp26", "rcp60", "rcp85","historical")
files_rainfall <- list.files(file.path(mydirection,"data/rainfall/"), pattern = "\\.nc4$", full.names = TRUE)
file_groups_rainfall <- split(files_rainfall, sapply(files_rainfall, get_model_rcp))

raster_groups_rainfall <- lapply(file_groups_rainfall, function(group_files) {
  r_list <- lapply(group_files, rast)
  names(r_list) <- NULL  
  r <- do.call(c, r_list) 
  return(r)})

rast_rainfall <- lapply(raster_groups_rainfall, function(r) {
  if (!same.crs(r, boundaries)) {boundaries <- project(boundaries, crs(r))}
  r <- crop(r, boundaries)
  r <- mask(r, boundaries)
  return(r)})

compute_mcwd <- function(r) {
  dates <- time(r)
  months <- month(dates)
  years <- year(dates)
  days_in_month <- days_in_month(dates)
  r_mm <- r * 86400 * days_in_month #kg/m²/s → mm/month
  hydro_years <- ifelse(months >= 10, years + 1, years)
  unique_hydro_years <- sort(unique(hydro_years))
  mcwd_stack <- rast()
  for (y in unique_hydro_years) {
    idx <- which(hydro_years == y)
    r_year <- r_mm[[idx]]
    n_months <- nlyr(r_year)
    cwd <- r_year[[1]]
    cwd <- ifel(cwd >= 100, 0, cwd - 100)
    mcwd=cwd
    for (i in 2:n_months) {
      rain <- r_year[[i]]
      deficit <- rain - 100
      cwd <- ifel(rain >= 100, 0, cwd + deficit)
      mcwd <- min(mcwd, cwd, na.rm = TRUE)}
    names(mcwd) <- paste0("MCWD_", y)
    mcwd_stack <- c(mcwd_stack, mcwd)}
  return(mcwd_stack)}

rain_hist <- rast_rainfall[grep("historical", names(rast_rainfall))]
rain_rcp <- rast_rainfall[!grepl("historical", names(rast_rainfall))]

#clean temporal dimension
start_date <- as.Date("1861-01-01")
rain_hist <- lapply(rain_hist, function(r) {
  n <- nlyr(r)
  dates <- seq(from = start_date, by = "1 month", length.out = n)
  terra::time(r) <- dates
  return(r)})

start_date <- as.Date("2006-01-01")
rain_rcp <- lapply(rain_rcp, function(r) {
  n <- nlyr(r)
  dates <- seq(from = start_date, by = "1 month", length.out = n)
  terra::time(r) <- dates
  return(r)})

mcwd_hist <- lapply(rain_hist, compute_mcwd)
mcwd_rcp <- lapply(rain_rcp, compute_mcwd)

baseline_means <- lapply(mcwd_hist, function(mcwd) {
  years <- as.integer(gsub("MCWD_", "", names(mcwd)))
  idx <- which(years >= 1981 & years <= 2004)
  mean(mcwd[[idx]], na.rm = TRUE)})

mcwd_anomaly <- mapply(function(mcwd_proj, baseline) {
  mcwd_proj - baseline  
}, mcwd_rcp, baseline_means, SIMPLIFY = FALSE)

output_dir <- file.path(mydirection, "output/rainfall")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

lapply(names(mcwd_anomaly), function(name) {
  r <- mcwd_anomaly[[name]]
  out_path <- file.path(output_dir, paste0(name, ".tif"))
  writeRaster(r, out_path, overwrite = TRUE)})


#prepare data
#match annual mean temp and annual MCWD anomaly data
#create regional fixed effects
#scale carbon losses by carbon heterogeneity in amazon (weights in regression)
MODEL <- c("gfdl", "hadgem","ipsl","miroc")
RCP <- c("rcp26", "rcp60", "rcp85")

get_raster_info <- function(r, model, rcp, variable_name, column_name) {
  xy <- as.data.frame(r, xy = TRUE)  
  xy <- xy %>%
    pivot_longer(cols = starts_with(column_name), names_to = "year", values_to = variable_name)
  xy$year=substr(xy$year,6,9)
  xy$model <- model
  xy$rcp <- rcp
  return(xy)}

df <- data.frame()
for (model in MODEL) {
  df_update_model<-data.frame()
  for (rcp in RCP) {
    df_update_rcp<- data.frame()
    mcwd_anomaly = rast(file.path(mydirection, paste0("/output/rainfall/",model,"_",rcp,".tif")))
    temp = rast(file.path(mydirection, paste0("/output/temperature/",model,"_",rcp,".tif")))
    mcwd_anomaly <- get_raster_info(mcwd_anomaly, model, rcp, "mcwd_anomaly", "MCWD_")
    temp <- get_raster_info(temp, model, rcp, "temperature", "year_")  
    df_update_rcp <- left_join(mcwd_anomaly, temp, by = c("x", "y", "model", "rcp", "year"))
    df_update_model=rbind(df_update_model,df_update_rcp)}
    df=rbind(df,df_update_model)}

carbon_stock_above=rast(file.path(mydirection, "/data/earthdata_carbon2010/aboveground_biomass_carbon_2010.tif"))
carbon_stock_below=rast(file.path(mydirection, "/data/earthdata_carbon2010/belowground_biomass_carbon_2010.tif"))
carbon_weight=carbon_stock_above+carbon_stock_below
carbon_weight=aggregate(carbon_weight,fact=180,fun=sum)
carbon_weight <- crop(carbon_weight, boundaries)
carbon_weight <- mask(carbon_weight, boundaries)
carbon_weight=carbon_weight/1000000
boundaries <- project(boundaries, carbon_weight)
raster_subregion <- rasterize(boundaries, carbon_weight, field = "Subregion", touches=TRUE)
points <- vect(df[, c("x", "y")], geom = c("x", "y"), crs = crs(carbon_weight))
extracted_vals <- terra::extract(carbon_weight, points)
df$carbon <- extracted_vals[, 2]
extracted_vals <- terra::extract(raster_subregion, points)
df$Subregion <- extracted_vals[,2]
#write.csv(df, file.path(mydirection,"output/final_dataframe.csv"), row.names = FALSE)

#df=read.csv(file.path(mydirection,"output/final_dataframe.csv"))
df=as.data.frame(df)

#some descriptive statistics
summary_stats_general <- df %>%
  summarise(
    mean_temp = mean(temperature, na.rm = TRUE),
    sd_temp = sd(temperature, na.rm = TRUE),
    min_temp = min(temperature, na.rm = TRUE),
    max_temp = max(temperature, na.rm = TRUE),
    mean_mcwd_anomaly = mean(mcwd_anomaly, na.rm = TRUE),
    sd_mcwd_anomaly = sd(mcwd_anomaly, na.rm = TRUE),
    min_mcwd_anomaly = min(mcwd_anomaly, na.rm = TRUE),
    max_mcwd_anomaly = max(mcwd_anomaly, na.rm = TRUE),
    n = n()
  )%>%
  mutate(model = "all models", rcp="all rcp")
summary_stats_model <- df %>%
  group_by(model) %>%
  summarise(
    mean_temp = mean(temperature, na.rm = TRUE),
    sd_temp = sd(temperature, na.rm = TRUE),
    min_temp = min(temperature, na.rm = TRUE),
    max_temp = max(temperature, na.rm = TRUE),
    mean_mcwd_anomaly = mean(mcwd_anomaly, na.rm = TRUE),
    sd_mcwd_anomaly = sd(mcwd_anomaly, na.rm = TRUE),
    min_mcwd_anomaly = min(mcwd_anomaly, na.rm = TRUE),
    max_mcwd_anomaly = max(mcwd_anomaly, na.rm = TRUE),
    n = n()
  ) %>%
  mutate(rcp="all rcp")
summary_stats_rcp <- df %>%
  group_by(rcp) %>%
  summarise(
    mean_temp = mean(temperature, na.rm = TRUE),
    sd_temp = sd(temperature, na.rm = TRUE),
    min_temp = min(temperature, na.rm = TRUE),
    max_temp = max(temperature, na.rm = TRUE),
    mean_mcwd_anomaly = mean(mcwd_anomaly, na.rm = TRUE),
    sd_mcwd_anomaly = sd(mcwd_anomaly, na.rm = TRUE),
    min_mcwd_anomaly = min(mcwd_anomaly, na.rm = TRUE),
    max_mcwd_anomaly = max(mcwd_anomaly, na.rm = TRUE),
    n = n()
  )  %>%
  mutate(model="all models")
summary_stats_combined <- bind_rows(summary_stats_general, summary_stats_model, summary_stats_rcp)
summary_stats_combined <- summary_stats_combined %>%
  mutate(across(where(is.numeric), ~ round(.x, 1)))
summary_stats_combined <- summary_stats_combined %>%
  select(model, rcp, n, mean_temp, sd_temp, min_temp, max_temp, mean_mcwd_anomaly, sd_mcwd_anomaly, min_mcwd_anomaly, max_mcwd_anomaly)

summary_stats_combined=summary_stats_combined[,c("model","rcp","n","mean_temp","sd_temp","mean_mcwd_anomaly","sd_mcwd_anomaly")]
summary_stats_combined$n=round(summary_stats_combined$n,0)/1000
summary_stats_combined$mean_temp=paste0(summary_stats_combined$mean_temp, " (", summary_stats_combined$sd_temp, ")")
summary_stats_combined$mean_mcwd=paste0(summary_stats_combined$mean_mcwd_anomaly, " (", summary_stats_combined$sd_mcwd_anomaly, ")")
summary_stats_combined=summary_stats_combined[,c("model","rcp","n","mean_temp","mean_mcwd")]
colnames(summary_stats_combined)=c("Model","RCP","n (thousands)", "Temperature average (sd)", "MWCD anomaly average (sd)")
summary_stats_combined_sorted <- summary_stats_combined %>%
  arrange(Model, RCP) %>%
  group_by(Model) %>%
  mutate(Model = if_else(row_number() == 1, Model, "")) %>%
  ungroup()

latex_table <- xtable(summary_stats_combined_sorted,
                      label = "tab:summary_stats",
                      align = c("l", "l","l", "r","r","r"))

print(latex_table,
      include.rownames = FALSE,
      file = file.path(mydirection,"output/summary_stat.tex"),
      floating = FALSE,
      sanitize.text.function = identity,
      overwrite = TRUE)

summary_stats_clean <- summary_stats_combined_sorted %>%
  rename(
    Model = Model,
    RCP = RCP,
    Observations = `n (thousands)`,
    Temperature = `Temperature average (sd)`,
    MCWD = `MWCD anomaly average (sd)`
  )

kable_table <- summary_stats_clean %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    digits = 2,
    align = c("l", "l", "r", "r", "r"),
    escape = FALSE
  ) %>%
  kable_styling(
    latex_options = c("hold_position", "scale_down"),
    font_size = 10
  )

save_kable(kable_table, file.path(mydirection,"output/summary_stat.tex"))

df_diff <- df %>%
  mutate(
    year = as.numeric(year),  
    rcp = as.character(rcp),
    model = as.factor(model)) %>%
  filter(year >= 2010 & year <= 2090) %>%
  mutate(period = case_when(
    year >= 2050 & year < 2090 ~ "future",
    year >= 2010 & year < 2050 ~ "past",
  TRUE ~ NA_character_ )) %>%
  filter(!is.na(period)) %>%
  group_by(model, rcp, period, x, y) %>%
  summarise(mean_mcwd = mean(mcwd_anomaly, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = mean_mcwd) %>%
  mutate(mean_diff = future - past)

model_levels <- unique(df_diff$model)
x_limits <- range(df_diff$mean_diff, na.rm = TRUE)

plot_cdf_by_rcp <- function(data, rcp_value, hide_y_axis = FALSE, hide_x_axis = FALSE) {
  p <- data %>%
    filter(rcp == rcp_value) %>%
    ggplot(aes(x = mean_diff, color = model)) +
    stat_ecdf(geom = "step", size = 1) +
    scale_color_brewer(palette = "Dark2", name = "Climate Model") +
    scale_x_continuous(limits = x_limits) +  # Appliquer la même limite x
    labs(
      title = paste("RCP -", rcp_value),
      x = "Difference in Maximum Cumulative Water Deficit (2090–2050 vs 2050–2010)",
      y = if (hide_y_axis) NULL else "Cumulative distribution") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title.y = element_text(size = 18),
      axis.title.x = element_text(size = 18),
      axis.text.y = element_text(size = 16),
      legend.position = "right",
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16))
  if (hide_y_axis) {
    p <- p + theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank())}
    if (hide_x_axis) {
    p <- p + theme(axis.title.x = element_blank())}
  return(p)}

df_diff <- df_diff %>%
  mutate(rcp = case_when(
    rcp == "rcp26" ~ "2.6",
    rcp == "rcp60" ~ "6.0",
    rcp == "rcp85" ~ "8.5",
    TRUE ~ rcp))

rcps <- unique(df_diff$rcp)

plots <- lapply(seq_along(rcps), function(i) {plot_cdf_by_rcp(df_diff, rcps[i], hide_y_axis = i != 1,  hide_x_axis = i != 2)})

final_plot <- patchwork::wrap_plots(plots, nrow = 1, guides = "collect") &
              ggplot2::theme(legend.position = "bottom")

ggsave(file.path(mydirection,"figure/mcwd_anomaly_models.pdf"), plot = final_plot, width = 14, height = 6)

#prepare panel data
#transform annual MCWD anomaly to carbon losses
df$carbon_mcwd_anomaly=df$mcwd_anomaly*mcwd_to_carbon

#mcwd_anomaly is mm
#mcwd to carbon maps carbon losses (MgC=1t) per mm per year per ha
#there are 5 years per period
#total area of subregions boundaries
area_m2 = expanse(boundaries, unit = "m")
total_ha <- sum(area_m2 / 10000, na.rm = TRUE)

#total carbon possibly lost is 75GtC (total dieback, Armstrong McKay, 2022)
df$carbon_mcwd_anomaly=total_ha*5*df$carbon_mcwd_anomaly/75/1e9

df <- df %>%
  mutate(id = paste0("cell_", x, "_", y, "_", rcp))
df <- df %>%
  mutate(id_rcp = paste0("cell_", x, "_", y))
df <- df %>%
  arrange(id, year, model)

df_gfdl <- pdata.frame(df[df$model=="gfdl",], index = c("id", "year"))
df_ipsl <- pdata.frame(df[df$model=="ipsl",], index = c("id", "year"))
df_hadgem <- pdata.frame(df[df$model=="hadgem",], index = c("id", "year"))
df_miroc <- pdata.frame(df[df$model=="miroc",], index = c("id", "year"))

#run some tests
run_panel_diagnostics <- function(df_p, dep_var, indep_var) {
  fmla <- as.formula(paste(dep_var, "~", indep_var))  
  pooled <- plm(fmla, data = df_p, model = "pooling")
  fe <- plm(fmla, data = df_p, model = "within")
  fe_time <- plm(update(fmla, . ~ . + factor(year)), data = df_p, model = "within")
  re <- plm(fmla, data = df_p, model = "random")
  f_test <- tryCatch(pFtest(fe, pooled), error = function(e) NA) #Ftest for FE
  lm_test <- tryCatch(plmtest(pooled, type = "bp"), error = function(e) NA) #LM for RE
  hausman <- tryCatch(phtest(fe, re), error = function(e) NA) #Hausman test
  f_time <- tryCatch(pFtest(fe_time, fe), error = function(e) NA) #F test for time FE
  bp <- tryCatch(bptest(fmla, data = df_p), error = function(e) NA) #Breusch-Pagan heteroskedasticity
  serial <- tryCatch(pwartest(fe), error = function(e) NA) #wooldridge serial correlation
  pesaran <- tryCatch(pcdtest(fe, test = "cd"), error = function(e) NA) #pesaran cross sectional dependence
  result <- tibble::tibble(
  Test = c("F test: FE vs Pooled",
           "LM test: RE vs Pooled",
           "Hausman: FE vs RE",
           "F test: Time FE",
           "Breusch-Pagan: Heteroskedasticity",
           "Wooldridge: Serial Correlation",
           "Pesaran: Cross-sectional Dependence"),
  `Test Statistic` = c(
    ifelse(!is.null(f_test), f_test$statistic, NA),
    ifelse(!is.null(lm_test), lm_test$statistic, NA),
    ifelse(!is.null(hausman), hausman$statistic, NA),
    ifelse(!is.null(f_time), f_time$statistic, NA),
    ifelse(!is.null(bp), bp$statistic, NA),
    ifelse(!is.null(serial), serial$statistic, NA),
    ifelse(!is.null(pesaran), pesaran$statistic, NA)),
  `p-value` = c(
    ifelse(!is.null(f_test), f_test$p.value, NA),
    ifelse(!is.null(lm_test), lm_test$p.value, NA),
    ifelse(!is.null(hausman), hausman$p.value, NA),
    ifelse(!is.null(f_time), f_time$p.value, NA),
    ifelse(!is.null(bp), bp$p.value, NA),
    ifelse(!is.null(serial), serial$p.value, NA),
    ifelse(!is.null(pesaran), pesaran$p.value, NA)))
  return(result)}

diagnostics_gfdl <- run_panel_diagnostics(df_gfdl, dep_var = "carbon_mcwd_anomaly", indep_var = "temperature")
diagnostics_ipsl <- run_panel_diagnostics(df_ipsl, dep_var = "carbon_mcwd_anomaly", indep_var = "temperature")
diagnostics_miroc <- run_panel_diagnostics(df_miroc, dep_var = "carbon_mcwd_anomaly", indep_var = "temperature")
diagnostics_hadgem <- run_panel_diagnostics(df_hadgem, dep_var = "carbon_mcwd_anomaly", indep_var = "temperature")

diagnostics_gfdl$model <- "GFDL"
diagnostics_ipsl$model <- "IPSL"
diagnostics_miroc$model <- "MIROC"
diagnostics_hadgem$model <- "HadGEM"

all_diagnostics <- bind_rows(diagnostics_gfdl, diagnostics_ipsl, diagnostics_miroc, diagnostics_hadgem)

all_diagnostics <- all_diagnostics %>%
  dplyr::select(Model = model, Test, `Test Statistic`, `p-value`) %>%
  mutate(
    Statistic = round(`Test Statistic`, 5),
    pvalue = round(`p-value`, 5))

results <- all_diagnostics %>%
  group_by(Model) %>%
  mutate(Model_disp = if_else(row_number() == 1, Model, "")) %>%
  ungroup() %>%
  dplyr::select(Model = Model_disp, Test,Statistic,pvalue)

xtab <- xtable(results,
               caption = "Panel data tests",
               label = "tab:fe_specs",
               align = c("l","l", "l", "r", "r"))
print(xtab,
      include.rownames = FALSE,
      file = file.path(mydirection,"output/test_fe_table.tex"),
      floating = FALSE,
      sanitize.text.function = identity,
      overwrite = TRUE)

xtab_aer <- xtable(xtab,
                   caption = "Panel data tests",
                   label = "tab:fe_specs",
                   align = c("l","l","l","r","r"))

print(xtab_aer,
      include.rownames = FALSE,
      floating = FALSE,      
      booktabs = TRUE,
      sanitize.text.function = identity,
      file = file.path(mydirection,"output/test_fe_table.tex"))

panel_spec_tests <- all_diagnostics %>%
  filter(Test %in% c("F test: FE vs Pooled", "Hausman: FE vs RE", "F test: Time FE")) %>%
  select(Model, Test, Statistic, pvalue) %>%
  mutate(stat_p = sprintf("%.2f (%.2f)", Statistic, pvalue)) %>%
  select(Model, Test, stat_p) %>%
  pivot_wider(names_from = Test,
              values_from = stat_p) %>%
  dplyr::select(Model, `F test: FE vs Pooled`, `Hausman: FE vs RE`, `F test: Time FE`)

error_quality_tests <- all_diagnostics %>%
  filter(Test %in% c("Breusch-Pagan: Heteroskedasticity", 
                     "Wooldridge: Serial Correlation", 
                     "Pesaran: Cross-sectional Dependence")) %>%
  select(Model, Test, Statistic, pvalue) %>%
  mutate(stat_p = sprintf("%.2f (%.2f)", Statistic, pvalue)) %>%
  select(Model, Test, stat_p) %>%
  pivot_wider(
    names_from = Test,
    values_from = stat_p
  ) %>%
  rename(
    `Heteroskedasticity` = `Breusch-Pagan: Heteroskedasticity`,
    `Serial Correlation` = `Wooldridge: Serial Correlation`,
    `Cross-sectional Dependence` = `Pesaran: Cross-sectional Dependence`
  ) %>%
  dplyr::select(Model,
         `Heteroskedasticity`,
         `Serial Correlation`,
         `Cross-sectional Dependence`)

panel_spec_tests
error_quality_tests

tab <- xtable(panel_spec_tests,
              caption = "Panel Specification Tests",
              label = "tab:panel_spec_tests",
              digits = 2)

cat("\\renewcommand{\\theadfont}{\\normalsize\\bfseries}\n")
cat("\\setlength{\\tabcolsep}{6pt}\n")

print(tab,
      file = file.path(mydirection, "output/panel_spec_tests.tex"),
      type = "latex",
      include.rownames = FALSE,
      sanitize.text.function = identity,
      floating = TRUE,
      table.placement = "ht",
      caption.placement = "top")

tab <- xtable(error_quality_tests,
              caption = "Panel Specification Tests",
              label = "tab:error_quality_tests",
              digits = 2)
print(tab,
      file = file.path(mydirection, "output/error_quality_tests.tex"),
      type = "latex",
      include.rownames = FALSE,
      sanitize.text.function = identity,
      floating = TRUE,
      table.placement = "ht",
      caption.placement = "top")


#run preferred regression #temperature anomaly in log #drop 27
model_list <- list(GFDL = df_gfdl, IPSL = df_ipsl, MIROC = df_miroc, HadGEM = df_hadgem)
results <- data.frame()

temperature_optimal=28
results <- data.frame()
for (name in names(model_list)) {
  df_full <- model_list[[name]]
  #5y per period
  #75 million ha
  #75 GtC that might be loss by tipping
  df_full$log_temperature=df_full$temperature-temperature_optimal
  df_full=df_full[df_full$log_temperature >0, ]
  df_full$log_temperature=log(df_full$log_temperature+1)
  specs <- list(
    "FE" = feols(carbon_mcwd_anomaly ~ log_temperature | id, data = df_full),
    "FE + Time FE" = feols(carbon_mcwd_anomaly ~ log_temperature | id + year, data = df_full),
    "FE + Reg FE" = feols(carbon_mcwd_anomaly ~ log_temperature | id + Subregion, data = df_full),
    "FE + Time + Reg FE" = feols(carbon_mcwd_anomaly ~ log_temperature | id + year + Subregion, data = df_full),
    "FE + Time + Reg FE + weights" = feols(carbon_mcwd_anomaly ~ log_temperature | id + year + Subregion, data = df_full, weights = ~carbon))
  for (spec_name in names(specs)) {
    model <- specs[[spec_name]]
    coef_est <- coef(model)["log_temperature"]
    se_classic <- se(model)["log_temperature"]
    se_dk <- se(model, cluster = c("year","id"))["log_temperature"]  # DK-like clustered SEs
    t_val <- coef_est / se_dk
    results <- rbind(results, data.frame(
      Model = name,
      RCP = "ALL",
      Specification = spec_name,
      Coefficient = round(coef_est, 4),
      SE = round(se_classic, 4),
      DKSE = round(se_dk, 4),
      tvalue = round(t_val, 2)))}

  for (rcp_level in unique(df_full$rcp)) {
    df <- subset(df_full, rcp == rcp_level)
    model <- feols(carbon_mcwd_anomaly ~ log_temperature | id_rcp + year + Subregion, data = df, weights = ~carbon)
    coef_est <- coef(model)["log_temperature"]
    se_classic <- se(model)["log_temperature"]
    se_dk <- se(model, cluster = c("id","year"))["log_temperature"]
    t_val <- coef_est / se_dk
    results <- rbind(results, data.frame(
      Model = name,
      RCP = rcp_level,
      Specification = "FE + Time + Reg FE + weights",
      Coefficient = round(coef_est, 4),
      SE = round(se_classic, 4),
      DKSE = round(se_dk, 4),
      tvalue = round(t_val, 2)))}}

#transform and store final coefficient
#coefficient for the whole rainforest  
#rainforest is 750million ha, coefficient is tC/ha/y, 5y per period
#total possible carbon losses in the forest because of a dieback are 75Gt (10e-9)
#mean coefficient is epsilon
results=results[results$Specification=="FE + Time + Reg FE + weights",]
coefficients=results[results$RCP=="ALL",]
rownames(coefficients) <- NULL

results_tex <- xtable(results,
                      caption = "FE Panel Regression Results for All Specifications and RCPs",
                      label = "tab:results_full",
                      digits = c(0, 0, 0, 0, 4, 4, 4, 2))  # chiffres significatifs

output_file <- file.path(mydirection, "output/results_full.tex")

print(results_tex,
      file = output_file,
      type = "latex",
      include.rownames = FALSE,    # pas de row names
      floating = TRUE,
      table.placement = "ht",
      caption.placement = "top")

#keep across all RCP
table1_data <- results %>%
  filter(RCP == "ALL") %>%
  dplyr::select(Model, Specification, Coefficient, DKSE) %>%
  mutate(coef_dk = sprintf("%.3f\n(%.3f)",Coefficient,DKSE),
         Specification = recode(Specification,
                                "FE" = "Cell",
                                "FE + Time FE" = "Cell + Year",
                                "FE + Reg FE" = "Cell + Reg",
                                "FE + Time + Reg FE" = "Cell + Year + Reg",
                                "FE + Time + Reg FE + weights" = "Cell + Year + Reg + Weights")) %>%
  select(Model, Specification, coef_dk)

table1_wide <- table1_data %>%
  pivot_wider(names_from = Specification, values_from = coef_dk) %>%
  select(Model,
         "Cell",
         "Cell + Year",
         "Cell + Year + Reg",
         "Cell + Year + Reg + Weights")

tab1_tex <- xtable(table1_wide,
                   caption = "Panel Regression Results (RCP = ALL). Coefficient (DKSE) in parentheses.",
                   label = "tab:main_results")
print(tab1_tex,
      file = file.path(mydirection, "output/table1_main_results.tex"),
      type = "latex",
      include.rownames = FALSE,
      sanitize.text.function = identity,
      floating = TRUE)

table2_data <- results %>%
  filter(Specification == "FE + Time + Reg FE + weights") %>%
  dplyr::select(Model, RCP, Coefficient, DKSE)

table2_data <- table2_data %>%
  mutate(coef_dk = sprintf("%.3f\n(%.3f)", Coefficient, DKSE)) %>%
  dplyr::select(Model, RCP, coef_dk)

table2_wide <- table2_data %>%
  pivot_wider(names_from = RCP,
              values_from = coef_dk)

tab2_tex <- xtable(table2_wide,
                   caption = "Preferred Specification Across RCP Scenarios. Driscoll-Kraay SE in parentheses.",
                   label = "tab:rcp_results")

print(tab2_tex,
      file = file.path(mydirection, "output/table2_rcp_results.tex"),
      type = "latex",
      include.rownames = FALSE,
      sanitize.text.function = identity,
      floating = TRUE)

coefficients$Coefficient=abs(coefficients$Coefficient)
epsilon_max=mean(coefficients$Coefficient)
epsilon_max
epsilons=c(epsilon_max = epsilon_max, epsilon_max1=coefficients$Coefficient[1], epsilon_max2=coefficients$Coefficient[2], epsilon_max3 = coefficients$Coefficient[3], epsilon_max4=coefficients$Coefficient[4])

write.table(
  data.frame(name = names(epsilons), value = epsilons),
  file = file.path(mydirection, "output/epsilons.csv"),
  sep = ";",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE)

results_print <- results %>%
  group_by(Model) %>%
  mutate(Model_disp = if_else(row_number() == 1, Model, "")) %>%
  group_by(Model, RCP, .add = TRUE) %>%
  mutate(RCP_disp = if_else(row_number() == 1, RCP, "")) %>%
  ungroup() %>%
  group_by(Model, Specification, .add = TRUE) %>% 
  mutate(Specification_disp = if_else(row_number() == 1, Specification, "")) %>%
  ungroup() %>%
  select(
    Model = Model_disp,
    RCP = RCP_disp,
    Specification = Specification_disp,
    Coefficient, SE, DKSE, tvalue)

xtab <- xtable(results_print,
               caption = "Fixed Effects Panel Estimations by Specification and RCP",
               label = "tab:fe_specs",
               align = c("l","l", "l", "l", "r", "r", "r", "r"))
print(xtab,
      include.rownames = FALSE,
      file = file.path(mydirection,"output/fe_table.tex"),
      floating = FALSE, # ⬅️ évite \begin{table}
      sanitize.text.function = identity,
      overwrite=TRUE)


#approximate temperature corridors from Kriegler with SSP
#use ERF temperature pathways from IPCC (AR6 WG1 7SM)
#arbitraly choose the one that approximate temperature corridors
#scale them to reproduce the central path of Kriegler corridor
Times <- seq(1, 37, by = 1)

ssp_scenarios <- c("119","126","245","370","434","460","534","585")
ERF_list <- lapply(ssp_scenarios, function(ssp) {
  file_path <- file.path(mydirection, paste0("data/temperature_ERF/ERF_ssp", ssp, "_1750-2500.csv"))
    if (ssp == "534") {file_path <- file.path(mydirection, "data/temperature_ERF/ERF_ssp534-over_1750-2500.csv")}
  df <- read.csv(file_path)
  df$total })
names(ERF_list) <- paste0("ERF_ssp", ssp_scenarios)
list2env(ERF_list, envir = .GlobalEnv)

Times_temp <- seq(1750, 2500, by = 1)
temperature2 <- function(time, state, parms) {
  with(as.list(c(state, parms)), {
      temperature0 <- ((ERF_ssp119[time-1748]-ERF_ssp119[time-1749]) + a*x0 - e*g*(x0-y0))/C
      temperature_ocean0<- g*(x0-y0)/Cd
      temperature1 <- ((ERF_ssp126[time-1748]-ERF_ssp126[time-1749]) + a*x1 - e*g*(x1-y1))/C
      temperature_ocean1<- g*(x1-y1)/Cd
      temperature2 <- ((ERF_ssp434[time-1748]-ERF_ssp434[time-1749]) + a*x2 - e*g*(x2-y2))/C
      temperature_ocean2<- g*(x2-y2)/Cd
      temperature3 <- ((ERF_ssp534[time-1748]-ERF_ssp534[time-1749]) + a*x3 - e*g*(x3-y3))/C
      temperature_ocean3<- g*(x3-y3)/Cd
      temperature4 <- ((ERF_ssp245[time-1748]-ERF_ssp245[time-1749]) + a*x4 - e*g*(x4-y4))/C
      temperature_ocean4<- g*(x4-y4)/Cd
      temperature5 <- ((ERF_ssp460[time-1748]-ERF_ssp460[time-1749]) + a*x5 - e*g*(x5-y5))/C
      temperature_ocean5<- g*(x5-y5)/Cd
      temperature6 <- ((ERF_ssp370[time-1748]-ERF_ssp370[time-1749]) + a*x6 - e*g*(x6-y6))/C
      temperature_ocean6<- g*(x6-y6)/Cd
      temperature7 <- ((ERF_ssp585[time-1748]-ERF_ssp585[time-1749]) + a*x7 - e*g*(x7-y7))/C
      temperature_ocean7<- g*(x7-y7)/Cd
    return(list(c(temperature0,temperature_ocean0,temperature1,temperature_ocean1,temperature2,temperature_ocean2,temperature3,temperature_ocean3,temperature4,temperature_ocean4,temperature5,temperature_ocean5,temperature6,temperature_ocean6,temperature7,temperature_ocean7))) })}
yini  <- c(x0=0,y0=0,x1=0,y1=0,x2=0,y2=0,x3=0,y3=0,x4=0,y4=0,x5=0,y5=0,x6=0,y6=0,x7=0,y7=0)
pars  <- c(a=-1.33,C=8.1, Cd=110,e=1.34,g=0.62) 
results2<- ode(y=yini, times=Times_temp, func=temperature2, parms = pars)
cols <- c(2,4,6,8,10,12,14,16)
ssp_names <- c("119","126","434","534","245","460","370","585")
temp_list <- lapply(cols, function(col) cumsum(results2[, col]))
names(temp_list) <- paste0("temp_ssp", ssp_names)
list2env(temp_list, envir = .GlobalEnv)
years <- seq(2000, 2500, by = 1)
temp=cbind(years, temp_ssp119[250:750], temp_ssp126[250:750], temp_ssp434[250:750],temp_ssp534[250:750],temp_ssp245[250:750],temp_ssp460[250:750],temp_ssp370[250:750],temp_ssp585[250:750])
temp=as.data.frame(temp)
temp_2200=temp[1:200,]
colnames(temp_2200)=c('years',"SSP1-1.9","SSP1-2.6","SSP4-3.4","SSP5-3.4","SSP2-4.5","SSP4-6.0","SSP3-7.0","SSP5-8.5")
colnames(temp)=c('years',"SSP1-1.9","SSP1-2.6","SSP4-3.4","SSP5-3.4","SSP2-4.5","SSP4-6.0","SSP3-7.0","SSP5-8.5")
#write.csv(temp, file = "temp_ssp_ipcc_complete.xlsx", row.names = FALSE)

temp2500=temp[17:501,] #SSP434, SSP460, SSP585
row=as.data.frame(rep(1:((nrow(temp2500)/5)), each=5))
temp2500=cbind(temp2500,row)
temp2500=as.data.frame(temp2500)
temp2500$label <- +(!duplicated(temp2500[,10]))
temp2500=as.data.frame(temp2500)
temp2500<-temp2500[temp2500$label == 1,]
temp2500=as.data.frame(temp2500)
temp2500<-temp2500[,1:9]
#write.csv(temp2500, file = file.path(mydirection,"data/temp_ssp_ipcc2500.csv"), row.names = FALSE)
temp_2200=temp2500[1:37,]
write.csv(temp_2200, file = file.path(mydirection,"data/temp_ssp_ipcc2200.csv"), row.names = FALSE)
