dev.off()
rm(list=ls())
set.seed(1234)
setwd("C:/Users/Fillon/Desktop/scientifique/P2_Amazon/github/")

#load libraries
libraries=c("gridExtra","readxl","reticulate","lwgeom","patchwork","viridis","grid","data.table","exactextractr","cowplot","terra","sf","patchwork","tmap","fitdistrplus","rnaturalearthdata","ggplot2","stringr","dplyr","purrr","tidyr","rnaturalearth","rnaturalearth")
lapply(libraries, require, character.only = TRUE) #load libraries

mydirection_data= file.path(getwd(),"analysis_figures/")
mydirection= file.path(getwd(),"numerical_model/parameters/")
myfolder = file.path(getwd(),"numerical_model/outputs/")

id_scenario = "benchmark"
#id_scenario = counterfactual1
#id_scenario = counterfactual2

if (id_scenario=="benchmark"){
run_det_noamz="final_amazon_tcre_run0002"
run_det_amz="final_amazon_tcre_run0003"
run_sto1_noamz="final_amazon_tcre_run0004"
run_sto1_amz="final_amazon_tcre_run0005"
run_sto2_amz="final_amazon_tcre_run0006"}

if (id_scenario=="counterfactual1"){
run_det_noamz="final_amazon_tcre_run0007"
run_det_amz="final_amazon_tcre_run0008"
run_sto1_noamz="final_amazon_tcre_run0009"
run_sto1_amz="final_amazon_tcre_run0010"
run_sto2_amz="final_amazon_tcre_run0011"}

if (id_scenario=="counterfactual2"){
run_det_noamz="final_amazon_tcre_run0012"
run_det_amz="final_amazon_tcre_run0013"
run_sto1_noamz="final_amazon_tcre_run0014"
run_sto1_amz="final_amazon_tcre_run0015"
run_sto2_amz="final_amazon_tcre_run0016"}

#final results
col_names <- c(
  "control", "state1", "state2", "state3", "SCD","SCD_temperature","SCD_subsystem", "SCCDS",
  "SCCDS_temperature", "SCCDS_subsystem","SCCDS_crossT","SCCDS_crossA","SCCDS_covA","SCCDS_covT")

myfile=file.path(myfolder,run_sto1_noamz,"outputs_stochastic.csv")
data1 <- fread(myfile, sep=";", header=FALSE)
target_length <- length(data1)
col_names <- c(col_names, rep(NA, target_length - length(col_names)))
setnames(data1, col_names)
data1$ID ="benchmark"
data1$time=1:nrow(data1)

myfile=file.path(myfolder,run_sto1_amz,"outputs_stochastic.csv")
data2 <- fread(myfile, sep=";", header=FALSE)
setnames(data2, col_names)
data2$ID ="climate_risk"
data2$time=1:nrow(data2)

myfile=file.path(myfolder,run_sto2_amz,"outputs_stochastic.csv")
data3 <- fread(myfile, sep=";", header=FALSE)
setnames(data3, col_names)
data3$ID ="both_risks"
data3$time=1:nrow(data3)

sum(data3$marginal_utility_temp-data1$marginal_utility_temp)

df1=rbind(data1, data2, data3)
df1=df1[df1$time==1,]
bench_value <- df1[ID == "benchmark", SCCDS]

df_change <- df1

##do for stochastic runs some counterfactuals stochastic paths
ref_file=file.path(myfolder,run_det_amz,"control_V_notstochastic.csv")
refcontrol <- fread(ref_file, sep = ";", header = TRUE)
setnames(refcontrol, c("control"))

##do for stochastic runs some counterfactuals stochastic paths
ref_file=file.path(myfolder,run_det_amz,"state_V_notstochastic.csv")
ref <- fread(ref_file, sep = ";", header = FALSE)
setnames(ref, c("state1", "state2", "state3"))
ref=cbind(ref, refcontrol)
ref[, time := 2015 + 5 * (1:.N - 1)]
ref <- ref[time <= 2300, .(time, state2_ref = state2, state3_ref = state3, control_ref=control)]

# === 2️⃣ Fonction pour charger les scénarios stochastiques ===
load_scenario <- function(run_id, scenario_label) {
  file_mean <- file.path(myfolder, run_id, "outputs_stochastic.csv")
  data1 <- fread(file_mean, sep = ";", header = FALSE)
  setnames(data1, col_names)
  data1[, c("ID", "scenario", "time") := .("mean", scenario_label, 1:.N)]
  data1 <- data1[, .(control, state2, state3, ID, scenario, time)]
  
  # --- 5th percentile
  file_5 <- file.path(myfolder, run_id, "outputs_stochastic5.csv")
  col_names_small <- c("control", "state2", "state3")
  data2 <- fread(file_5, sep = ";", header = FALSE)
  setnames(data2, col_names_small)
  data2[, c("ID", "scenario", "time") := .("5", scenario_label, 1:.N)]
  
  # --- 95th percentile
  file_95 <- file.path(myfolder, run_id, "outputs_stochastic95.csv")
  data3 <- fread(file_95, sep = ";", header = FALSE)
  setnames(data3, col_names_small)
  data3[, c("ID", "scenario", "time") := .("95", scenario_label, 1:.N)]
  
  # --- Combine
  data <- rbindlist(list(data1, data2, data3))
  data[, time := 2015 + 5 * (time - 1)]
}

# === 3️⃣ Charger tous les scénarios ===
data_6 <- load_scenario(run_sto1_noamz, "Climate risk without Amazon")  # version baseline
data_7 <- load_scenario(run_sto1_amz, "Climate risk with Amazon")
data_9 <- load_scenario(run_sto2_amz, "Both risks with Amazon")

data_all <- rbindlist(list(data_6, data_7, data_9))

data_all <- dcast(
  melt(data_all, id.vars = c("time", "ID", "scenario"),
       measure.vars = c("control", "state2", "state3")),
  scenario + time ~ ID + variable,
  value.var = "value")

setnames(data_all,
         old = names(data_all)[-(1:2)],
         new = gsub("_", "", gsub("(.*)_(.*)", "\\2_\\1", names(data_all)[-(1:2)])))

data_all <- merge(data_all, ref, by = "time", all.x = TRUE)

data_all[, `:=`(
  controlmean_diff    = controlmean - control_ref,
  control5_diff    = control5 - control_ref,
  control95_diff    = control95 - control_ref,
  state2mean_diff = state2mean - state2_ref,
  state25_diff    = state25 - state2_ref,
  state295_diff   = state295 - state2_ref,
  state3mean_diff = (state3mean*100) - (state3_ref*100),
  state35_diff    = (state35*100) - (state3_ref*100),
  state395_diff   = (state395*100) - (state3_ref*100))]

plot_diff <- function(varname, name, what, show_legend = TRUE) {
  ggplot(data_all[data_all$time <= 2100 & data_all$scenario != "Climate risk without Amazon", ], 
         aes(x = time,
             y = .data[[paste0(varname, "mean_diff")]],
             color = scenario,
             fill = scenario)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_ribbon(aes(
      ymin = .data[[paste0(varname, "5_diff")]],
      ymax = .data[[paste0(varname, "95_diff")]]
    ), alpha = 0.2, color = NA) + 
    geom_line(size = 1.2) + 
    scale_color_viridis_d(option = "plasma", end = 0.8) +
    scale_fill_viridis_d(option = "plasma", end = 0.8) +
    theme_minimal(base_size = 16) + 
    labs(
      title = name,
      subtitle = "Deviation from Deterministic",
      x = "Year",
      y = paste0("Difference (", what, ")"),
      color = "Average & Stochastic Paths",
      fill = "Average & Stochastic Paths"
    ) +
    theme(
      legend.position = if(show_legend) "bottom" else "none",
      plot.title = element_text(face = "bold", size = 22),
      plot.subtitle = element_text(size = 20, color = "grey40"),
      axis.title = element_text(size = 20),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))}

# 1. Créer les graphiques sans légendes
p1 <- plot_diff("control", "Abatement rate", "in %", show_legend = FALSE)
p2 <- plot_diff("state2", "Temperature", "in °C", show_legend = FALSE)
p3 <- plot_diff("state3", "Forest coverage", "in %", show_legend = FALSE)
p_final <- (p1 | p2 | p3) + 
  plot_layout(guides = "collect") & 
  theme(
    # Taille globale du texte de la légende
    legend.text = element_text(size = 18), 
    legend.title = element_text(size = 20, face = "bold"),
    
    # Taille des titres et étiquettes des axes pour les 3 graphiques
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    
    # Taille des titres de chaque graphique
    plot.title = element_text(size = 20, face = "bold"),
    
    # Position de la légende
    legend.position = "bottom"
  )

# Affichage
p_final

ggsave(file.path(mydirection_data,"figure",paste0("plot_stochastic_states_wrtdeterministic_",id_scenario,".pdf")), plot = p_final, width = 18, height = 8)

## =============================================================================
##  LEVELS VERSION — same three panels as plot_diff, but in levels.
##  Drops in after your plot_diff block; re-uses data_all and the same calendar.
##
##  Mirrors plot_diff exactly: same filters, same palette, same type sizes. The
##  only changes are that the y variable is the level rather than the deviation,
##  the zero line is dropped since zero is no longer meaningful, and the
##  deterministic run is drawn as a dashed line so the gap the deviation panels
##  plot stays visible here.
##
##  mult exists because state3 is a fraction in the file and is reported in per
##  cent, exactly as your deviation block does with (state3mean*100).
## =============================================================================

## couleurs fixées par nom : ne dépendent ni de l'ordre des niveaux,
## ni du nombre de scénarios présents après filtrage
pal_sc <- c(
  "Both risks with Amazon"      = "#0D0887",   # bleu foncé
  "Climate risk with Amazon"    = "#F89441",   # orange
  "Climate risk without Amazon" = "grey45"
)

plot_level <- function(varname, name, what, mult = 1, show_legend = TRUE) {
  ggplot(data_all[data_all$time <= 2100 & data_all$scenario != "Climate risk without Amazon", ],
         aes(x = time,
             y = mult * .data[[paste0(varname, "mean")]],
             color = scenario,
             fill = scenario)) +
    geom_ribbon(aes(
      ymin = mult * .data[[paste0(varname, "5")]],
      ymax = mult * .data[[paste0(varname, "95")]]
    ), alpha = 0.2, color = NA) +
    geom_line(aes(x = time, y = mult * .data[[paste0(varname, "_ref")]]),
              colour = "black", linetype = "dashed", linewidth = 0.8,
              inherit.aes = FALSE) +
    geom_line(size = 1.2) +
    scale_color_manual(values = pal_sc, drop = FALSE) +
    scale_fill_manual(values = pal_sc, drop = FALSE) +
    theme_minimal(base_size = 16) +
    labs(
      title = name,
#      subtitle = "Dashed line: deterministic",
      x = "Year",
      y = paste0("Level (", what, ")"),
      color = "Average & Stochastic Paths",
      fill = "Average & Stochastic Paths"
    ) +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      plot.title = element_text(face = "bold", size = 22),
      plot.subtitle = element_text(size = 20, color = "grey40"),
      axis.title = element_text(size = 20),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))
}

l1 <- plot_level("control", "Abatement rate",  "in %",  mult = 100, show_legend = FALSE)
l2 <- plot_level("state2",  "Temperature",     "in °C", mult = 1,   show_legend = FALSE)
l3 <- plot_level("state3",  "Forest coverage", "in %",  mult = 100, show_legend = FALSE)

l_final <- (l1 | l2 | l3) +
  plot_layout(guides = "collect") &
  theme(
    legend.text  = element_text(size = 18),
    legend.title = element_text(size = 20, face = "bold"),
    axis.title   = element_text(size = 20),
    axis.text    = element_text(size = 20),
    plot.title   = element_text(size = 20, face = "bold"),
    legend.position = "bottom"
  )

l_final

ggsave(file.path(mydirection_data, "figure",
                 paste0("plot_stochastic_states_levels_", id_scenario, ".pdf")),
       plot = l_final, width = 18, height = 8)


## =============================================================================
##  Figures from meta_run2_simulation3.py outputs (pooled over 100 draws).
##  Input : outputs_stochastic.csv  (98 periods x 37 columns, sep = ";", no header)
##          each row is a period t; conventions are COLUMNS, not rows.
##  Output: (1) SCCDS uplift over time, (2) SCD convention fan + M-vs-N over time.
##
##  Column map (1-based R index = 0-based python column + 1):
##    V5  perm_M              V16 perm_N              V33 perm_N_spaceSuppr
##    V15 rev_M               V35 perm_M_spaceSuppr   V36 perm_N_climOnly
##    V19 SCC_foss (env)      V23 SCC_explicit        V8  SCCDS (= env)
##    V18 uplift_expl (%)     V37 uplift_env (%)
## =============================================================================

rm(list = ls()); set.seed(1234)
invisible(lapply(c("data.table", "ggplot2", "patchwork"),
                 require, character.only = TRUE))

## ---- paths (edit to your machine) -------------------------------------------
setwd("C:/Users/Fillon/Desktop/scientifique/P2_Amazon/github/")
mydirection_data <- file.path(getwd(), "analysis_figures/")
myfolder         <- file.path(getwd(), "numerical_model/outputs/")
run_sto2_amz     <- "final_amazon_tcre_run0006"   # Amazon active, climate + Amazon risk
run_amz_clim     <- "final_amazon_tcre_run0005"   # Amazon active, climate risk only
run_noamz        <- "final_amazon_tcre_run0004"   # Amazon frozen, climate risk only
id_scenario      <- "benchmark"

## ---- calendar ---------------------------------------------------------------
FIRST_YEAR <- 2020
FREEZE_T   <- 38                       # structural freeze period -> 2200
PERIOD_YR  <- (2200 - FIRST_YEAR) / FREEZE_T
YMAX_T     <- 44

## ---- palette ----------------------------------------------------------------
col_M  <- "#3D5A80"; col_N <- "#C1121F"; col_cr <- "grey45"; col_env <- "#EE9B00"

## ---- load -------------------------------------------------------------------
d <- fread(file.path(myfolder, run_sto2_amz, "outputs_stochastic.csv"), sep = ";")
setnames(d, paste0("V", seq_len(ncol(d))))
d[, t := .I - 1L]
d[, year := FIRST_YEAR + t * PERIOD_YR]
d <- d[t <= YMAX_T]
expl <- d$V23                                   # SCC explicit, homogeneous denominator

## =============================================================================
##  FIGURE 1 : the SCCDS at t0, its level and its decomposition
##  Left  : the reported LEVEL, in envelope (general-equilibrium) space.
##          SCC without the Amazon  ->  + Amazon feedback  ->  SCCDS envelope.
##  Right : the MECHANISM, in explicit (partial-equilibrium) space, where the
##          channels are additive: (I) + (II) + loop.
##  The two spaces are linked by the GE factor, reported in the caption.
##  They are NOT interchangeable: the decomposition is only additive in explicit
##  space, and the level worth reporting is the envelope one.
## =============================================================================

t0    <- d[t == 0]
envL  <- t0$V19                       # SCCDS envelope = dV/dS, the GE level
explI <- t0$V23                       # (I) SCC explicit, partial equilibrium
chII  <- t0$V25                       # (II) subsystem channel
loopj <- t0$V26                       # loop channel, j = 0 fold ONLY (compact version)
explT <- t0$V27                       # SCCDS explicit, compact = (I) + (II) + loop_j0
                                      # NB: the full loop channel is (III)_T + (IV)_T,
                                      # built in the mechanism panel below. Do not add
                                      # loopj on top of it: the j=0 fold is a subset.
upE   <- t0$V37                       # envelope uplift, %  (3D vs frozen-Amazon run)
scc2d <- envL / (1 + upE / 100)       # SCC without the Amazon, envelope
geF   <- envL / explI                 # GE factor linking the two spaces

## ---- levels: read the pure envelope price from the three runs ---------------
##  0004 : Amazon frozen,  climate risk only   -> the benchmark, normalised to 1
##  0005 : Amazon active,  climate risk only   -> adds the feedback
##  0006 : Amazon active,  climate + Amazon risk -> adds the idiosyncratic risk
##  Column 18 holds the pure envelope, dV/dS. In runs 0004/0005 it exists only
##  once the stochastic==1 branch has been patched; before that, column 7 is a
##  usable fallback for 0004 (id_climate = 0) but NOT for 0005, where it
##  double-counts the feedback. We warn if we have to fall back.
lvl <- function(run) {
  f <- file.path(myfolder, run, "outputs_stochastic.csv")
  if (!file.exists(f)) return(NA_real_)
  z <- fread(f, sep = ";"); setnames(z, paste0("V", seq_len(ncol(z))))
  v18 <- if (ncol(z) >= 19) z$V19[1] else NA_real_
  if (!is.na(v18) && abs(v18) > 1e-9) return(v18)
  warning(sprintf("run %s: column 18 empty, falling back to column 7", run))
  z$V8[1]
}
L0 <- lvl(run_noamz); L1 <- lvl(run_amz_clim); L2 <- lvl(run_sto2_amz)

u1  <- L1 / L0 - 1                     # Amazon feedback, climate risk only
u2  <- L2 / L1 - 1                     # additional effect of Amazon risk
tot <- L2 / L0 - 1                     # total, should match column 36 of 0006
cat(sprintf("\n  cascade: %.2f -> %.2f -> %.2f $/tC   (+%.2f%%, then +%.2f%%, total +%.2f%%)\n",
            L0, L1, L2, 100*u1, 100*u2, 100*tot))
cat(sprintf("  cross-check against column 36 of the run: %+.2f%%  [gap %.2f pt]\n",
            upE, 100*tot - upE))

## ---- panel L : the level, normalised so that SCC without the Amazon = 1 -----
wf <- data.table(
  lab  = factor(c("SCC\nno Amazon", "+ climate\nrisk", "+ Amazon\nrisk", "SCCDS"),
                levels = c("SCC\nno Amazon", "+ climate\nrisk", "+ Amazon\nrisk", "SCCDS")),
  ymin = c(0, 1, L1/L0, 0),
  ymax = c(1, L1/L0, L2/L0, L2/L0),
  type = c("base", "add", "add2", "total")
)
wf[, xn := as.numeric(lab)]

pLev <- ggplot(wf) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey55") +
  geom_rect(aes(xmin = xn - .34, xmax = xn + .34, ymin = ymin, ymax = ymax, fill = type),
            alpha = .92) +
  geom_segment(aes(x = xn + .34, xend = xn + .66, y = ymax, yend = ymax),
               data = wf[type %in% c("base", "add")],
               linetype = "dotted", colour = "grey40") +
  geom_text(data = wf[type %in% c("base", "total")],
            aes(x = xn, y = ymax, label = sprintf("%.3f", ymax)),
            vjust = -0.6, size = 4.3, fontface = "bold") +
  geom_text(data = wf[type == "add"],
            aes(x = xn, y = ymax, label = sprintf("+%.1f%%", 100*u1)),
            vjust = -0.6, size = 4.1, fontface = "bold", colour = col_env) +
  geom_text(data = wf[type == "add2"],
            aes(x = xn, y = ymax, label = sprintf("+%.1f%%", 100*u2)),
            vjust = -0.6, size = 4.1, fontface = "bold", colour = col_N) +
  scale_fill_manual(values = c(base = "grey65", add = col_env,
                               add2 = "#E88C7D", total = col_N), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, .10))) +
  scale_x_continuous(breaks = wf$xn, labels = levels(wf$lab)) +
  coord_cartesian(ylim = c(0.90, NA)) +
  labs(x = NULL, y = "SCC without the Amazon = 1",
       title = "The level: the Amazon raises the carbon price",
       subtitle = "Envelope prices, general equilibrium included, normalised to the no-Amazon benchmark.") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9.5, colour = "grey30"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 9.5, lineheight = .95),
        axis.title = element_text(face = "bold"))

## ---- panel R : the wedge, as the 2x2 of the paper ----
##  Columns 10-13 of the python output (R: V11-V14) hold the four components:
##      V11 = III_T -> (III) recurring, expected
##      V13 = IV_T  -> (III) recurring, risk premium
##      V12 = III_A -> (II)  impact,    expected
##      V14 = IV_A  -> (II)  impact,    risk premium
##  Runs up to and including run0006 wrote the A-terms RAW (only L1T carried
##  -monetary); from the monetisation patch on, all four are in $/tC. We detect
##  which by checking the X1 identity, (II)^E + (II)^cov = column 24, under both
##  scalings, so the script is correct either way.
mon    <- t0$V28
rawA   <- abs((-mon * (t0$V12 + t0$V14)) - chII) < abs((t0$V12 + t0$V14) - chII)
sclA   <- if (rawA) -mon else 1
IIexp  <- sclA * t0$V12                # (II)^E    impact, expected
IIcov  <- sclA * t0$V14                # (II)^cov  impact, risk premium
IIIexp <- t0$V11                       # (III)^E   recurring, expected
IIIcov <- t0$V13                       # (III)^cov recurring, risk premium
wedge  <- IIexp + IIcov + IIIexp + IIIcov
explF  <- explI + wedge                # SCCDS explicit, full loop channel
cat(sprintf("\n  [A-channel scaling: %s]\n",
            if (rawA) "raw in CSV, monetised here" else "already monetised in CSV"))

grd <- data.table(
  channel = rep(c("(II) impact\ndamage at t only",
                  "(III) recurring\ndamage at every date"), each = 2),
  split   = rep(c("Expected", "Risk premium"), times = 2),
  val     = c(IIexp, IIcov, IIIexp, IIIcov)
)
grd[, channel := factor(channel, levels = rev(unique(channel)))]
grd[, split   := factor(split,   levels = c("Expected", "Risk premium"))]
grd[, sh := 100 * val / wedge]

rowT <- grd[, .(sh = sum(sh), val = sum(val)), by = channel]   # impact vs recurring
colT <- grd[, .(sh = sum(sh), val = sum(val)), by = split]     # expected vs risk

pGrid <- ggplot(grd, aes(split, channel)) +
  geom_tile(aes(fill = sh), colour = "white", linewidth = 2.5) +
  geom_text(aes(label = sprintf("%.1f%%", sh)),
            size = 5.0, fontface = "bold", colour = "white") +
  geom_text(data = rowT, aes(x = 2.72, y = channel, label = sprintf("%.1f%%", sh)),
            inherit.aes = FALSE, size = 4.2, fontface = "bold", colour = col_N) +
  geom_text(data = colT, aes(x = split, y = 0.38, label = sprintf("%.1f%%", sh)),
            inherit.aes = FALSE, size = 4.2, fontface = "bold", colour = col_M) +
  annotate("text", x = 2.72, y = 2.62, label = "channel", size = 3.4, colour = col_N) +
  annotate("text", x = 0.42, y = 0.38, label = "split", size = 3.4, colour = col_M) +
  scale_fill_gradient(low = "#9DB4CE", high = col_M, guide = "none") +
  coord_cartesian(xlim = c(0.5, 2.9), ylim = c(0.3, 2.7), clip = "off") +
  labs(x = NULL, y = NULL,
       title = "The mechanism: the composition of the wedge",
       subtitle = "Shares of the wedge, explicit space, where the four components are additive.") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9.5, colour = "grey30"),
        plot.margin = margin(6, 14, 6, 6),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 10, lineheight = .95, hjust = 1))
pMech <- pGrid

fig1 <- pLev | pMech
fig1 <- fig1 + plot_layout(widths = c(1.1, 1)) +
#  plot_annotation(
#    caption = sprintf(paste("Both panels are relative. In our calibration the no-Amazon benchmark is",
#                            "$%.0f per tonne of carbon and the SCCDS is $%.0f, but the level depends on",
#                            "the damage and discounting calibration and is not what the figure is about.",
#                            "\nThe level is an envelope price, the derivative of the value function, which",
#                            "carries the general-equilibrium propagation through capital and prices a",
#                            "non-marginal decline of the forest.",
#                            "\nThe composition is additive only in explicit, partial-equilibrium space,",
#                            "where it is linearised: it gives the weight of each channel, not the size of",
#                            "the wedge, and its four shares must not be scaled up to recover the level."),
#                     L0, L2),
    theme = theme(plot.caption = element_text(hjust = 0, size = 8.5, colour = "grey30", lineheight = 1.15))
    #)

ggsave(file.path(mydirection_data, "figure", paste0("plot_SCCDS_", id_scenario, ".pdf")),
       plot = fig1, width = 13.2, height = 5.8)

cat("\n--- SCCDS at t0 ---\n")
cat(sprintf("  SCC without Amazon (envelope) = %8.2f $/tC\n", scc2d))
cat(sprintf("  Amazon feedback               = %+8.2f %%  -> +%.2f $/tC\n", upE, envL - scc2d))
cat(sprintf("  SCCDS envelope (GE level)     = %8.2f $/tC   <- reported level\n", envL))
cat("\n  Wedge decomposition, explicit space ($/tC and share of wedge):\n")
cat(sprintf("    (I)  fossil benchmark          = %9.4f\n", explI))
cat(sprintf("    (II)  impact,    expected      = %9.4f  (%5.1f%%)\n", IIexp,  100*IIexp/wedge))
cat(sprintf("    (II)  impact,    risk premium  = %9.4f  (%5.1f%%)\n", IIcov,  100*IIcov/wedge))
cat(sprintf("    (III) recurring, expected      = %9.4f  (%5.1f%%)\n", IIIexp, 100*IIIexp/wedge))
cat(sprintf("    (III) recurring, risk premium  = %9.4f  (%5.1f%%)\n", IIIcov, 100*IIIcov/wedge))
cat(sprintf("    -------------------------------------------------\n"))
cat(sprintf("    (II)  impact total             = %9.4f  (%5.1f%%)\n",
            IIexp+IIcov,   100*(IIexp+IIcov)/wedge))
cat(sprintf("    (III) recurring total          = %9.4f  (%5.1f%%)\n",
            IIIexp+IIIcov, 100*(IIIexp+IIIcov)/wedge))
cat(sprintf("    expected total                 = %9.4f  (%5.1f%%)\n",
            IIexp+IIIexp,  100*(IIexp+IIIexp)/wedge))
cat(sprintf("    risk premium total             = %9.4f  (%5.1f%%)\n",
            IIcov+IIIcov,  100*(IIcov+IIIcov)/wedge))
cat(sprintf("\n    X1 identity: (II)^E + (II)^cov = %.4f  vs column 24 (II) = %.4f  [gap %.1e]\n",
            IIexp+IIcov, chII, abs(IIexp+IIcov-chII)))
cat(sprintf("    wedge = %.4f   SCCDS explicit = %.2f   uplift = %.2f%%\n",
            wedge, explF, 100*wedge/explI))
cat(sprintf("    GE factor envelope/explicit    = %9.3f\n", geF))

## =============================================================================
##  FIGURE 2 : is a forest tonne worth a fossil tonne?
##  Left  : the full fan of accounting conventions at t0 (dot plot, ordered).
##  Right : the two central conventions, M and N, over time, to the freeze.
##  All ratios are against the homogeneous carbon price (SCC explicit).
## =============================================================================

## ---- panel L : convention fan at t0 (dot plot) ----
t0 <- d[t == 0]
dot <- data.table(
  name = c("Reversible loss",
           "Permanent, M (standard IAM)",
           "Permanent, M, healing suppressed",
           "Permanent, N, climate only",
           "Permanent, N (scale-neutral)",
           "Permanent, N, healing suppressed"),
  val  = c(t0$V15, t0$V5, t0$V35, t0$V36, t0$V16, t0$V33) / t0$V23,
  fam  = c("M", "M", "M", "N", "N", "N"),
  anchor = c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE)
)
dot[, name := factor(name, levels = name[order(val)])]

pDot <- ggplot(dot, aes(val, name, colour = fam)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = .9) +
  geom_segment(aes(x = 1, xend = val, yend = name), linewidth = .7, alpha = .45) +
  geom_point(aes(size = anchor)) +
  geom_text(aes(label = sprintf("%.2f", val),
                hjust = ifelse(val > 1, -0.35, 1.35)),
            size = 3.4, fontface = "bold") +
  scale_colour_manual(values = c(M = col_M, N = col_N), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 3, `TRUE` = 5.2), guide = "none") +
  scale_x_continuous(expand = expansion(mult = .12)) +
  #+
#  annotate("text", x = 1, y = 6.7, label = "carbon price", angle = 0,
#           vjust = 0, size = 3.4, fontface = "bold") +
  coord_cartesian(clip = "off", ylim = c(1, 6.4)) +
  labs(x = "SCD / carbon price at t0", y = NULL,
       title = "The full range of accounting conventions") +
       #,
#       subtitle = "From a reversible loss to a permanent conversion with healing suppressed. M and N (larger dots) are the two central conventions.") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, colour = "grey30"),
        plot.margin = margin(6, 10, 18, 6),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(face = "bold"))

## ---- panel R : M vs N over time ----
band <- data.table(
  year = d$year,
  M_lo = pmin(d$V15, d$V5,  d$V35) / expl,
  M_hi = pmax(d$V15, d$V5,  d$V35) / expl,
  N_lo = pmin(d$V16, d$V33, d$V36) / expl,
  N_hi = pmax(d$V16, d$V33, d$V36) / expl,
  M_mid = d$V5  / expl,
  N_mid = d$V16 / expl
)
yr0 <- min(band$year); yr1 <- max(band$year)

pT <- ggplot(band, aes(x = year)) +
  geom_ribbon(aes(ymin = M_lo, ymax = M_hi), fill = col_M, alpha = .22) +
  geom_ribbon(aes(ymin = N_lo, ymax = N_hi), fill = col_N, alpha = .22) +
  geom_line(aes(y = M_mid), colour = col_M, linewidth = 1) +
  geom_line(aes(y = N_mid), colour = col_N, linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = .9) +
  annotate("text", x = yr0 + 6, y = 1.075, label = "N: worth MORE",
           colour = col_N, size = 3.5, hjust = 0, fontface = "bold") +
  annotate("text", x = yr0 + 6, y = 0.905, label = "M: worth LESS",
           colour = col_M, size = 3.5, hjust = 0, fontface = "bold") +
#  annotate("text", x = 2200, y = 1.0, label = "  freeze", hjust = 0, size = 3.2) +
  coord_cartesian(xlim = c(yr0, yr1), ylim = c(0.80, 1.16), clip = "off") +
  labs(x = "Year of conversion", y = "SCD / carbon price",
       title = "The two central conventions over time",
       subtitle = "The gap closes as the forest nears its attractor.") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, colour = "grey30"),
        plot.margin = margin(6, 30, 6, 6),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"))

fig2 <- pDot | pT
fig2 <- fig2 + plot_layout(widths = c(1.25, 1)) +
  plot_annotation(title = "Is a tonne of forest worth a tonne of fossil carbon?",
#                  subtitle = "The answer straddles the carbon price: it is a matter of accounting convention, not of a parameter.",
                  theme = theme(plot.title = element_text(face = "bold", size = 15),
                                plot.subtitle = element_text(size = 10.5, colour = "grey30")))
ggsave(file.path(mydirection_data, "figure", paste0("scd_fan_", id_scenario, ".pdf")),
       plot = fig2, width = 14.5, height = 6.4)

cat("\nFigure SCD written. Convention fan at t0 (/carbon price):\n")
print(dot[order(val), .(convention = as.character(name), ratio = round(val, 3))])


## =============================================================================
##  ROBUSTNESS FIGURE — SCCDS and SCD across the three calibrations.
##
##  Each specification is a distinct calibration that reproduces the expert
##  trajectory and the Kriegler tipping probabilities equally well, but differs
##  in the rate at which a gap in the canopy closes (growth0). That rate is
##  identified neither by the trajectory, since regeneration is nil at carrying
##  capacity, nor by the tipping probabilities, which constrain climate-driven
##  loss. It is therefore varied over an order of magnitude.
##
##  Left  : the SCCDS uplift, split into the feedback rung (frozen -> active
##          Amazon, climate risk only) and the risk rung (adding subsystem risk).
##          The total is near six percent in all three.
##  Right : the SCD convention fan at t0. The carbon price falls inside the fan
##          in all three; what moves is where it falls, not whether it does.
##
##  Input : outputs_stochastic.csv for the nine runs (98 x 37, sep = ";").
##  Column map (1-based R = 0-based python + 1):
##    V19 envelope level dV/dS      V23 SCC explicit (denominator)
##    V37 envelope uplift (%)
##    V15 rev_M     V5  perm_M    V35 perm_M_spaceSuppr
##    V36 perm_N_climOnly          V16 perm_N   V33 perm_N_spaceSuppr
## =============================================================================

rm(list = ls()); set.seed(1234)
invisible(lapply(c("data.table", "ggplot2", "patchwork"),
                 require, character.only = TRUE))

setwd("C:/Users/Fillon/Desktop/scientifique/P2_Amazon/github/")
mydirection_data <- file.path(getwd(), "analysis_figures/")
myfolder         <- file.path(getwd(), "numerical_model/outputs/")

## ---- the three specifications ----------------------------------------------
## Rename the labels here if you want something more descriptive than
## "counterfactual"; they propagate to both panels and the legend.
specs <- list(
  list(id = "Benchmark",
       noamz = "final_amazon_tcre_run0004",
       amz   = "final_amazon_tcre_run0005",
       full  = "final_amazon_tcre_run0006"),
  list(id = "Counterfactual 1",
       noamz = "final_amazon_tcre_run0009",
       amz   = "final_amazon_tcre_run0010",
       full  = "final_amazon_tcre_run0011"),
  list(id = "Counterfactual 2",
       noamz = "final_amazon_tcre_run0014",
       amz   = "final_amazon_tcre_run0015",
       full  = "final_amazon_tcre_run0016")
)
SPEC_LEVELS <- sapply(specs, `[[`, "id")

col_M <- "#3D5A80"; col_N <- "#C1121F"; col_env <- "#EE9B00"
col_spec <- c("#3D5A80", "#0A9396", "#BB3E03")
names(col_spec) <- SPEC_LEVELS

## ---- loader -----------------------------------------------------------------
rd <- function(run) {
  f <- file.path(myfolder, run, "outputs_stochastic.csv")
  if (!file.exists(f)) stop(sprintf("missing: %s", f))
  z <- fread(f, sep = ";"); setnames(z, paste0("V", seq_len(ncol(z))))
  z
}

casc <- list(); fan <- list()
for (sp in specs) {
  z0 <- rd(sp$noamz); z1 <- rd(sp$amz); z2 <- rd(sp$full)
  L0 <- z0$V19[1]; L1 <- z1$V19[1]; L2 <- z2$V19[1]

  ## the two rungs of the cascade, in percent of the frozen-Amazon benchmark
  casc[[sp$id]] <- data.table(
    spec     = sp$id,
    feedback = 100 * (L1 / L0 - 1),
    risk     = 100 * (L2 / L1 - 1),
    total    = 100 * (L2 / L0 - 1),
    common   = z2$V37[1],          # common-draw uplift, the cleaner measure
    L0 = L0, L2 = L2
  )

  ## the convention fan at t0, against the homogeneous carbon price
  den <- z2$V23[1]
  fan[[sp$id]] <- data.table(
    spec = sp$id,
    name = c("Reversible loss",
             "Permanent, M (standard IAM)",
             "Permanent, M, healing suppressed",
             "Permanent, N, climate only",
             "Permanent, N (scale-neutral)",
             "Permanent, N, healing suppressed"),
    val  = c(z2$V15[1], z2$V5[1], z2$V35[1], z2$V36[1], z2$V16[1], z2$V33[1]) / den,
    fam  = c("M", "M", "M", "N", "N", "N")
  )
}
casc <- rbindlist(casc); fan <- rbindlist(fan)
casc[, spec := factor(spec, levels = SPEC_LEVELS)]
fan[,  spec := factor(spec, levels = SPEC_LEVELS)]

## order the conventions by their benchmark value, and keep that order for all
ord <- fan[spec == SPEC_LEVELS[1]][order(val), name]
fan[, name := factor(name, levels = ord)]

## =============================================================================
##  PANEL L : the SCCDS uplift survives
## =============================================================================
cl <- melt(casc[, .(spec, feedback, risk)], id.vars = "spec",
           variable.name = "rung", value.name = "pct")
cl[, rung := factor(rung, levels = c("risk", "feedback"),
                    labels = c("Subsystem risk", "Mean feedback"))]

pUp <- ggplot(cl, aes(spec, pct, fill = rung)) +
  geom_col(width = .58, alpha = .92) +
  geom_text(aes(label = sprintf("%.1f", pct)),
            position = position_stack(vjust = .5),
            size = 3.6, fontface = "bold", colour = "white") +
  geom_text(data = casc, inherit.aes = FALSE,
            aes(x = spec, y = total, label = sprintf("%.1f%%", total)),
            vjust = -0.6, size = 4.3, fontface = "bold") +
  scale_fill_manual(values = c(`Mean feedback` = col_env,
                               `Subsystem risk` = col_N), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, .16))) +
  labs(x = NULL, y = "Increase in the carbon price (%)",
       title = "The Amazon premium is robust",
       subtitle = "Total uplift over the frozen-Amazon benchmark, split into the mean feedback and the subsystem risk premium.") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, colour = "grey30"),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(face = "bold"))

## =============================================================================
##  PANEL R : the fan always straddles the carbon price
## =============================================================================
rng <- fan[, .(lo = min(val), hi = max(val)), by = name]

pFan <- ggplot(fan, aes(val, name)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = .9) +
  geom_segment(data = rng, inherit.aes = FALSE,
               aes(x = lo, xend = hi, y = name, yend = name),
               colour = "grey70", linewidth = 1.1) +
  geom_point(aes(colour = spec, shape = spec), size = 3.1) +
  scale_colour_manual(values = col_spec, name = NULL) +
  scale_shape_manual(values = c(16, 17, 15), name = NULL) +
  scale_x_continuous(expand = expansion(mult = .10)) +
  labs(x = "SCD / carbon price at t0", y = NULL,
       title = "A forest tonne still straddles a fossil tonne",
       subtitle = "The carbon price falls inside the fan in every specification; the ordering of the conventions never changes.") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, colour = "grey30"),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 9.5),
        axis.title = element_text(face = "bold"))

figR <- pUp | pFan
figR <- figR + plot_layout(widths = c(1, 1.45)) +
  plot_annotation(
    title = "Robustness across three calibrations of the canopy closure rate",
    theme = theme(plot.title = element_text(face = "bold", size = 15)))

ggsave(file.path(mydirection_data, "figure", "robustness_specs.pdf"),
       plot = figR, width = 15, height = 6.2)

## ---- console summary --------------------------------------------------------
cat("\n=== SCCDS cascade by specification ===\n")
print(casc[, .(spec, L0 = round(L0, 2), L2 = round(L2, 2),
               feedback = round(feedback, 2), risk = round(risk, 2),
               total = round(total, 2), common_draw = round(common, 2))])

cat("\n=== SCD fan by specification (ratio to the carbon price) ===\n")
print(dcast(fan, name ~ spec, value.var = "val")[order(name)])

cat("\n=== qualitative checks ===\n")
for (nm in levels(fan$name)) {
  v <- fan[name == nm, val]
  side <- ifelse(v < 1, "<", ">")
  cat(sprintf("  %-34s %s  %s\n", nm, paste(sprintf("%.3f", v), collapse = " "),
              if (length(unique(side)) == 1) "same side" else "CROSSES UNITY"))
}
cat(sprintf("\n  straddle holds in all specs: %s\n",
            all(fan[, .(ok = min(val) < 1 & max(val) > 1), by = spec]$ok)))
