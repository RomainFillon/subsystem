#Figures for numerical illustration with flood maps
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

df_change <- df1[ID %in% c("climate_risk", "both_risks"),
                 .(ID, SCCDS_change = 100 * (SCCDS - bench_value) / abs(bench_value))]

state_colors <- c("climate_risk" = "black", "both_risks" = "black")

# Labels personnalisés pour l'axe X
x_labels <- c("climate_risk" = "Climate risk",
              "both_risks" = "Climate + \n Amazon risks")

plot_left <- ggplot(df_change, aes(x =  factor(ID, levels = c("climate_risk", "both_risks")), y = SCCDS_change, fill = ID)) +
  geom_col(width = 0.6, color = "black") +
  geom_text(aes(label = sprintf("%.2f%%", SCCDS_change)),
            vjust = -0.5, color = "black", size = 6, fontface = "bold") +  # taille augmentée
  scale_fill_manual(values = state_colors) +
  scale_x_discrete(labels = x_labels) +
  theme_minimal(base_size = 16) +  # augmente la taille de base
  labs(title = "Increase in social cost of carbon \n when accounting for Amazon feedback",
       x = NULL, y = "% increase relative to baseline model") +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),   # titre plus grand et bold
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = 14),  # labels X plus grands et bold
    axis.text.y = element_text(face = "bold", size = 14),  # labels Y plus grands et bold
    axis.title.y = element_text(face = "bold", size = 16)  # titre Y plus grand et bold
  ) +
  ylim(0, max(df_change$SCCDS_change)*1.2)

plot_left

df_comp <- df1[ID %in% c("climate_risk", "both_risks"),
               .(standardT = SCCDS_temperature,
                 standardA = SCCDS_subsystem,
                 cross = sum(SCCDS_crossT + SCCDS_crossA),
                 cov = sum(SCCDS_covA + SCCDS_covT)),
               by = ID]

# Mise au format long
df_comp_long <- melt(df_comp, id.vars = "ID",
                     variable.name = "Composante", value.name = "Valeur")

df_comp_long$Valeur[df_comp_long$Composante == "standardT"]=df_comp_long$Valeur[df_comp_long$Composante == "standardT"]-data1$SCCDS_temperature[1]

df_comp_long[, Pourcentage := 100 * abs(Valeur) / sum(abs(Valeur)), by = ID]
x_labels <- c("climate_risk" = "Climate risk",
              "both_risks" = "Climate + \n Amazon risks")

# Palette simplifiée
colors_simple <- c("standardT" = "black",  # bleu clair
                   "standardA" = "grey",
                   "cross"     = "grey",
                   "cov"       = "grey")

# Créer une nouvelle colonne pour la légende simplifiée
df_comp_long$Composante_simple <- df_comp_long$Composante
df_comp_long$Composante_simple <- as.character(df_comp_long$Composante_simple)
df_comp_long$Composante_simple[df_comp_long$Composante_simple == "standardT"] <- "Temperature Effect"
df_comp_long$Composante_simple[df_comp_long$Composante_simple == "standardA"] <- "Direct Amazon Effect"
df_comp_long$Composante_simple[df_comp_long$Composante_simple == "cross"] <- "Feedback Amplification"
df_comp_long$Composante_simple[df_comp_long$Composante_simple == "cov"] <- "Interaction Effects"
df_comp_long$Composante_simple <- factor(df_comp_long$Composante_simple)

palette_gray4 <- c(
  "Temperature Effect"        = "#E5E5E5",
  "Direct Amazon Effect"      = "#9E9E9E",
  "Feedback Amplification"    = "#4D4D4D",
  "Interaction Effects"       = "#000000"
)

plot_middle <- ggplot(df_comp_long, aes(x =  factor(ID, levels = c("climate_risk", "both_risks")), y = Pourcentage, fill = Composante_simple)) +
  geom_col(width = 0.6, color = NA)+
  scale_fill_manual(values = palette_gray4)+
  scale_x_discrete(labels = x_labels) +
  theme_minimal(base_size = 16) +
  labs(title = "What drives this increase?",
       x = NULL, y = "Contribution to total increase (%)") +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 16),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 14)) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
combined <- plot_grid(plot_left, plot_middle, nrow = 1, align = "hv")

# Add global title and subtitle
final_plot <- ggdraw() +
  draw_label("Accounting for Amazon Feedback Increases the Social Cost of Carbon",
             x = 0.5, y = 0.97,
             hjust = 0.5,
             fontface = "bold",
             size = 18) +
  draw_label("Decomposition of the increase relative to a stochastic model without Amazon dynamics",
             x = 0.5, y = 0.93,
             hjust = 0.5,
             size = 14) +
  draw_plot(combined, y = 0, height = 0.90)

final_plot
ggsave(file.path(mydirection_data,"figure",paste0("plot_SCCDS_",id_scenario,".pdf")), plot = final_plot, width = 16, height = 8)

df_change <- df1
state_colors <- c("benchmark"="black", "climate_risk" = "black", "both_risks" = "black")

x_labels <- c("benchmark" = "Deterministic",
            "climate_risk" = "Climate risk",
              "both_risks" = "Climate + \n Amazon risks")

df_change$ID <- factor(df_change$ID, 
                       levels = c("benchmark", "climate_risk", "both_risks"))

#ref_value <- df_change$SCD[df_change$ID == "benchmark"]
#df_change$ratio_rel <- df_change$SCD / ref_value
SCC_reference = df_change$SCCDS[df_change$ID=="benchmark"]
SCC_both = df_change$SCCDS[df_change$ID=="both_risks"]

#df_comp_long
df_comp <- df1[,c("ID","SCD_subsystem","SCD_temperature")]
df_comp_long <- melt(df_comp, id.vars = "ID",
                     variable.name = "Composante", value.name = "Valeur")

df_comp_long[, Pourcentage := 100 - 100 * Valeur / sum(abs(Valeur)), by = ID]
x_labels <- c("climate_risk" = "Climate risk only",
              "both_risks" = "Climate + \n Amazon risks")

# Palette simplifiée
colors_simple <- c("SCD_subsystem" = "black",  # bleu clair
                   "SCD_temperature" = "grey")

# Créer une nouvelle colonne pour la légende simplifiée
df_comp_long$Composante <- as.character(df_comp_long$Composante)
df_comp_long$Composante[df_comp_long$Composante=="SCD_subsystem"]= "Propagation via A"
df_comp_long$Composante[df_comp_long$Composante=="SCD_temperature"]= "Propagation via T"
df_comp_long$Composante <- factor(df_comp_long$Composante)

x_labels <- c("benchmark" = "Deterministic",
            "climate_risk" = "Climate risk",
              "both_risks" = "Climate + \n Amazon risks")

# Graphique
df_comp_long=df_comp_long[df_comp_long$ID!="benchmark",]

# reconstruire niveaux à partir des parts
df_stack <- merge(
  df_change[, c("ID", "SCD")],
  df_comp_long,
  by = "ID"
)

df_stack$value <- df_stack$SCD * df_stack$Pourcentage / 100

SCC_reference=df_change$SCCDS[df_change$ID=="benchmark"]
df_change = df_change[df_change$ID=="both_risks",]
df_change$SCD_ratio_to_SCCSD =100* (df_change$SCD-df_change$SCCDS)/df_change$SCCDS
df_change$SCD_ratio_to_SCC = 100*(df_change$SCD-SCC_reference)/SCC_reference

df_long <- df_change %>%
  # Sélectionne seulement les colonnes utiles
  select(ID, SCD_ratio_to_SCCSD, SCD_ratio_to_SCC) %>%
  pivot_longer(
    cols = c(SCD_ratio_to_SCCSD, SCD_ratio_to_SCC),
    names_to = "ratio_type",
    values_to = "ratio"
  ) %>%
  # Remplacer le nom de ID pour SCC et SCCDS selon ratio_type si besoin
  mutate(ID = case_when(
    ratio_type == "SCD_ratio_to_SCCSD" ~ ifelse(ID == "both_risks", "SCCDS", ID),
    ratio_type == "SCD_ratio_to_SCC"    ~ ifelse(ID == "both_risks", "SCC", ID),
    TRUE ~ ID
  )) %>%
  select(ID, ratio)

state_colors <- c("SCC" = "black", "SCCDS" = "black")

# Labels personnalisés pour l'axe X
x_labels <- c("SCC" = "Compared to setting\nwithout Amazon dynamics",
              "SCCDS" = "Compared to setting\nwith Amazon dynamics")

plot_left <- ggplot(df_long, aes(x =  factor(ID, levels = c("SCC", "SCCDS")), y = ratio, fill = ID)) +
  geom_col(width = 0.6, color = "black") +
  geom_text(aes(label = paste0(round(ratio, 0), "%")),
          vjust = -0.5, color = "black", size = 6, fontface = "bold") +
  scale_fill_manual(values = state_colors) +
  scale_x_discrete(labels = x_labels) +
  theme_minimal(base_size = 16) +  # augmente la taille de base
  labs(title = "How much more is a ton of carbon stored in the Amazon worth?\n In comparison with standard carbon emission", x = NULL, y = "Increase in value compared to a standard emission (%)") +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),   # titre plus grand et bold
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = 14),  # labels X plus grands et bold
    axis.text.y = element_text(face = "bold", size = 14),  # labels Y plus grands et bold
    axis.title.y = element_text(face = "bold", size = 16)  # titre Y plus grand et bold
  ) +
  ylim(0, max(df_long$ratio)*1.2)

plot_left

x_labels <- c("Propagation via A" = "Direct Amazon effect",
            "Propagation via T" = "Temperature effect")

palette_gray4 <- c("Direct Amazon effect" = "#9E9E9E",
                   "Temperature effect" = "#4D4D4D")

df_stack$Composante <- as.character(df_stack$Composante)
df_stack$Composante[df_stack$Composante == "Propagation via A"] <- "Direct Amazon effect"
df_stack$Composante[df_stack$Composante == "Propagation via T"] <- "Temperature effect"
df_stack$Composante <- factor(df_stack$Composante)

plot_right <- ggplot(df_stack[df_stack$ID=="both_risks",], aes(x =  factor(ID, levels = c("both_risks")), y = Pourcentage, fill = Composante)) +
  geom_col(width = 0.6, color = NA)+
  scale_fill_manual(values = palette_gray4)+
  theme_minimal(base_size = 16) +
  labs(title = "Contribution of each channel (in %)",
       x = NULL, y = "Share (in %)") +
    scale_x_discrete(labels = c(
    "both_risks"   = "Climate +\n Amazon risks"
  ))+
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 16),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 14)
  )

plot_SCD =plot_grid(plot_left, plot_right, nrow = 1, align = "hv")

# Add global title and subtitle
final_plot <- ggdraw() +
  draw_label("Valuing the carbon stored in the Amazon rainforest",
             x = 0.5, y = 0.97,
             hjust = 0.5,
             fontface = "bold",
             size = 18) +
  draw_label("Comparison with standard carbon value and what drives the increase",
             x = 0.5, y = 0.93,
             hjust = 0.5,
             size = 14) +
  draw_plot(plot_SCD, y = 0, height = 0.90)

ggsave(file.path(mydirection_data,"figure",paste0("plot_SCD_",id_scenario,".pdf")), plot = final_plot, width = 16, height = 8)

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

# check deterministic with and without the rainforest*
myfile=file.path(myfolder,run_det_noamz,"state_V_notstochastic.csv")
col_names <- c("state1", "state2", "state3")
data1 <- fread(myfile, sep=";", header=FALSE)
setnames(data1, col_names)
data1$ID ="Deterministic without Amazon"
data1$time=1:nrow(data1)

# check deterministic with and without the rainforest
myfile=file.path(myfolder,run_det_amz,"state_V_notstochastic.csv")
data2 <- fread(myfile, sep=";", header=FALSE)
setnames(data2, col_names)
data2$ID ="Deterministic with Amazon"
data2$time=1:nrow(data2)
data= rbind(data1, data2)

data$time=2015+5*(data$time-1)

#add stochastic
data_used=data_all[,c("time","scenario","state2mean","state3mean")]
colnames(data_used)= c("time","ID","state2","state3")
data=rbind(data[,c("time","ID","state2","state3")], data_used[data_used$ID=="Climate risk without Amazon"])
data_complete = data

state_colors <- c(
  "Deterministic without Amazon" = "black",
  "Deterministic with Amazon"    = "darkred",   # 
  "Climate risk without Amazon"  = "grey"   # dark grey
)

data$ID <- factor(
  data$ID,
  levels = c("Deterministic without Amazon","Climate risk without Amazon","Deterministic with Amazon")   # ordre souhaité dans la légende
)

data <- as.data.table(data)
ref <- data[ID == "Deterministic without Amazon", .(time, ref_temp = state2)]
data <- merge(data, ref, by = "time")
data[, diff_state2 := state2 - ref_temp]

data = data[data$time==2200,]
data[, Amazon := ifelse(grepl("with_Amazon", ID), "With Amazon", "Without Amazon")]

# Palette gris foncé/clair/noir
#palette_colors <- c("Deterministic_without_Amazon" = "grey60",
#                    "Deterministic_with_Amazon" = "grey30",
#                    "Climate_risk_without_Amazon" = "black")

ref_temp <- data$state2[1]
ref_forest <- data$state3[1] * 100  # En pourcentage pour le graph 2

# Couleurs AER (Noir, Gris moyen, Gris foncé)
aer_colors <- c("Deterministic without Amazon" = "black", 
                "Deterministic with Amazon"    = "grey30", 
                "Climate risk without Amazon"  = "grey60")

# 2. Thème AER avec textes agrandis
aer_theme_large <- theme_minimal() +
  theme(
    # Légende
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    # Axes
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 10)),
    # Grille et Titres
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    strip.text = element_text(size = 18)
  )

# 3. Graphique Température (Gauche)
p1 <- ggplot(data, aes(x = ID, y = state2, fill = ID)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = ref_temp, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_manual(values = aer_colors) +
  coord_cartesian(ylim = c(1.95, 1.99)) + 
  labs(title = "Temperature", 
       y = "Temperature (in °C)", 
       x = NULL) +
  aer_theme_large

# 4. Graphique Forêt (Droite)
# Note : on multiplie state3 par 100 pour avoir le pourcentage
p2 <- ggplot(data, aes(x = ID, y = state3 * 100, fill = ID)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = ref_forest, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_manual(values = aer_colors) +
  coord_cartesian(ylim = c(0, 100)) + 
  labs(title = "Amazon Forest Cover", 
       y = "Forest cover (%)", 
       x = NULL) +
  aer_theme_large

# 5. Assemblage avec Titre Global et Subtitle
combined_plot <- (p1 | p2) + 
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Stationary values for state variables (year 2200)\n",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 5)),
      legend.position = "bottom"
    )
  )

# Affichage
print(combined_plot)
ggsave(file.path(mydirection_data,"figure",paste0("plot_traj1_",id_scenario,".pdf")), plot = combined_plot, width = 16, height = 8)

data1 = data_complete[,c("time","ID","state2","state3")]
data2 = data_all[,c("time","scenario","state2mean","state3mean")]
colnames(data1)=c("time","scenario","state2mean","state3mean")
data=rbind(data1,data2)
#data = data[data$scenario!="Climate risk without Amazon"]

aer_colors <- c(
  "Deterministic with Amazon" = "black",
  "Climate risk with Amazon" = "grey60",
  "Both risks with Amazon"    = "grey30" 
)

data$scenario <- factor(
  data$scenario,
  levels = c("Deterministic with Amazon","Climate risk with Amazon","Both risks with Amazon")   # ordre souhaité dans la légende
)

data = data[data$time==2200,]

aer_theme_large <- theme_minimal() +
  theme(
    # Légende
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    # Axes
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 10)),
    # Grille et Titres
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    strip.text = element_text(size = 18)
  )


ref_temp <- data$state2mean[data$scenario=="Deterministic with Amazon"]
ref_forest <- data$state3mean[data$scenario=="Deterministic with Amazon"] * 100  # En pourcentage pour le graph 2

data=data[!is.na(data$scenario),]

p1 <- ggplot(data, aes(x = scenario, y = state2mean, fill = scenario)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = ref_temp, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_manual(values = aer_colors) +
  coord_cartesian(ylim = c(1.95, 1.99)) + 
  labs(title = "Temperature", 
       y = "Temperature (in °C)", 
       x = NULL) +
  aer_theme_large

# 4. Graphique Forêt (Droite)
# Note : on multiplie state3 par 100 pour avoir le pourcentage
p2 <- ggplot(data, aes(x = scenario, y = state3mean * 100, fill = scenario)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = ref_forest, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_manual(values = aer_colors) +
  coord_cartesian(ylim = c(50, 65)) + 
  labs(title = "Amazon Forest Cover", 
       y = "Forest cover (%)", 
       x = NULL) +
  aer_theme_large

# 5. Assemblage avec Titre Global et Subtitle
combined_plot <- (p1 | p2) + 
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Stationary values for state variables (year 2200)\n",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 5)),
      legend.position = "bottom"
    )
  )

# Affichage
print(combined_plot)

ggsave(file.path(mydirection_data,"figure",paste0("plot_traj2_",id_scenario,".pdf")), plot = combined_plot, width = 16, height = 8)

# 1. Run configuration
scenarios <- c("benchmark", "counterfactual1", "counterfactual2")
run_types <- c("run_det_amz", "run_sto1_amz", "run_sto2_amz")

run_ids <- list(
  "benchmark"       = c(run_det_amz="run0003", run_sto1_amz="run0005", run_sto2_amz="run0006"),
  "counterfactual1" = c(run_det_amz="run0008", run_sto1_amz="run0010", run_sto2_amz="run0011"),
  "counterfactual2" = c(run_det_amz="run0013", run_sto1_amz="run0015", run_sto2_amz="run0016")
)

prefix <- "final_amazon_tcre_"

col_names_base <- c(
  "control", "state1", "state2", "state3",
  "SCD","SCD_temperature","SCD_subsystem", 
  "SCCDS", "SCCDS_temperature", "SCCDS_subsystem",
  "SCCDS_crossT","SCCDS_crossA","SCCDS_covA","SCCDS_covT"
)

# 2. Data extraction function
fetch_data <- function(scen, type_name) {
  
  run_id <- run_ids[[scen]][type_name]
  run_folder <- paste0(prefix, run_id)
  myfile <- file.path(myfolder, run_folder, "outputs_stochastic.csv")
  
  if(!file.exists(myfile)) {
    warning(paste("File not found:", myfile))
    return(data.frame())
  }
  
  dt <- fread(myfile, sep=";", header=FALSE)
  
  # Assign column names
  target_length <- ncol(dt)
  cols <- c(col_names_base, rep(NA, target_length - length(col_names_base)))
  setnames(dt, cols)
  
  dt$time <- 1:nrow(dt)
  
  # Index for year 2200
  t_2200 <- round((2200 - 2015) / 5) + 1
  
  data.frame(
    Scenario = scen,
    RunType = type_name,
    SCD = dt[1, SCD],
    SCCDS = dt[1, SCCDS],
    state2 = dt[t_2200, state2],
    state3 = dt[t_2200, state3]
  )
}

# 3. Load all data
all_results <- expand.grid(Scenario = scenarios, RunType = run_types, stringsAsFactors = FALSE) %>%
  split(1:nrow(.)) %>%
  map_df(~fetch_data(.x$Scenario, .x$RunType))

if(nrow(all_results) == 0) {
  stop("No data could be loaded. Check 'myfolder' and folder names.")
}

# 4. Reshape + compute % differences
results_long <- all_results %>%
  pivot_longer(cols = c(SCD, SCCDS, state2, state3),
               names_to = "Variable",
               values_to = "Value")

bench_values <- results_long %>%
  filter(Scenario == "benchmark") %>%
  select(RunType, Variable, BenchValue = Value)

plot_data <- results_long %>%
  filter(Scenario != "benchmark") %>%
  left_join(bench_values, by = c("RunType", "Variable")) %>%
  mutate(PctDiff = (Value - BenchValue) / BenchValue * 100)

# 5. Clean labels + ordering
plot_data <- plot_data %>%
  mutate(
    RunType = recode(RunType,
      "run_det_amz" = "Stochastic 0",
      "run_sto1_amz" = "Stochastic 1",
      "run_sto2_amz" = "Stochastic 2"
    ),
    RunType = factor(RunType, levels = c("Stochastic 0", "Stochastic 1", "Stochastic 2")),
    Variable = factor(Variable, levels = c("SCD", "SCCDS", "state2", "state3"))
  )

# 6. Plot (independent y-axis per panel)
ggplot(plot_data, aes(x = Scenario, y = PctDiff, fill = Scenario)) +
  geom_col(color = "black", width = 0.7) +
  
  facet_wrap(
    ~ RunType + Variable,
    scales = "free_y",
    ncol = 4
  ) +
  
  scale_fill_manual(values = c(
    "counterfactual1" = "#56B4E9",
    "counterfactual2" = "#E69F00"
  )) +
  
  theme_minimal(base_size = 14) +
  labs(
    title = "Scenario Comparison: Percentage Difference from Benchmark",
    subtitle = "Rows: Stochastic specification | Columns: Variables",
    y = "% Difference vs Benchmark",
    x = ""
  ) +
  
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey90", color = NA),
    panel.spacing = unit(1.2, "lines"),
    panel.border = element_rect(colour = "grey80", fill = NA)
  )


ggsave(file.path(mydirection_data,"figure",paste0("plot_traj2_",id_scenario,".pdf")), plot = combined_plot, width = 16, height = 8)


#graph counterfactuals


scenarios <- c("benchmark", "counterfactual1", "counterfactual2")

run_types <- c("run_det_amz", "run_sto1_amz", "run_sto2_amz")

# Mapping des runs
run_ids <- list(
  "benchmark"       = c(run_det_amz="run0003", run_sto1_amz="run0005", run_sto2_amz="run0006"),
  "counterfactual1" = c(run_det_amz="run0008", run_sto1_amz="run0010", run_sto2_amz="run0011"),
  "counterfactual2" = c(run_det_amz="run0013", run_sto1_amz="run0015", run_sto2_amz="run0016")
)

prefix <- "final_amazon_tcre_"

col_names_base <- c(
  "control", "state1", "state2", "state3",
  "SCD","SCD_temperature","SCD_subsystem", 
  "SCCDS", "SCCDS_temperature", "SCCDS_subsystem",
  "SCCDS_crossT","SCCDS_crossA","SCCDS_covA","SCCDS_covT"
)

# =========================
# 2. Extraction function
# =========================

fetch_data <- function(scen, type_name) {
  
  run_id <- run_ids[[scen]][type_name]
  run_folder <- paste0(prefix, run_id)
  myfile <- file.path(myfolder, run_folder, "outputs_stochastic.csv")
  
  if(!file.exists(myfile)) {
    warning(paste("File not found:", myfile))
    return(data.frame())
  }
  
  dt <- data.table::fread(myfile, sep=";", header=FALSE)
  
  # Naming columns
  target_length <- ncol(dt)
  cols <- c(col_names_base, rep(NA, target_length - length(col_names_base)))
  data.table::setnames(dt, cols)
  
  dt$time <- 1:nrow(dt)
  
  # Index year 2200
  t_2200 <- round((2200 - 2015) / 5) + 1
  
  return(data.frame(
    Scenario = scen,
    RunType = type_name,
    SCD = dt[1, SCD],
    SCCDS = dt[1, SCCDS],
    state2 = dt[t_2200, state2],
    state3 = dt[t_2200, state3]
  ))
}

# =========================
# 3. Compile results
# =========================


all_results <- expand.grid(Scenario = scenarios, RunType = run_types, stringsAsFactors = FALSE) %>%
  split(1:nrow(.)) %>%
  map_df(~fetch_data(.x$Scenario, .x$RunType))

if(nrow(all_results) == 0) {
  stop("No data could be loaded. Check 'myfolder' and folder names.")
}
all_results <- all_results %>%
  mutate(RunType = recode(RunType,
    "run_det_amz" = "Deterministic",
    "run_sto1_amz" = "Climate Risk",
    "run_sto2_amz" = "Climate and Amazon risk"
  ))
# =========================
# 4. Reshape + rename variables
# =========================

results_long <- all_results %>%
  pivot_longer(cols = c(SCD, SCCDS, state2, state3),
               names_to = "Variable",
               values_to = "Value") %>%
  mutate(Variable = recode(Variable,
    "SCD" = "SCD",
    "SCCDS" = "SCCDS",
    "state2" = "Temperature",
    "state3" = "Amazon forest cover"
  ))

# Benchmark values
bench_values <- results_long %>%
  filter(Scenario == "benchmark") %>%
  select(RunType, Variable, BenchValue = Value)

# =========================
# 5. Compute differences + rename scenarios
# =========================

plot_data <- results_long %>%
  filter(Scenario != "benchmark") %>%
  left_join(bench_values, by = c("RunType", "Variable")) %>%
  mutate(
    PctDiff = (Value - BenchValue) / BenchValue * 100,
    Scenario = recode(Scenario,
      "counterfactual1" = "Counterfactual 1",
      "counterfactual2" = "Counterfactual 2"
    )
  )

# =========================
# 6. Plot
# =========================

plot <- ggplot(plot_data, aes(x = Scenario, y = PctDiff, fill = Scenario)) +
  geom_col(color = "black", width = 0.7) +
  
  facet_wrap(~ RunType + Variable, scales = "free_y", ncol = 4) +
  
  scale_fill_manual(values = c(
    "Counterfactual 1" = "#56B4E9",
    "Counterfactual 2" = "#E69F00"
  )) +
  
  theme_minimal(base_size = 14) +
  labs(
    title = "Scenario Comparison: % Difference from Benchmark",
    subtitle = "Columns: Variables | Rows: Model Specifications",
    y = "% Difference vs Benchmark",
    x = ""
  ) +
  
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1.2, "lines"),
    panel.border = element_rect(colour = "grey80", fill = NA)
  )

ggsave(file.path(mydirection_data,"figure",paste0("counterfactual.pdf")), plot = plot, width = 16, height = 8)

