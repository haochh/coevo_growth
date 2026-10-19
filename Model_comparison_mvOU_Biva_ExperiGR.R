load("/Users/user/Chunhui/Life_history/Code/Result/mvOU/mvOU_BAYOU_SyyDiagonal_niter100.RData")

library(tidyverse)
library(ggpubr)
library(ggh4x)
library(patchwork)


#----------------------------------------------------------#
#                       0 Preparation                      #
#----------------------------------------------------------#

### 1. Define predictor names ---
predictors <- c("Metabolic_capacity", "Metabolic_versatility", 
                "Metabolic_completeness", "GC_percentage", "Genome_size", "rRNA_operon_copy_number", 
                "Niche_breadth", "Oxygen", "Temperature", "Salinity", "pH")

### 2. Define Atype ---
Atype <- c("Diagonal", "UpperTri", "LowerTri", "Any")

### 3. Define regime type ---
Rtype <- c("null","regime_trait1", "regime_trait2", "regime_trait1_trait2")

### 4. Define error level ---
error <- c("low","mid","high")


#---------------------------------------------------------------------------------#
#                1 Model comparison; Syy = Diagonal; without M.error              #
#---------------------------------------------------------------------------------#

### 1. Load result: for each growth-trait pair, 100 times of model fitting were conducted for each combination of Atype and Regime type. 
load("mvOU_BAYOU_SyyDiagonal_niter100.RData")

### 2. Get the AIC, AICc, BIC etc for all iterations
AIC_mvOU_BAYOU_SyyDiagonal_niter100 <- data.frame()
for (var in predictors) {
  for (A in Atype) {
    for (R in Rtype) {
      for (i in 1:100) {
        aic.c <- mvOU_BAYOU_SyyDiagonal_niter100[[var]][[i]][["final"]][[A]][[R]][["model"]][["FinalFound"]][["ParamSummary"]][["aic.c"]]
        aic <- mvOU_BAYOU_SyyDiagonal_niter100[[var]][[i]][["final"]][[A]][[R]][["model"]][["FinalFound"]][["ParamSummary"]][["aic"]]
        bic <- mvOU_BAYOU_SyyDiagonal_niter100[[var]][[i]][["final"]][[A]][[R]][["model"]][["FinalFound"]][["ParamSummary"]][["bic"]]
        sic <- mvOU_BAYOU_SyyDiagonal_niter100[[var]][[i]][["final"]][[A]][[R]][["model"]][["FinalFound"]][["ParamSummary"]][["sic"]]
        LogLik <- mvOU_BAYOU_SyyDiagonal_niter100[[var]][[i]][["final"]][[A]][[R]][["model"]][["FinalFound"]][["ParamSummary"]][["LogLik"]]
        df <- data.frame(
          Predictor = var,
          Atype = A,
          Rtype = R,
          nIter = i,
          AICc = aic.c,
          AIC = aic,
          BIC = bic,
          SIC = sic,
          LogLik = LogLik
        )
        AIC_mvOU_BAYOU_SyyDiagonal_niter100 <- rbind(AIC_mvOU_BAYOU_SyyDiagonal_niter100, df)
      }
    }
  }
  
}

### 3. For each trait pair and each combination of Atype and Regime type, find the solution with the lowest AICc across all iterations
AIC_mvOU_BAYOU_SyyDiagonal_niter100_best.AICc <- AIC_mvOU_BAYOU_SyyDiagonal_niter100 %>%
  group_by(Predictor, Atype, Rtype) %>%
  slice_min(AICc, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Syytype = "Diagonal")


#-----------------------------------------------------------------------------#
#                2 Model comparison; Syy = Diagonal; with M.error             #
#-----------------------------------------------------------------------------#
### 1. Load result: for each growth-trait pair and each measurement error setting, 100 times of model fitting were conducted for each combination of Atype and Regime type. 
load("mvOU_BAYOU_SyyDiagonal_niter100_ME.RData")

### 2. Get the AIC, AICc, BIC etc for all iterations
AIC_mvOU_BAYOU_SyyDiagonal_niter100_ME <- data.frame()
for (err in error) {
  for (var in predictors) {
    for (A in Atype) {
      for (R in Rtype) {
        for (i in 1:100) {
          aic.c <- mvOU_BAYOU_SyyDiagonal_niter100_ME[[var]][[err]][[i]][["final"]][[A]][[R]][["model"]][["FinalFound"]][["ParamSummary"]][["aic.c"]]
          aic <- mvOU_BAYOU_SyyDiagonal_niter100_ME[[var]][[err]][[i]][["final"]][[A]][[R]][["model"]][["FinalFound"]][["ParamSummary"]][["aic"]]
          bic <- mvOU_BAYOU_SyyDiagonal_niter100_ME[[var]][[err]][[i]][["final"]][[A]][[R]][["model"]][["FinalFound"]][["ParamSummary"]][["bic"]]
          sic <- mvOU_BAYOU_SyyDiagonal_niter100_ME[[var]][[err]][[i]][["final"]][[A]][[R]][["model"]][["FinalFound"]][["ParamSummary"]][["sic"]]
          LogLik <- mvOU_BAYOU_SyyDiagonal_niter100_ME[[var]][[err]][[i]][["final"]][[A]][[R]][["model"]][["FinalFound"]][["ParamSummary"]][["LogLik"]]
          df <- data.frame(
            Predictor = var,
            Error = err,
            Atype = A,
            Rtype = R,
            nIter = i,
            AICc = aic.c,
            AIC = aic,
            BIC = bic,
            SIC = sic,
            LogLik = LogLik
          )
          AIC_mvOU_BAYOU_SyyDiagonal_niter100_ME <- rbind(AIC_mvOU_BAYOU_SyyDiagonal_niter100_ME, df)
        }
      }
    }
  }
}

### 3. Find the solution with the lowest AICc across all iterations
AIC_mvOU_BAYOU_SyyDiagonal_niter100_ME_best.AICc <- AIC_mvOU_BAYOU_SyyDiagonal_niter100_ME %>%
  group_by(Predictor, Error, Atype, Rtype) %>%
  slice_min(AICc, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Syytype = "Diagonal")


#----------------------------------------------------------------------------------------------#
#               3 Model comparison across different error settings; Syy = Diagonal             #
#----------------------------------------------------------------------------------------------#
### 1. Error = none 
# Create data frame for visualization
df_plot <- AIC_mvOU_BAYOU_SyyDiagonal_niter100_best.AICc %>%
  group_by(Predictor) %>%
  mutate(delta_AICc = AICc - min(AICc)) %>%
  ungroup()

# rename the predictors
df_plot <- df_plot %>%
  mutate(Predictor = case_when(
    Predictor == "Metabolic_capacity" ~ "Metabolic capacity",
    Predictor == "Genome_size" ~ "Genome size",
    Predictor == "Metabolic_versatility" ~ "Metabolic versatility",
    Predictor == "Metabolic_completeness" ~ "Metabolic completeness",
    Predictor == "GC_percentage" ~ "GC percentage",
    Predictor == "Niche_breadth" ~ "Niche breadth",
    Predictor == "Oxygen" ~ "Oxygen tolerance",
    Predictor == "Temperature" ~ "Opt.Temperature",
    Predictor == "Salinity" ~ "Opt.Salinity",
    Predictor == "pH" ~ "Opt.pH",
    Predictor == "rRNA_operon_copy_number" ~ "rrn Copy number"
  ))

# reorder the traits
df_plot$Predictor <- factor(df_plot$Predictor, levels = 
                              c("Metabolic completeness","Metabolic versatility", "Metabolic capacity",
                                "rrn Copy number","GC percentage", "Genome size", "Niche breadth", "Oxygen tolerance",
                                "Opt.Temperature", "Opt.Salinity", "Opt.pH"))

# reorder the regime
df_plot$Rtype <- factor(df_plot$Rtype, levels = 
                          rev(c("regime_trait1","regime_trait2","regime_trait1_trait2","null")))


### Figure
p.none <- ggplot(df_plot, aes(x = Atype, y = Rtype, fill = delta_AICc)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(delta_AICc, 1)), size = 4.5, color = "white") +
  scale_fill_viridis_c(option = "D", direction = -1) +
  geom_point(
    data = df_plot %>% filter(delta_AICc == 0),
    shape = 21, size = 4.5, stroke = 1.2, color = "black", fill = "white"
  ) +
  facet_wrap2(~ Predictor, scales = "fixed", axes = "margins", ncol = 4) +
  scale_y_discrete(labels = c(
    "null" = "No regime",
    "regime_trait1" = "Growth rate",
    "regime_trait2" = "Trait",
    "regime_trait1_trait2" = expression("Growth rate" * " × " * "Trait")
  )) +
  scale_x_discrete(labels = c(
    "Any" = "Bidirectional",
    "Diagonal" = "No coevolution",
    "LowerTri" = "Growth → Trait",
    "UpperTri" = "Trait → Growth"
  )) +
  labs(
    title = "No measurement error",
    x = "Mode of coevolution",
    y = "Regime definition",
    fill = expression(Delta*AIC[c])
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 18,       
      face = "bold",   
      hjust = 0.5      
    ),
    panel.grid = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 14   # increase this for larger panel titles
    ),
    axis.title.x = element_text(size = 14),  #  x-axis title size
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(
      angle = 30,
      hjust = 0.5,
      vjust = 0.7,
      size = 11
    ),
    axis.text.y = element_text(
      size = 13,
      hjust = 1,
    )
  )

p.none

### 2. Error = 0.01
# Create data frame for visualization
df_plot <- AIC_mvOU_BAYOU_SyyDiagonal_niter100_ME_best.AICc %>%
  filter(Error == "low") %>%
  group_by(Predictor) %>%
  mutate(delta_AICc = AICc - min(AICc)) %>%
  ungroup()

# rename the predictors
df_plot <- df_plot %>%
  mutate(Predictor = case_when(
    Predictor == "Metabolic_capacity" ~ "Metabolic capacity",
    Predictor == "Genome_size" ~ "Genome size",
    Predictor == "Metabolic_versatility" ~ "Metabolic versatility",
    Predictor == "Metabolic_completeness" ~ "Metabolic completeness",
    Predictor == "GC_percentage" ~ "GC percentage",
    Predictor == "Niche_breadth" ~ "Niche breadth",
    Predictor == "Oxygen" ~ "Oxygen tolerance",
    Predictor == "Temperature" ~ "Opt.Temperature",
    Predictor == "Salinity" ~ "Opt.Salinity",
    Predictor == "pH" ~ "Opt.pH",
    Predictor == "rRNA_operon_copy_number" ~ "rrn Copy number"
  ))

# reorder the trait
df_plot$Predictor <- factor(df_plot$Predictor, levels = 
                              c("Metabolic completeness","Metabolic versatility", "Metabolic capacity",
                                "rrn Copy number","GC percentage", "Genome size", "Niche breadth", "Oxygen tolerance",
                                "Opt.Temperature", "Opt.Salinity", "Opt.pH"))

# reorder the regime
df_plot$Rtype <- factor(df_plot$Rtype, levels = 
                          rev(c("regime_trait1","regime_trait2","regime_trait1_trait2","null")))


### Figure
p.low <- ggplot(df_plot, aes(x = Atype, y = Rtype, fill = delta_AICc)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(delta_AICc, 1)), size = 4.5, color = "white") +
  scale_fill_viridis_c(option = "D", direction = -1) +
  geom_point(
    data = df_plot %>% filter(delta_AICc == 0),
    shape = 21, size = 4.5, stroke = 1.2, color = "black", fill = "white"
  ) +
  facet_wrap2(~ Predictor, scales = "fixed", axes = "margins", ncol = 4) +
  scale_y_discrete(labels = c(
    "null" = "No regime",
    "regime_trait1" = "Growth rate",
    "regime_trait2" = "Trait",
    "regime_trait1_trait2" = expression("Growth rate" * " × " * "Trait")
  )) +
  scale_x_discrete(labels = c(
    "Any" = "Bidirectional",
    "Diagonal" = "No coevolution",
    "LowerTri" = "Growth → Trait",
    "UpperTri" = "Trait → Growth"
  )) +
  labs(
    title = "Measurement error = 0.01",
    x = "Mode of coevolution",
    y = "Regime definition",
    fill = expression(Delta*AIC[c])
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 18,       
      face = "bold",   
      hjust = 0.5      
    ),
    panel.grid = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 14   
    ),
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(
      angle = 30,
      hjust = 0.5,
      vjust = 0.7,
      size = 11
    ),
    axis.text.y = element_text(
      size = 13,
      hjust = 1,
    )
  )

p.low

### 3. Error = 0.05
# Create data frame for visualization
df_plot <- AIC_mvOU_BAYOU_SyyDiagonal_niter100_ME_best.AICc %>%
  filter(Error == "mid") %>%
  group_by(Predictor) %>%
  mutate(delta_AICc = AICc - min(AICc)) %>%
  ungroup()

# rename the predictors
df_plot <- df_plot %>%
  mutate(Predictor = case_when(
    Predictor == "Metabolic_capacity" ~ "Metabolic capacity",
    Predictor == "Genome_size" ~ "Genome size",
    Predictor == "Metabolic_versatility" ~ "Metabolic versatility",
    Predictor == "Metabolic_completeness" ~ "Metabolic completeness",
    Predictor == "GC_percentage" ~ "GC percentage",
    Predictor == "Niche_breadth" ~ "Niche breadth",
    Predictor == "Oxygen" ~ "Oxygen tolerance",
    Predictor == "Temperature" ~ "Opt.Temperature",
    Predictor == "Salinity" ~ "Opt.Salinity",
    Predictor == "pH" ~ "Opt.pH",
    Predictor == "rRNA_operon_copy_number" ~ "rrn Copy number"
  ))

# reorder the trait
df_plot$Predictor <- factor(df_plot$Predictor, levels = 
                              c("Metabolic completeness","Metabolic versatility", "Metabolic capacity",
                                "rrn Copy number","GC percentage", "Genome size", "Niche breadth", "Oxygen tolerance",
                                "Opt.Temperature", "Opt.Salinity", "Opt.pH"))

# reorder the regime
df_plot$Rtype <- factor(df_plot$Rtype, levels = 
                          rev(c("regime_trait1","regime_trait2","regime_trait1_trait2","null")))


### Figure
p.mid <- ggplot(df_plot, aes(x = Atype, y = Rtype, fill = delta_AICc)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(delta_AICc, 1)), size = 4.5, color = "white") +
  scale_fill_viridis_c(option = "D", direction = -1) +
  geom_point(
    data = df_plot %>% filter(delta_AICc == 0),
    shape = 21, size = 4.5, stroke = 1.2, color = "black", fill = "white"
  ) +
  facet_wrap2(~ Predictor, scales = "fixed", axes = "margins", ncol = 4) +
  scale_y_discrete(labels = c(
    "null" = "No regime",
    "regime_trait1" = "Growth rate",
    "regime_trait2" = "Trait",
    "regime_trait1_trait2" = expression("Growth rate" * " × " * "Trait")
  )) +
  scale_x_discrete(labels = c(
    "Any" = "Bidirectional",
    "Diagonal" = "No coevolution",
    "LowerTri" = "Growth → Trait",
    "UpperTri" = "Trait → Growth"
  )) +
  labs(
    title = "Measurement error = 0.05",
    x = "Mode of coevolution",
    y = "Regime definition",
    fill = expression(Delta*AIC[c])
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 18,       
      face = "bold",   
      hjust = 0.5      
    ),
    panel.grid = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 14   
    ),
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(
      angle = 30,
      hjust = 0.5,
      vjust = 0.7,
      size = 11
    ),
    axis.text.y = element_text(
      size = 13,
      hjust = 1,
    )
  )

p.mid


### 4. Error = 0.1 
# Create data frame for visualization
df_plot <- AIC_mvOU_BAYOU_SyyDiagonal_niter100_ME_best.AICc %>%
  filter(Error == "high") %>%
  group_by(Predictor) %>%
  mutate(delta_AICc = AICc - min(AICc)) %>%
  ungroup()

# rename the predictors
df_plot <- df_plot %>%
  mutate(Predictor = case_when(
    Predictor == "Metabolic_capacity" ~ "Metabolic capacity",
    Predictor == "Genome_size" ~ "Genome size",
    Predictor == "Metabolic_versatility" ~ "Metabolic versatility",
    Predictor == "Metabolic_completeness" ~ "Metabolic completeness",
    Predictor == "GC_percentage" ~ "GC percentage",
    Predictor == "Niche_breadth" ~ "Niche breadth",
    Predictor == "Oxygen" ~ "Oxygen tolerance",
    Predictor == "Temperature" ~ "Opt.Temperature",
    Predictor == "Salinity" ~ "Opt.Salinity",
    Predictor == "pH" ~ "Opt.pH",
    Predictor == "rRNA_operon_copy_number" ~ "rrn Copy number"
  ))

# reorder the trait
df_plot$Predictor <- factor(df_plot$Predictor, levels = 
                              c("Metabolic completeness","Metabolic versatility", "Metabolic capacity",
                                "rrn Copy number","GC percentage", "Genome size", "Niche breadth", "Oxygen tolerance",
                                "Opt.Temperature", "Opt.Salinity", "Opt.pH"))

# reorder the regime
df_plot$Rtype <- factor(df_plot$Rtype, levels = 
                          rev(c("regime_trait1","regime_trait2","regime_trait1_trait2","null")))


### Figure
p.high <- ggplot(df_plot, aes(x = Atype, y = Rtype, fill = delta_AICc)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(delta_AICc, 1)), size = 4.5, color = "white") +
  scale_fill_viridis_c(option = "D", direction = -1) +
  geom_point(
    data = df_plot %>% filter(delta_AICc == 0),
    shape = 21, size = 4.5, stroke = 1.2, color = "black", fill = "white"
  ) +
  facet_wrap2(~ Predictor, scales = "fixed", axes = "margins", ncol = 4) +
  scale_y_discrete(labels = c(
    "null" = "No regime",
    "regime_trait1" = "Growth rate",
    "regime_trait2" = "Trait",
    "regime_trait1_trait2" = expression("Growth rate" * " × " * "Trait")
  )) +
  scale_x_discrete(labels = c(
    "Any" = "Bidirectional",
    "Diagonal" = "No coevolution",
    "LowerTri" = "Growth → Trait",
    "UpperTri" = "Trait → Growth"
  )) +
  labs(
    title = "Measurement error = 0.1",
    x = "Mode of coevolution",
    y = "Regime definition",
    fill = expression(Delta*AIC[c])
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 18,       
      face = "bold",   
      hjust = 0.5      
    ),
    panel.grid = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 14   
    ),
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(
      angle = 30,
      hjust = 0.5,
      vjust = 0.7,
      size = 11
    ),
    axis.text.y = element_text(
      size = 13,
      hjust = 1,
    )
  )


##### Get the Figure S5 -----
p.error <- (p.none|p.low)/(p.mid|p.high)

#---------------------------------------------------------------------------------#
#                4 Model comparison; Syy = UpperTri; without M.error              #
#---------------------------------------------------------------------------------#
### 1. Load data
load(file = "mvOU_BAYOU_SyyUpperTri_niter100.RData")

### 2. Get the AIC, AICc, BIC etc for all iterations
AIC_mvOU_BAYOU_SyyUpperTri_niter100 <- data.frame()
for (var in predictors) {
  for (A in Atype) {
    for (R in Rtype) {
      for (i in 1:100) {
        aic.c <- mvOU_BAYOU_SyyUpperTri_niter100[[var]][[i]][["final"]][[A]][[R]][["model"]][["FinalFound"]][["ParamSummary"]][["aic.c"]]
        aic <- mvOU_BAYOU_SyyUpperTri_niter100[[var]][[i]][["final"]][[A]][[R]][["model"]][["FinalFound"]][["ParamSummary"]][["aic"]]
        bic <- mvOU_BAYOU_SyyUpperTri_niter100[[var]][[i]][["final"]][[A]][[R]][["model"]][["FinalFound"]][["ParamSummary"]][["bic"]]
        sic <- mvOU_BAYOU_SyyUpperTri_niter100[[var]][[i]][["final"]][[A]][[R]][["model"]][["FinalFound"]][["ParamSummary"]][["sic"]]
        LogLik <- mvOU_BAYOU_SyyUpperTri_niter100[[var]][[i]][["final"]][[A]][[R]][["model"]][["FinalFound"]][["ParamSummary"]][["LogLik"]]
        df <- data.frame(
          Predictor = var,
          Atype = A,
          Rtype = R,
          nIter = i,
          AICc = aic.c,
          AIC = aic,
          BIC = bic,
          SIC = sic,
          LogLik = LogLik
        )
        AIC_mvOU_BAYOU_SyyUpperTri_niter100 <- rbind(AIC_mvOU_BAYOU_SyyUpperTri_niter100, df)
      }
    }
  }
  
}
### 3. Find the solution with the lowest AICc across all iterations
AIC_mvOU_BAYOU_SyyUpperTri_niter100_best.AICc <- AIC_mvOU_BAYOU_SyyUpperTri_niter100 %>%
  group_by(Predictor, Atype, Rtype) %>%
  slice_min(AICc, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Syytype = "UpperTri")

#------------------------------------------------------------------------------#
#                 5 Model comparison: Syy = Diagonal vs. UpperTri              #
#------------------------------------------------------------------------------#
AIC_mvOU_BAYOU_niter100_best.AICc <- rbind(AIC_mvOU_BAYOU_SyyDiagonal_niter100_best.AICc, AIC_mvOU_BAYOU_SyyUpperTri_niter100_best.AICc)

p <- ggpaired(
  AIC_mvOU_BAYOU_niter100_best.AICc,
  x = "Syytype",
  y = "AICc",
  color = "Syytype",
  line.color = "gray70",
  line.size = 0.5,
  palette = "jco"
)

p2 <- p +
  stat_summary(fun = mean, geom = "point", size = 3, color = "black") +  # mean points
  stat_summary(fun = mean, geom = "line", aes(group = 1), 
               linetype = "dashed", color = "black") +                   # mean trend
  stat_compare_means(
    method = "wilcox.test",
    label.x.npc = "left",
    label.y.npc = "top"
  ) +
  labs(
    title = NULL,
    x = "Syy matrix",
    y = "AICc"
  ) +
  #theme_pubr(base_size = 14) +
  theme_bw(base_size = 14)+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.position = "none",
    axis.title.x = element_text(size = 14),  #  x-axis title size
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )

### Figure S8 ---
p2

#------------------------------------------------------------------------------------------------#
#               6 Model comparison across measurement error settings; Syy = Diagonal             #
#------------------------------------------------------------------------------------------------#

### 1. Combine model diagnosis outputs across 4 measurement error settings
AIC_mvOU_BAYOU_SyyDiagonal_niter100_best.AICc <- AIC_mvOU_BAYOU_SyyDiagonal_niter100_best.AICc %>%
  mutate(Error = "none")

AIC_mvOU_BAYOU_SyyDiagonal_niter100_best <- bind_rows(AIC_mvOU_BAYOU_SyyDiagonal_niter100_best.AICc,
                                                      AIC_mvOU_BAYOU_SyyDiagonal_niter100_ME_best.AICc)

### 2. Get the best model across measurement error settings 
AIC_mvOU_BAYOU_SyyDiagonal_niter100_best <- AIC_mvOU_BAYOU_SyyDiagonal_niter100_best %>%
  group_by(Predictor, Atype, Rtype) %>%
  slice_min(AICc, n = 1, with_ties = FALSE) %>%
  ungroup() 

### 3. Visualization
# Create data frame for visualization
df_plot <- AIC_mvOU_BAYOU_SyyDiagonal_niter100_best %>%
  group_by(Predictor) %>%
  mutate(delta_AICc = AICc - min(AICc)) %>%
  ungroup()

# rename the predictors
df_plot <- df_plot %>%
  mutate(Predictor = case_when(
    Predictor == "Metabolic_capacity" ~ "Metabolic capacity",
    Predictor == "Genome_size" ~ "Genome size",
    Predictor == "Metabolic_versatility" ~ "Metabolic versatility",
    Predictor == "Metabolic_completeness" ~ "Metabolic completeness",
    Predictor == "GC_percentage" ~ "GC percentage",
    Predictor == "Niche_breadth" ~ "Niche breadth",
    Predictor == "Oxygen" ~ "Oxygen tolerance",
    Predictor == "Temperature" ~ "Opt.Temperature",
    Predictor == "Salinity" ~ "Opt.Salinity",
    Predictor == "pH" ~ "Opt.pH",
    Predictor == "rRNA_operon_copy_number" ~ "rrn Copy number"
  ))

# reorder the predictor
df_plot$Predictor <- factor(df_plot$Predictor, levels = 
                              c("Metabolic completeness","Metabolic versatility", "Metabolic capacity",
                                "rrn Copy number","GC percentage", "Genome size", "Niche breadth", "Oxygen tolerance",
                                "Opt.Temperature", "Opt.Salinity", "Opt.pH"))

# reorder the regime
df_plot$Rtype <- factor(df_plot$Rtype, levels = 
                          rev(c("regime_trait1","regime_trait2","regime_trait1_trait2","null")))


### Figure
p.best <- ggplot(df_plot, aes(x = Atype, y = Rtype, fill = delta_AICc)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(delta_AICc, 1)), size = 4.5, color = "white") +
  scale_fill_viridis_c(option = "D", direction = -1) +
  geom_point(
    data = df_plot %>% filter(delta_AICc == 0),
    shape = 21, size = 4.5, stroke = 1.2, color = "black", fill = "white"
  ) +
  facet_wrap2(~ Predictor, scales = "fixed", axes = "margins", ncol = 4) +
  scale_y_discrete(labels = c(
    "null" = "No regime",
    "regime_trait1" = "Growth rate",
    "regime_trait2" = "Trait",
    "regime_trait1_trait2" = expression("Growth rate" * " × " * "Trait")
  )) +
  scale_x_discrete(labels = c(
    "Any" = "Bidirectional",
    "Diagonal" = "No coevolution",
    "LowerTri" = "Growth → Trait",
    "UpperTri" = "Trait → Growth"
  )) +
  labs(
    title = NULL,
    x = "Mode of coevolution",
    y = "Regime definition",
    fill = expression(Delta*AIC[c])
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 18,       
      face = "bold",   
      hjust = 0.5      
    ),
    panel.grid = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 14   # increase this for larger panel titles
    ),
    axis.title.x = element_text(size = 14),  #  x-axis title size
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(
      angle = 30,
      hjust = 0.5,
      vjust = 0.7,
      size = 11
    ),
    axis.text.y = element_text(
      size = 13,
      hjust = 1,
    )
  )

### Figure 4 ---
p.best

#------------------------------------------------------------------------------------------------#
#                     7 Model interpretation based on the best models; A matrix                  #
#------------------------------------------------------------------------------------------------#
### 1. Get the index of models with the smallest AICc 
AIC_mvOU_BAYOU_SyyDiagonal_niter100_best <- bind_rows(AIC_mvOU_BAYOU_SyyDiagonal_niter100_best.AICc,
                                                      AIC_mvOU_BAYOU_SyyDiagonal_niter100_ME_best.AICc)

AIC_best_mvOU_BAYOU_SyyDiagonal <- AIC_mvOU_BAYOU_SyyDiagonal_niter100_best %>%
  group_by(Predictor) %>%
  slice_min(AICc, n = 1, with_ties = FALSE) %>%
  ungroup() 

head(AIC_best_mvOU_BAYOU_SyyDiagonal)

### 2. Get the A matrix of the best model
Amatrix_mvOU_BAYOU_SyyDiagonal <- list()
for (i in seq_len(nrow(AIC_best_mvOU_BAYOU_SyyDiagonal))) {
  pre   <- as.character(AIC_best_mvOU_BAYOU_SyyDiagonal$Predictor[i])
  Atype <- as.character(AIC_best_mvOU_BAYOU_SyyDiagonal$Atype[i])
  Rtype <- as.character(AIC_best_mvOU_BAYOU_SyyDiagonal$Rtype[i])
  nIter <- as.numeric(AIC_best_mvOU_BAYOU_SyyDiagonal$nIter[i])
  Error <- as.character(AIC_best_mvOU_BAYOU_SyyDiagonal$Error[i])
  if (Error == "none") {
    A <- mvOU_BAYOU_SyyDiagonal_niter100[[pre]][[nIter]][["final"]][[Atype]][[Rtype]][["model"]][["FinalFound"]][["ParamsInModel"]][["A"]]
  } else {
    A <- mvOU_BAYOU_SyyDiagonal_niter100_ME[[pre]][[Error]][[nIter]][["final"]][[Atype]][[Rtype]][["model"]][["FinalFound"]][["ParamsInModel"]][["A"]]
  }
  Amatrix_mvOU_BAYOU_SyyDiagonal[[pre]] <- A
}

### 3. Specify the trait names
predictors <- c("Metabolic_completeness", "Metabolic_versatility", "Metabolic_capacity","rRNA_operon_copy_number",
                "GC_percentage", "Genome_size", "Niche_breadth", "Oxygen", "Temperature", "Salinity", "pH")


predictor2 <- c("Metabolic completeness","Metabolic versatility", "Metabolic capacity",
                "rrn Copy number","GC percentage", "Genome size", "Niche breadth", "Oxygen tolerance",
                "Opt.Temperature", "Opt.Salinity", "Opt.pH")

### 4. Function to generate the plot
plot_alpha <- function(A_list, predictors, predictors_pretty) {
  
  plot_list <- list()
  
  for (i in seq_along(predictors)) {
    
    pred <- predictors[i]
    pretty_name <- predictors_pretty[i]
    
    # Get matrix directly from list
    if (!pred %in% names(A_list)) next
    A <- A_list[[pred]]
    
    if (is.null(A)) next
    
    # Rename for plotting
    colnames(A) <- rownames(A) <- c("Growth rate", pretty_name)
    
    # Extract directional coefficients (column → row)
    val_GR_to_pred <- A[pretty_name, "Growth rate"]   # GR → Predictor
    val_pred_to_GR <- A["Growth rate", pretty_name]   # Predictor → GR
    
    # Colors
    col_GR_to_pred <- ifelse(val_GR_to_pred > 0, "red3",
                             ifelse(val_GR_to_pred < 0, "skyblue3", "grey60"))
    col_pred_to_GR <- ifelse(val_pred_to_GR > 0, "red3",
                             ifelse(val_pred_to_GR < 0, "skyblue3", "grey60"))
    
    # Layout
    x_GR <- 0
    x_pred <- 1
    y_label <- 0
    arrow_gap <- 0.02
    arrow_trim <- 0.3
    label_size <- 7
    value_size <- 6
    value_offset <- 0.04
    label_offset_GR <- 0.25
    label_offset_pred <- -0.25
    
    # Plot
    p <- ggplot() +
      # GR → Predictor
      geom_segment(
        aes(
          x = x_GR + arrow_trim, y = y_label + arrow_gap,
          xend = x_pred - arrow_trim, yend = y_label + arrow_gap
        ),
        arrow = arrow(length = unit(0.25, "cm"), ends = "last"),
        color = col_GR_to_pred, linewidth = 1
      ) +
      # Predictor → GR
      geom_segment(
        aes(
          x = x_pred - arrow_trim, y = y_label - arrow_gap,
          xend = x_GR + arrow_trim, yend = y_label - arrow_gap
        ),
        arrow = arrow(length = unit(0.25, "cm"), ends = "last"),
        color = col_pred_to_GR, linewidth = 1
      ) +
      # Labels
      annotate("text", x = x_GR + label_offset_GR, y = y_label,
               label = "Growth rate", size = label_size, hjust = 1) +
      annotate("text", x = x_pred + label_offset_pred, y = y_label,
               label = pretty_name, size = label_size, hjust = 0) +
      # Values
      annotate("text", x = 0.5, y = y_label + arrow_gap + value_offset,
               label = sprintf("%.2f", val_GR_to_pred),
               color = col_GR_to_pred, size = value_size) +
      annotate("text", x = 0.5, y = y_label - arrow_gap - value_offset,
               label = sprintf("%.2f", val_pred_to_GR),
               color = col_pred_to_GR, size = value_size) +
      coord_fixed(xlim = c(-0.1, 1.1), ylim = c(-0.15, 0.15), clip = "off") +
      theme_void() +
      theme(
        plot.margin = margin(t = 0, r = 25, b = 0, l = 25, unit = "mm")
      )
    
    plot_list[[pretty_name]] <- p
  }
  
  return(plot_list)
}

### 5. Plot and output the A matrix: Figure 5a
plots <- plot_alpha(
  Amatrix_mvOU_BAYOU_SyyDiagonal,
  predictors,
  predictor2
)

#---------------------------------------------------------------------------------------------------#
#                 8 Model interpretation based on the best models; Eigen-decomposition              #
#---------------------------------------------------------------------------------------------------#
### 1. Initialize the list of eigen-decompose
Eigen_mvOU_BAYOU_SyyDiagonal_best <- list()
for (i in seq_len(nrow(AIC_best_mvOU_BAYOU_SyyDiagonal))) {
  pre   <- as.character(AIC_best_mvOU_BAYOU_SyyDiagonal$Predictor[i])
  Atype <- as.character(AIC_best_mvOU_BAYOU_SyyDiagonal$Atype[i])
  Rtype <- as.character(AIC_best_mvOU_BAYOU_SyyDiagonal$Rtype[i])
  nIter <- as.numeric(AIC_best_mvOU_BAYOU_SyyDiagonal$nIter[i])
  Error <- as.character(AIC_best_mvOU_BAYOU_SyyDiagonal$Error[i])
  
  if (Error == "none") {
    half.t <- rbind(
      mvOU_BAYOU_SyyDiagonal_niter100[[pre]][[nIter]][["final"]][[Atype]][[Rtype]][["model"]][["FinalFound"]][["ParamSummary"]][["phyl.halflife"]][["directions"]],
      mvOU_BAYOU_SyyDiagonal_niter100[[pre]][[nIter]][["final"]][[Atype]][[Rtype]][["model"]][["FinalFound"]][["ParamSummary"]][["phyl.halflife"]][["halflives"]]
    )
  } else {
    half.t <- rbind(
      mvOU_BAYOU_SyyDiagonal_niter100_ME[[pre]][[Error]][[nIter]][["final"]][[Atype]][[Rtype]][["model"]][["FinalFound"]][["ParamSummary"]][["phyl.halflife"]][["directions"]],
      mvOU_BAYOU_SyyDiagonal_niter100_ME[[pre]][[Error]][[nIter]][["final"]][[Atype]][[Rtype]][["model"]][["FinalFound"]][["ParamSummary"]][["phyl.halflife"]][["halflives"]]
    )
  }
  
  Eigen_mvOU_BAYOU_SyyDiagonal_best[[pre]] <- half.t
}


### 3. Get the eigen-decompose list for visualization
eigen_matrices <- list()
for (i in 1:11) {
  pre <- predictors[i]
  A <- Amatrix_mvOU_BAYOU_SyyDiagonal[[pre]]
  eig <- eigen(A)
  
  lambda <- Re(eig$values)
  vectors <- Re(eig$vectors)
  
  # Half-life and scaling
  half_life <- log(2) / lambda
  half_life_percent <- half_life * 100
  
  predictor_name <- predictor2[i]
  
  # Build final matrix explicitly by columns
  final_matrix <- data.frame(
    `Growth rate` = vectors[1, ],
    `Factor` = vectors[2, ],
    `Eigenvalue (λ)` = lambda,
    `Half-life (t[1/2] = ln2/λ)` = half_life,
    `Tree Height (%)` = half_life_percent
  )
  
  # Rename columns for clarity
  colnames(final_matrix) <- c("Growth rate", "Factor", "Eigenvalue (λ)", 
                              "Half-life (t[1/2] = ln2/λ)", "Tree Height (%)")
  rownames(final_matrix) <- c(paste0("\u20D7e", 1:2, "_", predictor_name))
  
  eigen_matrices[[predictor_name]] <- t(final_matrix)
}

### 3. Combine all matrices horizontally to get the parameters of eigen-decomposition: Figure S10
combined_eigen_matrices <- do.call(cbind, eigen_matrices)

### 4. Visualization
# --- Extract info function ---
extract_info <- function(mat, pair_name) {
  m <- as.data.frame(mat, stringsAsFactors = FALSE)
  rn <- rownames(m)
  factor_row <- setdiff(rn, c("Growth rate", "Eigenvalue (λ)",
                              "Half-life (t[1/2] = ln2/λ)", "Tree Height (%)"))
  if (length(factor_row) > 1) {
    factor_row <- factor_row[grepl(pair_name, factor_row)]
    if (length(factor_row) == 0) factor_row <- setdiff(rn, c("Growth rate",
                                                             "Eigenvalue (λ)",
                                                             "Half-life (t[1/2] = ln2/λ)",
                                                             "Tree Height (%)"))[1]
  }
  
  data.frame(
    Pair = pair_name,
    eig = c("e1", "e2"),
    load_growth = as.numeric(m["Growth rate", 1:2]),
    load_factor = as.numeric(m[factor_row, 1:2]),
    lambda = as.numeric(m["Eigenvalue (λ)", 1:2]),
    half_life = as.numeric(m["Half-life (t[1/2] = ln2/λ)", 1:2]),
    stringsAsFactors = FALSE
  )
}

# --- Combine all matrices ---
all_info <- do.call(rbind,
                    lapply(names(eigen_matrices),
                           function(nm) extract_info(eigen_matrices[[nm]], nm))
)

# Compute endpoints (scaled by eigenvalue)
all_info <- all_info %>%
  mutate(xend = load_growth * lambda,
         yend = load_factor * lambda)

# --- Plot function (relaxed scaling) ---
plot_pair <- function(pair_name, df_pair) {
  hl_text <- paste0(
    "t½(e1) = ", round(df_pair$half_life[df_pair$eig == "e1"] * 100, 1), "%,  ",
    "t½(e2) = ", round(df_pair$half_life[df_pair$eig == "e2"] * 100, 1), "%"
  )
  
  
  ggplot(df_pair) +
    geom_segment(aes(x = 0, y = 0, xend = xend, yend = yend, color = eig),
                 arrow = arrow(length = unit(0.28, "cm")), size = 1.1) +
    geom_text(aes(x = xend * 1.08, y = yend * 1.08, label = eig),
              size = 5, fontface = "bold", color = "black") +
    scale_color_manual(values = c("e1" = "#1f78b4", "e2" = "#e31a1c"), guide = "none") +
    labs(
      title = pair_name,
      subtitle = hl_text,
      x = NULL,
      y = NULL
    ) +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14, color = "gray40"),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 13)
    )
}

# --- Generate all plots ---
plot_list <- lapply(predictor2, function(pn) {
  plot_pair(pn, filter(all_info, Pair == pn))
})

# Combine all panels (auto-scaled, no fixed coordinate limits)
combined <- wrap_plots(plot_list, ncol = 3)

### Figure 5b ---
combined


