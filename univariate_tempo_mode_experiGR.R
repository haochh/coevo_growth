library(tidyverse)
library(ape)
library(geiger)
library(gt)

#----------------------------------------------------------#
#                       0 Preparation                      #
#----------------------------------------------------------#

### 1. Define trait names ---
predictors <- c("Metabolic_capacity", "Metabolic_versatility", 
                "Metabolic_completeness", "GC_percentage", "Genome_size", "rRNA_operon_copy_number", 
                "Niche_breadth", "Oxygen", "Temperature", "Salinity", "pH")

predictors2 <- c("Metabolic capacity", "Metabolic versatility", 
                 "Metabolic completeness", "GC percentage", "Genome size", "rrn Copy number", 
                 "Niche breadth", "Oxygen tolerance", "Opt.Temperature", "Opt.Salinity", "Opt.pH")

variables <- c("Growth_rate_experiment", predictors)
variables2 <- c("Growth_rate_XU_withOGT", predictors)

### 2. Load data and tree ---
# Scaled data
load("GR_Madin_imputated_RRN_11_factor_scaled.RData")
# Tree (non-ultrametric)
tree_madin <- read.tree("tree_madin.tree")

### 3. Scale the tree to ensure height = 1  ---
tree_height <-mvSLOUCH::phyltree_paths(tree_madin)$tree_height
tree_height # = 2.21448
tree_madin_scaled <- tree_madin
tree_madin_scaled$edge.length <- tree_madin_scaled$edge.length/tree_height
# Check the height again, and this should = 1
mvSLOUCH::phyltree_paths(tree_madin_scaled)$tree_height # = 1

### 4.Prepare the data  ---
data_bivariate_mvOU_experiGR <- GR_Madin_imputated_RRN_11_factor_scaled %>%
  column_to_rownames(var = "gtdb_accession") %>%
  select(all_of(variables))

# Re-order the rows to ensure the order is the same with the tip labels in the tree
data_bivariate_mvOU_experiGR <- data_bivariate_mvOU_experiGR[tree_madin_scaled$tip.label, ]
# Check column orders and row orders are correct
# Column orders: should return all TRUE
# Row orders: confirm if this is done, should return TRUE
isTRUE(all.equal(tree_madin_scaled$tip.label,rownames(data_bivariate_mvOU_experiGR)))
# Row orders: should return all TRUE
table(row.names(data_bivariate_mvOU_experiGR) == tree_madin_scaled$tip.label)

### 5. Define model type ---

model <- c("BM","OU","EB")

#----------------------------------------------------------#
#                       1 Fit Continuous                   #
#----------------------------------------------------------#
Fit_Continuous_ExperiGR_11factors <- list()
for (var in variables) {
  
  trait <- setNames(
    data_bivariate_mvOU_experiGR[[var]],
    rownames(data_bivariate_mvOU_experiGR)
  )
  print(var)
  for (mod in model) {
    print(mod)
    fit <- fitContinuous(
      phy = tree_madin_scaled,
      dat = trait,
      model = mod,
      niter = 100
    )
    
    Fit_Continuous_ExperiGR_11factors[[var]][[mod]] <- fit
  }
}


#--------------------------------------------------------------#
#           2 Compare the goodness of fit using AICc           #
#--------------------------------------------------------------#

AICc_Fit_Continuous_ExperiGR_11factors <- data.frame()
for (var in variables) {
  for (mod in model) {
    aic.c <- Fit_Continuous_ExperiGR_11factors[[var]][[mod]][["opt"]][["aicc"]]
    df <- data.frame(
      Predictor = var,
      Model = mod,
      AICc = aic.c
    )
    
    AICc_Fit_Continuous_ExperiGR_11factors <- rbind(AICc_Fit_Continuous_ExperiGR_11factors, df)
  }
}


AICc_Fit_Continuous_ExperiGR_11factors <- AICc_Fit_Continuous_ExperiGR_11factors %>%
  mutate(Predictor = case_when(
    Predictor == "Growth_rate_experiment" ~ "Growth rate",
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


AICc_Fit_Continuous_ExperiGR_11factors$Predictor <- factor(AICc_Fit_Continuous_ExperiGR_11factors$Predictor, levels = 
                                                             c("Growth rate", "Metabolic completeness","Metabolic versatility", "Metabolic capacity",
                                                               "rrn Copy number","GC percentage", "Genome size", "Niche breadth", "Oxygen tolerance",
                                                               "Opt.Temperature", "Opt.Salinity", "Opt.pH"))


#------------------------------------------------------------------------------------------#
#           3 Visualization (Figure S1): Model comparisons for univariate models           #
#------------------------------------------------------------------------------------------#

df_plot <- AICc_Fit_Continuous_ExperiGR_11factors

# Step 1: Create a new column to identify the best model for each predictor
df_plot <- df_plot %>%
  group_by(Predictor) %>%
  mutate(best_model = Model[which.min(AICc)]) %>%
  ungroup()

# Step 2: Plot with color highlighting the best model
p<-ggplot(df_plot, aes(x = Model, y = AICc, color = Model == best_model)) +
  geom_point(position = position_dodge(width = 0.6), size = 4) +  
  theme_minimal() +  
  labs(x = "Model", y = "AICc", title = NULL) +
  scale_color_manual(values = c("TRUE" = "#1f77b4", "FALSE" = "#7f7f7f")) +  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.position = "none",
    axis.title.x = element_text(size = 14),  #  x-axis title size
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold")
  )+
  facet_wrap(~Predictor, scales = "free_y") 

### Figure S1 ---
p


#-----------------------------------------------------------------------#
#            4 Model parameter extraction + Phylogenetic signal         #
#-----------------------------------------------------------------------#
### 4.1 Phylogenetic signal ---
PhySignal_ExperiGR_11factors <- list()
for (var in variables) {
  trait <- setNames(
    data_bivariate_mvOU_experiGR[[var]],
    rownames(data_bivariate_mvOU_experiGR)
  )
  
  sig.k <- phylosig(tree_madin_scaled, trait, method="K", test=T)
  sig.lambda <- phylosig(tree_madin_scaled, trait, method="lambda", test=T)
  
  PhySignal_ExperiGR_11factors[[var]][["K"]] <- sig.k
  PhySignal_ExperiGR_11factors[[var]][["lambda"]] <- sig.lambda
}

# Extract parameters for phylogenetic signal analysis
PhySignal_ExperiGR_11factors_result <- data.frame()
for (var in variables) {
  K <- PhySignal_ExperiGR_11factors[[var]][["K"]][["K"]]
  K.p <- PhySignal_ExperiGR_11factors[[var]][["K"]][["P"]]
  lambda <- PhySignal_ExperiGR_11factors[[var]][["lambda"]][["lambda"]]
  lambda.p <- PhySignal_ExperiGR_11factors[[var]][["lambda"]][["P"]]
  
  df <- data.frame(
    Predictor = var,
    K = K,
    K_p = K.p,
    lambda = lambda,
    lambda_p = lambda.p
  )
  
  colnames(df) <- c("Predictor", "Blomberg's K", "p-value for K", "Pagel's λ", "p-value for λ")
  
  PhySignal_ExperiGR_11factors_result <- rbind(PhySignal_ExperiGR_11factors_result, df)
}


PhySignal_ExperiGR_11factors_result <- PhySignal_ExperiGR_11factors_result %>%
  mutate(Predictor = case_when(
    Predictor == "Growth_rate_experiment" ~ "Growth rate",
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


### 4.2 Model parameter extraction for the best models ---
# Get the best model
AICc_Fit_Continuous_ExperiGR_11factors_best <- AICc_Fit_Continuous_ExperiGR_11factors %>%
  group_by(Predictor) %>%
  slice_min(AICc, n = 1, with_ties = FALSE) %>%
  ungroup()

# Extract the parameters
Param_list <- list()
for (i in seq_len(nrow(AICc_Fit_Continuous_ExperiGR_11factors_best))) {
  
  var <- AICc_Fit_Continuous_ExperiGR_11factors_best$Predictor[i]
  mod <- AICc_Fit_Continuous_ExperiGR_11factors_best$Model[i]
  
  fit <- Fit_Continuous_ExperiGR_11factors[[var]][[mod]]
  opt <- fit$opt
  
  row <- list(
    Predictor = var,
    Model = mod,
    sigma.sq_BM = NA_real_,
    sigma.sq_OU = NA_real_,
    alpha_OU = NA_real_,
    sigma.sq_EB = NA_real_,
    a_EB = NA_real_
  )
  
  if (mod == "BM") {
    row$sigma.sq_BM <- opt$sigsq
  }
  
  if (mod == "OU") {
    row$sigma.sq_OU <- opt$sigsq
    row$alpha_OU    <- opt$alpha
  }
  
  if (mod == "EB") {
    row$sigma.sq_EB <- opt$sigsq
    
    # robust to naming differences
    if ("a" %in% names(opt)) {
      row$a_EB <- opt$a
    } else if ("r" %in% names(opt)) {
      row$a_EB <- opt$r
    }
  }
  
  Param_list[[i]] <- row
}

Param_Fit_Continuous <- bind_rows(Param_list)


# Compute half-life for OU models and rename columns
tree_depth <- max(node.depth.edgelength(tree_madin_scaled))

Param_Fit_Continuous <- Param_Fit_Continuous %>%
  # Compute half-lives for OU models
  mutate(
    `t(1/2)_OU` = ifelse(Model == "OU" & !is.na(alpha_OU), log(2) / alpha_OU, NA),
    `Tree height %` = ifelse(!is.na(`t(1/2)_OU`),
                             (`t(1/2)_OU` / tree_depth) * 100,
                             NA_real_),
    `Best model` = Model,
    # Rename columns for more professional notation
    `σ²_BM` = sigma.sq_BM,
    `σ²_OU` = sigma.sq_OU,
    `α_OU` = alpha_OU,
    `b_EB` = a_EB,
    `σ²_EB` = sigma.sq_EB
  ) %>%
  # Drop old columns (those with original names)
  select(Predictor, `Best model`, `σ²_BM`, `σ²_OU`, `α_OU`, `b_EB`, `σ²_EB`, `t(1/2)_OU`,`Tree height %`)

### Combine the results for univariate analyses
PhySignal_ExperiGR_11factors_result
# Check the results
Param_Fit_Continuous

Univariate_FitContinuous_final_output <- left_join(PhySignal_ExperiGR_11factors_result, Param_Fit_Continuous,
                                                   by = "Predictor")


#-------------------------------------------------------------------------------------#
#           5 Visualization (Table 1): Model output for univariate analyses           #
#-------------------------------------------------------------------------------------#

Univariate_FitContinuous_final_output$Predictor <- factor(Univariate_FitContinuous_final_output$Predictor, levels = 
                                                            c("Growth rate", "Metabolic completeness","Metabolic versatility", "Metabolic capacity",
                                                              "rrn Copy number","GC percentage", "Genome size", "Niche breadth", "Oxygen tolerance",
                                                              "Opt.Temperature", "Opt.Salinity", "Opt.pH"))

### rename the columns
table_pub <- Univariate_FitContinuous_final_output %>%
  rename(
    Trait = Predictor,
    `Blomberg's K` = `Blomberg's K`,
    `P (K)` = `p-value for K`,
    `Pagel's λ` = `Pagel's λ`,
    `P (λ)` = `p-value for λ`,
    `Best model` = `Best model`,
    `σ² (BM)` = `σ²_BM`,
    `σ² (OU)` = `σ²_OU`,
    `α (OU)` = `α_OU`,
    `b (EB)` = `b_EB`,
    `σ² (EB)` = `σ²_EB`,
    `t½ (OU)` = `t(1/2)_OU`,
    `t½ (% tree height)` = `Tree height %`
  )

### create the table
gt_table <- table_pub %>%
  gt() %>%
  
  # Format numbers
  fmt_number(
    columns = where(is.numeric),
    decimals = 3
  ) %>%
  
  # Scientific notation for p-values
  fmt_scientific(
    columns = c(`P (K)`, `P (λ)`),
    decimals = 2
  ) %>%
  
  # Replace NA with blank
  fmt_missing(
    everything(),
    missing_text = ""
  ) %>%
  
  # Align
  cols_align(
    align = "center",
    everything()
  ) %>%
  
  # Bold model column
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = `Best model`)
  ) %>%
  
  # Improve spacing
  tab_options(
    table.font.size = "small",
    data_row.padding = px(4)
  )

gt_table <- gt_table %>%
  tab_spanner(
    label = "Phylogenetic Signal",
    columns = c(`Blomberg's K`, `P (K)`, `Pagel's λ`, `P (λ)`)
  ) %>%
  tab_spanner(
    label = "Model Parameters",
    columns = c(`σ² (BM)`, `σ² (OU)`, `α (OU)`, `b (EB)`, `σ² (EB)`)
  ) %>%
  tab_spanner(
    label = "OU Dynamics",
    columns = c(`t½ (OU)`, `t½ (% tree height)`)
  )

### Table 1 ---
gt_table
