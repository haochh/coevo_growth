library(tidyverse)
library(bayou)
library(geiger)
library(patchwork)
library(mvSLOUCH)
library(PCMBaseCpp)


#----------------------------------------------------------#
#                       0 Preparation                      #
#----------------------------------------------------------#

### 1. Define predictor names ---
predictors <- c("Metabolic_capacity", "Metabolic_versatility", 
                "Metabolic_completeness", "GC_percentage", "Genome_size", "rRNA_operon_copy_number", 
                "Niche_breadth", "Oxygen", "Temperature", "Salinity", "pH")

predictors2 <- c("Metabolic capacity", "Metabolic versatility", 
                 "Metabolic completeness", "GC percentage", "Genome size", "rrn Copy number", 
                 "Niche breadth", "Oxygen tolerance", "Opt.Temperature", "Opt.Salinity", "Opt.pH")

variables <- c("Growth_rate_experiment", predictors)


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


#----------------------------------------------------------------#
#                         1 Fit OUM models                       #
#----------------------------------------------------------------#
Bayou_OUM_ExperiGR_11factors <- list()
for (var in variables) {
  
  trait <- setNames(
    data_bivariate_mvOU_experiGR[[var]],
    rownames(data_bivariate_mvOU_experiGR)
  )

  ## make the OUM prior distribution
  prior_OUM <- make.prior(tree_madin_scaled,
                        dists=list(dalpha="dhalfcauchy",
                                   dsig2="dhalfcauchy",
                                   dk="cdpois",
                                   dtheta="dnorm"),
                        param=list(dalpha=list(scale=0.1),
                                   dsig2=list(scale=0.1),
                                   dk=list(lambda=1, kmax=50),
                                   dsb=list(bmax=1, prob=1),
                                   dtheta=list(mean=mean(trait),
                                               sd=1.5*sd(trait))),
                        plot.prior=F)
  
  ## randomly select starting values for MCMC
  startpars <- priorSim(prior_OUM,tree_madin_scaled,plot=FALSE)$pars[[1]]
  
  ## set up the MCMC run
  mcmcOU <- bayou.makeMCMC(tree_madin_scaled,trait,
                           prior=prior_OUM,plot.freq=NULL,file.dir = tempdir(),
                           ticker.freq=1000000)
  
  ## Run the MCMC
  mcmcOU$run(100000)
  # Get the Chain 
  chainOUM <- mcmcOU$load()
  # Get the stepping stone for marginal likelihood
  ss_OUM <- mcmcOU$steppingstone(10000, chain = chainOUM, Bk = Bk,
                                 burnin = 0.3, plot = FALSE)
  # Get the results
  result <- list(
    MCMC = chainOUM,
    SS = ss_OUM
  )
  
  ## Store the results
  Bayou_OUM_ExperiGR_11factors[[var]] <- result
}


#---------------------------------------------------------------#
#                  2 Result interpretation (OUM)                #
#---------------------------------------------------------------#
### Load pre-computed Bayou results for OUM models
load("/Users/user/Chunhui/Life_history/Analysis_final/Univariate/bayou/Result/Bayou_OUM_ExperiGR_11factors.RData")

# Get the number of regimes for all traits: OUM
bayou_OUM_regime <- data.frame()
for (var in variables) {
  trait.rjMCMC <- Bayou_OUM_ExperiGR_11factors[[var]][["MCMC"]]
  trait.rjMCMC <- set.burnin(trait.rjMCMC,0.3)
  sum <- summary(trait.rjMCMC)
  k <- sum[["statistics"]][6,] %>% t() %>% as.data.frame()
  df <- k %>% mutate(Trait = var)
  bayou_OUM_regime <- rbind(bayou_OUM_regime, df)
}

#----------------------------------------------------------#
#                       3 Fit OU1 models                   #
#----------------------------------------------------------#
# Get the OU1 model 
Bayou_OU1_ExperiGR_11factors <- list()
for (var in variables) {
  ## Create a trait vector
  trait <- setNames(
    data_bivariate_mvOU_experiGR[[var]],
    rownames(data_bivariate_mvOU_experiGR)
  )
  
  ## make the OU1 prior distribution
  prior_OU1 <- make.prior(
    tree_madin_scaled,
    dists = list(
      dalpha ="dhalfcauchy",
      dsig2 ="dhalfcauchy",
      dk ="fixed",
      dsb ="fixed", 
      dtheta = "dnorm"
    ),
    param = list(
      dalpha=list(scale=0.1),
      dsig2=list(scale=0.1),
      dk="fixed",
      dsb ="fixed", 
      dtheta=list(mean=mean(trait), sd=1.5*sd(trait))
    ),
    fixed = list(k = 0, ntheta = 1, sb = 0),  
    model="OU",
    plot.prior=F
  )
  
  ## randomly select starting values for MCMC
  startpars <- priorSim(prior_OU1,tree_madin_scaled,plot=FALSE)$pars[[1]]
  
  ## set up the MCMC run
  mcmcOU <- bayou.makeMCMC(tree_madin_scaled,trait,
                           prior=prior_OU1,plot.freq=NULL,file.dir = tempdir(),
                           ticker.freq=1000000)
  
  
  ## Run the MCMC
  mcmcOU$run(100000)
  # Get the Chain 
  chainOU1 <- mcmcOU$load()
  # Get the stepping stone for marginal likelihood
  ss_OU1 <- mcmcOU$steppingstone(10000, chain = chainOU1, Bk = Bk,
                                 burnin = 0.3, plot = FALSE)
  # Get the results
  result <- list(
    MCMC = chainOU1,
    SS = ss_OU1
  )

  ## Store the results
  Bayou_OU1_ExperiGR_11factors[[var]] <- result
}


#----------------------------------------------------------#
#                        4 OU1 vs. OUM                     #
#----------------------------------------------------------#

load("Bayou_OU1_ExperiGR_11factors.RData")
load("Bayou_OUM_ExperiGR_11factors.RData")

Bayes_factor_OU1_OUM <- data.frame()
for (var in variables) {
  # Log marginal likelihoods
  lnr_OU1 <- Bayou_OU1_ExperiGR_11factors[[var]][["SS"]][["lnr"]]
  lnr_OUM <- Bayou_OUM_ExperiGR_11factors[[var]][["SS"]][["lnr"]]
  
  results <- data.frame(
    Trait = var,
    lnr_OU1 = lnr_OU1,
    lnr_OUM = lnr_OUM
  )
  # Log Bayes factor
  results$logBF_OUM_vs_OU1 <- results$lnr_OUM - results$lnr_OU1
  # Bayes factor
  results$BF_OUM_vs_OU1 <- exp(results$logBF_OUM_vs_OU1)
  # Add interpretation
  results$Support <- ifelse(
    results$logBF_OUM_vs_OU1 > 0,
    "OUM favoured",
    "OU1 favoured"
  )
  
  Bayes_factor_OU1_OUM <- rbind(Bayes_factor_OU1_OUM, results)
}


head(bayou_OUM_regime)

Bayou_results <- left_join(Bayes_factor_OU1_OUM, bayou_OUM_regime)

Bayou_results <- Bayou_results %>%
  mutate(Trait = case_when(
    Trait == "Growth_rate_experiment" ~ "Growth rate",
    Trait == "Metabolic_capacity" ~ "Metabolic capacity",
    Trait == "Genome_size" ~ "Genome size",
    Trait == "Metabolic_versatility" ~ "Metabolic versatility",
    Trait == "Metabolic_completeness" ~ "Metabolic completeness",
    Trait == "GC_percentage" ~ "GC percentage",
    Trait == "Niche_breadth" ~ "Niche breadth",
    Trait == "Oxygen" ~ "Oxygen tolerance",
    Trait == "Temperature" ~ "Opt.Temperature",
    Trait == "Salinity" ~ "Opt.Salinity",
    Trait == "pH" ~ "Opt.pH",
    Trait == "rRNA_operon_copy_number" ~ "rrn Copy number"
  ))


Bayou_results$Trait <- factor(Bayou_results$Trait, levels = 
                                rev(c("Growth rate", "Metabolic completeness","Metabolic versatility", "Metabolic capacity",
                                  "rrn Copy number","GC percentage", "Genome size", "Niche breadth", "Oxygen tolerance",
                                  "Opt.Temperature", "Opt.Salinity", "Opt.pH")))

### Result of bayou analysis ---
Bayou_results


#--------------------------------------------------------------------------#
#               5 Visualization: bayou results (Figure 3a)                 #
#--------------------------------------------------------------------------#

# ---- Panel 1: logBF bars ----
p1 <- ggplot(Bayou_results, aes(y = Trait, x = logBF_OUM_vs_OU1, fill = Support)) +
  geom_col(width = 0.45, alpha = 0.8) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c("OU1 favoured" = "skyblue3", "OUM favoured" = "red3")) +
  labs(x = NULL, y = NULL)+
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(size = 17),
    axis.title.x = element_text(size = 19),
    axis.title.y = element_text(size = 19),
    axis.text.y = element_text(size = 17),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13.5)
  )

# ---- Panel 2: forest plot  ----
p2 <- ggplot(Bayou_results, aes(y = Trait)) +
  geom_errorbarh(aes(xmin = HPD95Lower, xmax = HPD95Upper),
                 height = 0.25, na.rm = TRUE) +
  geom_point(aes(x = Mean, color = Support), size = 6, na.rm = TRUE) +
  scale_colour_manual(values = c("OU1 favoured" = "skyblue3", "OUM favoured" = "red3")) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 17),
    axis.title.x = element_text(size = 19),
    axis.title.y = element_text(size = 19),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13.5)
  )


# ---- Combine (aligned rows) ----
(p1 | p2) + plot_layout(widths = c(1.2, 1))


#------------------------------------------------------------------------------#
#               6 Visualize the location of regimes (Figure S2)                #
#------------------------------------------------------------------------------#
load("Bayou_OUM_ExperiGR_11factors.RData")

### Growth rate regime
trait.rjMCMC <- Bayou_OUM_ExperiGR_11factors[[1]][["MCMC"]]
trait.rjMCMC <- set.burnin(trait.rjMCMC,0.3)
sum <- summary(trait.rjMCMC)
sum$statistics

# Visualization

my.palette <- function(n) {
  c("#97a6c4","#f1a226")[1:n]
}

plotSimmap.mcmc(
  trait.rjMCMC,
  edge.type = "regimes",
  lwd = 2,
  pp.cutoff = 0.1,
  pal = my.palette,       
  show.tip.label = FALSE, 
  circles = TRUE,         # draw circles at nodes
  pp.labels = TRUE,       # show posterior probability numbers
  circle.cex.max = 5,     # max circle size
  pp.cex = 0.7
)

#---------------------------------------------------------------------------------------------------#
#               7 Visualize the location of regimes for associated traits (Figure S3)               #
#---------------------------------------------------------------------------------------------------#
### Metabolic_versatility
trait.rjMCMC <- Bayou_OUM_ExperiGR_11factors[["Metabolic_versatility"]][["MCMC"]]
trait.rjMCMC <- set.burnin(trait.rjMCMC,0.3)
sum <- summary(trait.rjMCMC)
sum$statistics

pp <- sum[["branch.posteriors"]]
# Get the theta of the regime 2
pp[which.max(pp$pp), ]

# Visualization
my.palette <- function(n) {
  c("#97a6c4","#f1a226")[1:n]
}


plotSimmap.mcmc(
  trait.rjMCMC,
  edge.type = "regimes",
  lwd = 2,
  pp.cutoff = 0.06,
  pal = my.palette,       
  show.tip.label = FALSE, 
  circles = F,         # draw circles at nodes
  pp.labels = F,       # show posterior probability numbers
  circle.cex.max = 5,     # max circle size
  type = "fan",
  pp.cex = 0.7
)

### Metabolic_capacity
trait.rjMCMC <- Bayou_OUM_ExperiGR_11factors[["Metabolic_capacity"]][["MCMC"]]
trait.rjMCMC <- set.burnin(trait.rjMCMC,0.3)
sum <- summary(trait.rjMCMC)
sum$statistics

pp <- sum[["branch.posteriors"]]
# Get the theta of the regime 2
pp[which.max(pp$pp), ]

# Visualization
my.palette <- function(n) {
  c("#97a6c4","#f1a226")[1:n]
}

plotSimmap.mcmc(
  trait.rjMCMC,
  edge.type = "regimes",
  lwd = 2,
  pp.cutoff = 0.06,
  pal = my.palette,       
  show.tip.label = FALSE, 
  circles = F,         # draw circles at nodes
  pp.labels = F,       # show posterior probability numbers
  circle.cex.max = 5,     # max circle size
  type = "fan",
  pp.cex = 0.7
)


### Genome_size
trait.rjMCMC <- Bayou_OUM_ExperiGR_11factors[["Genome_size"]][["MCMC"]]
trait.rjMCMC <- set.burnin(trait.rjMCMC,0.3)
sum <- summary(trait.rjMCMC)
sum$statistics

pp <- sum[["branch.posteriors"]]
# Get the theta of the regime 2
pp[which.max(pp$pp), ]

# Visualization
my.palette <- function(n) {
  c("#97a6c4","#f1a226")[1:n]
}


plotSimmap.mcmc(
  trait.rjMCMC,
  edge.type = "regimes",
  lwd = 2,
  pp.cutoff = 0.06,
  pal = my.palette,       
  show.tip.label = FALSE, 
  circles = F,         # draw circles at nodes
  pp.labels = F,       # show posterior probability numbers
  circle.cex.max = 5,     # max circle size
  type = "fan",
  pp.cex = 0.7
  
)


### Niche_breadth
trait.rjMCMC <- Bayou_OUM_ExperiGR_11factors[["Niche_breadth"]][["MCMC"]]
trait.rjMCMC <- set.burnin(trait.rjMCMC,0.3)
sum <- summary(trait.rjMCMC)
sum$statistics

pp <- sum[["branch.posteriors"]]
# Get the theta of the regime 2
pp[which.max(pp$pp), ]

# Visualization
my.palette <- function(n) {
  c("#97a6c4","#f1a226")[1:n]
}


plotSimmap.mcmc(
  trait.rjMCMC,
  edge.type = "regimes",
  lwd = 2,
  pp.cutoff = 0.1,
  pal = my.palette,       
  show.tip.label = FALSE, 
  circles = F,         # draw circles at nodes
  pp.labels = F,       # show posterior probability numbers
  circle.cex.max = 5,     # max circle size
  type = "fan",
  pp.cex = 0.7
  
)

### Oxygen
trait.rjMCMC <- Bayou_OUM_ExperiGR_11factors[["Oxygen"]][["MCMC"]]
trait.rjMCMC <- set.burnin(trait.rjMCMC,0.3)
sum <- summary(trait.rjMCMC)
sum$statistics

pp <- sum[["branch.posteriors"]]
# Get the theta of the regime 2
pp[which.max(pp$pp), ]

# Visualization
my.palette <- function(n) {
  c("#97a6c4","#f1a226")[1:n]
}

plotSimmap.mcmc(
  trait.rjMCMC,
  edge.type = "regimes",
  lwd = 2,
  pp.cutoff = 0.06,
  pal = my.palette,       
  show.tip.label = FALSE, 
  circles = F,         # draw circles at nodes
  pp.labels = F,       # show posterior probability numbers
  circle.cex.max = 5,     # max circle size
  type = "fan",
  pp.cex = 0.7
  
)


### Salinity
trait.rjMCMC <- Bayou_OUM_ExperiGR_11factors[["Salinity"]][["MCMC"]]
trait.rjMCMC <- set.burnin(trait.rjMCMC,0.3)
sum <- summary(trait.rjMCMC)
sum$statistics

pp <- sum[["branch.posteriors"]]
# Get the theta of the regime 2
pp[which.max(pp$pp), ]

# Visualization
my.palette <- function(n) {
  c("#97a6c4","#f1a226")[1:n]
}

plotSimmap.mcmc(
  trait.rjMCMC,
  edge.type = "regimes",
  lwd = 2,
  pp.cutoff = 0.06,
  pal = my.palette,       
  show.tip.label = FALSE, 
  circles = F,         # draw circles at nodes
  pp.labels = F,       # show posterior probability numbers
  circle.cex.max = 5,     # max circle size
  type = "fan",
  pp.cex = 0.7
  
)

#------------------------------------------------------------------------------#
#             8 Extract the location of regimes for mvOU analysis              #
#------------------------------------------------------------------------------#
### The index is post-order
tree_madin_scaled_postorder <- reorder(tree_madin_scaled, "postorder")

#---------------------------------------------
# Find all descendant edges of one focal edge
#---------------------------------------------
get_descendant_edges <- function(tree, focal_edge) {
  
  # the child node at the end of the focal edge
  child_node <- tree$edge[focal_edge, 2]
  
  # recursive function: collect all descendant nodes of a node
  get_desc_nodes <- function(node, edge_mat) {
    children <- edge_mat[edge_mat[, 1] == node, 2]
    if (length(children) == 0) {
      return(integer(0))
    } else {
      return(c(children, unlist(lapply(children, get_desc_nodes, edge_mat = edge_mat))))
    }
  }
  
  # all descendant nodes below the focal edge
  desc_nodes <- c(child_node, get_desc_nodes(child_node, tree$edge))
  
  # all edges whose parent node is in that descendant set
  desc_edges <- which(tree$edge[, 1] %in% desc_nodes)
  
  # include the focal edge itself
  sort(unique(c(focal_edge, desc_edges)))
}

#---------------------------------------------------------
# Build a 2-regime consensus from bayou branch.posteriors
#---------------------------------------------------------
bayou_consensus_highestpp <- function(tree, trait.rjMCMC,
                                      regime_desc = "regime1",
                                      regime_other = "regime2") {
  
  # bayou branch numbering is tied to postorder trees
  tree <- reorder(tree, "postorder")
  
  # summarize posterior
  smry <- summary(trait.rjMCMC)
  
  # extract the branch posterior table
  bp <- smry[["branch.posteriors"]]
  
  # safety check: rows of branch.posteriors should match rows of tree$edge
  if (nrow(bp) != nrow(tree$edge)) {
    stop("nrow(branch.posteriors) != nrow(tree$edge).")
  }
  
  # choose the branch with the highest posterior probability of a shift
  focal_edge <- which.max(bp$pp)
  
  # initialize all edges as the 'other' regime
  regimes <- rep(regime_other, nrow(tree$edge))
  
  # find focal edge + all descendant edges
  desc_edges <- get_descendant_edges(tree, focal_edge)
  
  # assign descendant clade to the focal regime
  regimes[desc_edges] <- regime_desc
  
  # return useful output
  list(
    tree = tree,
    regimes = regimes,
    focal_edge = focal_edge,
    focal_pp = bp$pp[focal_edge],
    focal_theta2 = bp$magnitude.of.theta2[focal_edge],
    focal_rel_location = bp$rel.location[focal_edge],
    branch.posteriors = bp,
    descendant_edges = desc_edges
  )
}

#------------------------------------------------------------------
# Create a list which store the regime info for each trait 
#------------------------------------------------------------------
### Extract regime for growth rate only 
trait.rjMCMC.1 <- Bayou_OUM_ExperiGR_11factors[[1]][["MCMC"]]
trait.rjMCMC.1 <- set.burnin(trait.rjMCMC.1,0.3)

cons1 <- bayou_consensus_highestpp(
  tree = tree_madin_scaled,
  trait.rjMCMC = trait.rjMCMC.1,
  regime_desc = "Trait1_regime1",
  regime_other = "Trait1_regime2"
)

### Extract regime for each trait
Regime_list <- list()
for (i in seq_along(predictors)) {
  pre <- predictors[i]
  trait.rjMCMC.2 <- Bayou_OUM_ExperiGR_11factors[[pre]][["MCMC"]]
  trait.rjMCMC.2 <- set.burnin(trait.rjMCMC.2,0.3)
  
  cons2 <- bayou_consensus_highestpp(
    tree = tree_madin_scaled,
    trait.rjMCMC = trait.rjMCMC.2,
    regime_desc = "Trait2_regime1",
    regime_other = "Trait2_regime2"
  )
  
  # Get the data frame of regime
  regime_df <- data.frame(
    edge_order = seq_len(length(cons1$regimes)),
    regime_trait1 = cons1$regimes, # Regime definition is only based on growth-rate regime
    regime_trait2 = cons2$regimes, # Regime definition is only based on trait-specific regime
    regime_trait1_trait2 = paste(cons1$regimes, cons2$regimes, sep = "__"), # Regime definition is based on growth*trait regime
    stringsAsFactors = FALSE
  )
  
  Regime_list[[pre]] <- regime_df
}

### Get the final regime annotation files ---
Regime_list








