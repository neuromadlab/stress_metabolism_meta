###########################################################################################
# Script: Meta-analysis on metabolic/hormonal modulators of cortisol stress reactivity
# Authors: Madeleine Kördel, Maria Meier, Anne Kühnel, Nils B. Kroemer
# Date: 2025-08-12
#
# Description:
# This R script performs a Bayesian meta-analysis to assess how glucose, progesterone,
# and estradiol influence cortisol reactivity to acute stress. Analyses include:
#   - Aggregating effect sizes
#   - Bayesian meta-analysis models
#   - Bayes factor robustness checks
#   - Posterior distribution visualization
#   - Funnel plots for publication bias
#   - Outlier detection and sensitivity analyses
###########################################################################################


# Clear global environment and load required packages
rm(list=ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(bayesmeta, compute.es, cowplot, data.table, dplyr, esc, forestplot, 
               ggplot2, knitr, MAd, metafor, reactable, readr, readxl, rmarkdown, 
               R.rsp, stringr, tidyr, scales, writexl, ggdist) 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Source helper functions
source("HelperFunctions.R")

# Load "df.RDa"
load(file = "df.RDa")

# Create separate data frames for each condition
df1 <- df %>% filter(Glucose==1) 
df2 <- df %>% filter(Progesterone==1 & `Sex diff` == 0)
df3 <- df %>% filter(Estradiol==1 & `Sex diff`==0)

df4 <- df %>% filter(Progesterone == 1)
df5 <- df %>% filter(Estradiol ==1)

### BAYES META-ANALYSIS: GLUCOSE ####
# This section:
# 1. Aggregates effect sizes for studies where glucose was the main variable
# 2. Fits a Bayesian meta-analysis model with a half-Cauchy prior on tau (scale = 0.5)
# 3. Uses a normal prior for the overall effect size mu (mean = 0, sd = 1.5)
# 4. Produces a forest plot of the model estimates

# Aggregate effect sizes for glucose
aggES1 <- agg(id     = study,
              es     = yi,
              var    = vi,
              data   = df1,
              cor = .5,
              method = "BHHR")

# Merge aggregated effect sizes with original data
MA1 <- merge(x = aggES1, y = df1, by.x = "id", by.y = "study")
MA1 <- unique(setDT(MA1) [sort.list(id)], by = "id")
MA1 <- with(MA1, MA1[order(MA1$es)])

# Generate bayesmeta-object "bma1", which stores all relevant results
bma1 <- bayesmeta(y = MA1$es,sigma = sqrt(MA1$var), labels = MA1$id, 
                  tau.prior = function(t) dhalfcauchy(t, scale = 0.5), FE =T,
                  mu.prior = c("mean" = 0, "sd" = 1.5))

# OUTPUT: Forest plot for glucose meta-analysis (main text, Figure 3)
forestplot(bma1)


### BAYES META-ANALYSIS: PROGESTERONE ####
# This section:
# 1. Aggregates effect sizes for studies where Progesterone was the main variable
# 2. Fits a Bayesian meta-analysis model with a half-Cauchy prior on tau (scale = 0.5)
# 3. Uses a normal prior for the overall effect size mu (mean = 0, sd = 1.5)
# 4. Produces a forest plot of the model estimates

# Aggregate effect sizes for progesterone
aggES2 <- agg(id     = study,
              es     = yi,
              var    = vi,
              data   = df2,
              cor = .5,
              method = "BHHR")
# Merge aggregated effect sizes with original data
MA2 <- merge(x = aggES2, y = df2, by.x = "id", by.y = "study")
MA2 <- unique(setDT(MA2) [sort.list(id)], by = "id")
MA2 <- with(MA2, MA2[order(MA2$es)])


# Generate bayesmeta-object "bma2", which stores all relevant results
bma2 <- bayesmeta(y = MA2$es,sigma = sqrt(MA2$var), labels = MA2$id, 
                  tau.prior = function(t) dhalfcauchy(t, scale = 0.5), FE =T,
                  mu.prior = c("mean" = 0, "sd" = 1.5))

# OUTPUT: Forest plot for glucose meta-analysis (main text, Figure 5)
forestplot(bma2)


### BAYES META-ANALYSIS: ESTRADIOL ####
# This section:
# 1. Aggregates effect sizes for studies where Estradiol was the main variable
# 2. Fits a Bayesian meta-analysis model with a half-Cauchy prior on tau (scale = 0.5)
# 3. Uses a normal prior for the overall effect size mu (mean = 0, sd = 1.5)
# 4. Produces a forest plot of the model estimates

# Aggregate effect sizes for estradiol
aggES3 <- agg(id     = study,
              es     = yi,
              var    = vi,
              data   = df3,
              cor = .5,
              method = "BHHR")

#  Merge aggregated effect sizes with original data
MA3 <- merge(x = aggES3, y = df3, by.x = "id", by.y = "study")
MA3 <- unique(setDT(MA3) [sort.list(id)], by = "id")
MA3 <- with(MA3, MA3[order(MA3$es)])

# Generate bayesmeta-object "bma3", which stores all relevant results
bma3 <- bayesmeta(y = MA3$es,sigma = sqrt(MA3$var), labels = MA3$id, 
                  tau.prior = function(t) dhalfcauchy(t, scale = 0.5), FE =T,
                  mu.prior = c("mean" = 0, "sd" = 1.5))

# OUTPUT: Forest plot for glucose meta-analysis (main text, Figure 5)
forestplot(bma3)



### Bayes Factor Robustness Check ####
# This section evaluates how robust the Bayes factor is for Glucose, Progesterone and Estradiol meta-analysis
# by varying the prior standard deviation for the effect size (mu).

#### Glucose ####

SD = 1.5  # Baseline prior SD for mu from bma1$mu.prior

# Bayes factor with a narrow prior (half the baseline SD)
narrow <- bayesmeta(
  y = MA1$es,               # observed effect sizes
  sigma = sqrt(MA1$var),    # standard errors from variance
  labels = MA1$id,           # study labels
  tau.prior = bma1$tau.prior,  # prior for heterogeneity
  mu.prior = c("mean" = 0, "sd" = (SD/2))  # narrow prior for mu
)

# Bayes factor with default prior (baseline SD)
default <- bayesmeta(y = MA1$es,sigma = sqrt(MA1$var), labels = MA1$id, 
                     tau.prior = bma1$tau.prior, 
                     mu.prior = c("mean" = 0, "sd" = SD)) # default prior

# Bayes factor with wide prior (SD + 1)
wide <- bayesmeta(y = MA1$es,sigma = sqrt(MA1$var), labels = MA1$id, 
                  tau.prior = bma1$tau.prior, 
                  mu.prior = c("mean" = 0, "sd" = SD+1))  # wider prior

# Bayes factor with ultra-wide prior (SD + 2)
ultrawide <- bayesmeta(y = MA1$es,sigma = sqrt(MA1$var), labels = MA1$id, 
                       tau.prior = bma1$tau.prior,
                       mu.prior = c("mean" = 0, "sd" = SD+2)) # ultra-wide prior

# Extract Bayes Factors (comparing alternative vs. null hypothesis)
BFS1 <- c(narrow$bayesfactor[1,2], default$bayesfactor[1,2], wide$bayesfactor[1,2], ultrawide$bayesfactor[1,2])

# Convert BF01 to BF10
BFS1_BF10 <- 1 / BFS1
print(BFS1_BF10)

# Combine Bayes Factors into a data frame
BFS1_BF10 <- data.frame(
  Condition = c("Glucose"),
  Narrow = c(BFS1_BF10[1]),
  Default = c(BFS1_BF10[2]),
  Wide = c(BFS1_BF10[3]),
  UltraWide = c(BFS1_BF10[4]))


#### Progesterone ####

SD=1.5 # Baseline prior SD for mu from bma2$mu.prior

# Bayes factor with a narrow prior (half the baseline SD)
narrow <- bayesmeta(y = MA2$es,sigma = sqrt(MA2$var), labels = MA2$id, 
                    tau.prior = bma2$tau.prior, #tauprior, 
                    mu.prior = c("mean" = 0, "sd" = (SD/2))) 
# Bayes factor with default prior (baseline SD)
default <- bayesmeta(y = MA2$es,sigma = sqrt(MA2$var), labels = MA2$id, 
                     tau.prior = bma2$tau.prior, #tauprior, 
                     mu.prior = c("mean" = 0, "sd" = SD))
# Bayes factor with wide prior (SD + 1)
wide <- bayesmeta(y = MA2$es,sigma = sqrt(MA2$var), labels = MA2$id, 
                  tau.prior = bma2$tau.prior, #tauprior, 
                  mu.prior = c("mean" = 0, "sd" = SD+1))
# Bayes factor with ultra-wide prior (SD + 2)
ultrawide <- bayesmeta(y = MA2$es,sigma = sqrt(MA2$var), labels = MA2$id, 
                       tau.prior = bma2$tau.prior, #tauprior,  
                       mu.prior = c("mean" = 0, "sd" = SD+2))

# Extract Bayes Factors (comparing alternative vs. null hypothesis)
BFS2 <- c(narrow$bayesfactor[1,2], default$bayesfactor[1,2], wide$bayesfactor[1,2], ultrawide$bayesfactor[1,2])


#### Estradiol ####

SD=1.5  # Baseline prior SD for mu from bma2$mu.prior

# Bayes factor with a narrow prior (half the baseline SD)
narrow <- bayesmeta(y = MA3$es,sigma = sqrt(MA3$var), labels = MA3$id, 
                    tau.prior = bma3$tau.prior, #tauprior, 
                    mu.prior = c("mean" = 0, "sd" = (SD/2))) 
# Bayes factor with default prior (baseline SD)
default <- bayesmeta(y = MA3$es,sigma = sqrt(MA3$var), labels = MA3$id, 
                     tau.prior = bma3$tau.prior, #tauprior, 
                     mu.prior = c("mean" = 0, "sd" = SD))
# Bayes factor with wide prior (SD + 1)
wide <- bayesmeta(y = MA3$es,sigma = sqrt(MA3$var), labels = MA3$id, 
                  tau.prior = bma3$tau.prior, #tauprior, 
                  mu.prior = c("mean" = 0, "sd" = SD+1))
# Bayes factor with ultra-wide prior (SD + 2)
ultrawide <- bayesmeta(y = MA3$es,sigma = sqrt(MA3$var), labels = MA3$id, 
                       tau.prior = bma3$tau.prior, #tauprior,  
                       mu.prior = c("mean" = 0, "sd" = SD+2))

# Extract Bayes Factors (comparing alternative vs. null hypothesis)
BFS3 <- c(narrow$bayesfactor[1,2], default$bayesfactor[1,2], wide$bayesfactor[1,2], ultrawide$bayesfactor[1,2])

# Combine all Bayes Factors into a data frame
bayes_factors <- data.frame(
  Condition = c("Glucose", "Progesterone", "Estradiol"),
  Narrow = c(BFS1[1], BFS2[1], BFS3[1]),
  Default = c(BFS1[2], BFS2[2], BFS3[2]),
  Wide = c(BFS1[3], BFS2[3], BFS3[3]),
  UltraWide = c(BFS1[4], BFS2[4], BFS3[4])
)


### plot bayes factors ####
# This section creates a plot of the Bayes Factors for Glucose, Progesterone, and Estradiol meta-analyses
# OUTPUT: bayes factors plot (main text, Figure 6)

# Convert Bayes Factors to a data frame for plotting
priors <- data.frame(
  prior_sd = c(0.75, 1.5, 2.5, 3.5))

# Combine Bayes Factors into a data frame
BFs <- cbind(priors, Glucose = BFS1, Progesterone = BFS2, Estradiol = BFS3)

# Convert to long format for ggplot
BFs_long <- BFs %>%
  pivot_longer(cols = -prior_sd, names_to = "Analysis", values_to = "BF01") %>%
  mutate(logBF01 = log(BF01))  # For natural log of Bayes Factors

# Create the plot
ggplot(BFs_long, aes(x = prior_sd, y = logBF01, color = Analysis)) +
  
  # Lines and points with better styling
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 3, alpha = 0.3, color = "black") +
  
  # Improved text labels with better positioning and styling
  geom_text(aes(label = round(logBF01, 2)), 
            vjust = -0.8, hjust = 0.5, size = 4, 
            fontface = "bold", show.legend = FALSE, color = "black") +
  geom_text(data = subset(BFs_long, Analysis == "Glucose"), 
            aes(label = c("narrow", "default", "wide", "ultrawide")), 
            vjust = 1.5, hjust = 0.5, size = 4, 
            fontface = "plain", show.legend = FALSE, color = "black") +
  
  # Reference lines with labels
  geom_hline(yintercept = 1, linetype = "dashed", 
             color = "darkgray", linewidth = 0.8, alpha = 0.7) +
  geom_hline(yintercept = -1, linetype = "dashed", 
             color = "darkgray", linewidth = 0.8, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "darkgray", linewidth = 0.8, alpha = 0.7) +
  
  # Add annotations for reference lines
  annotate("text", x = max(BFs_long$prior_sd) * 1, y = log(6.4) + 0.15, 
           label = "Moderate evidence for H0", size = 5, color = "darkgray") +
  annotate("text", x = max(BFs_long$prior_sd) * 1, y = log(2) + 0.15, 
           label = "Anecdotal evidence for H0", size = 5, color = "darkgray") +
  annotate("text", x = max(BFs_long$prior_sd) * 1, y = log(1/6.4) - 0.15, 
           label = "Moderate evidence for H1", size = 5, color = "darkgray") +
  annotate("text", x = max(BFs_long$prior_sd) * 1, y = log(1/2) - 0.15, 
           label = "Anecdotal evidence for H1", size = 5, color = "darkgray") +
  
  # Improved scales
  scale_x_continuous(
    name = "Standard deviation", 
    breaks =  BFs$prior_sd,
    expand = expansion(mult = c(0.05, 0.15))  # Add space for labels
  ) +
  scale_y_continuous(
    name = "ln(BF01)",
    breaks = seq(-2.0, 2.0, by = 1.0),
    expand = expansion(mult = c(0.3, 0.2))
  ) +
  
  # Better color palette
  scale_color_viridis_d(name = NULL, option = "plasma", end = 0.8) +
  
  # Enhanced theme
  theme_classic(base_size = 12) +
  theme(
    
    # Legend improvements
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),   # Anchor point: left, top
    legend.title = element_blank(),   # Remove legend title completely
    legend.text = element_text(size = 15, face = "bold"),
    
    # Axis improvements
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 15),
    axis.line = element_line(color = "black", linewidth = 0.7))


### posterior distributions ####
# This section extracts posterior samples from the Bayesian meta-analysis models
# and visualizes the posterior distributions for each condition (Glucose, Progesterone, Estradiol).
# OUTPUT: posterior distributions plot (main text, Figure 7)

# Extract posterior samples
samples1 <- data.frame(bma1$rposterior(5000))
samples1$analysis <- "Glucose"

samples2 <- data.frame(bma2$rposterior(5000))
samples2$analysis <- "Progesterone"

samples3 <- data.frame(bma3$rposterior(5000))
samples3$analysis <- "Estradiol"


# Combine all samples
posterior_all <- rbind(samples1, samples2, samples3)

# change order
posterior_all <- bind_rows(samples1, samples2, samples3) %>%
  mutate(analysis = factor(analysis, levels = c("Estradiol", "Progesterone", "Glucose")))

custom_colors <- c(
  "Glucose" = "#B12A90FF",
  "Progesterone" = "#FCA636FF",
  "Estradiol" = "#0D0887FF"
)

# Plot posterior distributions
ggplot(posterior_all, aes(x = mu, y = analysis, fill = analysis)) +
  # Reference line
  geom_vline(
    xintercept = 0, 
    linetype = "dashed", 
    color = "gray20", 
    linewidth = 0.8,
    alpha = 0.7
  ) +
  # Use stat_halfeye for better density visualization
  stat_halfeye(
    adjust = 1.3,  # Increased for smoother densities
    .width = c(0.5, 0.8, 0.95),  # Multiple credible intervals
    point_interval = median_qi,
    alpha = 0.85,
    normalize = "panels"  # Normalize within each panel
  ) +
  # Improved color palette
  scale_fill_manual(values = custom_colors) +
  # Enhanced axis formatting
  scale_x_continuous(
    limits = c(-0.75, 0.75),
    breaks = seq(-0.75, 0.75, by = 0.25),
    labels = function(x) sprintf("%.2f", x),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  # Improved labels
  labs(
    x = expression(bold("Posterior distribution of ") * bolditalic(mu)),
    y = NULL,
    title = NULL,
    caption = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(face = "bold", size = 16, color = "black"),
    axis.text = element_text(size = 15, color = "black")
  )


### funnel plot #####
# This section creates a funnel plot to visualize publication bias in the meta-analysis results.
# The funnel plot includes studies on glucose, progesterone, and estradiol, with different colors for each condition.
# OUTPUT: funnel plot (main text, Figure 4)

# Fit models
res1 <- rma(yi = es, vi = var, data = MA1, slab = MA1$id)
res2 <- rma(yi = es, vi = var, data = MA2, slab = MA2$id)
res3 <- rma(yi = es, vi = var, data = MA3, slab = MA3$id)

# Draw base funnel plot using one dataset (e.g., MA1)
par(mar = c(5, 5, 3, 3) + 0.1)

funnel(res1,
       level = c(90, 95, 99),
       shade = c("white", "gray80", "gray60"),
       xlab = "Observed outcome",
       ylab = "Standard error",
       main = "",
       legend = list(studies=FALSE),
       refline = 0,
       cex.lab = 1.7,
       cex.axis = 1.5,
       font.lab = 2, 
       label = "",
       back = "white")

# Add glucose studies
points(res1$yi, sqrt(res1$vi), pch = 16, col = "#B12A90FF")
to_label1 <- which(res1$yi > 1)
text(res1$yi[to_label1], sqrt(res1$vi[to_label1]),
     labels = res1$slab[to_label1],
     pos = 2, cex = 1.2, col = "#B12A90FF", font = 2)

# Add progesterone studies
points(res2$yi, sqrt(res2$vi), pch = 17, col = "#FCA636FF")
to_label2 <- which(res2$yi > 1 | res2$yi < -0.5)
text(res2$yi[to_label2], sqrt(res2$vi[to_label2]),
     labels = res2$slab[to_label2],
     pos = ifelse(res2$yi[to_label2] > 0, 4, 2),
     cex = 1.2, col = "#FCA636FF", font = 2)

# Add estradiol studies
points(res3$yi, sqrt(res3$vi), pch = 18, col = "#0D0887FF")
to_label3 <- which(res3$yi > 0.5 | res3$yi < -0.5)
text(res3$yi[to_label3], sqrt(res3$vi[to_label3]),
     labels = res3$slab[to_label3],
     pos = 4, cex = 1.2, col = "#0D0887FF", font = 2)

# Add a legend
legend("topleft",
       legend = c("Glucose", "Progesterone", "Estradiol"),
       col = c("#B12A90FF", "#FCA636FF", "#0D0887FF"),
       pch = c(16, 17, 18),
       cex = 1.5,
       bty = "n") 



### outlier plots ####
# This section identifies and visualizes outliers in the meta-analysis data for glucose, progesterone, and estradiol.
# The outlier detection is based on the interquartile range (IQR) method.

# Function to identify outliers based on IQR
is_outlier <- function(x) {
  return(x < quantile(x, 0.25, na.rm = TRUE) - 1.5 * IQR(x, na.rm = TRUE) |
           x > quantile(x, 0.75, na.rm = TRUE) + 1.5 * IQR(x, na.rm = TRUE))
}

#### plot outliers glucose ####

# create outlier dataframe for glucose (o = outlier)
MA1_o <- MA1 %>% tibble::rownames_to_column(var = "outlier") %>% mutate(is_outlier=ifelse(is_outlier(es), es, as.numeric(NA)))
MA1_o$id[which(is.na(MA1_o$is_outlier))] <- as.numeric(NA)

# plot outliers glucose
ggplot(MA1_o, aes(x = factor(0), es)) +
  geom_boxplot(outlier.size = 3.5, outlier.colour = "#D55E00", outlier.shape = 18, fill = "lightgrey") +
  geom_text(aes(label=id),na.rm = T, nudge_y = 0.02, nudge_x = 0.05) +
  stat_boxplot(geom="errorbar", width = 0.05) +
  scale_x_discrete(breaks = NULL) +
  xlab(NULL) + ylab("Hedges' g") +
  theme_minimal_hgrid(12)


#### plot outliers progesterone #### 

# create outlier dataframe for progesterone (o = outlier)
MA2_o <- MA2 %>% tibble::rownames_to_column(var = "outlier") %>% mutate(is_outlier=ifelse(is_outlier(es), es, as.numeric(NA)))
MA2_o$id[which(is.na(MA2_o$is_outlier))] <- as.numeric(NA)

# plot outliers progesterone
ggplot(MA2_o, aes(x = factor(0), es)) +
  geom_boxplot(outlier.size = 3.5, outlier.colour = "#D55E00", outlier.shape = 18, fill = "lightgrey") +
  geom_text(aes(label=id),na.rm = T, nudge_y = 0.02, nudge_x = 0.05) +
  stat_boxplot(geom="errorbar", width = 0.05) +
  scale_x_discrete(breaks = NULL) +
  xlab(NULL) + ylab("Hedges' g") +
  theme_minimal_hgrid(12)

#### plot outliers estradiol ####

# create outlier dataframe for estradiol (o = outlier)
MA3_o <- MA3 %>% tibble::rownames_to_column(var = "outlier") %>% mutate(is_outlier=ifelse(is_outlier(es), es, as.numeric(NA)))
MA3_o$id[which(is.na(MA3_o$is_outlier))] <- as.numeric(NA)

# plot outliers estradiol
ggplot(MA3_o, aes(x = factor(0), es)) +
  geom_boxplot(outlier.size = 3.5, outlier.colour = "#D55E00", outlier.shape = 18, fill = "lightgrey") +
  geom_text(aes(label=id),na.rm = T, nudge_y = 0.02, nudge_x = 0.05) +
  stat_boxplot(geom="errorbar", width = 0.05) +
  scale_x_discrete(breaks = NULL) +
  xlab(NULL) + ylab("Hedges' g") +
  theme_minimal_hgrid(12)



### re-run bayes meta-analysis without outliers #####
# Outlier exclusion criteria:
# Studies were flagged if effect sizes fell outside 1.5 × IQR from the 25th or 75th percentile
# The following studies were excluded based on this rule (possible measurement bias or extreme values):

#### Glucose ####

# exclude outliers from MA1 - Glucose ##
excluded_MA1 <- c("Gonzalez-Bono, 2002", "Kirschbaum, 1997")

# Exclude specific studies and run bayes meta (no_o = outlier are excluded)
bma1_no_o <- MA1 %>%
  filter(!id %in% excluded_MA1) %>%
  { bayesmeta(
    y = .$es,
    sigma = sqrt(.$var),
    labels = .$id,
    tau.prior = function(t) dhalfcauchy(t, scale = 0.5),
    FE = TRUE,
    mu.prior = c("mean" = 0, "sd" = 1.5)
  )
  }
forestplot(bma1_no_o)

#### Progesterone ####

# exclude outliers from MA2 - Progesterone
excluded_MA2 <- c("Boisseau, 2013")

# Exclude specific studies and run bayesmeta (no_o = outlier are excluded)
bma2_no_o <- MA2 %>%
  filter(!id %in% excluded_MA2) %>%
  { bayesmeta(
    y = .$es,
    sigma = sqrt(.$var),
    labels = .$id,
    tau.prior = function(t) dhalfcauchy(t, scale = 0.5),
    FE = TRUE,
    mu.prior = c("mean" = 0, "sd" = 1.5)
  )
  }
forestplot(bma2_no_o)


#### Estradiol ####

# exclude outliers from MA3 - Estradiol
excluded_MA3 <- c("Boisseau, 2013", "Bürger, 2025", "Barel, 2018")

# Exclude specific studies and run bayesmeta (no_o = outlier are excluded)
bma3_no_o <- MA3 %>%
  filter(!id %in% excluded_MA3) %>%
  { bayesmeta(
    y = .$es,
    sigma = sqrt(.$var),
    labels = .$id,
    tau.prior = function(t) dhalfcauchy(t, scale = 0.5),
    FE = TRUE,
    mu.prior = c("mean" = 0, "sd" = 1.5)
  )
  }
forestplot(bma3_no_o)


### re-run bayes factor without outliers ####
# This section evaluates how robust the Bayes factor is for Glucose, Progesterone and Estradiol meta-analysis
# by varying the prior standard deviation for the effect size (mu) after excluding outliers.

## Bayes factor robustness without outliers for 
#### Glucose ####

# Dataset glucose without outliers (no_o = outlier are excluded)
MA1_no_o <- MA1 %>% filter(!id %in% excluded_MA1)

SD=1.5 # from bma1_no_o$mu.prior (no outlier)

narrow <- bayesmeta(y = MA1_no_o$es,sigma = sqrt(MA1_no_o$var), labels = MA1_no_o$id, 
                    tau.prior = bma1_no_o$tau.prior, #tauprior, 
                    mu.prior = c("mean" = 0, "sd" = (SD/2))) 
default <- bayesmeta(y = MA1_no_o$es,sigma = sqrt(MA1_no_o$var), labels = MA1_no_o$id, 
                     tau.prior = bma1_no_o$tau.prior, #tauprior, 
                     mu.prior = c("mean" = 0, "sd" = SD))
wide <- bayesmeta(y = MA1_no_o$es,sigma = sqrt(MA1_no_o$var), labels = MA1_no_o$id, 
                  tau.prior = bma1_no_o$tau.prior, #tauprior, 
                  mu.prior = c("mean" = 0, "sd" = SD+1))
ultrawide <- bayesmeta(y = MA1_no_o$es,sigma = sqrt(MA1_no_o$var), labels = MA1_no_o$id, 
                       tau.prior = bma1_no_o$tau.prior, #tauprior,  
                       mu.prior = c("mean" = 0, "sd" = SD+2))

# Extract Bayes Factors (against null) for glucose
BFS1_no_o <- c(narrow$bayesfactor[1,2], default$bayesfactor[1,2], wide$bayesfactor[1,2], ultrawide$bayesfactor[1,2])
print(BFS1_no_o)

# Convert BFS1_no_o to BF10
BFS1_no_o_BF10 <- 1 / BFS1_no_o
print(BFS1_no_o_BF10)

# Combine Bayes Factors into a data frame
BFS1_no_o_BF10 <- data.frame(
  Condition = c("Glucose"),
  Narrow = c(BFS1_no_o_BF10[1]),
  Default = c(BFS1_no_o_BF10[2]),
  Wide = c(BFS1_no_o_BF10[3]),
  UltraWide = c(BFS1_no_o_BF10[4])
)

## Bayes factor robustness without outliers for
#### Progesterone ####

# Dataset progesterone without outliers (no_o = outlier are excluded)
MA2_no_o <- MA2 %>% filter(!id %in% excluded_MA2)

SD=1.5 # from bma2_no_o$mu.prior (no outliers)

narrow <- bayesmeta(y = MA2_no_o$es,sigma = sqrt(MA2_no_o$var), labels = MA2_no_o$id, 
                    tau.prior = bma2_no_o$tau.prior, #tauprior, 
                    mu.prior = c("mean" = 0, "sd" = (SD/2)))
default <- bayesmeta(y = MA2_no_o$es,sigma = sqrt(MA2_no_o$var), labels = MA2_no_o$id, 
                     tau.prior = bma2_no_o$tau.prior, #tauprior, 
                     mu.prior = c("mean" = 0, "sd" = SD))
wide <- bayesmeta(y = MA2_no_o$es,sigma = sqrt(MA2_no_o$var), labels = MA2_no_o$id,
                  tau.prior = bma2_no_o$tau.prior, #tauprior, 
                  mu.prior = c("mean" = 0, "sd" = SD+1))
ultrawide <- bayesmeta(y = MA2_no_o$es,sigma = sqrt(MA2_no_o$var), labels = MA2_no_o$id, 
                       tau.prior = bma2_no_o$tau.prior, #tauprior,  
                       mu.prior = c("mean" = 0, "sd" = SD+2))

# Extract Bayes Factors (against null) for progesterone
BFS2_no_o <- c(narrow$bayesfactor[1,2], default$bayesfactor[1,2], wide$bayesfactor[1,2], ultrawide$bayesfactor[1,2])
print(BFS2_no_o)


## Bayes factor robustness without outliers for
#### Estradiol ####

# Dataset estradiol without outliers (no_o = outlier are excluded)
MA3_no_o <- MA3 %>% filter(!id %in% excluded_MA3)

SD=1.5 # from bma3_no_o$mu.prior (no outliers)

narrow <- bayesmeta(y = MA3_no_o$es,sigma = sqrt(MA3_no_o$var), labels = MA3_no_o$id, 
                    tau.prior = bma3_no_o$tau.prior, #tauprior, 
                    mu.prior = c("mean" = 0, "sd" = (SD/2)))
default <- bayesmeta(y = MA3_no_o$es,sigma = sqrt(MA3_no_o$var), labels = MA3_no_o$id, 
                     tau.prior = bma3_no_o$tau.prior, #tauprior, 
                     mu.prior = c("mean" = 0, "sd" = SD))
wide <- bayesmeta(y = MA3_no_o$es,sigma = sqrt(MA3_no_o$var), labels = MA3_no_o$id,
                  tau.prior = bma3_no_o$tau.prior, #tauprior, 
                  mu.prior = c("mean" = 0, "sd" = SD+1))
ultrawide <- bayesmeta(y = MA3_no_o$es,sigma = sqrt(MA3_no_o$var), labels = MA3_no_o$id, 
                       tau.prior = bma3_no_o$tau.prior, #tauprior,  
                       mu.prior = c("mean" = 0, "sd" = SD+2))

# Extract Bayes Factors (against null) for estradiol
BFS3_no_o <- c(narrow$bayesfactor[1,2], default$bayesfactor[1,2], wide$bayesfactor[1,2], ultrawide$bayesfactor[1,2])
print(BFS3_no_o)


# Combine Bayes Factors into a data frame
bayes_factors_no_o <- data.frame(
  Condition = c("Glucose", "Progesterone", "Estradiol"),
  Narrow = c(BFS1_no_o[1], BFS2_no_o[1], BFS3_no_o[1]),
  Default = c(BFS1_no_o[2], BFS2_no_o[2], BFS3_no_o[2]),
  Wide = c(BFS1_no_o[3], BFS2_no_o[3], BFS3_no_o[3]),
  UltraWide = c(BFS1_no_o[4], BFS2_no_o[4], BFS3_no_o[4])
)