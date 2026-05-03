library(tidyverse)
library(lubridate)
library(spOccupancy)
library(readr)
library(ggplot2)
library(patchwork)

# Raw data of the trap dataset
trap <- read.csv("C:/Users/Zuhan Shahria/OneDrive/Documents/R Studio Codes/Project/TrapEffort_2014-2017_NJD.csv")
# Raw data of the mammal detection dataset
detect <- read.csv("C:/Users/Zuhan Shahria/OneDrive/Documents/R Studio Codes/Project/Mammal_DetectionData_2014-2017_NJD.csv")



# Converting into date class
detect$Date <- dmy(detect$Date)
trap$Date_On  <- dmy(trap$Date_On)
trap$Date_Off <- dmy(trap$Date_Off)

# Renaming site variable
detect <- detect %>% rename(Site_ID = Location)

# New table with only the single observations
detect_ind <- detect %>% filter(New_Contact == "New Contact")



# Reclassifying Random sites by their Block value before merging, so if a site has a random habitat feature but its block is continuous logged forest then that is also the habitat feature
trap <- trap %>%
  mutate(Habitat_Feature = if_else(
    str_trim(Habitat_Feature) == "Random",
    str_trim(Block),   # replacing Random with Block label
    str_trim(Habitat_Feature)
  ))



# Joining the two datasets, so now we have a large table 
det_matched <- detect_ind %>%
  inner_join(trap, by = "Site_ID") %>%
  filter(Date >= Date_On & Date <= Date_Off)




# Creating a new class called week which represents each ocasion of our deployment
det_matched <- det_matched %>%
  mutate(week = floor(as.numeric(Date - Date_On)/ 7) + 1)



# Merging into ecological categories, replacing the current habitat features with the catogorical ones
trap <- trap %>%
  mutate(Habitat_Feature = case_when(
    
    # Trail types (excluding riparian trail)
    Habitat_Feature %in% c("Trail",
                           "Trail (ridge)",
                           "Trail convergence",
                           "Random/trail")              ~ "Trail",
    
    # Riparian types (including riparian trail)
    Habitat_Feature %in% c("Riparian",
                           "Riparian logged forest",
                           "Riparian trail",
                           "Logged Forest Control River") ~ "Riparian",
    
    # Forest types (including oil palm riparian reserve)
    Habitat_Feature %in% c("Forest fragment",
                           "Forest Block",
                           "Oil palm",
                           "Logged Forest Control",
                           "Oil palm Riparian Reserve")  ~ "Forest",
    
    # Clearings
    Habitat_Feature %in% c("Clearing", "Forest clearing") ~ "Clearing",
    
    # Disturbed tracks
    Habitat_Feature %in% c("Logging road", "Skid trail") ~ "Disturbed Track",
    
    
    # Anything remaining from Block reclassification
    TRUE ~ Habitat_Feature
  )) %>%
  mutate(Habitat_Feature = factor(Habitat_Feature,
                                  levels = c("Trail",
                                             "Clearing",
                                             "Riparian",
                                             "Forest",
                                             "Disturbed Track"
                                             )))









# Getting site list and species list
species <- sort(unique(det_matched$Species))
sites   <- sort(unique(trap$Site_ID))
n.week  <- max(det_matched$week, na.rm = TRUE)
weeks   <- 1:n.week

n.sp  <- length(species)
n.site <- length(sites)


# This is the list of detection matricies for each species
Y <- array(0, dim = c(n.sp, n.site, n.week))

# This FOR loop will create each detection matrix for each species and store them in the list
for(k in 1:n.sp){
  
  sp_data <- det_matched %>%
    filter(Species == species[k]) %>%
    group_by(Site_ID, week) %>%
    summarise(.groups = "drop")
  
  for(r in 1:nrow(sp_data)){
    i <- which(sites == sp_data$Site_ID[r])
    j <- which(weeks == sp_data$week[r])
    
    if (length(i) > 0 && length(j) > 0)
      Y[k, i, j] <- 1
  }
}

# Checking
dim(Y)
as.data.frame(Y[7, , ])

# Site covariates
site_covs <- trap %>%
  group_by(Site_ID) %>%
  summarise(Habitat_Feature = first(Habitat_Feature),.groups = "drop") %>%
  filter(Site_ID %in% sites) %>%
  arrange(match(Site_ID, sites))

# Removing sites with missing or blank habitat
good_sites <- !(is.na(site_covs$Habitat_Feature) | str_trim(as.character(site_covs$Habitat_Feature)) == "")

# Reassigning the detection array so that only the good sites remain
Y           <- Y[, good_sites, ]
site_covs   <- site_covs[good_sites, ]
sites_clean <- sites[good_sites]

cat("The number of sites which we can use", sum(good_sites), "\n")

habitat <- site_covs$Habitat_Feature

# Occupancy covariate list
occ.covs <- list(Habitat = habitat)


# Building a detection covariate matrix of dimensions [n.site x n.week]
# Habitat is a site level covariate which does not change across ocasions so each row of this matrix will have the same labels
habitat_det_matrix <- matrix(
  rep(as.character(habitat), times = n.week),
  nrow  = length(sites_clean),
  ncol  = n.week,
  byrow = FALSE
)



habitat_det_matrix <- apply(habitat_det_matrix, 2, function(col)
  factor(col, levels = levels(habitat)))

# Converting into a data frame
habitat_det_matrix <- as.data.frame(habitat_det_matrix)

# This list be passed to the model
det.covs <- list(Habitat = habitat_det_matrix)



# Marking unsampled site x week cells as NA

effort <- matrix(0, nrow = length(sites_clean), ncol = n.week)

for (s in seq_along(sites_clean)) {
  site_traps <- trap %>% filter(Site_ID == sites_clean[s])
  for (r in 1:nrow(site_traps)) {
    if (!is.na(site_traps$Date_On[r]) && !is.na(site_traps$Date_Off[r])) {
      active_days  <- as.numeric(site_traps$Date_Off[r] - site_traps$Date_On[r])
      active_weeks <- floor(active_days / 7) + 1
      active_weeks <- pmin(active_weeks, n.week)
      effort[s, 1:active_weeks] <- 1
    }
  }
}

for (k in 1:n.sp) {
  Y[k, , ][effort == 0] <- NA
}
for (j in 1:n.week) {
  det.covs$Habitat[[j]][effort[, j] == 0] <- NA
}


# The Model


# Occupancy covariates remain as habitat
data.list <- list(
  y        = Y,
  occ.covs = occ.covs,
  det.covs = NULL        # no detection covariates, species-specific constant
)

set.seed(42)
out <- msPGOcc(
  occ.formula = ~ Habitat,   # occupancy varies by habitat
  det.formula  = ~ 1,        # detection is intercept only, constant per species
  data         = data.list,
  n.samples    = 10000,
  n.burn       = 5000,
  n.thin       = 5,
  n.chains     = 3,
  verbose      = TRUE
)

# Saving after fitting
saveRDS(out, "fitted_model_final.rds")

# To reload in a future session without refitting:
# out <- readRDS("fitted_model_final.rds")





# CONVERGENCE CHECK/DIAGNOSTICS
# Values are close to 1.1
rhat_vals <- out$rhat
cat("Max R-hat (occupancy betas):", max(rhat_vals$beta, na.rm = TRUE), "\n")
cat("Max R-hat (detection alphas):", max(rhat_vals$alpha, na.rm = TRUE), "\n")





# Extracting dimensions
n.post  <- dim(out$beta.samples)[1]
n.coef  <- 5   # intercept + 4 habitat coefficients for occupancy
n.alpha <- 1   # intercept only for detection

cat("Posterior draws:", n.post, "\n")
cat("Occupancy coefficients:", n.coef, "\n")
cat("Detection coefficients:", n.alpha, "\n")

# Reshaping beta samples 
beta_array  <- array(out$beta.samples,  dim = c(n.post, n.sp, n.coef))

# Reshaping alpha samples
alpha_array <- array(out$alpha.samples, dim = c(n.post, n.sp, n.alpha))




# =============================================================================
# OCCUPANCY PROBABILITY PER SPECIES PER HABITAT
# =============================================================================
hab_names <- c("Trail", "Clearing", "Riparian", "Forest", "Disturbed Track")

habitat_types <- data.frame(
  Habitat = factor(hab_names, levels = levels(habitat))
)
X_habitat <- model.matrix(~ Habitat, data = habitat_types)
n.hab     <- nrow(X_habitat)

psi_by_habitat <- array(NA, dim = c(n.post, n.sp, n.hab))

for (d in 1:n.post) {
  for (k in 1:n.sp) {
    for (h in 1:n.hab) {
      beta_vec  <- beta_array[d, k, ]
      logit_psi <- sum(beta_vec * X_habitat[h, ])
      psi_by_habitat[d, k, h] <- 1 / (1 + exp(-logit_psi))
    }
  }
}

# Posterior mean and 95% credible interval for occupancy
psi_hab_mean  <- apply(psi_by_habitat, c(2, 3), mean)
psi_hab_lower <- apply(psi_by_habitat, c(2, 3), quantile, 0.025)
psi_hab_upper <- apply(psi_by_habitat, c(2, 3), quantile, 0.975)

dimnames(psi_hab_mean)  <- list(species, hab_names)
dimnames(psi_hab_lower) <- list(species, hab_names)
dimnames(psi_hab_upper) <- list(species, hab_names)





# =============================================================================
# DETECTION PROBABILITY PER SPECIES
# =============================================================================
p_by_species <- matrix(NA, nrow = n.post, ncol = n.sp)

for (d in 1:n.post) {
  for (k in 1:n.sp) {
    logit_p        <- alpha_array[d, k, 1]
    p_by_species[d, k] <- 1 / (1 + exp(-logit_p))
  }
}

# Posterior mean and 95% credible interval for detection per species
p_mean  <- apply(p_by_species, 2, mean)
p_lower <- apply(p_by_species, 2, quantile, 0.025)
p_upper <- apply(p_by_species, 2, quantile, 0.975)

names(p_mean)  <- species
names(p_lower) <- species
names(p_upper) <- species



# =============================================================================
# BUILDING TIDY DATA FRAMES
# =============================================================================

# Occupancy tidy data frame — species x habitat
psi_df <- as.data.frame(psi_hab_mean) %>%
  rownames_to_column("Species") %>%
  pivot_longer(-Species, names_to = "Habitat", values_to = "Mean") %>%
  left_join(
    as.data.frame(psi_hab_lower) %>%
      rownames_to_column("Species") %>%
      pivot_longer(-Species, names_to = "Habitat", values_to = "Lower"),
    by = c("Species", "Habitat")
  ) %>%
  left_join(
    as.data.frame(psi_hab_upper) %>%
      rownames_to_column("Species") %>%
      pivot_longer(-Species, names_to = "Habitat", values_to = "Upper"),
    by = c("Species", "Habitat")
  ) %>%
  mutate(Habitat = factor(Habitat, levels = hab_names))

# Detection tidy data frame
p_df <- data.frame(
  Species = species,
  Mean    = p_mean,
  Lower   = p_lower,
  Upper   = p_upper
)

# =============================================================================
# PLOT 1 — Occupancy probability per species per habitat
# =============================================================================
# Spliting species into groups of 10 for separate plots
species_list <- unique(psi_df$Species)
n_species    <- length(species_list)
group_size   <- 10
n_groups     <- ceiling(n_species / group_size)

for (g in 1:n_groups) {
  
  # Geting the species for this group
  start_idx    <- (g - 1) * group_size + 1
  end_idx      <- min(g * group_size, n_species)
  species_group <- species_list[start_idx:end_idx]
  
  # Filtering data to only include species in this group
  psi_df_group <- psi_df %>% filter(Species %in% species_group)
  
  # Building the plot for this group
  plot_psi_group <- ggplot(psi_df_group,
                           aes(x = Habitat, y = Mean, colour = Habitat)) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper),
                  width = 0.3, linewidth = 0.6) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               colour = "grey60", linewidth = 0.4) +
    facet_wrap(~ Species, ncol = 5) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title    = paste0("Occupancy Probability (\u03c8) per Species per Habitat",
                        " — Group ", g, " of ", n_groups),
      subtitle = "Points = posterior mean | Error bars = 95% credible interval",
      x        = "Habitat type",
      y        = "Occupancy probability (\u03c8)"
    ) +
    theme_bw(base_size = 10) +   # larger base size since fewer panels
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "none",
      strip.text      = element_text(size = 9, face = "bold")
    )
  
  print(plot_psi_group)
  
  # Saving each group as a separate file
  ggsave(paste0("occupancy_species_group_", g, ".png"),
         plot_psi_group,
         width  = 14,
         height = 8,
         dpi    = 200)
}

# =============================================================================
# PLOT 2 — Detection probability per species
# =============================================================================
plot_p <- ggplot(p_df,
                 aes(x = reorder(Species, Mean), y = Mean)) +
  geom_point(size = 2.5, colour = "steelblue") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.3, linewidth = 0.6, colour = "steelblue") +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             colour = "grey60", linewidth = 0.4) +
  scale_y_continuous(limits = c(0, 0.5)) +   # standardised y-axis
  coord_flip() +
  labs(
    title    = "Detection Probability (p) per Species",
    subtitle = "Species-specific constant | Points = posterior mean | Error bars = 95% CRI",
    x        = NULL,
    y        = "Detection probability (p)"
  ) +
  theme_bw(base_size = 9)

print(plot_p)
ggsave("detection_per_species.png", plot_p,
       width = 10, height = 12, dpi = 200)

# =============================================================================
# PLOT 3 — Community level mean occupancy per habitat
# =============================================================================
psi_community <- psi_df %>%
  group_by(Habitat) %>%
  summarise(
    Community_Mean  = mean(Mean),
    Community_Lower = mean(Lower),
    Community_Upper = mean(Upper),
    .groups = "drop"
  )

plot_community_psi <- ggplot(psi_community,
                             aes(x = reorder(Habitat, Community_Mean),
                                 y = Community_Mean,
                                 fill = Habitat)) +
  geom_col(alpha = 0.8) +
  geom_errorbar(aes(ymin = Community_Lower,
                    ymax = Community_Upper),
                width = 0.3) +
  scale_y_continuous(limits = c(0, 1)) +   # standardised y-axis
  labs(
    title    = "Community-Level Mean Occupancy by Habitat Type",
    subtitle = "Averaged across all species | Error bars = mean 95% CRI",
    x        = "Habitat type",
    y        = "Mean occupancy probability (ψ)"
  ) +
  coord_flip() +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

print(plot_community_psi)
ggsave("community_occupancy_by_habitat.png", plot_community_psi,
       width = 10, height = 6, dpi = 200)

# =============================================================================
# PLOT 4 — Occupancy heatmap across species and habitats
# =============================================================================
# Computing mean occupancy across all habitats per species
# This is used to order species on the y-axis from highest to lowest
species_order <- psi_df %>%
  group_by(Species) %>%
  summarise(overall_mean = mean(Mean, na.rm = TRUE), .groups = "drop") %>%
  arrange(overall_mean) %>%   # ascending so darkest appears at top with coordinate system
  pull(Species)

plot_heatmap_psi <- ggplot(psi_df,
                           aes(x = Habitat,
                               y = factor(Species, levels = species_order),
                               fill = Mean)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low      = "white",
    mid      = "steelblue",
    high     = "darkblue",
    midpoint = 0.5,
    name     = "\u03c8",
    limits   = c(0, 1)
  ) +
  labs(
    title = "Occupancy Probability Heatmap — Species x Habitat",
    x     = "Habitat type",
    y     = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_heatmap_psi)
ggsave("occupancy_heatmap.png", plot_heatmap_psi,
       width = 10, height = 12, dpi = 200)

# =============================================================================
# saving results 
# =============================================================================
write.csv(
  psi_df %>% mutate(across(c(Mean, Lower, Upper), ~ round(., 3))),
  "occupancy_by_species_habitat.csv",
  row.names = FALSE
)

write.csv(
  p_df %>% mutate(across(c(Mean, Lower, Upper), ~ round(., 3))),
  "detection_by_species.csv",
  row.names = FALSE
)







