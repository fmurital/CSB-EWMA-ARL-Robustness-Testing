# CSB-EWMA ARL₁ Robustness Testing - Multiple Distributions
# Testing performance across different data-generating processes

library(parallel)
library(doParallel)
library(foreach)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stats)

# Define distributions to test
distributions <- c("normal", "laplace", "uniform", "exponential")
dist_names <- c("Normal", "Laplace", "Uniform", "Exponential")

# Define the shift magnitudes (matching the image format)
delta_vec <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50)

# Fixed parameters
p0 <- 0.5  # In-control proportion
k <- 10    # Number of streams

# Define optimal parameter combinations (See https://arxiv.org/pdf/2601.09968)
optimal_combinations <- data.frame(
  lambda = c(0.175, 0.2, 0.35, 0.4, 0.55, 0.65, 0.725, 0.8, 0.9,
             0.15, 0.2, 0.3, 0.45, 0.525, 0.65, 0.7, 0.85, 0.925),
  L = c(1.375, 1.4, 1.5, 1.525, 1.575, 1.625, 1.65, 1.675, 1.7,
        1.525, 1.575, 1.65, 1.725, 1.75, 1.8, 1.8, 1.875, 1.875),
  ARL0_target = c(rep(370, 9), rep(500, 9)),
  Group = c(rep("ARL0_370", 9), rep("ARL0_500", 9)),
  stringsAsFactors = FALSE
)

# Create a descriptive label for each combination 
optimal_combinations$combo_label <- sprintf("λ=%.3f, L=%.3f", 
                                            optimal_combinations$lambda, 
                                            optimal_combinations$L)

# Create shorter labels for table display (like "Job 80(4.1-170)" in image)
optimal_combinations$table_label <- sprintf("λ=%.3f, L=%.3f", 
                                            optimal_combinations$lambda, 
                                            optimal_combinations$L)

cat("Total optimal combinations:", nrow(optimal_combinations), "\n")

############################# DATA GENERATION FUNCTIONS #############################

# Helper function for Laplace distribution
rlaplace <- function(n, location = 0, scale = 1) {
  u <- runif(n, -0.5, 0.5)
  return(location - scale * sign(u) * log(1 - 2 * abs(u)))
}

# Function to generate continuous data from different distributions
generate_continuous_data <- function(distribution, n, shift = 0, p0 = 0.5) {
  if (distribution == "normal") {
    threshold_in <- qnorm(p0)
    p1 <- p0 + shift
    p1 <- min(max(p1, 0.01), 0.99)
    mean_shift <- threshold_in - qnorm(1 - p1)
    data <- rnorm(n, mean = mean_shift, sd = 1)
    
  } else if (distribution == "laplace") {
    threshold_in <- ifelse(p0 < 0.5, log(2 * p0), -log(2 * (1 - p0)))
    p1 <- p0 + shift
    p1 <- min(max(p1, 0.01), 0.99)
    location <- threshold_in + log(2 * (1 - p1))
    data <- rlaplace(n, location = location, scale = 1)
    
  } else if (distribution == "uniform") {
    threshold_in <- p0
    p1 <- p0 + shift
    p1 <- min(max(p1, 0.01), 0.99)
    shift_amount <- threshold_in + p1 - 1
    data <- runif(n, min = shift_amount, max = 1 + shift_amount)
    
  } else if (distribution == "exponential") {
    threshold_in <- qexp(p0, rate = 1)
    p1 <- p0 + shift
    p1 <- min(max(p1, 0.01), 0.99)
    rate <- -log(1 - p1) / threshold_in
    data <- rexp(n, rate = rate)
    
  } else {
    stop("Unknown distribution:", distribution)
  }
  
  return(data)
}

# Function to dichotomize continuous data to binary (0/1)
dichotomize_data <- function(data, distribution, p0 = 0.5) {
  if (distribution == "normal") {
    threshold <- qnorm(p0)
  } else if (distribution == "laplace") {
    threshold <- ifelse(p0 < 0.5, log(2 * p0), -log(2 * (1 - p0)))
  } else if (distribution == "uniform") {
    threshold <- p0
  } else if (distribution == "exponential") {
    threshold <- qexp(p0, rate = 1)
  } else {
    threshold <- quantile(data, p0)
  }
  
  binary_data <- as.integer(data > threshold)
  return(binary_data)
}

############################# ROBUST ARL₁ SIMULATION #############################

simulate_csb_ewma_arl1_robust <- function(lambda, L, distribution, delta,
                                          num_simulations = 50000, 
                                          max_time = 1000) {
  tryCatch({
    run_lengths <- numeric(num_simulations)
    mu0 <- k * p0
    sigma2_0 <- k * p0 * (1 - p0)
    
    p1 <- p0 + delta
    p1 <- min(max(p1, 0.01), 0.99)
    
    for (sim in 1:num_simulations) {
      t <- 1
      signal <- FALSE
      r_previous <- 0
      cumulative_sum <- 0
      
      while (!signal && t <= max_time) {
        continuous_data <- generate_continuous_data(distribution, k, shift = delta, p0 = p0)
        binary_data <- dichotomize_data(continuous_data, distribution, p0)
        C_t <- sum(binary_data)
        cumulative_sum <- cumulative_sum + C_t
        W_t <- (cumulative_sum - mu0 * t) / sqrt(t * sigma2_0)
        r_t <- lambda * W_t + (1 - lambda) * r_previous
        
        UCL <- L
        LCL <- -L
        
        if (r_t > UCL || r_t < LCL) {
          run_lengths[sim] <- t
          signal <- TRUE
        }
        
        r_previous <- r_t
        t <- t + 1
      }
      
      if (!signal) {
        run_lengths[sim] <- max_time
      }
    }
    
    mean_arl1 <- mean(run_lengths)
    sd_arl1 <- sd(run_lengths)
    ci_lower <- mean_arl1 - 1.96 * sd_arl1 / sqrt(num_simulations)
    ci_upper <- mean_arl1 + 1.96 * sd_arl1 / sqrt(num_simulations)
    
    return(list(
      ARL1 = mean_arl1,
      SD = sd_arl1,
      CI_lower = ci_lower,
      CI_upper = ci_upper,
      Simulations = num_simulations
    ))
  }, error = function(e) {
    return(list(
      ARL1 = NA,
      SD = NA,
      CI_lower = NA,
      CI_upper = NA,
      Simulations = 0
    ))
  })
}

############################# MAIN EXECUTION #############################

cat("\n=== CSB-EWMA ROBUST ARL₁ CALCULATION ===\n")

# Calculate for each distribution separately (for cleaner tables)
all_results <- data.frame()

for (dist in distributions) {
  cat("\nProcessing distribution:", dist, "...\n")
  
  # Create parameter grid for this distribution
  param_list <- list()
  idx <- 1
  
  for (i in 1:nrow(optimal_combinations)) {
    for (delta in delta_vec) {
      param_list[[idx]] <- data.frame(
        lambda = optimal_combinations$lambda[i],
        L = optimal_combinations$L[i],
        combo_label = optimal_combinations$combo_label[i],
        table_label = optimal_combinations$table_label[i],
        ARL0_target = optimal_combinations$ARL0_target[i],
        Group = optimal_combinations$Group[i],
        Distribution = dist,
        delta = delta,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }
  
  param_grid <- do.call(rbind, param_list)
  
  # Calculate ARL₁ for this distribution
  dist_results <- foreach(i = 1:nrow(param_grid), .combine = rbind) %do% {
    lambda <- param_grid$lambda[i]
    L <- param_grid$L[i]
    combo_label <- param_grid$combo_label[i]
    table_label <- param_grid$table_label[i]
    distribution <- param_grid$Distribution[i]
    delta <- param_grid$delta[i]
    ARL0_target <- param_grid$ARL0_target[i]
    Group <- param_grid$Group[i]
    
    sim_result <- simulate_csb_ewma_arl1_robust(
      lambda, L, distribution, delta,
      num_simulations = 50000, max_time = 1000
    )
    
    data.frame(
      lambda = lambda,
      L = L,
      combo_label = combo_label,
      table_label = table_label,
      Distribution = distribution,
      delta = delta,
      ARL0_target = ARL0_target,
      Group = Group,
      ARL1 = sim_result$ARL1,
      SD = sim_result$SD,
      Simulations = sim_result$Simulations,
      stringsAsFactors = FALSE
    )
  }
  
  all_results <- rbind(all_results, dist_results)
}

# Add distribution names
dist_name_map <- data.frame(
  Distribution = distributions,
  Dist_Name = dist_names,
  stringsAsFactors = FALSE
)
all_results <- all_results %>% left_join(dist_name_map, by = "Distribution")

############################# CREATE EXACT FORMAT TABLES #############################

cat("\n=== CREATING EXACT FORMAT TABLES ===\n")

# Function to create table in exact format from image 1
create_format_table <- function(results_df, group_filter, dist_filter) {
  # Filter data
  filtered <- results_df %>%
    filter(Group == group_filter, Dist_Name == dist_filter) %>%
    select(table_label, delta, ARL1) %>%
    arrange(table_label, delta)
  
  # Create wide format
  table_wide <- filtered %>%
    pivot_wider(
      names_from = delta,
      values_from = ARL1,
      names_prefix = "δ="
    )
  
  # Format ARL1 values (round to integers like in the image)
  for (col in colnames(table_wide)[-1]) {
    table_wide[[col]] <- round(table_wide[[col]])
  }
  
  return(table_wide)
}

# Create tables for each distribution and group
for (group in c("ARL0_370", "ARL0_500")) {
  group_name <- ifelse(group == "ARL0_370", "370", "500")
  
  for (dist in dist_names) {
    table_data <- create_format_table(all_results, group, dist)
    
    # Save to CSV
    filename <- sprintf("Table_ARL1_%s_Target%s.csv", dist, group_name)
    write.csv(table_data, filename, row.names = FALSE)
    cat(sprintf("Saved: %s\n", filename))
    
    # Also print to console
    cat(sprintf("\n=== ARL₁ Table for %s (Target ARL₀ = %s) ===\n", dist, group_name))
    print(table_data)
    cat("\n")
  }
}

############################# CREATE HEATMAP TABLE FORMAT #############################

cat("\n=== CREATING HEATMAP TABLE FORMAT ===\n")

# Create a combined table similar to image 1
create_heatmap_table <- function(results_df, dist_filter) {
  # Calculate average ARL1 across both groups for this distribution
  summary_table <- results_df %>%
    filter(Dist_Name == dist_filter) %>%
    group_by(table_label, delta) %>%
    summarise(
      Avg_ARL1 = mean(ARL1, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    pivot_wider(
      names_from = delta,
      values_from = Avg_ARL1,
      names_prefix = "δ="
    )
  
  # Round to integers
  for (col in colnames(summary_table)[-1]) {
    summary_table[[col]] <- round(summary_table[[col]])
  }
  
  return(summary_table)
}

# Create heatmap tables for each distribution
for (dist in dist_names) {
  heatmap_table <- create_heatmap_table(all_results, dist)
  filename <- sprintf("Heatmap_Table_%s.csv", dist)
  write.csv(heatmap_table, filename, row.names = FALSE)
  cat(sprintf("Saved heatmap table: %s\n", filename))
  
  # Print first few rows
  cat(sprintf("\nHeatmap Table for %s (first 5 rows):\n", dist))
  print(head(heatmap_table, 5))
  cat("\n")
}

############################# CREATE PLOTS IN EXACT FORMAT #############################

cat("\n=== CREATING EXACT FORMAT PLOTS ===\n")

# 1. Plot similar to image 2: ARL₁ across distributions for each target ARL₀
library(ggplot2)
library(RColorBrewer)

# Calculate average ARL1 for each distribution, delta, and group
plot_data <- all_results %>%
  group_by(Dist_Name, delta, Group, ARL0_target) %>%
  summarise(
    Mean_ARL1 = mean(ARL1, na.rm = TRUE),
    SD_ARL1 = sd(ARL1, na.rm = TRUE),
    .groups = 'drop'
  )

# Create plot for each target ARL₀
for (target in c(370, 500)) {
  group_name <- paste0("ARL0_", target)
  plot_subset <- plot_data %>% filter(ARL0_target == target)
  
  p <- ggplot(plot_subset, aes(x = delta, y = Mean_ARL1, color = Dist_Name, group = Dist_Name)) +
    geom_line(size = 1.5) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = Mean_ARL1 - SD_ARL1, ymax = Mean_ARL1 + SD_ARL1), 
                  width = 0.02, alpha = 0.5) +
    labs(
      title = sprintf("CSB-EWMA ARL₁ Robustness Across Distributions (Target ARL₀ = %d)", target),
      subtitle = "Performance comparison for different data-generating processes",
      x = "Shift Magnitude (δ)",
      y = "ARL₁",
      color = "Distribution"
    ) +
    scale_color_brewer(palette = "Set1") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) +
    scale_x_continuous(breaks = delta_vec, labels = sprintf("%.2f", delta_vec))
  
  filename <- sprintf("ARL1_Robustness_Target%d.png", target)
  ggsave(filename, p, width = 10, height = 7, dpi = 300)
  cat(sprintf("Saved plot: %s\n", filename))
}

# 2. Plot similar to image 3: Average Coefficient of Variation (CV)
cv_data <- all_results %>%
  group_by(Group, ARL0_target, delta, Dist_Name) %>%
  summarise(
    Mean_ARL1 = mean(ARL1, na.rm = TRUE),
    SD_ARL1 = sd(ARL1, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  group_by(Group, ARL0_target, delta) %>%
  summarise(
    Avg_ARL1 = mean(Mean_ARL1, na.rm = TRUE),
    Avg_SD = mean(SD_ARL1, na.rm = TRUE),
    CV = Avg_SD / Avg_ARL1,
    .groups = 'drop'
  )

# Create CV plot
p_cv <- ggplot(cv_data, aes(x = delta, y = CV, color = factor(ARL0_target), group = factor(ARL0_target))) +
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  labs(
    title = "Average Coefficient of Variation (CV) of ARL₁",
    subtitle = "CV = SD / Mean (lower CV indicates more consistent performance across distributions)",
    x = "Shift Magnitude (δ)",
    y = "Avg CV",
    color = "Target ARL₀"
  ) +
  scale_color_manual(values = c("370" = "blue", "500" = "red")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 12)
  ) +
  scale_x_continuous(breaks = delta_vec, labels = sprintf("%.2f", delta_vec)) +
  scale_y_continuous(limits = c(0, max(cv_data$CV) * 1.1))

ggsave("Average_CV_ARL1.png", p_cv, width = 10, height = 7, dpi = 300)
cat("Saved plot: Average_CV_ARL1.png\n")

############################# CREATE SUMMARY REPORT #############################

cat("\n=== SUMMARY REPORT ===\n")

# Identify best performing combinations (lowest ARL1 for each delta)
best_performers <- all_results %>%
  group_by(Dist_Name, delta, Group) %>%
  filter(ARL1 == min(ARL1, na.rm = TRUE)) %>%
  select(Dist_Name, delta, Group, table_label, ARL1) %>%
  arrange(Dist_Name, delta, Group)

cat("\nBest performing (λ, L) combinations for each distribution and delta:\n")
print(best_performers)

# Calculate overall performance metrics
overall_metrics <- all_results %>%
  group_by(Dist_Name) %>%
  summarise(
    Avg_ARL1 = mean(ARL1, na.rm = TRUE),
    Min_ARL1 = min(ARL1, na.rm = TRUE),
    Max_ARL1 = max(ARL1, na.rm = TRUE),
    SD_ARL1 = sd(ARL1, na.rm = TRUE),
    CV = SD_ARL1 / Avg_ARL1,
    .groups = 'drop'
  )

cat("\nOverall performance metrics by distribution:\n")
print(overall_metrics)

# Save summary statistics
write.csv(best_performers, "Best_Performers_Summary.csv", row.names = FALSE)
write.csv(overall_metrics, "Overall_Performance_Metrics.csv", row.names = FALSE)
write.csv(all_results, "Complete_ARL1_Results.csv", row.names = FALSE)

cat("\n=== FILES CREATED ===\n")
cat("1. Table_ARL1_[Distribution]_Target[370/500].csv - Individual tables for each distribution and target\n")
cat("2. Heatmap_Table_[Distribution].csv - Combined heatmap tables\n")
cat("3. ARL1_Robustness_Target[370/500].png - Performance comparison plots\n")
cat("4. Average_CV_ARL1.png - Coefficient of variation plot\n")
cat("5. Best_Performers_Summary.csv - Best performing combinations\n")
cat("6. Overall_Performance_Metrics.csv - Overall statistics\n")
cat("7. Complete_ARL1_Results.csv - All simulation results\n")

cat("\n=== EXECUTION COMPLETE ===\n")


############################# ADDITIONAL PLOTS BASED ON COMPLETE RESULTS #############################
library(tidyverse)
library(ggplot2)
library(gridExtra)

# Read the complete results
results <- read.csv("Complete_ARL1_Results.csv")

# Create cleaner labels with proper subscript formatting
results <- results %>%
  mutate(
    combo_label_clean = sprintf("(%.3f, %.3f)", lambda, L),
    ARL0_Target = ifelse(Group == "ARL0_370", "Target ARL₀ = 370", "Target ARL₀ = 500"),
    shift_label = factor(sprintf("δ = %.2f", delta), 
                         levels = sprintf("δ = %.2f", seq(0.05, 0.5, 0.05)))
  )

# Separate data for each target
results_370 <- results %>% filter(Group == "ARL0_370")
results_500 <- results %>% filter(Group == "ARL0_500")

# ------------------------------------------------------------
# OPTION 4: CONNECTED LINES (CLEANEST VISUALIZATION)
# ------------------------------------------------------------

# For ARL₀ = 370
p_370_lines <- ggplot(results_370, aes(x = reorder(combo_label_clean, lambda), 
                                       y = ARL1,
                                       color = Dist_Name,
                                       group = Dist_Name)) +
  geom_line(size = 0.8, alpha = 0.7) +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_wrap(~ shift_label, ncol = 5, scales = "free_x") +
  scale_color_manual(
    values = c(
      "Normal" = "#1f77b4", 
      "Laplace" = "#ff7f0e", 
      "Uniform" = "#2ca02c", 
      "Exponential" = "#d62728"
    )
  ) +
  labs(
    title = expression(bold("CSB-EWMA ARL"[1]~"Performance by Distribution (Target ARL"[0]~"= 370)")),
    subtitle = "Connected lines show performance trends across (λ, L) combinations",
    x = "(λ, L) Combinations (Ordered by λ)",
    y = expression("ARL"[1]),
    color = "Distribution"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    strip.text = element_text(face = "bold", size = 9),
    axis.title = element_text(size = 11),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.8, "lines")
  ) +
  guides(color = guide_legend(nrow = 1))

ggsave("ARL1_ConnectedLines_ARL0_370.png", p_370_lines, 
       width = 16, height = 8, dpi = 300, bg = "white")

# For ARL₀ = 500
p_500_lines <- ggplot(results_500, aes(x = reorder(combo_label_clean, lambda), 
                                       y = ARL1,
                                       color = Dist_Name,
                                       group = Dist_Name)) +
  geom_line(size = 0.8, alpha = 0.7) +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_wrap(~ shift_label, ncol = 5, scales = "free_x") +
  scale_color_manual(
    values = c(
      "Normal" = "#1f77b4", 
      "Laplace" = "#ff7f0e", 
      "Uniform" = "#2ca02c", 
      "Exponential" = "#d62728"
    )
  ) +
  labs(
    title = expression(bold("CSB-EWMA ARL"[1]~"Performance by Distribution (Target ARL"[0]~"= 500)")),
    subtitle = "Connected lines show performance trends across (λ, L) combinations",
    x = "(λ, L) Combinations (Ordered by λ)",
    y = expression("ARL"[1]),
    color = "Distribution"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    strip.text = element_text(face = "bold", size = 9),
    axis.title = element_text(size = 11),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.8, "lines")
  ) +
  guides(color = guide_legend(nrow = 1))

ggsave("ARL1_ConnectedLines_ARL0_500.png", p_500_lines, 
       width = 16, height = 8, dpi = 300, bg = "white")

# ------------------------------------------------------------
# OPTION 2: HEATMAP STYLE WITH VALUES 
# ------------------------------------------------------------

# Function to create heatmap with values for a specific ARL0 target
create_heatmap_with_values <- function(data, target_value, target_label) {
  # Calculate average ARL1 for each combination and shift
  heatmap_data <- data %>%
    group_by(combo_label_clean, delta, Dist_Name) %>%
    summarise(
      ARL1_avg = round(mean(ARL1, na.rm = TRUE)),
      .groups = 'drop'
    ) %>%
    mutate(
      shift_label = factor(sprintf("δ = %.2f", delta), 
                           levels = sprintf("δ = %.2f", seq(0.05, 0.5, 0.05)))
    ) %>%
    arrange(desc(lambda))  # Sort by lambda for ordering
  
  # Reorder the combinations by lambda value (extracted from the string)
  heatmap_data <- heatmap_data %>%
    mutate(
      lambda_val = as.numeric(str_extract(combo_label_clean, "(?<=\\()\\d+\\.\\d+")),
      combo_label_clean = reorder(combo_label_clean, lambda_val)
    )
  
  p <- ggplot(heatmap_data, aes(x = shift_label, 
                                y = reorder(combo_label_clean, desc(lambda_val)), 
                                fill = ARL1_avg)) +
    geom_tile(color = "white", size = 0.8, width = 0.95, height = 0.95) +
    geom_text(aes(label = ARL1_avg), color = "black", size = 3.5, fontface = "bold") +
    facet_wrap(~ Dist_Name, ncol = 4) +
    scale_fill_gradient2(
      low = "#006837", 
      mid = "#ffffbf", 
      high = "#a50026", 
      midpoint = median(heatmap_data$ARL1_avg),
      name = expression("ARL"[1]),
      breaks = seq(1, max(heatmap_data$ARL1_avg), by = 5)
    ) +
    labs(
      title = sprintf("CSB-EWMA ARL1 Heatmap: All (λ, L) Combinations (Target ARL0 = %s)", target_label),
      subtitle = "Darker green = lower ARL1 (better detection), Red = higher ARL1",
      x = "Shift Magnitude (δ)",
      y = "(λ, L) Combinations"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      strip.text = element_text(face = "bold", size = 10),
      axis.title = element_text(size = 11),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      panel.spacing = unit(1.2, "lines"),
      panel.grid = element_blank()
    ) +
    guides(fill = guide_colorbar(barwidth = 1, barheight = 15))
  
  return(p)
}

# Create heatmaps
p_370_heatmap <- create_heatmap_with_values(results_370, 370, "370")
ggsave("ARL1_Heatmap_with_Values_ARL0_370.png", p_370_heatmap, 
       width = 14, height = 8, dpi = 300, bg = "white")

p_500_heatmap <- create_heatmap_with_values(results_500, 500, "500")
ggsave("ARL1_Heatmap_with_Values_ARL0_500.png", p_500_heatmap, 
       width = 14, height = 8, dpi = 300, bg = "white")

# ------------------------------------------------------------
# OPTION 3: COEFFICIENT OF VARIATION PLOT
# ------------------------------------------------------------

# Calculate CV for each shift magnitude and ARL0 target
cv_data <- results %>%
  group_by(Group, delta, Dist_Name) %>%
  summarise(
    Mean_ARL1 = mean(ARL1, na.rm = TRUE),
    SD_ARL1 = sd(ARL1, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  group_by(Group, delta) %>%
  summarise(
    Avg_CV = mean(SD_ARL1 / Mean_ARL1, na.rm = TRUE),
    SD_CV = sd(SD_ARL1 / Mean_ARL1, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    ARL0_Label = ifelse(Group == "ARL0_370", "Target ARL₀ = 370", "Target ARL₀ = 500"),
    delta_label = factor(sprintf("%.2f", delta), 
                         levels = sprintf("%.2f", seq(0.05, 0.5, 0.05)))
  )

# Create CV plot exactly like the attachment
p_cv <- ggplot(cv_data, aes(x = delta, y = Avg_CV, 
                            color = ARL0_Label, 
                            group = ARL0_Label,
                            shape = ARL0_Label)) +
  geom_line(size = 1.2) +
  geom_point(size = 4) +
  geom_text(aes(label = sprintf("%.3f", Avg_CV)), 
            vjust = -0.8, size = 3.5, fontface = "bold") +
  scale_color_manual(
    values = c("Target ARL₀ = 370" = "#1f77b4", 
               "Target ARL₀ = 500" = "#ff7f0e")
  ) +
  scale_shape_manual(
    values = c("Target ARL₀ = 370" = 16, 
               "Target ARL₀ = 500" = 17)
  ) +
  labs(
    title = expression(bold("Average Coefficient of Variation (CV) of ARL"[1])),
    subtitle = "CV = SD / Mean (lower CV indicates more consistent performance across distributions)",
    x = expression("Shift Magnitude (δ)"),
    y = "Average CV",
    color = "Target ARL₀",
    shape = "Target ARL₀"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11, margin = margin(b = 15)),
    axis.title = element_text(size = 11),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  scale_x_continuous(
    breaks = seq(0.05, 0.5, by = 0.05),
    labels = c("0.05", "0.10", "0.15", "0.20", "0.25", 
               "0.30", "0.35", "0.40", "0.45", "0.50")
  ) +
  scale_y_continuous(
    limits = c(0, max(cv_data$Avg_CV) * 1.1),
    breaks = seq(0, 1.5, by = 0.3)
  ) +
  guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))

ggsave("Average_CV_ARL1.png", p_cv, 
       width = 10, height = 7, dpi = 300, bg = "white")

# Create a combined version with annotations
p_cv_annotated <- p_cv +
  annotate("text", x = 0.2, y = 0.9, 
           label = "Lower CV = More Consistent\nPerformance Across Distributions",
           size = 4, fontface = "italic", color = "darkred") +
  annotate("segment", x = 0.25, xend = 0.35, y = 0.85, yend = 0.6,
           arrow = arrow(type = "closed", length = unit(0.2, "cm")),
           color = "darkred", size = 0.8)

ggsave("Average_CV_ARL1_Annotated.png", p_cv_annotated, 
       width = 10, height = 7, dpi = 300, bg = "white")

# ------------------------------------------------------------
# CREATE A SUMMARY TABLE
# ------------------------------------------------------------

# Create a summary table of CV values
cv_summary_table <- cv_data %>%
  select(ARL0_Label, delta, Avg_CV) %>%
  pivot_wider(
    names_from = delta,
    values_from = Avg_CV,
    names_prefix = "δ = "
  ) %>%
  arrange(desc(ARL0_Label))

write.csv(cv_summary_table, "CV_Summary_Table.csv", row.names = FALSE)

cat("All plots created and saved:\n")
cat("1. ARL1_ConnectedLines_ARL0_370.png - Connected lines for ARL₀ = 370\n")
cat("2. ARL1_ConnectedLines_ARL0_500.png - Connected lines for ARL₀ = 500\n")
cat("3. ARL1_Heatmap_with_Values_ARL0_370.png - Heatmap with values for ARL₀ = 370\n")
cat("4. ARL1_Heatmap_with_Values_ARL0_500.png - Heatmap with values for ARL₀ = 500\n")
cat("5. Average_CV_ARL1.png - Coefficient of variation plot\n")
cat("6. Average_CV_ARL1_Annotated.png - CV plot with annotations\n")
cat("7. CV_Summary_Table.csv - Table of CV values for dissertation\n")

# Display the CV table
cat("\nCoefficient of Variation Summary Table:\n")

print(cv_summary_table)

