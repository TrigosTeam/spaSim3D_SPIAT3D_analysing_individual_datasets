# Code to plot SD1 analysis, showing only cell pair A/B (ignoring cell pairs A/A, B/A, B/B).
# Set seed ----
set.seed(678999821)
# Libraries ----
library(cowplot)
library(ggplot2)
library(S4Vectors)
library(stringr)
library(dplyr)

# Read data and set up ----
setwd("~/R/SD1_data")

metric_df_list <- readRDS("SD1_metric_df_list.RDS")

parameters_df <- readRDS("SD1_parameters_df.RDS")


# Functions ----

# Get plot_df by sampling a random slice
get_plot_df_for_random_slice <- function(metric_df_list, 
                                         metrics,
                                         parameters_df) {
  
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change slice column to numeric for safety
    metric_df$slice <- as.numeric(metric_df$slice)
    
    if (metric %in% c("EBSAC", "EBP_AUC")) {
      # For EBSAC and EBP_AUC, assume pair is the same as cell_types for consistency
      metric_df$pair <- gsub(',', '/', metric_df$cell_types)
    }
    else if (metric %in% c("ANE_AUC")) {
      # For ANE_AUC, assume pair is the same as target_cell_type for consistency (as target is of form A,B already)
      metric_df$pair <- gsub(',', '/', metric_df$target)
    }
    else {
      # Add reference-target column
      metric_df$pair <- paste(metric_df$reference, metric_df$target, sep = "/")
    }
    
    metric_df$metric <- metric  
    
    # Select 3D values, when slice == 0, then remove from data frame
    metric_values3D <- metric_df[[metric]][metric_df[["slice"]] == 0]
    metric_df <- metric_df[metric_df[["slice"]] != 0, ]
    
    # Select a random 2D slice
    metric_df <- metric_df %>%
      group_by(simulation, pair) %>%
      slice_sample(n = 1) %>%   # pick one random slice per simulation
      ungroup() %>%
      select(simulation, pair, metric, value2D = .data[[metric]])
    
    # Combine 3D and 2D values
    metric_df[["value3D"]] <- metric_values3D
    
    # Ensure simulation column is numeric
    metric_df$simulation <- as.numeric(metric_df$simulation)
    
    # Merge metric_df and parameters_df
    parameters_df$simulation <- seq(nrow(parameters_df))
    metric_df <- metric_df %>% left_join(parameters_df[, c("simulation", "arrangement", "shape")], by = "simulation")
    
    # Get structure column
    metric_df$structure <- paste(metric_df$arrangement, metric_df$shape, sep = "_")
    
    metric_df$metric <- metric  
    
    # Subset for A/B pair
    metric_df <- metric_df[metric_df$pair == "A/B", ]
    
    
    # Add metric_df to plot_df
    plot_df <- rbind(plot_df, metric_df[ , c("value3D", "value2D", "pair", "metric", "structure")])
    
  }
  return(plot_df)
}

# Plot 2D vs 3D for each metric and pair for a random slice
plot_2D_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot <- function(plot_df,
                                                                               metrics) {

  # Update metrics (remove _AUC from specific metrics)
  metrics <- sub("_AUC$", "", metrics)
  plot_df$metric <- sub("_AUC$", "", plot_df$metric)
  
  # For nicer tick labels
  sci_clean_threshold <- function(x) {
    
    sapply(x, function(v) {
      if (is.na(v)) return("")
      
      if (abs(v) < 1000) {
        return(as.character(v))
      }
      
      # scientific notation with 1 decimal place
      s <- formatC(v, format = "e", digits = 1)   # e.g. "1.5e+03"
      
      # remove "+" in exponent
      s <- gsub("e\\+", "e", s)
      
      # split mantissa and exponent
      parts <- strsplit(s, "e")[[1]]
      mant <- parts[1]
      exp  <- parts[2]
      
      # remove trailing .0 (so 1.0e3 → 1e3)
      mant <- sub("\\.0$", "", mant)
      
      # remove leading zeros in exponent
      exp <- sub("^0+", "", exp)
      
      paste0(mant, "e", exp)
    })
  }
  
  # Factor metrics
  plot_df$metric <- factor(plot_df$metric, metrics)
  
  # Get correlation
  corr_df <- plot_df %>% 
    group_by(metric) %>% 
    summarise(
      corr = cor(value3D, value2D, method = "spearman", use = "complete.obs"), 
      .groups = "drop")
  
  write.csv(corr_df, "~/R/values_from_figures/random_slice_corr_df.csv")
  
  fig <- ggplot(plot_df, aes(x = value3D, y = value2D)) +
    geom_point(alpha = 0.25, color = "#0062c5", size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                color = "#bb0036", linewidth = 1) +
    # labs(
    #   title = "Scatterplots showing 2D vs 3D, for each metric, a random slice and cell pair A/B",
    #   x = "3D value",
    #   y = "2D value"
    # ) +
    labs(
      x = "",
      y = ""
    ) +

    facet_wrap(~ interaction(metric), scales = "free", ncol = length(metrics)) +
    
    scale_x_continuous(n.breaks = 3, labels = sci_clean_threshold) + 
    scale_y_continuous(n.breaks = 3, labels = sci_clean_threshold) +
    
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      
      # Axis tick labels
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      
      # Axis titles
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      
      # Plot title
      plot.title = element_text(size = 18),
      
      # Facet strip titles
      strip.text = element_text(size = 15)
    ) +
    geom_text(
      data = corr_df,
      aes(x = -Inf, y = Inf, label = paste0("r: ", round(corr, 3))),
      inherit.aes = FALSE,
      hjust = -0.1, 
      vjust = 1.4,
      size = 5,   # r-value font size
      color = "black"
    )
  
  
  return(fig)
}

# Plot percentage difference vs 3D for each metric and pair for a random slice
plot_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot <- function(plot_df,
                                                                                                  metrics) {
  
  # Update metrics (remove _AUC from specific metrics)
  metrics <- sub("_AUC$", "", metrics)
  plot_df$metric <- sub("_AUC$", "", plot_df$metric)
  
  # For nicer tick labels
  sci_clean_threshold <- function(x) {

    x[!(x %in% range(x, na.rm = T))] <- NA

    sapply(x, function(v) {
      if (is.na(v)) return("")
      
      if (abs(v) < 1000) {
        return(as.character(v))
      }
      
      # scientific notation with 1 decimal place
      s <- formatC(v, format = "e", digits = 1)   # e.g. "1.5e+03"
      
      # remove "+" in exponent
      s <- gsub("e\\+", "e", s)
      
      # split mantissa and exponent
      parts <- strsplit(s, "e")[[1]]
      mant <- parts[1]
      exp  <- parts[2]
      
      # remove trailing .0 (so 1.0e3 → 1e3)
      mant <- sub("\\.0$", "", mant)
      
      # remove leading zeros in exponent
      exp <- sub("^0+", "", exp)
      
      paste0(mant, "e", exp)
    })
  }
  
  # Factor for metric
  plot_df$metric <- factor(plot_df$metric, metrics)
  
  # Get error column
  plot_df$error <- (plot_df$value2D - plot_df$value3D) / plot_df$value3D * 100
  
  fig <- ggplot(plot_df, aes(x = value3D, y = error)) +
    geom_point(alpha = 0.25, color = "#0062c5", size = 1) +
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) +
    labs(
      x = "3D value",
      y = "Percentage difference (%)"
    ) +
    facet_wrap(~ interaction(metric), scales = "free", ncol = length(metrics)) +
    scale_x_continuous(n.breaks = 5, labels = sci_clean_threshold) +
    scale_y_continuous(limits = c(-1000, 1000), n.breaks = 5, labels = sci_clean_threshold) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x  = element_text(size = 11),
      axis.text.y  = element_text(size = 11),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      strip.text   = element_text(size = 11)
    )
  
  
  return(fig)
}

# Plot percentage difference vs metric for each pair for a random slice
plot_percentage_difference_vs_metric_by_pair_A_B_for_random_slice_box_plot <- function(plot_df,
                                                                                       metrics) {
  
  # Update metrics (remove _AUC from specific metrics)
  metrics <- sub("_AUC$", "", metrics)
  plot_df$metric <- sub("_AUC$", "", plot_df$metric)
  
  # Get error column
  plot_df$error <- (plot_df$value2D - plot_df$value3D) / plot_df$value3D * 100
  
  # Get correlation
  corr_df <- plot_df %>% 
    group_by(metric) %>% 
    summarise(
      corr = cor(value3D, value2D, method = "spearman", use = "complete.obs"), 
      .groups = "drop")
  
  # map corr ∈ [-1,1] → error ∈ [y_min, y_max]
  y_min <- -1000
  y_max <- 1000
  corr_to_error <- function(c) scales::rescale(c, to = c(y_min, y_max), from = c(0, 1))
  
  # inverse transform for the secondary axis
  error_to_corr <- function(e) scales::rescale(e, to = c(0, 1), from = c(y_min, y_max))
  
  corr_df$y_trans <- corr_to_error(corr_df$corr)
  
  # For nicer tick labels
  sci_clean_threshold <- function(x) {
    sapply(x, function(v) {
      if (is.na(v)) return("")
      
      if (abs(v) <= 1000) {
        return(as.character(v))
      }
      
      # scientific notation with 1 decimal place
      s <- formatC(v, format = "e", digits = 1)   # e.g. "1.5e+03"
      
      # remove "+" in exponent
      s <- gsub("e\\+", "e", s)
      
      # split mantissa and exponent
      parts <- strsplit(s, "e")[[1]]
      mant <- parts[1]
      exp  <- parts[2]
      
      # remove trailing .0 (so 1.0e3 → 1e3)
      mant <- sub("\\.0$", "", mant)
      
      # remove leading zeros in exponent
      exp <- sub("^0+", "", exp)
      
      paste0(mant, "e", exp)
    })
  }
  
  # Factor metrics
  plot_df$metric <- factor(plot_df$metric, metrics)
  
  # Get summary of error to save
  summary_df <- plot_df %>%
    group_by(metric) %>%
    summarise(
      min     = min(error, na.rm = TRUE),
      max     = max(error, na.rm = TRUE),
      range   = max - min,
      median  = median(error, na.rm = TRUE),
      average = mean(error, na.rm = TRUE)
    )
  write.csv(summary_df, paste("~/R/values_from_figures/random_slice_summary_df.csv", sep = ""))
  
  # Get extreme values to save
  extreme_values_df <- plot_df %>%
    filter(error > 1000 | error < -1000) %>%
    arrange(metric, error)
  write.csv(extreme_values_df, paste("~/R/values_from_figures/random_slice_extreme_values_df.csv", sep = ""))
  
  fig <- ggplot(plot_df, aes(x = metric, y = error)) +
    geom_boxplot(fill = "lightgray") +
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) +
    
    # median labels
    stat_summary(
      fun = function(z) median(z, na.rm = TRUE),
      geom = "text",
      aes(label = after_stat(sprintf("%.1f", y))),
      vjust = -0.5,
      size = 5,
      color = "black"
    ) +
    
    # ⭐ correlation stars
    geom_point(data = corr_df, 
               aes(x = metric, y = y_trans), 
               shape = 8, # star 
               size = 5, 
               color = "#0062c5") +
    
    labs(
      # title = "Box plots showing percentage difference between 2D and 3D metrics, for a random slice and cell pair A/B",
      x = "Metric",
      y = "Percentage difference (%)"
    ) +
    
    scale_y_continuous(n.breaks = 5, labels = sci_clean_threshold,
                       sec.axis = sec_axis(~ error_to_corr(.), name = "Spearman Correlation")) +
    
    coord_cartesian(ylim = c(y_min, y_max)) +
    
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      
      # Font sizes
      plot.title      = element_text(size = 16),
      axis.title.x    = element_text(size = 16),
      axis.title.y    = element_text(size = 16),
      axis.text.x     = element_text(size = 16, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y     = element_text(size = 16)
    )
  
  
  
  return(fig)
}



# Plot 2D vs 3D for each metric and pair for averaged slice
plot_2D_vs_3D_by_metric_and_pair_A_B_for_averaged_slice_scatter_plot <- function(metric_df_list,
                                                                                 metrics) {
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change slice column to numeric for safety
    metric_df$slice <- as.numeric(metric_df$slice)
    
    # Add pair column
    if (metric %in% c("EBSAC", "EBP_AUC")) {
      # For EBSAC and EBP_AUC, assume pair is the same as cell_types for consistency
      metric_df$pair <- gsub(',', '/', metric_df$cell_types)
    }
    else if (metric %in% c("ANE_AUC")) {
      # For ANE_AUC, assume pair is the same as target_cell_type for consistency (as target is of form A,B already)
      metric_df$pair <- gsub(',', '/', metric_df$target)
    }
    else {
      # Add reference-target column
      metric_df$pair <- paste(metric_df$reference, metric_df$target, sep = "/")
    }
    
    # Subset for A/B pair
    metric_df <- metric_df[metric_df$pair == "A/B", ]
    
    # Modify the metric_df so that there is a column of 3D values, and a column of 2D values
    colnames(metric_df)[colnames(metric_df) == metric] <- "metric_value"
    metric_df$sim_pair <- paste(metric_df$simulation, metric_df$pair, sep = "_")
    
    metric_df3D <- metric_df[metric_df[["slice"]] == 0, ]
    colnames(metric_df3D)[colnames(metric_df3D) == "metric_value"] <- "value3D"
    
    # Take the average of all slices for each simulation and pair to get 2D values
    metric_df2D <- metric_df[metric_df[["slice"]] != 0, ] %>%
      group_by(sim_pair) %>%
      summarise(value2D = mean(metric_value, na.rm = T))
    
    # Merge metric_df3D and metric_df2D
    metric_df <- merge(metric_df3D, metric_df2D, by = "sim_pair")
    
    metric_df$metric <- metric  
    
    # Add metric_df to plot_df
    plot_df <- rbind(plot_df, metric_df[ , c("value3D", "value2D", "pair", "metric")])
  }
  
  # Update metrics (remove _AUC from specific metrics)
  metrics <- sub("_AUC$", "", metrics)
  plot_df$metric <- sub("_AUC$", "", plot_df$metric)
  
  # For nicer tick labels
  sci_clean_threshold <- function(x) {
    sapply(x, function(v) {
      if (is.na(v)) return("")
      
      if (abs(v) < 1000) {
        return(as.character(v))
      }
      
      # scientific notation with 1 decimal place
      s <- formatC(v, format = "e", digits = 1)   # e.g. "1.5e+03"
      
      # remove "+" in exponent
      s <- gsub("e\\+", "e", s)
      
      # split mantissa and exponent
      parts <- strsplit(s, "e")[[1]]
      mant <- parts[1]
      exp  <- parts[2]
      
      # remove trailing .0 (so 1.0e3 → 1e3)
      mant <- sub("\\.0$", "", mant)
      
      # remove leading zeros in exponent
      exp <- sub("^0+", "", exp)
      
      paste0(mant, "e", exp)
    })
  }
  
  # Factor for metric
  plot_df$metric <- factor(plot_df$metric, metrics)
  
  # Get correlation
  corr_df <- plot_df %>% 
    group_by(metric) %>% 
    summarise(
      corr = cor(value3D, value2D, method = "spearman", use = "complete.obs"), 
      .groups = "drop")
  
  write.csv(corr_df, "~/R/values_from_figures/averaged_slice_corr_df.csv")
  
  fig <- ggplot(plot_df, aes(x = value3D, y = value2D)) +
    geom_point(alpha = 0.25, color = "#0062c5", size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "#bb0036", linewidth = 1) + # Dotted line of the equation y = x
    labs(title = "Scatter plots showing 2D vs 3D with spearman correlation, for each metric, for averaged slices and cell pair A/B",
         x = "3D value",
         y = "2D value") +
    facet_wrap(~ interaction(metric), scales = "free", ncol = length(metrics)) +
    
    scale_x_continuous(n.breaks = 3, labels = sci_clean_threshold) + 
    scale_y_continuous(n.breaks = 3, labels = sci_clean_threshold) +
    
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      
      # Axis tick labels
      axis.text.x = element_text(size = 4),
      axis.text.y = element_text(size = 4),
      
      # Axis titles
      axis.title.x = element_text(size = 4),
      axis.title.y = element_text(size = 4),
      
      # Plot title
      plot.title = element_text(size = 6),
      
      # Facet strip titles
      strip.text = element_text(size = 4)
    ) +
    
    # Add r-value text
    geom_text(
      data = corr_df,
      aes(x = -Inf, y = Inf, label = paste0("r: ", round(corr, 3))),
      inherit.aes = FALSE,
      hjust = -0.1, 
      vjust = 1.4,
      size = 4,   # r-value font size
      color = "black"
    )
  
  return(fig)
}

# Plot percentage difference vs 3D for each metric and pair for averaged slice
plot_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_averaged_slice_scatter_plot <- function(metric_df_list,
                                                                                                    metrics) {
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change slice column to numeric for safety
    metric_df$slice <- as.numeric(metric_df$slice)
    
    # Add pair column
    if (metric %in% c("EBSAC", "EBP_AUC")) {
      # For EBSAC and EBP_AUC, assume pair is the same as cell_types for consistency
      metric_df$pair <- gsub(',', '/', metric_df$cell_types)
    }
    else if (metric %in% c("ANE_AUC")) {
      # For ANE_AUC, assume pair is the same as target_cell_type for consistency (as target is of form A,B already)
      metric_df$pair <- gsub(',', '/', metric_df$target)
    }
    else {
      # Add reference-target column
      metric_df$pair <- paste(metric_df$reference, metric_df$target, sep = "/")
    }
    
    # Subset for A/B pair
    metric_df <- metric_df[metric_df$pair == "A/B", ]
    
    # Modify the metric_df so that there is a column of 3D values, and a column of 2D values
    colnames(metric_df)[colnames(metric_df) == metric] <- "metric_value"
    metric_df$sim_pair <- paste(metric_df$simulation, metric_df$pair, sep = "_")
    
    metric_df3D <- metric_df[metric_df[["slice"]] == 0, ]
    colnames(metric_df3D)[colnames(metric_df3D) == "metric_value"] <- "value3D"
    
    # Take the average of all slices for each simulation and pair to get 2D values
    metric_df2D <- metric_df[metric_df[["slice"]] != 0, ] %>%
      group_by(sim_pair) %>%
      summarise(value2D = mean(metric_value, na.rm = T))
    
    # Merge metric_df3D and metric_df2D
    metric_df <- merge(metric_df3D, metric_df2D, by = "sim_pair")
    
    metric_df$metric <- metric  
    
    # Add metric_df to plot_df
    plot_df <- rbind(plot_df, metric_df[ , c("value3D", "value2D", "pair", "metric")])
    
  }
  
  # Update metrics (remove _AUC from specific metrics)
  metrics <- sub("_AUC$", "", metrics)
  plot_df$metric <- sub("_AUC$", "", plot_df$metric)
  
  # For nicer tick labels
  sci_clean_threshold <- function(x) {
    
    x[!(x %in% range(x, na.rm = T))] <- NA
    
    sapply(x, function(v) {
      if (is.na(v)) return("")
      
      if (abs(v) < 1000) {
        return(as.character(v))
      }
      
      # scientific notation with 1 decimal place
      s <- formatC(v, format = "e", digits = 1)   # e.g. "1.5e+03"
      
      # remove "+" in exponent
      s <- gsub("e\\+", "e", s)
      
      # split mantissa and exponent
      parts <- strsplit(s, "e")[[1]]
      mant <- parts[1]
      exp  <- parts[2]
      
      # remove trailing .0 (so 1.0e3 → 1e3)
      mant <- sub("\\.0$", "", mant)
      
      # remove leading zeros in exponent
      exp <- sub("^0+", "", exp)
      
      paste0(mant, "e", exp)
    })
  }
  
  # Factor for metric
  plot_df$metric <- factor(plot_df$metric, metrics)
  
  # Get error column
  plot_df$error <- (plot_df$value2D - plot_df$value3D) / plot_df$value3D * 100
  
  fig <- ggplot(plot_df, aes(x = value3D, y = error)) +
    geom_point(alpha = 0.25, color = "#0062c5", size = 1) +
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    labs(
      # title = "Scatter plots showing percentage difference between 2D and 3D metrics vs 3D, for each metric, for averaged slices and cell pair A/B",
      x = "3D value",
      y = "Percentage difference (%)"
    ) +
    facet_wrap(~ interaction(metric), scales = "free", ncol = length(metrics)) +  
    scale_x_continuous(n.breaks = 3, labels = sci_clean_threshold) + 
    scale_y_continuous(limits = c(-1000, 1000), n.breaks = 3, labels = sci_clean_threshold) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x  = element_text(size = 11),
      axis.text.y  = element_text(size = 11),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      strip.text   = element_text(size = 11)
    )
  
  return(fig)
}

# Plot percentage difference vs metric for each pair for averaged slice
plot_percentage_difference_vs_metric_by_pair_A_B_for_averaged_slice_box_plot <- function(metric_df_list,
                                                                                         metrics) {
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change slice column to numeric for safety
    metric_df$slice <- as.numeric(metric_df$slice)
    
    # Add pair column
    if (metric %in% c("EBSAC", "EBP_AUC")) {
      # For EBSAC and EBP_AUC, assume pair is the same as cell_types for consistency
      metric_df$pair <- gsub(',', '/', metric_df$cell_types)
    }
    else if (metric %in% c("ANE_AUC")) {
      # For ANE_AUC, assume pair is the same as target_cell_type for consistency (as target is of form A,B already)
      metric_df$pair <- gsub(',', '/', metric_df$target)
    }
    else {
      # Add reference-target column
      metric_df$pair <- paste(metric_df$reference, metric_df$target, sep = "/")
    }
    
    # Subset for A/B pair
    metric_df <- metric_df[metric_df$pair == "A/B", ]
    
    # Modify the metric_df so that there is a column of 3D values, and a column of 2D values
    colnames(metric_df)[colnames(metric_df) == metric] <- "metric_value"
    metric_df$sim_pair <- paste(metric_df$simulation, metric_df$pair, sep = "_")
    
    metric_df3D <- metric_df[metric_df[["slice"]] == 0, ]
    colnames(metric_df3D)[colnames(metric_df3D) == "metric_value"] <- "value3D"
    
    # Take the average of all slices for each simulation and pair to get 2D values
    metric_df2D <- metric_df[metric_df[["slice"]] != 0, ] %>%
      group_by(sim_pair) %>%
      summarise(value2D = mean(metric_value, na.rm = T))
    
    # Merge metric_df3D and metric_df2D
    metric_df <- merge(metric_df3D, metric_df2D, by = "sim_pair")
    
    metric_df$metric <- metric  
    
    # Add metric_df to plot_df
    plot_df <- rbind(plot_df, metric_df[ , c("value3D", "value2D", "pair", "metric")])
    
  }
  
  # Update metrics (remove _AUC from specific metrics)
  metrics <- sub("_AUC$", "", metrics)
  plot_df$metric <- sub("_AUC$", "", plot_df$metric)
  
  # Get error column
  plot_df$error <- (plot_df$value2D - plot_df$value3D) / plot_df$value3D * 100
  
  # Get correlation
  corr_df <- plot_df %>% 
    group_by(metric) %>% 
    summarise(
      corr = cor(value3D, value2D, method = "spearman", use = "complete.obs"), 
      .groups = "drop")
  
  # map corr ∈ [-1,1] → error ∈ [-1000,1000]
  y_min <- -1000
  y_max <- 1000
  corr_to_error <- function(c) scales::rescale(c, to = c(y_min, y_max), from = c(0, 1))
  
  # inverse transform for the secondary axis
  error_to_corr <- function(e) scales::rescale(e, to = c(0, 1), from = c(y_min, y_max))
  
  corr_df$y_trans <- corr_to_error(corr_df$corr)
  
  # For nicer tick labels
  sci_clean_threshold <- function(x) {
    sapply(x, function(v) {
      if (is.na(v)) return("")
      
      if (abs(v) <= 1000) {
        return(as.character(v))
      }
      
      # scientific notation with 1 decimal place
      s <- formatC(v, format = "e", digits = 1)   # e.g. "1.5e+03"
      
      # remove "+" in exponent
      s <- gsub("e\\+", "e", s)
      
      # split mantissa and exponent
      parts <- strsplit(s, "e")[[1]]
      mant <- parts[1]
      exp  <- parts[2]
      
      # remove trailing .0 (so 1.0e3 → 1e3)
      mant <- sub("\\.0$", "", mant)
      
      # remove leading zeros in exponent
      exp <- sub("^0+", "", exp)
      
      paste0(mant, "e", exp)
    })
  }
  
  # Factor for metric
  plot_df$metric <- factor(plot_df$metric, metrics)
  
  # Get summary of error to save
  summary_df <- plot_df %>%
    group_by(metric) %>%
    summarise(
      min     = min(error, na.rm = TRUE),
      max     = max(error, na.rm = TRUE),
      range   = max - min,
      median  = median(error, na.rm = TRUE),
      average = mean(error, na.rm = TRUE)
    )
  write.csv(summary_df, paste("~/R/values_from_figures/average_slice_summary_df.csv", sep = ""))
  
  # Get extreme values to save
  extreme_values_df <- plot_df %>%
    filter(error > 1000 | error < -1000) %>%
    arrange(metric, error)
  write.csv(extreme_values_df, paste("~/R/values_from_figures/average_slice_extreme_values_df.csv", sep = ""))
  
  fig <- ggplot(plot_df, aes(x = metric, y = error)) +
    geom_boxplot(fill = "lightgray") +
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    
    # # Median labels
    # stat_summary(
    #   fun = function(z) median(z, na.rm = TRUE),
    #   geom = "text",
    #   aes(label = after_stat(sprintf("%.1f", y))),
    #   vjust = -0.5,
    #   size = 5,
    #   color = "black"
    # ) +
    
    # ⭐ correlation stars
    geom_point(data = corr_df, 
               aes(x = metric, y = y_trans), 
               shape = 8, # star 
               size = 5, 
               color = "#0062c5") +
    
    labs(
      # title = "Box plots showing percentage difference between 2D and 3D metrics with spearman correlation, for averaged slices for cell pair A/B",
       x = "Metric",
       y = "Percentage difference (%)"
      ) +
    
    scale_y_continuous(n.breaks = 5, labels = sci_clean_threshold,
                       sec.axis = sec_axis(~ error_to_corr(.), name = "Spearman Correlation")) +
    
    coord_cartesian(ylim = c(y_min, y_max)) +
    
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      
      # Font sizes
      plot.title      = element_text(size = 16),
      axis.title.x    = element_text(size = 16),
      axis.title.y    = element_text(size = 16),
      axis.text.x     = element_text(size = 16, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y     = element_text(size = 16)
    )
  
  return(fig)
}



# Plot 2D vs 3D for each metric and pair for 3 slices
plot_2D_vs_3D_by_metric_and_pair_A_B_for_three_slices_scatter_plot <- function(metric_df_list,
                                                                               metrics) {
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change slice column to numeric for safety
    metric_df$slice <- as.numeric(metric_df$slice)
    
    # Add pair column
    if (metric %in% c("EBSAC", "EBP_AUC")) {
      # For EBSAC and EBP_AUC, assume pair is the same as cell_types for consistency
      metric_df$pair <- gsub(',', '/', metric_df$cell_types)
    }
    else if (metric %in% c("ANE_AUC")) {
      # For ANE_AUC, assume pair is the same as target_cell_type for consistency (as target is of form A,B already)
      metric_df$pair <- gsub(',', '/', metric_df$target)
    }
    else {
      # Add reference-target column
      metric_df$pair <- paste(metric_df$reference, metric_df$target, sep = "/")
    }
    
    # Subset for A/B pair
    metric_df <- metric_df[metric_df$pair == "A/B", ]
    
    # Modify the metric_df so that there is a column of 3D values, and a column of 2D values
    colnames(metric_df)[colnames(metric_df) == metric] <- "metric_value"
    metric_df$sim_pair <- paste(metric_df$simulation, metric_df$pair, sep = "_")
    
    metric_df3D <- metric_df[metric_df[["slice"]] == 0, ]
    colnames(metric_df3D)[colnames(metric_df3D) == "metric_value"] <- "value3D"
    
    # Take 3 slices for each simulation and pair to get 2D values
    # Slice indices 7, 10, 13 correspond with 145, 175, 205 z coord
    metric_df2D <- metric_df[metric_df[["slice"]] %in% c(7, 10, 13), ] 
    
    # Merge metric_df3D and metric_df2D
    metric_df <- metric_df2D %>% left_join(metric_df3D[, c("sim_pair", "value3D")], by = "sim_pair")
    colnames(metric_df)[colnames(metric_df) == "metric_value"] <- "value2D"
    
    metric_df$metric <- metric  
    
    # Add metric_df to plot_df
    plot_df <- rbind(plot_df, metric_df[ , c("value3D", "value2D", "pair", "metric", "slice")])
    
  }
  
  # Update metrics (remove _AUC from specific metrics)
  metrics <- sub("_AUC$", "", metrics)
  plot_df$metric <- sub("_AUC$", "", plot_df$metric)
  
  # For nicer tick labels
  sci_clean_threshold <- function(x) {
    sapply(x, function(v) {
      if (is.na(v)) return("")
      
      if (abs(v) < 1000) {
        return(as.character(v))
      }
      
      # scientific notation with 1 decimal place
      s <- formatC(v, format = "e", digits = 1)   # e.g. "1.5e+03"
      
      # remove "+" in exponent
      s <- gsub("e\\+", "e", s)
      
      # split mantissa and exponent
      parts <- strsplit(s, "e")[[1]]
      mant <- parts[1]
      exp  <- parts[2]
      
      # remove trailing .0 (so 1.0e3 → 1e3)
      mant <- sub("\\.0$", "", mant)
      
      # remove leading zeros in exponent
      exp <- sub("^0+", "", exp)
      
      paste0(mant, "e", exp)
    })
  }
  
  # Factor for metric
  plot_df$metric <- factor(plot_df$metric, metrics)
  
  # Change slice column back to character for plotting
  plot_df$slice <- as.character(metric_df$slice)
  
  plot_df$slice <- factor(plot_df$slice, levels = c('7', '10', '13'))
  
  fig <- ggplot(plot_df, aes(x = value3D, y = value2D, color = slice)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "#bb0036", linewidth = 1) + # Dotted line of the equation y = x
    labs(title = "Scatter plots showing 2D vs 3D with spearman correlation, for each metric, for 3 slices and for cell pair A/B",
         x = "3D value",
         y = "2D value") +
    facet_wrap(~ interaction(metric), scales = "free", ncol = length(metrics)) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(size = 4),  # make x-axis text smaller
      axis.text.y = element_text(size = 4)   # make y-axis text smaller
    ) +
    scale_x_continuous(n.breaks = 3, labels = sci_clean_threshold) +
    scale_y_continuous(n.breaks = 3, labels = sci_clean_threshold) +
    scale_color_manual(
      values = c(
        "7" = "#9437a8",
        "10" = "#007128",
        "13" = "#b8db50"
        
      )
    )
  
  return(fig)
}

# Plot percentage difference vs 3D for each metric and pair for 3 slices
plot_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_three_slices_scatter_plot <- function(metric_df_list,
                                                                                                  metrics) {
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change slice column to numeric for safety
    metric_df$slice <- as.numeric(metric_df$slice)
    
    # Add pair column
    if (metric %in% c("EBSAC", "EBP_AUC")) {
      # For EBSAC and EBP_AUC, assume pair is the same as cell_types for consistency
      metric_df$pair <- gsub(',', '/', metric_df$cell_types)
    }
    else if (metric %in% c("ANE_AUC")) {
      # For ANE_AUC, assume pair is the same as target_cell_type for consistency (as target is of form A,B already)
      metric_df$pair <- gsub(',', '/', metric_df$target)
    }
    else {
      # Add reference-target column
      metric_df$pair <- paste(metric_df$reference, metric_df$target, sep = "/")
    }
    
    # Subset for A/B pair
    metric_df <- metric_df[metric_df$pair == "A/B", ]
    
    # Modify the metric_df so that there is a column of 3D values, and a column of 2D values
    colnames(metric_df)[colnames(metric_df) == metric] <- "metric_value"
    metric_df$sim_pair <- paste(metric_df$simulation, metric_df$pair, sep = "_")
    
    metric_df3D <- metric_df[metric_df[["slice"]] == 0, ]
    colnames(metric_df3D)[colnames(metric_df3D) == "metric_value"] <- "value3D"
    
    # Take 3 slices for each simulation and pair to get 2D values
    # Slice indices 7, 10, 13 correspond with 145, 175, 205 z coord
    metric_df2D <- metric_df[metric_df[["slice"]] %in% c(7, 10, 13), ] 
    
    # Merge metric_df3D and metric_df2D
    metric_df <- metric_df2D %>% left_join(metric_df3D[, c("sim_pair", "value3D")], by = "sim_pair")
    colnames(metric_df)[colnames(metric_df) == "metric_value"] <- "value2D"
    
    metric_df$metric <- metric  
    
    # Add metric_df to plot_df
    plot_df <- rbind(plot_df, metric_df[ , c("value3D", "value2D", "pair", "metric", "slice")])
    
  }
  
  # Update metrics (remove _AUC from specific metrics)
  metrics <- sub("_AUC$", "", metrics)
  plot_df$metric <- sub("_AUC$", "", plot_df$metric)
  
  # Change slice column back to character for plotting
  plot_df$slice <- as.character(metric_df$slice)
  
  plot_df$slice <- factor(plot_df$slice, levels = c('7', '10', '13'))
  
  # For nicer tick labels
  sci_clean_threshold <- function(x) {
    
    x[!(x %in% range(x, na.rm = T))] <- NA
    
    sapply(x, function(v) {
      if (is.na(v)) return("")
      
      if (abs(v) < 1000) {
        return(as.character(v))
      }
      
      # scientific notation with 1 decimal place
      s <- formatC(v, format = "e", digits = 1)   # e.g. "1.5e+03"
      
      # remove "+" in exponent
      s <- gsub("e\\+", "e", s)
      
      # split mantissa and exponent
      parts <- strsplit(s, "e")[[1]]
      mant <- parts[1]
      exp  <- parts[2]
      
      # remove trailing .0 (so 1.0e3 → 1e3)
      mant <- sub("\\.0$", "", mant)
      
      # remove leading zeros in exponent
      exp <- sub("^0+", "", exp)
      
      paste0(mant, "e", exp)
    })
  }
  
  # Factor for metric
  plot_df$metric <- factor(plot_df$metric, metrics)
  
  # Get error column
  plot_df$error <- (plot_df$value2D - plot_df$value3D) / plot_df$value3D * 100
  
  fig <- ggplot(plot_df, aes(x = value3D, y = error, color = slice)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    labs(
      # title = "Scatter plots showing percentage difference between 2D and 3D metrics vs 3D, for each metric, for three slices and cell pair A/B",
      x = "3D value",
      y = "Percentage difference (%)"
    ) +
    facet_wrap(~ interaction(metric), scales = "free", ncol = length(metrics)) +  
    scale_x_continuous(n.breaks = 3, labels = sci_clean_threshold) +
    scale_y_continuous(limits = c(-1000, 1000), n.breaks = 3, labels = sci_clean_threshold) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x  = element_text(size = 11),
      axis.text.y  = element_text(size = 11),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      strip.text   = element_text(size = 11),
      legend.position = "none"
    ) +
    scale_color_manual(
      values = c(
        "7" = "#9437a8",
        "10" = "#007128",
        "13" = "#b8db50"
        
      )
    )
  
  return(fig)
}

# Plot percentage difference vs metric for each pair for 3 slices
plot_percentage_difference_vs_metric_by_pair_A_B_for_three_slice_box_plot <- function(metric_df_list,
                                                                                      metrics) {
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change slice column to numeric for safety
    metric_df$slice <- as.numeric(metric_df$slice)
    
    # Add pair column
    if (metric %in% c("EBSAC", "EBP_AUC")) {
      # For EBSAC and EBP_AUC, assume pair is the same as cell_types for consistency
      metric_df$pair <- gsub(',', '/', metric_df$cell_types)
    }
    else if (metric %in% c("ANE_AUC")) {
      # For ANE_AUC, assume pair is the same as target_cell_type for consistency (as target is of form A,B already)
      metric_df$pair <- gsub(',', '/', metric_df$target)
    }
    else {
      # Add reference-target column
      metric_df$pair <- paste(metric_df$reference, metric_df$target, sep = "/")
    }
    
    # Subset for A/B pair
    metric_df <- metric_df[metric_df$pair == "A/B", ]
    
    # Modify the metric_df so that there is a column of 3D values, and a column of 2D values
    colnames(metric_df)[colnames(metric_df) == metric] <- "metric_value"
    metric_df$sim_pair <- paste(metric_df$simulation, metric_df$pair, sep = "_")
    
    metric_df3D <- metric_df[metric_df[["slice"]] == 0, ]
    colnames(metric_df3D)[colnames(metric_df3D) == "metric_value"] <- "value3D"
    
    # Take 3 slices for each simulation and pair to get 2D values
    # Slice indices 7, 10, 13 correspond with 145, 175, 205 z coord
    metric_df2D <- metric_df[metric_df[["slice"]] %in% c(7, 10, 13), ] 
    
    # Merge metric_df3D and metric_df2D
    metric_df <- metric_df2D %>% left_join(metric_df3D[, c("sim_pair", "value3D")], by = "sim_pair")
    colnames(metric_df)[colnames(metric_df) == "metric_value"] <- "value2D"
    
    metric_df$metric <- metric  
    
    # Add metric_df to plot_df
    plot_df <- rbind(plot_df, metric_df[ , c("value3D", "value2D", "pair", "metric", "slice")])
    
  }
  
  # Update metrics (remove _AUC from specific metrics)
  metrics <- sub("_AUC$", "", metrics)
  plot_df$metric <- sub("_AUC$", "", plot_df$metric)
  
  # Change slice column back to character for plotting
  plot_df$slice <- as.character(metric_df$slice)
  
  plot_df$slice <- factor(plot_df$slice, levels = c('7', '10', '13'))
  
  # Get error column
  plot_df$error <- (plot_df$value2D - plot_df$value3D) / plot_df$value3D * 100
  
  # Factor for metric
  plot_df$metric <- factor(plot_df$metric, metrics)
  
  # Get correlation
  corr_df <- plot_df %>% 
    group_by(metric, slice) %>% 
    summarise(
      corr = cor(value3D, value2D, method = "spearman", use = "complete.obs"), 
      .groups = "drop")
  
  write.csv(corr_df, "~/R/values_from_figures/three_slice_corr_df.csv")
  
  # map corr ∈ [-1,1] → error ∈ [-1000,1000]
  y_min <- -1000
  y_max <- 1000
  corr_to_error <- function(c) scales::rescale(c, to = c(y_min, y_max), from = c(0, 1))
  
  # inverse transform for the secondary axis
  error_to_corr <- function(e) scales::rescale(e, to = c(0, 1), from = c(y_min, y_max))
  
  corr_df$y_trans <- corr_to_error(corr_df$corr)
  
  # For nicer tick labels
  sci_clean_threshold <- function(x) {
    sapply(x, function(v) {
      if (is.na(v)) return("")
      
      if (abs(v) <= 1000) {
        return(as.character(v))
      }
      
      # scientific notation with 1 decimal place
      s <- formatC(v, format = "e", digits = 1)   # e.g. "1.5e+03"
      
      # remove "+" in exponent
      s <- gsub("e\\+", "e", s)
      
      # split mantissa and exponent
      parts <- strsplit(s, "e")[[1]]
      mant <- parts[1]
      exp  <- parts[2]
      
      # remove trailing .0 (so 1.0e3 → 1e3)
      mant <- sub("\\.0$", "", mant)
      
      # remove leading zeros in exponent
      exp <- sub("^0+", "", exp)
      
      paste0(mant, "e", exp)
    })
  }
  
  # Define slice colors
  slice_cols <- c(
    "7"  = "#9437a8",
    "10" = "#007128",
    "13" = "#b8db50"
  )
  
  # Get extreme values to save
  extreme_values_df <- plot_df %>%
    filter(error > 1000 | error < -1000) %>%
    arrange(metric, error)
  write.csv(extreme_values_df, paste("~/R/values_from_figures/three_slice_extreme_values_df.csv", sep = ""))
  
  
  fig <- ggplot(plot_df, aes(x = metric, y = error, fill = slice)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = slice_cols) +
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    
    # # Median labels
    # stat_summary(
    #   fun = function(z) median(z, na.rm = TRUE),
    #   geom = "text",
    #   aes(label = after_stat(sprintf("%.1f", y)), group = slice),
    #   position = position_dodge(width = 0.8),
    #   vjust = -0.5,
    #   size = 3.5,
    #   color = "black"
    # ) +
    
    
    # ⭐ correlation stars
    geom_point(
      data = corr_df,
      aes(x = metric, y = y_trans, group = slice),
      shape = 8,
      size = 5,
      color = "#0062c5",
      position = position_dodge(width = 0.8)
    ) +
    
    
    labs(
      # title = "Box plots showing percentage difference between 2D and 3D metrics with spearman correlation, for three slices and cell pair A/B",
       x = "Metric",
       y = "Percentage difference (%)"
      ) +
    
    scale_y_continuous(n.breaks = 3, labels = sci_clean_threshold,
                       sec.axis = sec_axis(~ error_to_corr(.), name = "Spearman Correlation")) +
    
    coord_cartesian(ylim = c(y_min, y_max)) +
    
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      
      # Font sizes
      plot.title      = element_text(size = 16),
      axis.title.x    = element_text(size = 16),
      axis.title.y    = element_text(size = 16),
      axis.text.x     = element_text(size = 16, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y     = element_text(size = 16),
      
      legend.position = "none"
      
    )
  
  return(fig)
}



# Plot 2D vs 3D for each metric and pair for a random slice, annotating for tissue structure
plot_2D_vs_3D_by_metric_and_pair_A_B_for_random_slice_showing_structure_scatter_plot <- function(plot_df,
                                                                                                 metrics,
                                                                                                 parameters_df) {

  # Update metrics (remove _AUC from specific metrics)
  metrics <- sub("_AUC$", "", metrics)
  plot_df$metric <- sub("_AUC$", "", plot_df$metric)
  
  # For nicer tick labels
  sci_clean_threshold <- function(x) {
    sapply(x, function(v) {
      if (is.na(v)) return("")
      
      if (abs(v) < 1000) {
        return(as.character(v))
      }
      
      # scientific notation with 1 decimal place
      s <- formatC(v, format = "e", digits = 1)   # e.g. "1.5e+03"
      
      # remove "+" in exponent
      s <- gsub("e\\+", "e", s)
      
      # split mantissa and exponent
      parts <- strsplit(s, "e")[[1]]
      mant <- parts[1]
      exp  <- parts[2]
      
      # remove trailing .0 (so 1.0e3 → 1e3)
      mant <- sub("\\.0$", "", mant)
      
      # remove leading zeros in exponent
      exp <- sub("^0+", "", exp)
      
      paste0(mant, "e", exp)
    })
  }
  
  # Factor for metric
  plot_df$metric <- factor(plot_df$metric, metrics)
  
  fig <- ggplot(plot_df, aes(x = value3D, y = value2D, color = structure)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "#bb0036", linewidth = 1) + # Dotted line of the equation y = x
    labs(title = "Scatterplots showing 2D vs 3D, for each metric, for a random slice and cell pair A/B, showing structure",
         x = "3D value",
         y = "2D value") +
    facet_wrap(~ interaction(metric), scales = "free", ncol = length(metrics)) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(size = 4),  # make x-axis text smaller
      axis.text.y = element_text(size = 4)   # make y-axis text smaller
    ) +
    scale_x_continuous(n.breaks = 3, labels = sci_clean_threshold) +
    scale_y_continuous(n.breaks = 3, labels = sci_clean_threshold) +
    scale_color_manual(
      values = c(
        "mixed_ellipsoid" = "#007128",
        "mixed_network" = "#b8db50",
        "ringed_ellipsoid" = "#9437a8",
        "ringed_network" = "#d99dff",
        "separated_ellipsoid" = "#770026",
        "separated_network" = "#bb0036"
      )
    )
  
  return(fig)
}

# Plot percentage difference vs 3D for each metric and pair for a random slice, annotating for tissue structure
plot_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_random_slice_showing_structure_scatter_plot <- function(plot_df,
                                                                                                                    metrics,
                                                                                                                    parameters_df) {
  # Update metrics (remove _AUC from specific metrics)
  metrics <- sub("_AUC$", "", metrics)
  plot_df$metric <- sub("_AUC$", "", plot_df$metric)
  
  # For nicer tick labels
  sci_clean_threshold <- function(x) {
    
    x[!(x %in% range(x, na.rm = T))] <- NA
    
    sapply(x, function(v) {
      if (is.na(v)) return("")
      
      if (abs(v) < 1000) {
        return(as.character(v))
      }
      
      # scientific notation with 1 decimal place
      s <- formatC(v, format = "e", digits = 1)   # e.g. "1.5e+03"
      
      # remove "+" in exponent
      s <- gsub("e\\+", "e", s)
      
      # split mantissa and exponent
      parts <- strsplit(s, "e")[[1]]
      mant <- parts[1]
      exp  <- parts[2]
      
      # remove trailing .0 (so 1.0e3 → 1e3)
      mant <- sub("\\.0$", "", mant)
      
      # remove leading zeros in exponent
      exp <- sub("^0+", "", exp)
      
      paste0(mant, "e", exp)
    })
  }
  
  # Factor for metric
  plot_df$metric <- factor(plot_df$metric, metrics)
  
  # Get error column
  plot_df$error <- (plot_df$value2D - plot_df$value3D) / plot_df$value3D * 100
  
  fig <- ggplot(plot_df, aes(x = value3D, y = error, color = structure)) +
    geom_point(alpha = 1, size = 1) +
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    labs(
      # title = "Scatterplots showing percentage difference between 2D and 3D metrics vs 3D, for each metric, for a random slice and cell pair A/B, showing structure",
      x = "3D value",
      y = "Percentage difference (%)"
      ) +
    facet_wrap(~ interaction(metric), scales = "free", ncol = length(metrics)) +  
    scale_x_continuous(n.breaks = 3, labels = sci_clean_threshold) +
    scale_y_continuous(limits = c(-1000, 1000), n.breaks = 3, labels = sci_clean_threshold) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x  = element_text(size = 11),
      axis.text.y  = element_text(size = 11),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      strip.text   = element_text(size = 11),
      legend.position = "none"
    ) +
    scale_color_manual(
      values = c(
        "mixed_ellipsoid" = "#007128",
        "mixed_network" = "#b8db50",
        "ringed_ellipsoid" = "#9437a8",
        "ringed_network" = "#d99dff",
        "separated_ellipsoid" = "#770026",
        "separated_network" = "#bb0036"
      )
    )
  
  return(fig)
}

# Plot 2D vs 3D correlation vs tissue structure for each metric and pair for a random slice
plot_2D_vs_3D_correlation_vs_structure_by_metric_and_pair_A_B_for_random_slice_bar_plot <- function(plot_df,
                                                                                                    metrics,
                                                                                                    parameters_df) {
  
  # Update metrics (remove _AUC from specific metrics)
  metrics <- sub("_AUC$", "", metrics)
  plot_df$metric <- sub("_AUC$", "", plot_df$metric)
  
  # Compute spearman correlation
  corr_df <- plot_df %>%
    group_by(pair, metric, structure) %>%
    summarise(
      corr = cor(value3D, value2D, method = "spearman", use = "complete.obs"),
      .groups = "drop"
    )
  
  # Factor for metric
  plot_df$metric <- factor(plot_df$metric, metrics)
  
  fig <- ggplot(corr_df, aes(x = structure, y = corr, fill = structure)) + 
    geom_col(alpha = 0.8) + 
    labs(title = "Bar plots showing spearman correlation vs structure, for each metric, for a random slice and cell pair A/B", x = "Structure", y = "Spearman Correlation") + 
    facet_wrap(~ interaction(metric), scales = "free", ncol = length(metrics)) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal() + 
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
          axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_manual(
      values = c(
        "mixed_ellipsoid" = "#007128",
        "mixed_network" = "#b8db50",
        "ringed_ellipsoid" = "#9437a8",
        "ringed_network" = "#d99dff",
        "separated_ellipsoid" = "#770026",
        "separated_network" = "#bb0036"
      )
    )
  
  return(fig)
}

# plot correlation vs structure by metric and pair A/B for random slice
plot_correlation_vs_structure_by_metric_and_pair_A_B_for_random_slice_dot_plot <- function(plot_df,
                                                                                           metrics) {
  
  # Update metrics (remove _AUC from specific metrics)
  metrics <- sub("_AUC$", "", metrics)
  plot_df$metric <- sub("_AUC$", "", plot_df$metric)
  
  # Compute spearman correlation
  corr_df <- plot_df %>%
    group_by(pair, metric, structure) %>%
    summarise(
      corr = cor(value3D, value2D, method = "spearman", use = "complete.obs"),
      .groups = "drop"
    )
  
  fig <- ggplot(corr_df, aes(x = structure, y = corr)) + 
    geom_point() +
    labs(title = "Dot plots showing spearman correlation vs structure, for each metric, for a random slice and cell pair A/B", x = "Structure", y = "Spearman Correlation") + 
    scale_y_continuous(limits = c(-1, 1)) +
    theme_minimal() + 
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(fig)
}

# Plot percentage difference vs metric for each pair for a random slice, annotating for tissue structure
plot_percentage_difference_vs_metric_by_pair_A_B_for_random_slice_showing_structure_box_plot <- function(plot_df,
                                                                                                         metrics,
                                                                                                         parameters_df) {
  
  # Update metrics (remove _AUC from specific metrics)
  metrics <- sub("_AUC$", "", metrics)
  plot_df$metric <- sub("_AUC$", "", plot_df$metric)
  
  # Get error column
  plot_df$error <- (plot_df$value2D - plot_df$value3D) / plot_df$value3D * 100
  
  # Factor for metric
  plot_df$metric <- factor(plot_df$metric, metrics)
  
  # Compute spearman correlation
  corr_df <- plot_df %>%
    group_by(pair, metric, structure) %>%
    summarise(
      corr = cor(value3D, value2D, method = "spearman", use = "complete.obs"),
      .groups = "drop"
    )
  
  write.csv(corr_df, "~/R/values_from_figures/showing_structure_corr_df.csv")
  
  # map corr ∈ [-1,1] → error ∈ [-1000,1000]
  y_min <- -1000
  y_max <- 1000
  corr_to_error <- function(c) scales::rescale(c, to = c(y_min, y_max), from = c(0, 1))
  
  # inverse transform for the secondary axis
  error_to_corr <- function(e) scales::rescale(e, to = c(0, 1), from = c(y_min, y_max))
  
  corr_df$y_trans <- corr_to_error(corr_df$corr)
  
  # For nicer tick labels
  sci_clean_threshold <- function(x) {
    sapply(x, function(v) {
      if (is.na(v)) return("")
      
      if (abs(v) < 1000) {
        return(as.character(v))
      }
      
      # scientific notation with 1 decimal place
      s <- formatC(v, format = "e", digits = 1)   # e.g. "1.5e+03"
      
      # remove "+" in exponent
      s <- gsub("e\\+", "e", s)
      
      # split mantissa and exponent
      parts <- strsplit(s, "e")[[1]]
      mant <- parts[1]
      exp  <- parts[2]
      
      # remove trailing .0 (so 1.0e3 → 1e3)
      mant <- sub("\\.0$", "", mant)
      
      # remove leading zeros in exponent
      exp <- sub("^0+", "", exp)
      
      paste0(mant, "e", exp)
    })
  }
  
  # Structure colors
  structure_cols <- c(
    "mixed_ellipsoid" = "#007128",
    "mixed_network" = "#b8db50",
    "ringed_ellipsoid" = "#9437a8",
    "ringed_network" = "#d99dff",
    "separated_ellipsoid" = "#770026",
    "separated_network" = "#bb0036"
  )
  
  fig <- ggplot(plot_df, aes(x = metric, y = error, fill = structure)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = structure_cols) +
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    
    # # Median labels
    # stat_summary(
    #   fun = function(z) median(z, na.rm = TRUE),
    #   geom = "text",
    #   aes(label = after_stat(sprintf("%.1f", y)), group = structure),
    #   position = position_dodge(width = 0.8),
    #   vjust = -0.5,
    #   size = 5,
    #   color = "black"
    # ) +
    
    
    # ⭐ correlation stars
    geom_point(
      data = corr_df,
      aes(x = metric, y = y_trans, group = structure),
      shape = 8,
      size = 3,
      color = "#0062c5",
      position = position_dodge(width = 0.8)
    ) +
    
    
    labs(
      # title = "Box plots showing percentage difference between 2D and 3D metrics with spearman correlation, for a random slice for cell pair A/B, showing tissue structure",
       x = "Metric",
       y = "Percentage difference (%)"
      ) +
    
    scale_y_continuous(n.breaks = 3, labels = sci_clean_threshold,
                       sec.axis = sec_axis(~ error_to_corr(.), name = "Spearman Correlation")) +
    
    coord_cartesian(ylim = c(y_min, y_max)) +
    
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      
      # Font sizes
      plot.title      = element_text(size = 16),
      axis.title.x    = element_text(size = 16),
      axis.title.y    = element_text(size = 16),
      axis.text.x     = element_text(size = 16, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y     = element_text(size = 16),
      
      legend.position = "none"
    )
  
  return(fig)
}




# Running the functions ----

metrics <- c("AMD",
             "ANC_AUC", "ACIN_AUC", "ANE_AUC",
             "MS_AUC", "NMS_AUC",
             "CK_AUC", "CL_AUC", "CG_AUC",
             "COO_AUC",
             "PBP_AUC", "EBP_AUC", "PBSAC", "EBSAC")



plot_df <- get_plot_df_for_random_slice(metric_df_list,
                                        metrics,
                                        parameters_df)

setwd("~/R/plots/SD1")
fig_2D_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot <- plot_2D_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot(plot_df,
                                                                                                                                        metrics)
pdf("fig_2D_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot.pdf", width = 30, height = 4)
print(fig_2D_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot)
dev.off()


setwd("~/R/plots/SD1")
fig_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot <- plot_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot(plot_df,
                                                                                                                                                                              metrics)
pdf("fig_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot.pdf", width = 17.5, height = 2.5)
print(fig_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_random_slice_scatter_plot)
dev.off()


setwd("~/R/plots/SD1")
fig_percentage_difference_vs_metric_by_pair_A_B_for_random_slice_box_plot <- plot_percentage_difference_vs_metric_by_pair_A_B_for_random_slice_box_plot(plot_df,
                                                                                                                                                        metrics)
pdf("fig_percentage_difference_vs_metric_by_pair_A_B_for_random_slice_box_plot.pdf", width = 9, height = 6)
print(fig_percentage_difference_vs_metric_by_pair_A_B_for_random_slice_box_plot)
dev.off()






setwd("~/R/plots/SD1")
fig_2D_vs_3D_by_metric_and_pair_A_B_for_averaged_slice_scatter_plot <- plot_2D_vs_3D_by_metric_and_pair_A_B_for_averaged_slice_scatter_plot(metric_df_list,
                                                                                                                                            metrics)
pdf("fig_2D_vs_3D_by_metric_and_pair_A_B_for_averaged_slice_scatter_plot.pdf", width = 24, height = 10)
print(fig_2D_vs_3D_by_metric_and_pair_A_B_for_averaged_slice_scatter_plot)
dev.off()


setwd("~/R/plots/SD1")
fig_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_averaged_slice_scatter_plot <- plot_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_averaged_slice_scatter_plot(metric_df_list,
                                                                                                                                                                                  metrics)
pdf("fig_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_averaged_slice_scatter_plot.pdf", width = 17.5, height = 2.5)
print(fig_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_averaged_slice_scatter_plot)
dev.off()


setwd("~/R/plots/SD1")
fig_percentage_difference_vs_metric_by_pair_A_B_for_averaged_slice_box_plot <- plot_percentage_difference_vs_metric_by_pair_A_B_for_averaged_slice_box_plot(metric_df_list,
                                                                                                                                                            metrics)
pdf("fig_percentage_difference_vs_metric_by_pair_A_B_for_averaged_slice_box_plot.pdf", width = 15, height = 4)
print(fig_percentage_difference_vs_metric_by_pair_A_B_for_averaged_slice_box_plot)
dev.off()





setwd("~/R/plots/SD1")
fig_2D_vs_3D_by_metric_and_pair_A_B_for_three_slices_scatter_plot <- plot_2D_vs_3D_by_metric_and_pair_A_B_for_three_slices_scatter_plot(metric_df_list,
                                                                                                                                        metrics)
pdf("fig_2D_vs_3D_by_metric_and_pair_A_B_for_three_slices_scatter_plot.pdf", width = 24, height = 10)
print(fig_2D_vs_3D_by_metric_and_pair_A_B_for_three_slices_scatter_plot)
dev.off()


setwd("~/R/plots/SD1")
fig_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_three_slices_scatter_plot <- plot_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_three_slices_scatter_plot(metric_df_list,
                                                                                                                                                                              metrics)
pdf("fig_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_three_slices_scatter_plot.pdf", width = 17.5, height = 2.5)
print(fig_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_three_slices_scatter_plot)
dev.off()

setwd("~/R/plots/SD1")
fig_percentage_difference_vs_metric_by_pair_A_B_for_three_slice_box_plot <- plot_percentage_difference_vs_metric_by_pair_A_B_for_three_slice_box_plot(metric_df_list,
                                                                                                                                                      metrics)
pdf("fig_percentage_difference_vs_metric_by_pair_A_B_for_three_slice_box_plot.pdf", width = 15, height = 4)
print(fig_percentage_difference_vs_metric_by_pair_A_B_for_three_slice_box_plot)
dev.off()








setwd("~/R/plots/SD1")
fig_2D_vs_3D_by_metric_and_pair_A_B_for_random_slice_showing_structure_scatter_plot <- plot_2D_vs_3D_by_metric_and_pair_A_B_for_random_slice_showing_structure_scatter_plot(plot_df,
                                                                                                                                                                            metrics,
                                                                                                                                                                            parameters_df)
pdf("fig_2D_vs_3D_by_metric_and_pair_A_B_for_random_slice_showing_structure_scatter_plot.pdf", width = 24, height = 10)
print(fig_2D_vs_3D_by_metric_and_pair_A_B_for_random_slice_showing_structure_scatter_plot)
dev.off()


setwd("~/R/plots/SD1")
fig_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_random_slice_showing_structure_scatter_plot <- plot_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_random_slice_showing_structure_scatter_plot(plot_df,
                                                                                                                                                                                                                  metrics,
                                                                                                                                                                                                                  parameters_df)
pdf("fig_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_random_slice_showing_structure_scatter_plot.pdf", width = 17.5, height = 2.5)
print(fig_percentage_difference_vs_3D_by_metric_and_pair_A_B_for_random_slice_showing_structure_scatter_plot)
dev.off()


setwd("~/R/plots/SD1")
fig_2D_vs_3D_correlation_vs_structure_by_metric_and_pair_A_B_for_random_slice_bar_plot <- plot_2D_vs_3D_correlation_vs_structure_by_metric_and_pair_A_B_for_random_slice_bar_plot(plot_df,
                                                                                                                                                                                  metrics,
                                                                                                                                                                                  parameters_df)
pdf("fig_2D_vs_3D_correlation_vs_structure_by_metric_and_pair_A_B_for_random_slice_bar_plot.pdf", width = 24, height = 10)
print(fig_2D_vs_3D_correlation_vs_structure_by_metric_and_pair_A_B_for_random_slice_bar_plot)
dev.off()


setwd("~/R/plots/SD1")
fig_correlation_vs_structure_by_metric_and_pair_A_B_for_random_slice_dot_plot <- plot_correlation_vs_structure_by_metric_and_pair_A_B_for_random_slice_dot_plot(plot_df,
                                                                                                                                                                metrics)
pdf("fig_correlation_vs_structure_by_metric_and_pair_A_B_for_random_slice_dot_plot.pdf", width = 15, height = 4)
print(fig_correlation_vs_structure_by_metric_and_pair_A_B_for_random_slice_dot_plot)
dev.off()



setwd("~/R/plots/SD1")
fig_percentage_difference_vs_metric_by_pair_A_B_for_random_slice_showing_structure_box_plot <- plot_percentage_difference_vs_metric_by_pair_A_B_for_random_slice_showing_structure_box_plot(plot_df,
                                                                                                                                                                                            metrics,
                                                                                                                                                                                            parameters_df)
pdf("fig_percentage_difference_vs_metric_by_pair_A_B_for_random_slice_showing_structure_box_plot.pdf", width = 15, height = 4)
print(fig_percentage_difference_vs_metric_by_pair_A_B_for_random_slice_showing_structure_box_plot)
dev.off()


