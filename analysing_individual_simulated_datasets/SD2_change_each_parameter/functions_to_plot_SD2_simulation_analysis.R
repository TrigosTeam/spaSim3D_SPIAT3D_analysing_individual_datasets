# Code to plot SD2 analysis.
# Libraries ----
library(cowplot)
library(ggplot2)
library(S4Vectors)
library(stringr)
library(dplyr)
library(scales)
library(DescTools)

# Read data and set up ----
setwd("~/R/SD2_data")
metric_df_list <- readRDS("SD2_metric_df_list.RDS")
parameters_df <- readRDS("SD2_parameters_df.RDS")


# Functions for 4 pair metrics-----
# For metrics which use all cell pairs (A/A, A/B, B/A, B/B)
# Subset metric_df_list and parameters_df
add_pair_to_metric_df <- function(metric_df, metric) {
  
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
  return(metric_df)
}

get_parameters <- function(arrangement, shape) {
  parameters <- c()
  
  # Add other parameters specific to arrangement and shape
  if (arrangement == "mixed") {
    parameters <- c(parameters, "cluster_prop_A")
  }
  else if (arrangement == "ringed") {
    parameters <- c(parameters, "ring_width_factor")        
  }
  else if (arrangement == "separated") {
    parameters <- c(parameters, "distance")
  }
  
  if (shape == "ellipsoid") {
    parameters <- c(parameters, "E_volume")
  }
  else if (shape == "network") {
    parameters <- c(parameters, "N_width")
  }
  
  # All plots include bg_prop_A and bg_prop_B as parameters
  parameters <- c(parameters, "bg_prop_A", "bg_prop_B")
  
  return(parameters)
}

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

combine_metric_and_parameters_df <- function(metric_df, parameters_df) {
  pairs <- unique(metric_df$pair)
  metric_and_parameters_df <- data.frame()
  for (slice in unique(metric_df$slice)) {
    metric_and_parameters_df <- rbind(metric_and_parameters_df, 
                                      cbind(metric_df[metric_df$slice == slice, ], 
                                            parameters_df[rep(1:nrow(parameters_df), each = length(pairs)), ]))
  }
  return(metric_and_parameters_df)
}

plot_3D_vs_parameters_for_non_gradient_metrics_scatter_plot <- function(metric_df_list, parameters_df, metric) {
  
  fig_list <- list()
  
  # Define parameters
  arrangements <- c("mixed", "ringed", "separated")
  shapes <- c("ellipsoid", "network")
  
  # Subset the metric_df_list
  metric_df <- metric_df_list[[metric]]
  
  # Remove Inf rows
  metric_df <- metric_df[!is.infinite(metric_df[[metric]]), ]
  
  # Add 'pair' column to metric_df
  metric_df <- add_pair_to_metric_df(metric_df, metric)
  
  pairs <- unique(metric_df$pair)
  
  # Add to parameters_df
  parameters_df$distance <- 450 - parameters_df$cluster1_x_coord # 450 is the x-coordinate of cluster2
  parameters_df$E_volume <- (4 / 3) * pi * parameters_df$E_radius_x * parameters_df$E_radius_y * parameters_df$E_radius_z
  
  # Update parameters_df
  parameters_df$variable_parameter[parameters_df$variable_parameter == "E_radius_x"] <- "E_volume"
  parameters_df$variable_parameter[parameters_df$variable_parameter == "cluster1_x_coord"] <- "distance"
  
  # Combine metric_df with parameters_df
  metric_df <- combine_metric_and_parameters_df(metric_df, parameters_df)
  
  # Subset for 3D values only
  metric_df <- metric_df[metric_df$slice == 0, ]
  
  # Update metric (remove _AUC from specific metrics)
  colnames(metric_df)[colnames(metric_df) == metric] <- sub("_AUC$", "", metric)
  metric <- sub("_AUC$", "", metric)
  
  for (arrangement in arrangements) {
    for (shape in shapes) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      fig_list[[arrangement_shape]] <- list()
      
      # Subset metric_df for arrangement and shape
      metric_arrangement_shape_df <- metric_df[metric_df$arrangement == arrangement & metric_df$shape == shape, ]
      
      # Get parameters
      parameters <- get_parameters(arrangement, shape)
      
      for (pair in pairs) {
        # Subset for pair
        metric_arrangement_shape_pair_df <- metric_arrangement_shape_df[metric_arrangement_shape_df$pair == pair, ]
        fig_list[[arrangement_shape]][[pair]] <- list()
        
        for (parameter in parameters) {
          
          # Subset for parameter
          plot_df <- metric_arrangement_shape_pair_df[metric_arrangement_shape_pair_df$variable_parameter == parameter, ]
          
          fig <- ggplot(plot_df, aes_string(parameter, metric)) +
            geom_point(size = 0.5) +
            theme_minimal() +
            theme(
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              plot.title = element_text(size = 10),
              axis.text.x = element_text(size = 8),  # make x-axis text smaller
              axis.text.y = element_text(size = 8)   # make y-axis text smaller
            ) +
            scale_x_continuous(breaks = pretty_breaks(n = 3), labels = sci_clean_threshold) + 
            scale_y_continuous(breaks = pretty_breaks(n = 3), labels = sci_clean_threshold) 
          
          # Get correlation and p-value       
          p_value <- "N/A"
          correlation <- "N/A"
          if (sum(!is.na(plot_df[[metric]])) >= 2) {
            
            correlation_test <- cor.test(plot_df[[parameter]], plot_df[[metric]], method = "spearman")
            correlation <- round(correlation_test$estimate, 3)
            p_value <- correlation_test$p.value
            
            if (!is.na(p_value)) {
              # Format correlation and p-value
              formatCustomSci <- function(x) {
                x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
                alpha <- round(as.numeric(x_sci[ , 1]), 1)
                power <- as.integer(x_sci[ , 2])
                paste(alpha, power, sep = "e")
              }
              if (p_value == 0) p_value <- 2.2e-308
              if (0 < p_value && p_value < 1e-3)  {
                p_value <- formatCustomSci(p_value)
              }
              else {
                p_value  <- round(p_value, 3)
              }
            }
          }
          fig <- fig + ggtitle(paste("r:", correlation, ", p:", p_value))
          fig_list[[arrangement_shape]][[pair]][[parameter]] <- fig
        }
      }
    }
  }
  
  arrangement_shape_figs <- list()
  
  for (shape in shapes) {  
    for (arrangement in arrangements) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      
      pair_figs <- list()
      
      for (pair in pairs) {
        
        parameters <- get_parameters(arrangement, shape)
        
        pair_fig <- plot_grid(plotlist = fig_list[[arrangement_shape]][[pair]],
                              ncol = length(fig_list[[arrangement_shape]][[pair]]))
        title <- ggdraw() + draw_label(paste("pair:", pair))
        pair_fig <- plot_grid(title, pair_fig, ncol = 1, rel_heights = c(0.1, 1))
        
        pair_figs[[pair]] <- pair_fig
      }
      arrangement_shape_fig <- plot_grid(plotlist = pair_figs, ncol = 1)
      
      title <- ggdraw() + draw_label(paste(arrangement, "-", shape, sep = ""), fontface = "bold")
      arrangement_shape_fig <- plot_grid(title, arrangement_shape_fig, ncol = 1, rel_heights = c(0.02, 1))  + 
        theme(plot.margin = margin(10, 10, 10, 10),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1))  
      
      arrangement_shape_figs[[arrangement_shape]] <- arrangement_shape_fig
    }
  }
  
  fig <- plot_grid(plotlist = arrangement_shape_figs,
                   nrow = length(shapes),
                   ncol = length(arrangements))
  
  return(fig)
}

plot_3D_and_2D_vs_parameters_for_non_gradient_metrics_scatter_plot <- function(metric_df_list, parameters_df, metric) {
  
  fig_list <- list()
  
  # Define parameters
  arrangements <- c("mixed", "ringed", "separated")
  shapes <- c("ellipsoid", "network")
  
  # Subset the metric_df_list
  metric_df <- metric_df_list[[metric]]
  
  # Remove Inf rows
  metric_df <- metric_df[!is.infinite(metric_df[[metric]]), ]
  
  # Add 'pair' column to metric_df
  metric_df <- add_pair_to_metric_df(metric_df, metric)
  
  pairs <- unique(metric_df$pair)
  
  # Make 'slice' column categorical
  metric_df$slice <- as.character(metric_df$slice)
  
  # Add to parameters_df
  parameters_df$distance <- 450 - parameters_df$cluster1_x_coord # 450 is the x-coordinate of cluster2
  parameters_df$E_volume <- (4 / 3) * pi * parameters_df$E_radius_x * parameters_df$E_radius_y * parameters_df$E_radius_z
  
  # Update parameters_df
  parameters_df$variable_parameter[parameters_df$variable_parameter == "E_radius_x"] <- "E_volume"
  parameters_df$variable_parameter[parameters_df$variable_parameter == "cluster1_x_coord"] <- "distance"
  
  # Combine metric_df with parameters_df
  metric_df <- combine_metric_and_parameters_df(metric_df, parameters_df)
  
  # Update metric (remove _AUC from specific metrics)
  colnames(metric_df)[colnames(metric_df) == metric] <- sub("_AUC$", "", metric)
  metric <- sub("_AUC$", "", metric)
  
  for (arrangement in arrangements) {
    for (shape in shapes) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      fig_list[[arrangement_shape]] <- list()
      
      # Subset metric_df for arrangement and shape
      metric_arrangement_shape_df <- metric_df[metric_df$arrangement == arrangement & metric_df$shape == shape, ]
      
      # Get parameters
      parameters <- get_parameters(arrangement, shape)
      
      for (pair in pairs) {
        # Subset for pair
        metric_arrangement_shape_pair_df <- metric_arrangement_shape_df[metric_arrangement_shape_df$pair == pair, ]
        fig_list[[arrangement_shape]][[pair]] <- list()
        
        for (parameter in parameters) {
          # Subset for parameter
          plot_df <- metric_arrangement_shape_pair_df[metric_arrangement_shape_pair_df$variable_parameter == parameter, ]
          
          fig <- ggplot(plot_df, aes_string(parameter, metric)) +
            
            # 1) First draw all *non‑zero* slices
            geom_point(
              data = subset(plot_df, slice != "0"),
              aes(color = slice),
              size = 0.5,
              alpha = 0.5
            ) +
            
            # 2) Then draw slice 0 on top
            geom_point(
              data = subset(plot_df, slice == "0"),
              color = "black",
              size = 0.5,
              alpha = 0.5
            ) +
            
            theme_minimal() +
            theme(
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              plot.title = element_text(size = 10),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              legend.position = "none"
            ) +
            scale_x_continuous(breaks = pretty_breaks(n = 3), labels = sci_clean_threshold) +
            scale_y_continuous(breaks = pretty_breaks(n = 3), labels = sci_clean_threshold) +
            scale_color_manual(
              values = c(
                "1" = "#9437a8",
                "2" = "#007128",
                "3" = "#b8db50"
                # no need to include "0" here since we draw it manually
              )
            )
          
          
          fig_list[[arrangement_shape]][[pair]][[parameter]] <- fig
        }
      }
    }
  }
  
  arrangement_shape_figs <- list()
  
  for (shape in shapes) {  
    for (arrangement in arrangements) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      
      pair_figs <- list()
      
      for (pair in pairs) {
        
        pair_fig <- plot_grid(plotlist = fig_list[[arrangement_shape]][[pair]],
                              ncol = length(fig_list[[arrangement_shape]][[pair]]))
        title <- ggdraw() + draw_label(paste("pair:", pair))
        pair_fig <- plot_grid(title, pair_fig, ncol = 1, rel_heights = c(0.1, 1))
        
        pair_figs[[pair]] <- pair_fig
      }
      arrangement_shape_fig <- plot_grid(plotlist = pair_figs, ncol = 1)
      
      title <- ggdraw() + draw_label(paste(arrangement, "-", shape, sep = ""), fontface = "bold")
      arrangement_shape_fig <- plot_grid(title, arrangement_shape_fig, ncol = 1, rel_heights = c(0.02, 1))  + 
        theme(plot.margin = margin(10, 10, 10, 10),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1))  
      
      arrangement_shape_figs[[arrangement_shape]] <- arrangement_shape_fig
    }
  }
  
  fig <- plot_grid(plotlist = arrangement_shape_figs,
                   nrow = length(shapes),
                   ncol = length(arrangements))
  
  return(fig) 
}

plot_percentage_difference_vs_parameters_for_non_gradient_metrics_scatter_plot <- function(metric_df_list, parameters_df, metric) {
  
  fig_list <- list()
  
  # Define parameters
  arrangements <- c("mixed", "ringed", "separated")
  shapes <- c("ellipsoid", "network")
  
  # Subset the metric_df_list
  metric_df <- metric_df_list[[metric]]
  
  # Remove Inf rows
  metric_df <- metric_df[!is.infinite(metric_df[[metric]]), ]
  
  # Add 'pair' column to metric_df
  metric_df <- add_pair_to_metric_df(metric_df, metric)
  
  pairs <- unique(metric_df$pair)
  
  # Make 'slice' column categorical
  metric_df$slice <- as.character(metric_df$slice)
  
  # Add to parameters_df
  parameters_df$distance <- 450 - parameters_df$cluster1_x_coord # 450 is the x-coordinate of cluster2
  parameters_df$E_volume <- (4 / 3) * pi * parameters_df$E_radius_x * parameters_df$E_radius_y * parameters_df$E_radius_z
  
  # Update parameters_df
  parameters_df$variable_parameter[parameters_df$variable_parameter == "E_radius_x"] <- "E_volume"
  parameters_df$variable_parameter[parameters_df$variable_parameter == "cluster1_x_coord"] <- "distance"
  
  # Combine metric_df with parameters_df
  metric_df <- combine_metric_and_parameters_df(metric_df, parameters_df)
  
  # Add error column to metric_df
  metric_df$error <- ((metric_df[[metric]] -  rep(metric_df[[metric]][metric_df$slice == "0"], 4)) / rep(metric_df[[metric]][metric_df$slice == "0"], 4)) * 100
  
  # Update metric (remove _AUC from specific metrics)
  colnames(metric_df)[colnames(metric_df) == metric] <- sub("_AUC$", "", metric)
  metric <- sub("_AUC$", "", metric)
  
  for (arrangement in arrangements) {
    for (shape in shapes) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      fig_list[[arrangement_shape]] <- list()
      
      # Subset metric_df for arrangement and shape
      metric_arrangement_shape_df <- metric_df[metric_df$arrangement == arrangement & metric_df$shape == shape, ]
      
      # Get parameters
      parameters <- get_parameters(arrangement, shape)
      
      for (pair in pairs) {
        # Subset for pair
        metric_arrangement_shape_pair_df <- metric_arrangement_shape_df[metric_arrangement_shape_df$pair == pair, ]
        fig_list[[arrangement_shape]][[pair]] <- list()
        
        for (parameter in parameters) {
          # Subset for parameter
          plot_df <- metric_arrangement_shape_pair_df[metric_arrangement_shape_pair_df$variable_parameter == parameter, ]
          
          # Remove 3D values
          plot_df <- plot_df[plot_df$slice != "0", ]
          
          fig <- ggplot(plot_df, aes_string(parameter, "error", color = "slice")) +
            geom_point(size = 0.5) +
            geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
            theme_minimal() +
            labs(y = paste(metric, "percentage diff (%)")) +
            theme(
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              plot.title = element_text(size = 10),
              axis.text.x = element_text(size = 8),  # make x-axis text smaller
              axis.text.y = element_text(size = 8),   # make y-axis text smaller
              legend.position = "none"
            ) +
            scale_x_continuous(breaks = pretty_breaks(n = 3), labels = sci_clean_threshold) + 
            scale_y_continuous(breaks = pretty_breaks(n = 3), labels = sci_clean_threshold) +
            scale_color_manual(
              values = c(
                "1" = "#9437a8",
                "2" = "#007128",
                "3" = "#b8db50"
              )
            )
          
          fig_list[[arrangement_shape]][[pair]][[parameter]] <- fig
        }
      }
    }
  }
  
  arrangement_shape_figs <- list()
  
  for (shape in shapes) {  
    for (arrangement in arrangements) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      
      pair_figs <- list()
      
      for (pair in pairs) {
        
        pair_fig <- plot_grid(plotlist = fig_list[[arrangement_shape]][[pair]],
                              ncol = length(fig_list[[arrangement_shape]][[pair]]))
        title <- ggdraw() + draw_label(paste("pair:", pair))
        pair_fig <- plot_grid(title, pair_fig, ncol = 1, rel_heights = c(0.1, 1))
        
        pair_figs[[pair]] <- pair_fig
      }
      arrangement_shape_fig <- plot_grid(plotlist = pair_figs, ncol = 1)
      
      title <- ggdraw() + draw_label(paste(arrangement, "-", shape, sep = ""), fontface = "bold")
      arrangement_shape_fig <- plot_grid(title, arrangement_shape_fig, ncol = 1, rel_heights = c(0.02, 1))  + 
        theme(plot.margin = margin(10, 10, 10, 10),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1))  
      
      arrangement_shape_figs[[arrangement_shape]] <- arrangement_shape_fig
    }
  }
  
  fig <- plot_grid(plotlist = arrangement_shape_figs,
                   nrow = length(shapes),
                   ncol = length(arrangements))
  
  return(fig) 
}

plot_2D_vs_slice_for_non_gradient_metrics_violin_plot <- function(metric_df_list, parameters_df, metric) {
  
  fig_list <- list()
  
  # Define parameters
  arrangements <- c("mixed", "ringed", "separated")
  shapes <- c("ellipsoid", "network")
  
  # Subset the metric_df_list
  metric_df <- metric_df_list[[metric]]
  
  # Remove Inf rows
  metric_df <- metric_df[!is.infinite(metric_df[[metric]]), ]
  
  # Add 'pair' column to metric_df
  metric_df <- add_pair_to_metric_df(metric_df, metric)
  
  pairs <- unique(metric_df$pair)
  
  # Make 'slice' column categorical
  metric_df$slice <- as.character(metric_df$slice)
  
  # Add to parameters_df
  parameters_df$distance <- 450 - parameters_df$cluster1_x_coord # 450 is the x-coordinate of cluster2
  parameters_df$E_volume <- (4 / 3) * pi * parameters_df$E_radius_x * parameters_df$E_radius_y * parameters_df$E_radius_z
  
  # Update parameters_df
  parameters_df$variable_parameter[parameters_df$variable_parameter == "E_radius_x"] <- "E_volume"
  parameters_df$variable_parameter[parameters_df$variable_parameter == "cluster1_x_coord"] <- "distance"
  
  # Combine metric_df with parameters_df
  metric_df <- combine_metric_and_parameters_df(metric_df, parameters_df)
  
  # Update metric (remove _AUC from specific metrics)
  colnames(metric_df)[colnames(metric_df) == metric] <- sub("_AUC$", "", metric)
  metric <- sub("_AUC$", "", metric)
  
  for (arrangement in arrangements) {
    for (shape in shapes) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      fig_list[[arrangement_shape]] <- list()
      
      # Subset metric_df for arrangement and shape
      metric_arrangement_shape_df <- metric_df[metric_df$arrangement == arrangement & metric_df$shape == shape, ]
      
      # Get parameters
      parameters <- get_parameters(arrangement, shape)
      
      for (pair in pairs) {
        # Subset for pair
        metric_arrangement_shape_pair_df <- metric_arrangement_shape_df[metric_arrangement_shape_df$pair == pair, ]
        fig_list[[arrangement_shape]][[pair]] <- list()
        
        for (parameter in parameters) {
          # Subset for parameter
          plot_df <- metric_arrangement_shape_pair_df[metric_arrangement_shape_pair_df$variable_parameter == parameter, ]
          
          # Remove 3D values
          plot_df <- plot_df[plot_df$slice != "0", ]
          
          # Violin plot
          fig <- ggplot(plot_df, aes_string(x = "slice", y = metric, fill = "slice")) +
            geom_violin() +
            theme_minimal() +
            theme(
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              plot.title = element_text(size = 10),
              axis.text.x = element_text(size = 8),  # make x-axis text smaller
              axis.text.y = element_text(size = 8),   # make y-axis text smaller
              legend.position = "none"
            ) +
            scale_y_continuous(breaks = pretty_breaks(n = 3), labels = sci_clean_threshold) +
            scale_fill_manual(
              values = c(
                "1" = "#9437a8",
                "2" = "#007128",
                "3" = "#b8db50"
              )
            )
          
          # Get JT p-value if possible
          p_value <- "N/A"
          if (sum(!is.na(plot_df[[metric]][plot_df$slice == "1"])) > 0 &&
              sum(!is.na(plot_df[[metric]][plot_df$slice == "2"])) > 0 &&
              sum(!is.na(plot_df[[metric]][plot_df$slice == "3"])) > 0) {
            
            # Factor slice column
            plot_df$slice <- factor(plot_df$slice, c("1", "2", "3"), ordered = TRUE)
            
            # Get JT value
            slice1_mean_value <- mean(plot_df[[metric]][plot_df$slice == "1"], na.rm = TRUE)
            slice3_mean_value <- mean(plot_df[[metric]][plot_df$slice == "3"], na.rm = TRUE)
            
            if (slice3_mean_value > slice1_mean_value) {
              trend_direction <- "increasing"
            }
            else {
              trend_direction <- "decreasing"
            }
            jt_test <- JonckheereTerpstraTest(plot_df[[metric]], plot_df$slice, alternative = trend_direction)
            
            p_value <- jt_test$p.value
            
            if (!is.na(p_value)) {
              # Format correlation and p-value
              formatCustomSci <- function(x) {
                x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
                alpha <- round(as.numeric(x_sci[ , 1]), 1)
                power <- as.integer(x_sci[ , 2])
                paste(alpha, power, sep = "e")
              }
              if (p_value == 0) p_value <- 2.2e-308
              if (0 < p_value && p_value < 1e-3)  {
                p_value <- formatCustomSci(p_value)
              }
              else {
                p_value  <- round(p_value, 3)
              }
            }
          }
          fig <- fig + ggtitle(paste("p:", p_value)) 
          fig_list[[arrangement_shape]][[pair]][[parameter]] <- fig
        }
      }
    }
  }
  
  arrangement_shape_figs <- list()
  
  for (shape in shapes) {  
    for (arrangement in arrangements) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      
      pair_figs <- list()
      
      for (pair in pairs) {
        
        # If it is the final pair, add parameters to bottom
        if (pair == pairs[length(pairs)]) {
          parameters <- get_parameters(arrangement, shape) 
          for (parameter in parameters) {
            subtitle <- ggdraw() + draw_label(parameter)
            fig_list[[arrangement_shape]][[pair]][[parameter]] <- plot_grid(fig_list[[arrangement_shape]][[pair]][[parameter]], 
                                                                            subtitle, ncol = 1, rel_heights = c(1, 0.1))
          }
        }
        
        pair_fig <- plot_grid(plotlist = fig_list[[arrangement_shape]][[pair]],
                              ncol = length(fig_list[[arrangement_shape]][[pair]]))
        title <- ggdraw() + draw_label(paste("pair:", pair))
        pair_fig <- plot_grid(title, pair_fig, ncol = 1, rel_heights = c(0.1, 1))
        pair_figs[[pair]] <- pair_fig 
      }
      
      arrangement_shape_fig <- plot_grid(plotlist = pair_figs, ncol = 1)
      
      title <- ggdraw() + draw_label(paste(arrangement, "-", shape, sep = ""), fontface = "bold")
      arrangement_shape_fig <- plot_grid(title, arrangement_shape_fig, ncol = 1, rel_heights = c(0.02, 1))  + 
        theme(plot.margin = margin(10, 10, 10, 10),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1))  
      
      arrangement_shape_figs[[arrangement_shape]] <- arrangement_shape_fig
    }
  }
  
  fig <- plot_grid(plotlist = arrangement_shape_figs,
                   nrow = length(shapes),
                   ncol = length(arrangements))
  
  return(fig) 
}

plot_3D_and_2D_vs_slice_for_non_gradient_metrics_violin_plot <- function(metric_df_list, parameters_df, metric) {
  
  fig_list <- list()
  
  # Define parameters
  arrangements <- c("mixed", "ringed", "separated")
  shapes <- c("ellipsoid", "network")
  
  # Subset the metric_df_list
  metric_df <- metric_df_list[[metric]]
  
  # Remove Inf rows
  metric_df <- metric_df[!is.infinite(metric_df[[metric]]), ]
  
  # Add 'pair' column to metric_df
  metric_df <- add_pair_to_metric_df(metric_df, metric)
  
  pairs <- unique(metric_df$pair)
  
  # Make 'slice' column categorical
  metric_df$slice <- as.character(metric_df$slice)
  
  # Add to parameters_df
  parameters_df$distance <- 450 - parameters_df$cluster1_x_coord # 450 is the x-coordinate of cluster2
  parameters_df$E_volume <- (4 / 3) * pi * parameters_df$E_radius_x * parameters_df$E_radius_y * parameters_df$E_radius_z
  
  # Update parameters_df
  parameters_df$variable_parameter[parameters_df$variable_parameter == "E_radius_x"] <- "E_volume"
  parameters_df$variable_parameter[parameters_df$variable_parameter == "cluster1_x_coord"] <- "distance"
  
  # Combine metric_df with parameters_df
  metric_df <- combine_metric_and_parameters_df(metric_df, parameters_df)
  
  # Update metric (remove _AUC from specific metrics)
  colnames(metric_df)[colnames(metric_df) == metric] <- sub("_AUC$", "", metric)
  metric <- sub("_AUC$", "", metric)
  
  for (arrangement in arrangements) {
    for (shape in shapes) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      fig_list[[arrangement_shape]] <- list()
      
      # Subset metric_df for arrangement and shape
      metric_arrangement_shape_df <- metric_df[metric_df$arrangement == arrangement & metric_df$shape == shape, ]
      
      # Get parameters
      parameters <- get_parameters(arrangement, shape)
      
      for (pair in pairs) {
        # Subset for pair
        metric_arrangement_shape_pair_df <- metric_arrangement_shape_df[metric_arrangement_shape_df$pair == pair, ]
        fig_list[[arrangement_shape]][[pair]] <- list()
        
        for (parameter in parameters) {
          # Subset for parameter
          plot_df <- metric_arrangement_shape_pair_df[metric_arrangement_shape_pair_df$variable_parameter == parameter, ]
          
          # Change slice index for 3D value to "3D"
          plot_df$slice[plot_df$slice == "0"] <- "3D"
          
          # Violin plot
          fig <- ggplot(plot_df, aes_string(x = "slice", y = metric, fill = "slice")) +
            geom_violin() +
            theme_minimal() +
            theme(
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              plot.title = element_text(size = 10),
              axis.text.x = element_text(size = 8),  # make x-axis text smaller
              axis.text.y = element_text(size = 8),   # make y-axis text smaller
              legend.position = "none"
            ) +
            scale_y_continuous(breaks = pretty_breaks(n = 3), labels = sci_clean_threshold) +
            scale_fill_manual(
              values = c(
                "1" = "#9437a8",
                "2" = "#007128",
                "3" = "#b8db50",
                "3D" = "white"
              )
            )
          
          # Get JT p-value if possible
          p_value <- "N/A"
          if (sum(!is.na(plot_df[[metric]][plot_df$slice == "1"])) > 0 &&
              sum(!is.na(plot_df[[metric]][plot_df$slice == "2"])) > 0 &&
              sum(!is.na(plot_df[[metric]][plot_df$slice == "3"])) > 0) {
            
            # Create separate df with only 2D values
            df2D <- plot_df[plot_df$slice != "3D", ]
            
            # Factor slice column
            df2D$slice <- factor(df2D$slice, c("1", "2", "3"), ordered = TRUE)
            
            # Get JT value
            slice1_mean_value <- mean(df2D[[metric]][df2D$slice == "1"], na.rm = TRUE)
            slice3_mean_value <- mean(df2D[[metric]][df2D$slice == "3"], na.rm = TRUE)
            
            if (slice3_mean_value > slice1_mean_value) {
              trend_direction <- "increasing"
            }
            else {
              trend_direction <- "decreasing"
            }
            jt_test <- JonckheereTerpstraTest(df2D[[metric]], df2D$slice, alternative = trend_direction)
            
            p_value <- jt_test$p.value
            
            if (!is.na(p_value)) {
              # Format correlation and p-value
              formatCustomSci <- function(x) {
                x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
                alpha <- round(as.numeric(x_sci[ , 1]), 1)
                power <- as.integer(x_sci[ , 2])
                paste(alpha, power, sep = "e")
              }
              if (p_value == 0) p_value <- 2.2e-308
              if (0 < p_value && p_value < 1e-3)  {
                p_value <- formatCustomSci(p_value)
              }
              else {
                p_value  <- round(p_value, 3)
              }
            }
          }
          fig <- fig + ggtitle(paste("p:", p_value)) 
          fig_list[[arrangement_shape]][[pair]][[parameter]] <- fig
        }
      }
    }
  }
  
  arrangement_shape_figs <- list()
  
  for (shape in shapes) {  
    for (arrangement in arrangements) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      
      pair_figs <- list()
      
      for (pair in pairs) {
        
        # If it is the final pair, add parameters to bottom
        if (pair == pairs[length(pairs)]) {
          parameters <- get_parameters(arrangement, shape) 
          for (parameter in parameters) {
            subtitle <- ggdraw() + draw_label(parameter)
            fig_list[[arrangement_shape]][[pair]][[parameter]] <- plot_grid(fig_list[[arrangement_shape]][[pair]][[parameter]], 
                                                                            subtitle, ncol = 1, rel_heights = c(1, 0.1))
          }
        }
        
        pair_fig <- plot_grid(plotlist = fig_list[[arrangement_shape]][[pair]],
                              ncol = length(fig_list[[arrangement_shape]][[pair]]))
        title <- ggdraw() + draw_label(paste("pair:", pair))
        pair_fig <- plot_grid(title, pair_fig, ncol = 1, rel_heights = c(0.1, 1))
        pair_figs[[pair]] <- pair_fig 
      }
      
      arrangement_shape_fig <- plot_grid(plotlist = pair_figs, ncol = 1)
      
      title <- ggdraw() + draw_label(paste(arrangement, "-", shape, sep = ""), fontface = "bold")
      arrangement_shape_fig <- plot_grid(title, arrangement_shape_fig, ncol = 1, rel_heights = c(0.02, 1))  + 
        theme(plot.margin = margin(10, 10, 10, 10),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1))  
      
      arrangement_shape_figs[[arrangement_shape]] <- arrangement_shape_fig
    }
  }
  
  fig <- plot_grid(plotlist = arrangement_shape_figs,
                   nrow = length(shapes),
                   ncol = length(arrangements))
  
  return(fig) 
}



# Running the functions for 4 pair metrics ----
# For metrics which use all cell pairs (A/A, A/B, B/A, B/B)
metrics_with_4_pairs <- c("AMD",
                          "ANC_AUC",
                          "CK_AUC", "CL_AUC", "CG_AUC",
                          "FCO_AUC")


setwd("~/R/plots/SD2/metrics_with_4_cell_pairs")
pdf("fig_3D_vs_parameters_for_non_gradient_metrics_4_pairs_scatter_plot.pdf", width = 26, height = 25)

for (metric in metrics_with_4_pairs) {
  
  fig_3D_vs_parameters_for_non_gradient_metrics_scatter_plot <- plot_3D_vs_parameters_for_non_gradient_metrics_scatter_plot(metric_df_list,
                                                                                                                            parameters_df,
                                                                                                                            metric)
  print(fig_3D_vs_parameters_for_non_gradient_metrics_scatter_plot)
  
}
dev.off()







setwd("~/R/plots/SD2/metrics_with_4_cell_pairs")
pdf("fig_3D_and_2D_vs_parameters_for_non_gradient_metrics_4_pairs_scatter_plot.pdf", width = 26, height = 25)

for (metric in metrics_with_4_pairs) {
  
  fig_3D_and_2D_vs_parameters_for_non_gradient_metrics_scatter_plot <- plot_3D_and_2D_vs_parameters_for_non_gradient_metrics_scatter_plot(metric_df_list,
                                                                                                                                          parameters_df,
                                                                                                                                          metric)
  print(fig_3D_and_2D_vs_parameters_for_non_gradient_metrics_scatter_plot)
  
}
dev.off()





setwd("~/R/plots/SD2/metrics_with_4_cell_pairs")
pdf("fig_percentage_difference_vs_parameters_for_non_gradient_metrics_4_pairs_scatter_plot.pdf", width = 26, height = 25)

for (metric in metrics_with_4_pairs) {
  
  fig_percentage_difference_vs_parameters_for_non_gradient_metrics_scatter_plot <- plot_percentage_difference_vs_parameters_for_non_gradient_metrics_scatter_plot(metric_df_list,
                                                                                                                                                                  parameters_df,
                                                                                                                                                                  metric)
  print(fig_percentage_difference_vs_parameters_for_non_gradient_metrics_scatter_plot)
  
}
dev.off()







setwd("~/R/plots/SD2/metrics_with_4_cell_pairs")
pdf("fig_2D_vs_slice_for_non_gradient_metrics_4_pairs_violin_plot.pdf", width = 26, height = 25)

for (metric in metrics_with_4_pairs) {
  
  fig_2D_vs_slice_for_non_gradient_metrics_violin_plot <- plot_2D_vs_slice_for_non_gradient_metrics_violin_plot(metric_df_list,
                                                                                                                parameters_df,
                                                                                                                metric)
  print(fig_2D_vs_slice_for_non_gradient_metrics_violin_plot)
  
}
dev.off()






setwd("~/R/plots/SD2/metrics_with_4_cell_pairs")
pdf("fig_3D_and_2D_vs_slice_for_non_gradient_metrics_4_pairs_violin_plot.pdf", width = 26, height = 25)

for (metric in metrics_with_4_pairs) {
  
  fig_3D_and_2D_vs_slice_for_non_gradient_metrics_violin_plot <- plot_3D_and_2D_vs_slice_for_non_gradient_metrics_violin_plot(metric_df_list,
                                                                                                                              parameters_df,
                                                                                                                              metric)
  print(fig_3D_and_2D_vs_slice_for_non_gradient_metrics_violin_plot)
  
}
dev.off()



# Functions for 2 pair metrics-----
# For metrics which use only 2 cell pairs (A/B, B/A)
# The only change from functions for 4 pair metrics:
# arrangement_shape_fig <- plot_grid(title, arrangement_shape_fig, ncol = 1, rel_heights = c(0.02, 1))  BECOMES
# arrangement_shape_fig <- plot_grid(title, arrangement_shape_fig, ncol = 1, rel_heights = c(0.04, 1))    0.02 --> 0.04 

# Subset metric_df_list and parameters_df
add_pair_to_metric_df <- function(metric_df, metric) {
  
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
  return(metric_df)
}

get_parameters <- function(arrangement, shape) {
  parameters <- c()
  
  # Add other parameters specific to arrangement and shape
  if (arrangement == "mixed") {
    parameters <- c(parameters, "cluster_prop_A")
  }
  else if (arrangement == "ringed") {
    parameters <- c(parameters, "ring_width_factor")        
  }
  else if (arrangement == "separated") {
    parameters <- c(parameters, "distance")
  }
  
  if (shape == "ellipsoid") {
    parameters <- c(parameters, "E_volume")
  }
  else if (shape == "network") {
    parameters <- c(parameters, "N_width")
  }
  
  # All plots include bg_prop_A and bg_prop_B as parameters
  parameters <- c(parameters, "bg_prop_A", "bg_prop_B")
  
  return(parameters)
}

# For nicer tick labels
sci_clean_threshold <- function(x) {
  # x[!(x %in% range(x, na.rm = T))] <- NA
  
  sapply(x, function(v) {
    if (is.na(v)) {
      return('')
    }
    if (abs(v) < 1000) {
      return(as.character(v))   # keep normal numbers
    }
    # scientific notation
    s <- format(v, scientific = TRUE)   # e.g. "1e+03"
    s <- gsub("\\+", "", s)             # remove "+"
    s <- gsub("e0+", "e", s)            # remove leading zeros in exponent
    s
  })
}

combine_metric_and_parameters_df <- function(metric_df, parameters_df) {
  pairs <- unique(metric_df$pair)
  metric_and_parameters_df <- data.frame()
  for (slice in unique(metric_df$slice)) {
    metric_and_parameters_df <- rbind(metric_and_parameters_df, 
                                      cbind(metric_df[metric_df$slice == slice, ], 
                                            parameters_df[rep(1:nrow(parameters_df), each = length(pairs)), ]))
  }
  return(metric_and_parameters_df)
}

plot_3D_vs_parameters_for_non_gradient_metrics_scatter_plot <- function(metric_df_list, parameters_df, metric) {
  
  
  fig_list <- list()
  
  # Define parameters
  arrangements <- c("mixed", "ringed", "separated")
  shapes <- c("ellipsoid", "network")
  
  # Subset the metric_df_list
  metric_df <- metric_df_list[[metric]]
  
  # Remove Inf rows
  metric_df <- metric_df[!is.infinite(metric_df[[metric]]), ]
  
  # Add 'pair' column to metric_df
  metric_df <- add_pair_to_metric_df(metric_df, metric)
  
  pairs <- unique(metric_df$pair)
  
  # Add to parameters_df
  parameters_df$distance <- 450 - parameters_df$cluster1_x_coord # 450 is the x-coordinate of cluster2
  parameters_df$E_volume <- (4 / 3) * pi * parameters_df$E_radius_x * parameters_df$E_radius_y * parameters_df$E_radius_z
  
  # Update parameters_df
  parameters_df$variable_parameter[parameters_df$variable_parameter == "E_radius_x"] <- "E_volume"
  parameters_df$variable_parameter[parameters_df$variable_parameter == "cluster1_x_coord"] <- "distance"
  
  # Combine metric_df with parameters_df
  metric_df <- combine_metric_and_parameters_df(metric_df, parameters_df)
  
  # Subset for 3D values only
  metric_df <- metric_df[metric_df$slice == 0, ]
  
  # Update metric (remove _AUC from specific metrics)
  colnames(metric_df)[colnames(metric_df) == metric] <- sub("_AUC$", "", metric)
  metric <- sub("_AUC$", "", metric)
  
  for (arrangement in arrangements) {
    for (shape in shapes) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      fig_list[[arrangement_shape]] <- list()
      
      # Subset metric_df for arrangement and shape
      metric_arrangement_shape_df <- metric_df[metric_df$arrangement == arrangement & metric_df$shape == shape, ]
      
      # Get parameters
      parameters <- get_parameters(arrangement, shape)
      
      for (pair in pairs) {
        # Subset for pair
        metric_arrangement_shape_pair_df <- metric_arrangement_shape_df[metric_arrangement_shape_df$pair == pair, ]
        fig_list[[arrangement_shape]][[pair]] <- list()
        
        for (parameter in parameters) {
          
          # Subset for parameter
          plot_df <- metric_arrangement_shape_pair_df[metric_arrangement_shape_pair_df$variable_parameter == parameter, ]
          
          fig <- ggplot(plot_df, aes_string(parameter, metric)) +
            geom_point(size = 0.5) +
            theme_minimal() +
            theme(
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              plot.title = element_text(size = 10),
              axis.text.x = element_text(size = 8),  # make x-axis text smaller
              axis.text.y = element_text(size = 8)   # make y-axis text smaller
            ) +
            scale_x_continuous(breaks = pretty_breaks(n = 3), labels = sci_clean_threshold) + 
            scale_y_continuous(breaks = pretty_breaks(n = 3), labels = sci_clean_threshold) 
          
          # Get correlation and p-value       
          p_value <- "N/A"
          correlation <- "N/A"
          if (sum(!is.na(plot_df[[metric]])) >= 2) {
            
            correlation_test <- cor.test(plot_df[[parameter]], plot_df[[metric]], method = "spearman")
            correlation <- round(correlation_test$estimate, 3)
            p_value <- correlation_test$p.value
            
            if (!is.na(p_value)) {
              # Format correlation and p-value
              formatCustomSci <- function(x) {
                x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
                alpha <- round(as.numeric(x_sci[ , 1]), 1)
                power <- as.integer(x_sci[ , 2])
                paste(alpha, power, sep = "e")
              }
              if (p_value == 0) p_value <- 2.2e-308
              if (0 < p_value && p_value < 1e-3)  {
                p_value <- formatCustomSci(p_value)
              }
              else {
                p_value  <- round(p_value, 3)
              }
            }
          }
          fig <- fig + ggtitle(paste("r:", correlation, ", p:", p_value))
          fig_list[[arrangement_shape]][[pair]][[parameter]] <- fig
        }
      }
    }
  }
  
  arrangement_shape_figs <- list()
  
  for (shape in shapes) {  
    for (arrangement in arrangements) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      
      pair_figs <- list()
      
      for (pair in pairs) {
        
        parameters <- get_parameters(arrangement, shape)
        
        pair_fig <- plot_grid(plotlist = fig_list[[arrangement_shape]][[pair]],
                              ncol = length(fig_list[[arrangement_shape]][[pair]]))
        title <- ggdraw() + draw_label(paste("pair:", pair))
        pair_fig <- plot_grid(title, pair_fig, ncol = 1, rel_heights = c(0.1, 1))
        
        pair_figs[[pair]] <- pair_fig
      }
      arrangement_shape_fig <- plot_grid(plotlist = pair_figs, ncol = 1)
      
      title <- ggdraw() + draw_label(paste(arrangement, "-", shape, sep = ""), fontface = "bold")
      arrangement_shape_fig <- plot_grid(title, arrangement_shape_fig, ncol = 1, rel_heights = c(0.04, 1))  + 
        theme(plot.margin = margin(10, 10, 10, 10),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1))  
      
      arrangement_shape_figs[[arrangement_shape]] <- arrangement_shape_fig
    }
  }
  
  fig <- plot_grid(plotlist = arrangement_shape_figs,
                   nrow = length(shapes),
                   ncol = length(arrangements))
  
  return(fig)
}

plot_3D_and_2D_vs_parameters_for_non_gradient_metrics_scatter_plot <- function(metric_df_list, parameters_df, metric) {
  
  fig_list <- list()
  
  # Define parameters
  arrangements <- c("mixed", "ringed", "separated")
  shapes <- c("ellipsoid", "network")
  
  # Subset the metric_df_list
  metric_df <- metric_df_list[[metric]]
  
  # Remove Inf rows
  metric_df <- metric_df[!is.infinite(metric_df[[metric]]), ]
  
  # Add 'pair' column to metric_df
  metric_df <- add_pair_to_metric_df(metric_df, metric)
  
  pairs <- unique(metric_df$pair)
  
  # Make 'slice' column categorical
  metric_df$slice <- as.character(metric_df$slice)
  
  # Add to parameters_df
  parameters_df$distance <- 450 - parameters_df$cluster1_x_coord # 450 is the x-coordinate of cluster2
  parameters_df$E_volume <- (4 / 3) * pi * parameters_df$E_radius_x * parameters_df$E_radius_y * parameters_df$E_radius_z
  
  # Update parameters_df
  parameters_df$variable_parameter[parameters_df$variable_parameter == "E_radius_x"] <- "E_volume"
  parameters_df$variable_parameter[parameters_df$variable_parameter == "cluster1_x_coord"] <- "distance"
  
  # Combine metric_df with parameters_df
  metric_df <- combine_metric_and_parameters_df(metric_df, parameters_df)
  
  # Update metric (remove _AUC from specific metrics)
  colnames(metric_df)[colnames(metric_df) == metric] <- sub("_AUC$", "", metric)
  metric <- sub("_AUC$", "", metric)
  
  for (arrangement in arrangements) {
    for (shape in shapes) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      fig_list[[arrangement_shape]] <- list()
      
      # Subset metric_df for arrangement and shape
      metric_arrangement_shape_df <- metric_df[metric_df$arrangement == arrangement & metric_df$shape == shape, ]
      
      # Get parameters
      parameters <- get_parameters(arrangement, shape)
      
      for (pair in pairs) {
        # Subset for pair
        metric_arrangement_shape_pair_df <- metric_arrangement_shape_df[metric_arrangement_shape_df$pair == pair, ]
        fig_list[[arrangement_shape]][[pair]] <- list()
        
        for (parameter in parameters) {
          # Subset for parameter
          plot_df <- metric_arrangement_shape_pair_df[metric_arrangement_shape_pair_df$variable_parameter == parameter, ]
          
          fig <- ggplot(plot_df, aes_string(parameter, metric)) +
            
            # 1) First draw all *non‑zero* slices
            geom_point(
              data = subset(plot_df, slice != "0"),
              aes(color = slice),
              size = 0.5,
              alpha = 0.5
            ) +
            
            # 2) Then draw slice 0 on top
            geom_point(
              data = subset(plot_df, slice == "0"),
              color = "black",
              size = 0.5,
              alpha = 0.5
            ) +
            
            theme_minimal() +
            theme(
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              plot.title = element_text(size = 10),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              legend.position = "none"
            ) +
            scale_x_continuous(breaks = pretty_breaks(n = 3), labels = sci_clean_threshold) +
            scale_y_continuous(breaks = pretty_breaks(n = 3), labels = sci_clean_threshold) +
            scale_color_manual(
              values = c(
                "1" = "#9437a8",
                "2" = "#007128",
                "3" = "#b8db50"
                # no need to include "0" here since we draw it manually
              )
            )
          
          
          fig_list[[arrangement_shape]][[pair]][[parameter]] <- fig
        }
      }
    }
  }
  
  arrangement_shape_figs <- list()
  
  for (shape in shapes) {  
    for (arrangement in arrangements) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      
      pair_figs <- list()
      
      for (pair in pairs) {
        
        pair_fig <- plot_grid(plotlist = fig_list[[arrangement_shape]][[pair]],
                              ncol = length(fig_list[[arrangement_shape]][[pair]]))
        title <- ggdraw() + draw_label(paste("pair:", pair))
        pair_fig <- plot_grid(title, pair_fig, ncol = 1, rel_heights = c(0.1, 1))
        
        pair_figs[[pair]] <- pair_fig
      }
      arrangement_shape_fig <- plot_grid(plotlist = pair_figs, ncol = 1)
      
      title <- ggdraw() + draw_label(paste(arrangement, "-", shape, sep = ""), fontface = "bold")
      arrangement_shape_fig <- plot_grid(title, arrangement_shape_fig, ncol = 1, rel_heights = c(0.04, 1))  + 
        theme(plot.margin = margin(10, 10, 10, 10),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1))  
      
      arrangement_shape_figs[[arrangement_shape]] <- arrangement_shape_fig
    }
  }
  
  fig <- plot_grid(plotlist = arrangement_shape_figs,
                   nrow = length(shapes),
                   ncol = length(arrangements))
  
  return(fig) 
}

plot_percentage_difference_vs_parameters_for_non_gradient_metrics_scatter_plot <- function(metric_df_list, parameters_df, metric) {
  
  fig_list <- list()
  
  # Define parameters
  arrangements <- c("mixed", "ringed", "separated")
  shapes <- c("ellipsoid", "network")
  
  # Subset the metric_df_list
  metric_df <- metric_df_list[[metric]]
  
  # Remove Inf rows
  metric_df <- metric_df[!is.infinite(metric_df[[metric]]), ]
  
  # Add 'pair' column to metric_df
  metric_df <- add_pair_to_metric_df(metric_df, metric)
  
  pairs <- unique(metric_df$pair)
  
  # Make 'slice' column categorical
  metric_df$slice <- as.character(metric_df$slice)
  
  # Add to parameters_df
  parameters_df$distance <- 450 - parameters_df$cluster1_x_coord # 450 is the x-coordinate of cluster2
  parameters_df$E_volume <- (4 / 3) * pi * parameters_df$E_radius_x * parameters_df$E_radius_y * parameters_df$E_radius_z
  
  # Update parameters_df
  parameters_df$variable_parameter[parameters_df$variable_parameter == "E_radius_x"] <- "E_volume"
  parameters_df$variable_parameter[parameters_df$variable_parameter == "cluster1_x_coord"] <- "distance"
  
  # Combine metric_df with parameters_df
  metric_df <- combine_metric_and_parameters_df(metric_df, parameters_df)
  
  # Add error column to metric_df
  metric_df$error <- ((metric_df[[metric]] -  rep(metric_df[[metric]][metric_df$slice == "0"], 4)) / rep(metric_df[[metric]][metric_df$slice == "0"], 4)) * 100
  
  # Update metric (remove _AUC from specific metrics)
  colnames(metric_df)[colnames(metric_df) == metric] <- sub("_AUC$", "", metric)
  metric <- sub("_AUC$", "", metric)
  
  for (arrangement in arrangements) {
    for (shape in shapes) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      fig_list[[arrangement_shape]] <- list()
      
      # Subset metric_df for arrangement and shape
      metric_arrangement_shape_df <- metric_df[metric_df$arrangement == arrangement & metric_df$shape == shape, ]
      
      # Get parameters
      parameters <- get_parameters(arrangement, shape)
      
      for (pair in pairs) {
        # Subset for pair
        metric_arrangement_shape_pair_df <- metric_arrangement_shape_df[metric_arrangement_shape_df$pair == pair, ]
        fig_list[[arrangement_shape]][[pair]] <- list()
        
        for (parameter in parameters) {
          # Subset for parameter
          plot_df <- metric_arrangement_shape_pair_df[metric_arrangement_shape_pair_df$variable_parameter == parameter, ]
          
          # Remove 3D values
          plot_df <- plot_df[plot_df$slice != "0", ]
          
          fig <- ggplot(plot_df, aes_string(parameter, "error", color = "slice")) +
            geom_point(size = 0.5) +
            geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
            theme_minimal() +
            labs(y = paste(metric, "percentage diff (%)")) +
            theme(
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              plot.title = element_text(size = 10),
              axis.text.x = element_text(size = 8),  # make x-axis text smaller
              axis.text.y = element_text(size = 8),   # make y-axis text smaller
              legend.position = "none"
            ) +
            scale_x_continuous(breaks = pretty_breaks(n = 3), labels = sci_clean_threshold) + 
            scale_y_continuous(breaks = pretty_breaks(n = 3), labels = sci_clean_threshold) +
            scale_color_manual(
              values = c(
                "1" = "#9437a8",
                "2" = "#007128",
                "3" = "#b8db50"
              )
            )
          
          fig_list[[arrangement_shape]][[pair]][[parameter]] <- fig
        }
      }
    }
  }
  
  arrangement_shape_figs <- list()
  
  for (shape in shapes) {  
    for (arrangement in arrangements) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      
      pair_figs <- list()
      
      for (pair in pairs) {
        
        pair_fig <- plot_grid(plotlist = fig_list[[arrangement_shape]][[pair]],
                              ncol = length(fig_list[[arrangement_shape]][[pair]]))
        title <- ggdraw() + draw_label(paste("pair:", pair))
        pair_fig <- plot_grid(title, pair_fig, ncol = 1, rel_heights = c(0.1, 1))
        
        pair_figs[[pair]] <- pair_fig
      }
      arrangement_shape_fig <- plot_grid(plotlist = pair_figs, ncol = 1)
      
      title <- ggdraw() + draw_label(paste(arrangement, "-", shape, sep = ""), fontface = "bold")
      arrangement_shape_fig <- plot_grid(title, arrangement_shape_fig, ncol = 1, rel_heights = c(0.04, 1))  + 
        theme(plot.margin = margin(10, 10, 10, 10),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1))  
      
      arrangement_shape_figs[[arrangement_shape]] <- arrangement_shape_fig
    }
  }
  
  fig <- plot_grid(plotlist = arrangement_shape_figs,
                   nrow = length(shapes),
                   ncol = length(arrangements))
  
  return(fig) 
}

plot_2D_vs_slice_for_non_gradient_metrics_violin_plot <- function(metric_df_list, parameters_df, metric) {
  
  fig_list <- list()
  
  # Define parameters
  arrangements <- c("mixed", "ringed", "separated")
  shapes <- c("ellipsoid", "network")
  
  # Subset the metric_df_list
  metric_df <- metric_df_list[[metric]]
  
  # Remove Inf rows
  metric_df <- metric_df[!is.infinite(metric_df[[metric]]), ]
  
  # Add 'pair' column to metric_df
  metric_df <- add_pair_to_metric_df(metric_df, metric)
  
  pairs <- unique(metric_df$pair)
  
  # Make 'slice' column categorical
  metric_df$slice <- as.character(metric_df$slice)
  
  # Add to parameters_df
  parameters_df$distance <- 450 - parameters_df$cluster1_x_coord # 450 is the x-coordinate of cluster2
  parameters_df$E_volume <- (4 / 3) * pi * parameters_df$E_radius_x * parameters_df$E_radius_y * parameters_df$E_radius_z
  
  # Update parameters_df
  parameters_df$variable_parameter[parameters_df$variable_parameter == "E_radius_x"] <- "E_volume"
  parameters_df$variable_parameter[parameters_df$variable_parameter == "cluster1_x_coord"] <- "distance"
  
  # Combine metric_df with parameters_df
  metric_df <- combine_metric_and_parameters_df(metric_df, parameters_df)
  
  # Update metric (remove _AUC from specific metrics)
  colnames(metric_df)[colnames(metric_df) == metric] <- sub("_AUC$", "", metric)
  metric <- sub("_AUC$", "", metric)
  
  for (arrangement in arrangements) {
    for (shape in shapes) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      fig_list[[arrangement_shape]] <- list()
      
      # Subset metric_df for arrangement and shape
      metric_arrangement_shape_df <- metric_df[metric_df$arrangement == arrangement & metric_df$shape == shape, ]
      
      # Get parameters
      parameters <- get_parameters(arrangement, shape)
      
      for (pair in pairs) {
        # Subset for pair
        metric_arrangement_shape_pair_df <- metric_arrangement_shape_df[metric_arrangement_shape_df$pair == pair, ]
        fig_list[[arrangement_shape]][[pair]] <- list()
        
        for (parameter in parameters) {
          # Subset for parameter
          plot_df <- metric_arrangement_shape_pair_df[metric_arrangement_shape_pair_df$variable_parameter == parameter, ]
          
          # Remove 3D values
          plot_df <- plot_df[plot_df$slice != "0", ]
          
          # Violin plot
          fig <- ggplot(plot_df, aes_string(x = "slice", y = metric, fill = "slice")) +
            geom_violin() +
            theme_minimal() +
            theme(
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              plot.title = element_text(size = 10),
              axis.text.x = element_text(size = 8),  # make x-axis text smaller
              axis.text.y = element_text(size = 8),   # make y-axis text smaller
              legend.position = "none"
            ) +
            scale_y_continuous(breaks = pretty_breaks(n = 3), labels = sci_clean_threshold) +
            scale_fill_manual(
              values = c(
                "1" = "#9437a8",
                "2" = "#007128",
                "3" = "#b8db50"
              )
            )
          
          # Get JT p-value if possible
          p_value <- "N/A"
          if (sum(!is.na(plot_df[[metric]][plot_df$slice == "1"])) > 0 &&
              sum(!is.na(plot_df[[metric]][plot_df$slice == "2"])) > 0 &&
              sum(!is.na(plot_df[[metric]][plot_df$slice == "3"])) > 0) {
            
            # Factor slice column
            plot_df$slice <- factor(plot_df$slice, c("1", "2", "3"), ordered = TRUE)
            
            # Get JT value
            slice1_mean_value <- mean(plot_df[[metric]][plot_df$slice == "1"], na.rm = TRUE)
            slice3_mean_value <- mean(plot_df[[metric]][plot_df$slice == "3"], na.rm = TRUE)
            
            if (slice3_mean_value > slice1_mean_value) {
              trend_direction <- "increasing"
            }
            else {
              trend_direction <- "decreasing"
            }
            jt_test <- JonckheereTerpstraTest(plot_df[[metric]], plot_df$slice, alternative = trend_direction)
            
            p_value <- jt_test$p.value
            
            if (!is.na(p_value)) {
              # Format correlation and p-value
              formatCustomSci <- function(x) {
                x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
                alpha <- round(as.numeric(x_sci[ , 1]), 1)
                power <- as.integer(x_sci[ , 2])
                paste(alpha, power, sep = "e")
              }
              if (p_value == 0) p_value <- 2.2e-308
              if (0 < p_value && p_value < 1e-3)  {
                p_value <- formatCustomSci(p_value)
              }
              else {
                p_value  <- round(p_value, 3)
              }
            }
          }
          fig <- fig + ggtitle(paste("p:", p_value)) 
          fig_list[[arrangement_shape]][[pair]][[parameter]] <- fig
        }
      }
    }
  }
  
  arrangement_shape_figs <- list()
  
  for (shape in shapes) {  
    for (arrangement in arrangements) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      
      pair_figs <- list()
      
      for (pair in pairs) {
        
        # If it is the final pair, add parameters to bottom
        if (pair == pairs[length(pairs)]) {
          parameters <- get_parameters(arrangement, shape) 
          for (parameter in parameters) {
            subtitle <- ggdraw() + draw_label(parameter)
            fig_list[[arrangement_shape]][[pair]][[parameter]] <- plot_grid(fig_list[[arrangement_shape]][[pair]][[parameter]], 
                                                                            subtitle, ncol = 1, rel_heights = c(1, 0.1))
          }
        }
        
        pair_fig <- plot_grid(plotlist = fig_list[[arrangement_shape]][[pair]],
                              ncol = length(fig_list[[arrangement_shape]][[pair]]))
        title <- ggdraw() + draw_label(paste("pair:", pair))
        pair_fig <- plot_grid(title, pair_fig, ncol = 1, rel_heights = c(0.1, 1))
        pair_figs[[pair]] <- pair_fig 
      }
      
      arrangement_shape_fig <- plot_grid(plotlist = pair_figs, ncol = 1)
      
      title <- ggdraw() + draw_label(paste(arrangement, "-", shape, sep = ""), fontface = "bold")
      arrangement_shape_fig <- plot_grid(title, arrangement_shape_fig, ncol = 1, rel_heights = c(0.04, 1))  + 
        theme(plot.margin = margin(10, 10, 10, 10),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1))  
      
      arrangement_shape_figs[[arrangement_shape]] <- arrangement_shape_fig
    }
  }
  
  fig <- plot_grid(plotlist = arrangement_shape_figs,
                   nrow = length(shapes),
                   ncol = length(arrangements))
  
  return(fig) 
}

plot_3D_and_2D_vs_slice_for_non_gradient_metrics_violin_plot <- function(metric_df_list, parameters_df, metric) {
  
  fig_list <- list()
  
  # Define parameters
  arrangements <- c("mixed", "ringed", "separated")
  shapes <- c("ellipsoid", "network")
  
  # Subset the metric_df_list
  metric_df <- metric_df_list[[metric]]
  
  # Remove Inf rows
  metric_df <- metric_df[!is.infinite(metric_df[[metric]]), ]
  
  # Add 'pair' column to metric_df
  metric_df <- add_pair_to_metric_df(metric_df, metric)
  
  pairs <- unique(metric_df$pair)
  
  # Make 'slice' column categorical
  metric_df$slice <- as.character(metric_df$slice)
  
  # Add to parameters_df
  parameters_df$distance <- 450 - parameters_df$cluster1_x_coord # 450 is the x-coordinate of cluster2
  parameters_df$E_volume <- (4 / 3) * pi * parameters_df$E_radius_x * parameters_df$E_radius_y * parameters_df$E_radius_z
  
  # Update parameters_df
  parameters_df$variable_parameter[parameters_df$variable_parameter == "E_radius_x"] <- "E_volume"
  parameters_df$variable_parameter[parameters_df$variable_parameter == "cluster1_x_coord"] <- "distance"
  
  # Combine metric_df with parameters_df
  metric_df <- combine_metric_and_parameters_df(metric_df, parameters_df)
  
  # Update metric (remove _AUC from specific metrics)
  colnames(metric_df)[colnames(metric_df) == metric] <- sub("_AUC$", "", metric)
  metric <- sub("_AUC$", "", metric)
  
  for (arrangement in arrangements) {
    for (shape in shapes) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      fig_list[[arrangement_shape]] <- list()
      
      # Subset metric_df for arrangement and shape
      metric_arrangement_shape_df <- metric_df[metric_df$arrangement == arrangement & metric_df$shape == shape, ]
      
      # Get parameters
      parameters <- get_parameters(arrangement, shape)
      
      for (pair in pairs) {
        # Subset for pair
        metric_arrangement_shape_pair_df <- metric_arrangement_shape_df[metric_arrangement_shape_df$pair == pair, ]
        fig_list[[arrangement_shape]][[pair]] <- list()
        
        for (parameter in parameters) {
          # Subset for parameter
          plot_df <- metric_arrangement_shape_pair_df[metric_arrangement_shape_pair_df$variable_parameter == parameter, ]
          
          # Change slice index for 3D value to "3D"
          plot_df$slice[plot_df$slice == "0"] <- "3D"
          
          # Violin plot
          fig <- ggplot(plot_df, aes_string(x = "slice", y = metric, fill = "slice")) +
            geom_violin() +
            theme_minimal() +
            theme(
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              plot.title = element_text(size = 10),
              axis.text.x = element_text(size = 8),  # make x-axis text smaller
              axis.text.y = element_text(size = 8),   # make y-axis text smaller
              legend.position = "none"
            ) +
            scale_y_continuous(breaks = pretty_breaks(n = 3), labels = sci_clean_threshold) +
            scale_fill_manual(
              values = c(
                "1" = "#9437a8",
                "2" = "#007128",
                "3" = "#b8db50",
                "3D" = "white"
              )
            )
          
          # Get JT p-value if possible
          p_value <- "N/A"
          if (sum(!is.na(plot_df[[metric]][plot_df$slice == "1"])) > 0 &&
              sum(!is.na(plot_df[[metric]][plot_df$slice == "2"])) > 0 &&
              sum(!is.na(plot_df[[metric]][plot_df$slice == "3"])) > 0) {
            
            # Create separate df with only 2D values
            df2D <- plot_df[plot_df$slice != "3D", ]
            
            # Factor slice column
            df2D$slice <- factor(df2D$slice, c("1", "2", "3"), ordered = TRUE)
            
            # Get JT value
            slice1_mean_value <- mean(df2D[[metric]][df2D$slice == "1"], na.rm = TRUE)
            slice3_mean_value <- mean(df2D[[metric]][df2D$slice == "3"], na.rm = TRUE)
            
            if (slice3_mean_value > slice1_mean_value) {
              trend_direction <- "increasing"
            }
            else {
              trend_direction <- "decreasing"
            }
            jt_test <- JonckheereTerpstraTest(df2D[[metric]], df2D$slice, alternative = trend_direction)
            
            p_value <- jt_test$p.value
            
            if (!is.na(p_value)) {
              # Format correlation and p-value
              formatCustomSci <- function(x) {
                x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
                alpha <- round(as.numeric(x_sci[ , 1]), 1)
                power <- as.integer(x_sci[ , 2])
                paste(alpha, power, sep = "e")
              }
              if (p_value == 0) p_value <- 2.2e-308
              if (0 < p_value && p_value < 1e-3)  {
                p_value <- formatCustomSci(p_value)
              }
              else {
                p_value  <- round(p_value, 3)
              }
            }
          }
          fig <- fig + ggtitle(paste("p:", p_value)) 
          fig_list[[arrangement_shape]][[pair]][[parameter]] <- fig
        }
      }
    }
  }
  
  arrangement_shape_figs <- list()
  
  for (shape in shapes) {  
    for (arrangement in arrangements) {
      
      # Get a merged arrangement shape variable
      arrangement_shape <- paste(arrangement, shape, sep = "_")
      
      pair_figs <- list()
      
      for (pair in pairs) {
        
        # If it is the final pair, add parameters to bottom
        if (pair == pairs[length(pairs)]) {
          parameters <- get_parameters(arrangement, shape) 
          for (parameter in parameters) {
            subtitle <- ggdraw() + draw_label(parameter)
            fig_list[[arrangement_shape]][[pair]][[parameter]] <- plot_grid(fig_list[[arrangement_shape]][[pair]][[parameter]], 
                                                                            subtitle, ncol = 1, rel_heights = c(1, 0.1))
          }
        }
        
        pair_fig <- plot_grid(plotlist = fig_list[[arrangement_shape]][[pair]],
                              ncol = length(fig_list[[arrangement_shape]][[pair]]))
        title <- ggdraw() + draw_label(paste("pair:", pair))
        pair_fig <- plot_grid(title, pair_fig, ncol = 1, rel_heights = c(0.1, 1))
        pair_figs[[pair]] <- pair_fig 
      }
      
      arrangement_shape_fig <- plot_grid(plotlist = pair_figs, ncol = 1)
      
      title <- ggdraw() + draw_label(paste(arrangement, "-", shape, sep = ""), fontface = "bold")
      arrangement_shape_fig <- plot_grid(title, arrangement_shape_fig, ncol = 1, rel_heights = c(0.04, 1))  + 
        theme(plot.margin = margin(10, 10, 10, 10),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1))  
      
      arrangement_shape_figs[[arrangement_shape]] <- arrangement_shape_fig
    }
  }
  
  fig <- plot_grid(plotlist = arrangement_shape_figs,
                   nrow = length(shapes),
                   ncol = length(arrangements))
  
  return(fig) 
}



# Running the functions for 2 pair metrics ----
# For metrics which use only 2 cell pairs (A/B, B/A)
metrics_with_2_pairs <- c("ACIN_AUC", "ANE_AUC",
                          "MS_AUC", "NMS_AUC",
                          "PBP_AUC", "EBP_AUC", "PBSAC", "EBSAC")


setwd("~/R/plots/SD2/metrics_with_2_cell_pairs")
pdf("fig_3D_vs_parameters_for_non_gradient_metrics_2_pairs_scatter_plot.pdf", width = 26, height = 13)

for (metric in metrics_with_2_pairs) {
  
  fig_3D_vs_parameters_for_non_gradient_metrics_scatter_plot <- plot_3D_vs_parameters_for_non_gradient_metrics_scatter_plot(metric_df_list,
                                                                                                                            parameters_df,
                                                                                                                            metric)
  print(fig_3D_vs_parameters_for_non_gradient_metrics_scatter_plot)
  
}
dev.off()







setwd("~/R/plots/SD2/metrics_with_2_cell_pairs")
pdf("fig_3D_and_2D_vs_parameters_for_non_gradient_metrics_2_pairs_scatter_plot.pdf", width = 26, height = 13)

for (metric in metrics_with_2_pairs) {
  
  fig_3D_and_2D_vs_parameters_for_non_gradient_metrics_scatter_plot <- plot_3D_and_2D_vs_parameters_for_non_gradient_metrics_scatter_plot(metric_df_list,
                                                                                                                                          parameters_df,
                                                                                                                                          metric)
  print(fig_3D_and_2D_vs_parameters_for_non_gradient_metrics_scatter_plot)
  
}
dev.off()





setwd("~/R/plots/SD2/metrics_with_2_cell_pairs")
pdf("fig_percentage_difference_vs_parameters_for_non_gradient_metrics_2_pairs_scatter_plot.pdf", width = 26, height = 13)

for (metric in metrics_with_2_pairs) {
  
  fig_percentage_difference_vs_parameters_for_non_gradient_metrics_scatter_plot <- plot_percentage_difference_vs_parameters_for_non_gradient_metrics_scatter_plot(metric_df_list,
                                                                                                                                                                  parameters_df,
                                                                                                                                                                  metric)
  print(fig_percentage_difference_vs_parameters_for_non_gradient_metrics_scatter_plot)
  
}
dev.off()









setwd("~/R/plots/SD2/metrics_with_2_cell_pairs")
pdf("fig_2D_vs_slice_for_non_gradient_metrics_2_pairs_violin_plot.pdf", width = 26, height = 13)

for (metric in metrics_with_2_pairs) {
  
  fig_2D_vs_slice_for_non_gradient_metrics_violin_plot <- plot_2D_vs_slice_for_non_gradient_metrics_violin_plot(metric_df_list,
                                                                                                                parameters_df,
                                                                                                                metric)
  print(fig_2D_vs_slice_for_non_gradient_metrics_violin_plot)
  
}
dev.off()







setwd("~/R/plots/SD2/metrics_with_2_cell_pairs")
pdf("fig_3D_and_2D_vs_slice_for_non_gradient_metrics_2_pairs_violin_plot.pdf", width = 26, height = 13)

for (metric in metrics_with_2_pairs) {
  
  fig_3D_and_2D_vs_slice_for_non_gradient_metrics_violin_plot <- plot_3D_and_2D_vs_slice_for_non_gradient_metrics_violin_plot(metric_df_list,
                                                                                                                              parameters_df,
                                                                                                                              metric)
  print(fig_3D_and_2D_vs_slice_for_non_gradient_metrics_violin_plot)
  
}
dev.off()


