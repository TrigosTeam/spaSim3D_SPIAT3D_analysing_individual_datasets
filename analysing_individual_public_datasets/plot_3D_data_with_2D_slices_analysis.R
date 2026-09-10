# Code to plot public 3D dataset analysis.

### Read data and set file name and save directory (THE ONLY PART YOU NEED TO CHANGE) ------
metric_df_list <- readRDS("***your metric_df_list.RDS***") # e.g. "~/R/public_data_analysis/CyCIF_metric_df_list.RDS"

file_name_prefix <- "*** prefix for file name ***" # e.g. "colorectal_cancer"
save_directory <- "*** directory to save plots ***" # e.g. "~/R/plots/public_data/colorectal_cancer"

### Libraries -----
library(SpatialExperiment)
library(dbscan)
library(alphashape3d)
library(apcluster)
library(plotly)
library(dplyr)
library(reshape2)
library(gtools)
library(cowplot)
library(Hmisc)
library(alphashape3d)

### Format data -----

get_gradient <- function(metric) {
  if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "COO", "CK", "CL", "CG")) {
    return("radius")
  }
  else if (metric %in% c("PBP", "EBP")) {
    return("threshold")  
  }
  else {
    stop("Invalid metric. Must be gradient-based")
  }
}


## Turn gradient radii metrics into AUC and add to metric_df list
get_AUC_for_radii_gradient_metrics <- function(y) {
  x <- radii
  h <- diff(x)[1]
  n <- length(x)
  
  AUC <- (h / 2) * (y[1] + 2 * sum(y[2:(n - 1)]) + y[n])
  
  return(AUC)
}


radii <- seq(20, 100, 10)
radii_colnames <- paste("r", radii, sep = "")

gradient_radii_metrics <- c("MS", "NMS", "ACIN", "ANE", "ANC", "COO", "CK", "CL", "CG")


for (metric in gradient_radii_metrics) {
  metric_AUC_name <- paste(metric, "AUC", sep = "_")
  
  if (metric %in% c("MS", "NMS", "ANC", "COO", "CK", "CL", "CG")) {
    subset_colnames <- c("slice", "reference", "target", metric_AUC_name)
  }
  else {
    subset_colnames <- c("slice", "reference", metric_AUC_name)
  }
  
  df <- metric_df_list[[metric]]
  df[[metric_AUC_name]] <- apply(df[ , radii_colnames], 1, get_AUC_for_radii_gradient_metrics)
  metric_df_list[[metric_AUC_name]] <- df
  
}

## Turn threshold radii metrics into AUC and add to metric_df list
thresholds <- seq(0.01, 1, 0.01)
thresholds_colnames <- paste("t", thresholds, sep = "")

# PBP_AUC 3D
PBP_df <- metric_df_list[["PBP"]]
PBP_df$PBP_AUC <- apply(PBP_df[ , thresholds_colnames], 1, function(y) {
  sum(diff(thresholds) * (head(y, -1) + tail(y, -1)) / 2)
})
PBP_AUC_df <- PBP_df[ , c("slice", "reference", "target", "PBP_AUC")]
metric_df_list[["PBP_AUC"]] <- PBP_AUC_df

# EBP_AUC 3D
EBP_df <- metric_df_list[["EBP"]]
EBP_df$EBP_AUC <- apply(EBP_df[ , thresholds_colnames], 1, function(y) {
  sum(diff(thresholds) * (head(y, -1) + tail(y, -1)) / 2)
})
EBP_AUC_df <- EBP_df[ , c("slice", "cell_types", "EBP_AUC")]
metric_df_list[["EBP_AUC"]] <- EBP_AUC_df


### Plotting functions -----

## Functions to plot

# Utility function to add pair column
add_pair_column_to_metric_df <- function(metric_df,
                                         metric) {
  
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
  
  return(metric_df)
}

plot_3D_and_2D_metric_vs_pair_box_plot <- function(metric_df_list,
                                                   metric) {
  
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
  
  # Get metric_df for current metric
  metric_df <- metric_df_list[[metric]]
  
  # Change and further subset columns of metric_df_subset
  colnames(metric_df)[colnames(metric_df) == metric] <- "value"
  
  # Add pair column
  metric_df <- add_pair_column_to_metric_df(metric_df, metric)
  
  # Update metric (remove _AUC from specific metrics)
  metric <- sub("_AUC$", "", metric)
  
  fig <- ggplot(metric_df, aes(x = pair, y = value)) +
    geom_boxplot(data = metric_df[metric_df$slice != as.character(max(as.numeric(metric_df$slice))), ],
                 outlier.shape = NA, fill = "lightgray") +
    geom_jitter(data = metric_df[metric_df$slice != as.character(max(as.numeric(metric_df$slice))), ],
                width = 0.2, alpha = 0.5, color = "#0062c5") +
    geom_point(data = metric_df[metric_df$slice == as.character(max(as.numeric(metric_df$slice))), ],
               shape = 8, color = "#bb0036", size = 3) +  # Red stars for 3D value
    labs(title = paste("Boxplots showing 3D and 2D ", metric, " vs cell pairs", sep = ""), x = "", y = metric) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)  # Rotate x-axis labels vertically
    ) +
    scale_y_continuous(labels = sci_clean_threshold)
  
  
  return(fig)
}

plot_percentage_difference_vs_pair_box_plot <- function(metric_df_list,
                                                        metric) {
  
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
  
  # Get metric_df for current metric
  metric_df <- metric_df_list[[metric]]
  
  # Change and further subset columns of metric_df
  colnames(metric_df)[colnames(metric_df) == metric] <- "value"
  
  # Obtain 3D value
  value_3D <- metric_df[["value"]][metric_df[["slice"]] == as.character(max(as.numeric(metric_df$slice)))]
  
  # Calculate error
  metric_df[["error"]] <- ((metric_df[["value"]] - value_3D) / value_3D) * 100
  
  # Remove value column (only using error now)
  metric_df["value"] <- NULL
  
  # Remove 3D data (integrated into error)
  metric_df <- metric_df[metric_df[["slice"]] != max(as.integer(metric_df$slice)), ]
  
  # Add pair column
  metric_df <- add_pair_column_to_metric_df(metric_df, metric)
  
  # Update metric (remove _AUC from specific metrics)
  metric <- sub("_AUC$", "", metric)
  
  fig <- ggplot(metric_df, aes(x = pair, y = error)) +
    geom_boxplot(outlier.shape = NA, fill = "lightgray") +  # Hide default outliers to avoid duplication
    geom_jitter(width = 0.2, alpha = 0.5, color = "#0062c5") +  # Add dots with slight horizontal jitter
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    labs(title = paste("Boxplots showing percentage difference between 2D and 3D ", metric, " vs cell pairs", sep = ""),
         x = "Pair",
         y = "Percentage difference (%)") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)  # Rotate x-axis labels vertically
    ) +
    scale_y_continuous(labels = sci_clean_threshold)
  
  return(fig)
}

plot_percentage_difference_vs_slice_box_plot <- function(metric_df_list,
                                                         metric) {
  
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
  
  # Get metric_df for current metric
  metric_df <- metric_df_list[[metric]]
  
  # Change and further subset columns of metric_df
  colnames(metric_df)[colnames(metric_df) == metric] <- "value"
  
  # Obtain 3D value
  value_3D <- metric_df[["value"]][metric_df[["slice"]] == as.character(max(as.numeric(metric_df$slice)))]
  
  # Calculate error
  metric_df[["error"]] <- ((metric_df[["value"]] - value_3D) / value_3D) * 100
  
  # Remove value column (only using error now)
  metric_df["value"] <- NULL
  
  # Remove 3D data (integrated into error)
  metric_df <- metric_df[metric_df[["slice"]] != max(as.integer(metric_df$slice)), ]
  
  # Factor slice column so it is from 1, 2, 3, ...
  metric_df$slice <- factor(metric_df$slice, as.character(sort(as.integer(unique(metric_df$slice)))))
  
  # Update metric (remove _AUC from specific metrics)
  metric <- sub("_AUC$", "", metric)
  
  fig <- ggplot(metric_df, aes(x = slice, y = error)) +
    geom_boxplot(outlier.shape = NA, fill = "lightgray") +  # Hide default outliers to avoid duplication
    geom_jitter(width = 0.2, alpha = 0.5, color = "#0062c5") +  # Add dots with slight horizontal jitter
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    labs(title = paste("Boxplots showing percentage difference between 2D and 3D ", metric, " vs slice index", sep = ""),
         x = "Slice index",
         y = "Percentage difference (%)") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)  # Rotate x-axis labels vertically
    ) +
    scale_y_continuous(labels = sci_clean_threshold)
  
  return(fig)
}


plot_median_percentage_difference_of_each_pair_vs_metrics_box_plot <- function(metric_df_list,
                                                                               metrics) {
  
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
  
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change and further subset columns of metric_df
    colnames(metric_df)[colnames(metric_df) == metric] <- "value"
    
    # Obtain 3D value
    value_3D <- metric_df[["value"]][metric_df[["slice"]] == as.character(max(as.numeric(metric_df$slice)))]
    
    # Calculate error
    metric_df[["error"]] <- ((metric_df[["value"]] - value_3D) / value_3D) * 100
    
    # Remove value column (only using error now)
    metric_df["value"] <- NULL
    
    # Remove 3D data (integrated into error)
    metric_df <- metric_df[metric_df[["slice"]] != max(as.integer(metric_df$slice)), ]
    
    # Add pair column
    metric_df <- add_pair_column_to_metric_df(metric_df, metric)
    
    # Extract median values for error for each pair
    median_df <- metric_df %>%
      group_by(pair) %>%
      dplyr::summarize(median_error = median(error, na.rm = TRUE), .groups = "drop")
    
    median_df$metric <- metric  
    
    # Add median_df to plot_df
    plot_df <- rbind(plot_df, median_df)
  }
  
  # Update metrics (remove _AUC from specific metrics)
  metrics <- sub("_AUC$", "", metrics)
  plot_df$metric <- sub("_AUC$", "", plot_df$metric)
  
  # Factor for metric
  plot_df$metric <- factor(plot_df$metric, metrics)
  
  fig <- ggplot(plot_df, aes(x = metric, y = median_error)) +
    geom_boxplot(outlier.shape = NA, fill = "lightgray") +  # Hide default outliers to avoid duplication
    geom_jitter(width = 0.2, alpha = 0.5, color = "#0062c5") +  # Add dots with slight horizontal jitter
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    labs(title = "Box plots showing median percentage difference between 2D and 3D metrics of each pair",
         x = "Metric",
         y = "Percentage difference (%)") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    scale_y_continuous(labels = sci_clean_threshold)
  
  return(fig)
}

plot_median_percentage_difference_of_each_slice_vs_metrics_box_plot <- function(metric_df_list,
                                                                                metrics) {
  
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
  
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change and further subset columns of metric_df
    colnames(metric_df)[colnames(metric_df) == metric] <- "value"
    
    # Obtain 3D value
    value_3D <- metric_df[["value"]][metric_df[["slice"]] == as.character(max(as.numeric(metric_df$slice)))]
    
    # Calculate error
    metric_df[["error"]] <- ((metric_df[["value"]] - value_3D) / value_3D) * 100
    
    # Remove value column (only using error now)
    metric_df["value"] <- NULL
    
    # Remove 3D data (integrated into error)
    metric_df <- metric_df[metric_df[["slice"]] != max(as.integer(metric_df$slice)), ]
    
    # Factor slice column so it is from 1, 2, 3, ...
    metric_df$slice <- factor(metric_df$slice, as.character(sort(as.integer(unique(metric_df$slice)))))
    
    
    median_df <- metric_df %>%
      group_by(slice) %>%
      dplyr::summarize(median_error = median(error, na.rm = TRUE), .groups = "drop")
    
    median_df$metric <- metric  
    
    # Add median_df to plot_df
    plot_df <- rbind(plot_df, median_df)
  }
  
  # Update metrics (remove _AUC from specific metrics)
  metrics <- sub("_AUC$", "", metrics)
  plot_df$metric <- sub("_AUC$", "", plot_df$metric)
  
  # Factor metrics
  plot_df$metric <- factor(plot_df$metric, metrics)
  
  summary_df <- plot_df %>%
    group_by(metric) %>%
    summarise(
      min     = min(median_error, na.rm = TRUE),
      max     = max(median_error, na.rm = TRUE),
      range   = max - min,
      median  = median(median_error, na.rm = TRUE),
      average = mean(median_error, na.rm = TRUE)
    )
  write.csv(summary_df, paste("~/R/values_from_figures/", file_name_prefix, "_summary_df.csv", sep = ""))
  
  fig <- ggplot(plot_df, aes(x = metric, y = median_error)) +
    geom_boxplot(outlier.shape = NA, fill = "lightgray") +  # Hide default outliers to avoid duplication
    geom_jitter(width = 0.2, alpha = 0.5, color = "#0062c5") +  # Add dots with slight horizontal jitter
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    labs(
      # title = "Box plots showing median percentage difference between 2D and 3D metrics of each slice",
       x = "Metric",
       y = "Percentage difference (%)"
      ) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      
      # Axis tick labels 
      axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 15), 
      
      # Axis titles 
      axis.title.x = element_text(size = 15), 
      axis.title.y = element_text(size = 15), 
      
      # Plot title 
      plot.title = element_text(size = 15)
    ) +
    scale_y_continuous(labels = sci_clean_threshold)
  
  return(fig)
}


plot_percentage_difference_vs_metrics_for_all_pairs_and_slices_box_plot <- function(metric_df_list,
                                                                                    metrics) {
  
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
  
  plot_df <- data.frame()
  
  for (metric in metrics) {
    
    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]
    
    # Change and further subset columns of metric_df
    colnames(metric_df)[colnames(metric_df) == metric] <- "value"
    
    # Obtain 3D value
    value_3D <- metric_df[["value"]][metric_df[["slice"]] == as.character(max(as.numeric(metric_df$slice)))]
    
    # Calculate error
    metric_df[["error"]] <- ((metric_df[["value"]] - value_3D) / value_3D) * 100
    
    # Remove value column (only using error now)
    metric_df["value"] <- NULL
    
    # Remove 3D data (integrated into error)
    metric_df <- metric_df[metric_df[["slice"]] != max(as.integer(metric_df$slice)), ]
    
    # Add metric column
    metric_df$metric <- metric
    
    # Keep error and metric column
    metric_df <- metric_df[ , c("error", "metric")]
    
    # Add median_df to plot_df
    plot_df <- rbind(plot_df, metric_df)
  }
  
  # Update metrics (remove _AUC from specific metrics)
  metrics <- sub("_AUC$", "", metrics)
  plot_df$metric <- sub("_AUC$", "", plot_df$metric)
  
  # Factor for metric
  plot_df$metric <- factor(plot_df$metric, metrics)
  
  fig <- ggplot(plot_df, aes(x = metric, y = error)) +
    geom_boxplot(outlier.shape = NA, fill = "lightgray") +  # Hide default outliers to avoid duplication
    geom_hline(yintercept = 0, color = "#bb0036", linetype = "dotted", linewidth = 1) + # Red dotted line at y = 0
    labs(title = "Box plots showing percentage difference between 2D and 3D metrics for all cell pairs and slices",
         x = "Metric",
         y = "Percentage difference (%)") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    scale_y_continuous(labels = sci_clean_threshold)
  
  return(fig)
}


## Get plot
metric <- "AMD"
metrics <- c("AMD",
             "ANC_AUC", "ACIN_AUC", "ANE_AUC",
             "MS_AUC", "NMS_AUC",
             "CK_AUC", "CL_AUC", "CG_AUC",
             "COO_AUC",
             "PBP_AUC", "EBP_AUC", "PBSAC", "EBSAC")


### Plotting and upload ------
setwd(save_directory)
# This is for a SINGLE metric
fig_3D_and_2D_metric_vs_pair_box_plot <- plot_3D_and_2D_metric_vs_pair_box_plot(metric_df_list,
                                                                                metric)
pdf(paste(file_name_prefix, "_fig_3D_and_2D_metric_vs_pair_box_plot.pdf", sep = ""), width = 9, height = 6)
print(fig_3D_and_2D_metric_vs_pair_box_plot)
dev.off()



setwd(save_directory)
# This is for a SINGLE metric
fig_percentage_difference_vs_pair_box_plot <- plot_percentage_difference_vs_pair_box_plot(metric_df_list,
                                                                                          metric)
pdf(paste(file_name_prefix, "_fig_percentage_difference_vs_pair_box_plot.pdf", sep = ""), width = 9, height = 6)
print(fig_percentage_difference_vs_pair_box_plot)
dev.off()



setwd(save_directory)
# This is for a SINGLE metric
fig_percentage_difference_vs_slice_box_plot <- plot_percentage_difference_vs_slice_box_plot(metric_df_list,
                                                                                            metric)
pdf(paste(file_name_prefix, "_fig_percentage_difference_vs_slice_box_plot.pdf", sep = ""), width = 9, height = 6)
print(fig_percentage_difference_vs_slice_box_plot)
dev.off()



setwd(save_directory)
fig_median_percentage_difference_of_each_pair_vs_metrics_box_plot <- 
  plot_median_percentage_difference_of_each_pair_vs_metrics_box_plot(metric_df_list,
                                                                     metrics)
pdf(paste(file_name_prefix, "_fig_median_percentage_difference_of_each_pair_vs_metrics_box_plot.pdf", sep = ""), width = 9, height = 6)
print(fig_median_percentage_difference_of_each_pair_vs_metrics_box_plot)
dev.off()




setwd(save_directory)
fig_median_percentage_difference_of_each_slice_vs_metrics_box_plot <- 
  plot_median_percentage_difference_of_each_slice_vs_metrics_box_plot(metric_df_list,
                                                                      metrics)
pdf(paste(file_name_prefix, "_fig_median_percentage_difference_of_each_slice_vs_metrics_box_plot.pdf", sep = ""), width = 9, height = 5)
print(fig_median_percentage_difference_of_each_slice_vs_metrics_box_plot)
dev.off()



setwd(save_directory)
fig_percentage_difference_vs_metrics_for_all_pairs_and_slices_box_plot <- 
  plot_percentage_difference_vs_metrics_for_all_pairs_and_slices_box_plot(metric_df_list,
                                                                          metrics)
pdf(paste(file_name_prefix, "_fig_percentage_difference_vs_metrics_for_all_pairs_and_slices_box_plot.pdf", sep = ""), width = 9, height = 6)
print(fig_percentage_difference_vs_metrics_for_all_pairs_and_slices_box_plot)
dev.off()



