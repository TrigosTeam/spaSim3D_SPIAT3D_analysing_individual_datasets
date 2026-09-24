# Code to analyse a public 3D dataset made up of mutliple slices.

# ******** READ THIS *********
# This function assumes the following:
# data3D is a data frame with column names: "Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position", ""Cell.Type"
# "Cell.Z.Position" contains discrete values.
# cell_types consists of only cells in the "Cell.Type" column, and at least two cells.

analyse_3D_data_with_2D_slices <- function(
    data3D,
    cell_types,
    radii = seq(20, 100, 10),
    n_splits = 10,
    thresholds = seq(0.01, 1, 0.01)
) {
  
  slice_z_coords <- unique(data3D$Cell.Z.Position)
  n_slices <- length(unique(data3D$Cell.Z.Position))
  n_cell_type_combinations <- length(cell_types)^2
  radii_colnames <- paste("r", radii, sep = "")
  thresholds_colnames <- paste("t", thresholds, sep = "")
  
  create_empty_metric_df_list <- function(
    cell_types,
    n_slices,
    radii_colnames,
    thresholds_colnames
  ) {
    n_cell_type_combinations <- length(cell_types)^2
    
    # Define AMD data frames as well as constants
    AMD_df_colnames <- c("slice", "reference", "target", "AMD")
    AMD_df <- data.frame(matrix(nrow = (n_slices + 1) * n_cell_type_combinations, ncol = length(AMD_df_colnames)))
    colnames(AMD_df) <- AMD_df_colnames
    
    # Define MS, NMS, ANC, ACIN, COO, ANE data frames as well as constants
    radii_colnames <- paste("r", radii, sep = "")
    
    MS_df_colnames <- c("slice", "reference", "target", radii_colnames)
    MS_df <- data.frame(matrix(nrow = (n_slices + 1) * n_cell_type_combinations, ncol = length(MS_df_colnames)))
    colnames(MS_df) <- MS_df_colnames
    
    NMS_df <- ANC_df <- ANE_df <- ACIN_df <- COO_df <- CK_df <- CL_df <- CG_df <- MS_df
    
    # Define SAC and prevalence data frames as well as constants
    thresholds_colnames <- paste("t", thresholds, sep = "")
    
    PBSAC_df_colnames <- c("slice", "reference", "target", "PBSAC")
    PBSAC_df <- data.frame(matrix(nrow = (n_slices + 1) * n_cell_type_combinations, ncol = length(PBSAC_df_colnames)))
    colnames(PBSAC_df) <- PBSAC_df_colnames
    
    PBP_df_colnames <- c("slice", "reference", "target", thresholds_colnames)
    PBP_df <- data.frame(matrix(nrow = (n_slices + 1) * n_cell_type_combinations, ncol = length(PBP_df_colnames)))
    colnames(PBP_df) <- PBP_df_colnames
    
    EBSAC_df_colnames <- c("slice", "cell_types", "EBSAC")
    EBSAC_df <- data.frame(matrix(nrow = (n_slices + 1) * n_cell_type_combinations, ncol = length(EBSAC_df_colnames)))
    colnames(EBSAC_df) <- EBSAC_df_colnames
    
    EBP_df_colnames <- c("slice", "cell_types", thresholds_colnames)
    EBP_df <- data.frame(matrix(nrow = (n_slices + 1) * n_cell_type_combinations, ncol = length(EBP_df_colnames)))
    colnames(EBP_df) <- EBP_df_colnames
    
    
    # Add all to list:
    metric_df_list <- list(AMD = AMD_df,
                           MS = MS_df,
                           NMS = NMS_df,
                           ACIN = ACIN_df,
                           ANE = ANE_df,
                           ANC = ANC_df,
                           COO = COO_df,
                           CK = CK_df,
                           CL = CL_df,
                           CG = CG_df,
                           PBSAC = PBSAC_df,
                           PBP = PBP_df,
                           EBSAC = EBSAC_df,
                           EBP = EBP_df)
    
    return(metric_df_list)
  }
  
  metric_df_list <- create_empty_metric_df_list(
    cell_types,
    n_slices,
    radii_colnames,
    thresholds_colnames
  )
  
  for (i in seq(n_slices + 1, 1)) {
    print(i)
    
    # i represents the current slice index
    if (i == n_slices + 1) {
      df <- data3D
    }
    # if i is the last index, analyse in 3D instead
    else {
      df <- data3D[data3D$Cell.Z.Position == slice_z_coords[i], ]
    }
    
    ### 3D analysis -----------------------------
    if (i == n_slices + 1) {
      
      index <- n_cell_type_combinations * (i - 1) + 1 
      
      minimum_distance_data <- calculate_minimum_distances_between_cell_types3D(df,
                                                                                cell_types,
                                                                                show_summary = F,
                                                                                plot_image = F,
                                                                                feature_colname = "Cell.Type")
      
      minimum_distance_data_summary <- summarise_distances_between_cell_types3D(minimum_distance_data)
      
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "slice"] <- i
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "reference"] <- minimum_distance_data_summary$reference
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "target"] <- minimum_distance_data_summary$target
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "AMD"] <- minimum_distance_data_summary$mean
      
      # Need a new index for gradient-based data which increments after each target cell type
      pair_index <- n_cell_type_combinations * (i - 1) + 1 
      
      for (reference_cell_type in cell_types) {
        gradient_data <- calculate_all_gradient_cc_metrics3D(df,
                                                             reference_cell_type,
                                                             cell_types,
                                                             radii,
                                                             feature_colname = "Cell.Type")
        
        for (target_cell_type in cell_types) {
          print(paste(reference_cell_type, target_cell_type, sep = "/"))
          metric_df_list[["ANC"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["ACIN"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["COO"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["CK"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["CL"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["CG"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["MS"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["NMS"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["ANE"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, 
                                                                                         paste(reference_cell_type, target_cell_type, sep = ","))
          
          
          if (is.null(gradient_data)) {
            metric_df_list[["ANC"]][pair_index, radii_colnames] <- NA
            metric_df_list[["COO"]][pair_index, radii_colnames] <- NA
            metric_df_list[["CK"]][pair_index, radii_colnames] <- NA
            metric_df_list[["CL"]][pair_index, radii_colnames] <- NA
            metric_df_list[["CG"]][pair_index, radii_colnames] <- NA
            metric_df_list[["ACIN"]][pair_index, radii_colnames] <- NA
            metric_df_list[["MS"]][pair_index, radii_colnames] <- NA
            metric_df_list[["NMS"]][pair_index, radii_colnames] <- NA
            metric_df_list[["ANE"]][pair_index, radii_colnames] <- NA
          }
          else {
            metric_df_list[["ANC"]][pair_index, radii_colnames] <- gradient_data[["neighbourhood_counts"]][[target_cell_type]]
            metric_df_list[["COO"]][pair_index, radii_colnames] <- gradient_data[["co_occurrence"]][[target_cell_type]]
            metric_df_list[["CK"]][pair_index, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]] - gradient_data[["cross_K"]][["expected"]]
            metric_df_list[["CL"]][pair_index, radii_colnames] <- gradient_data[["cross_L"]][[target_cell_type]] - gradient_data[["cross_L"]][["expected"]]
            metric_df_list[["CG"]][pair_index, radii_colnames] <- gradient_data[["cross_G"]][[target_cell_type]][["observed_cross_G"]] - gradient_data[["cross_G"]][[target_cell_type]][["expected_cross_G"]]
            
            if (reference_cell_type != target_cell_type) {
              metric_df_list[["ACIN"]][pair_index, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
              metric_df_list[["MS"]][pair_index, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$mixing_score
              metric_df_list[["NMS"]][pair_index, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$normalised_mixing_score
              metric_df_list[["ANE"]][pair_index, radii_colnames] <- gradient_data[["neighbourhood_entropy"]][[target_cell_type]]
            }
          }
          if (reference_cell_type == target_cell_type) {
            metric_df_list[["ACIN"]][pair_index, radii_colnames] <- Inf
            metric_df_list[["MS"]][pair_index, radii_colnames] <- Inf
            metric_df_list[["NMS"]][pair_index, radii_colnames] <- Inf
            metric_df_list[["ANE"]][pair_index, radii_colnames] <- Inf
          }
          
          # Spatial heterogeneity metrics
          metric_df_list[["PBSAC"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["PBP"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          
          metric_df_list[["EBSAC"]][pair_index, c("slice", "cell_types")] <- c(i, paste(reference_cell_type, target_cell_type, sep = ","))
          metric_df_list[["EBP"]][pair_index, c("slice", "cell_types")] <- c(i, paste(reference_cell_type, target_cell_type, sep = ","))
          
          if (reference_cell_type != target_cell_type) {
            proportion_grid_metrics <- calculate_cell_proportion_grid_metrics3D(df, 
                                                                                n_splits,
                                                                                reference_cell_type, 
                                                                                target_cell_type,
                                                                                plot_image = F,
                                                                                feature_colname = "Cell.Type")
            
            if (is.null(proportion_grid_metrics)) {
              metric_df_list[["PBSAC"]][pair_index, "PBSAC"] <- NA
              metric_df_list[["PBP"]][pair_index, thresholds_colnames] <- NA
            }
            else {
              PBSAC <- calculate_spatial_autocorrelation3D(proportion_grid_metrics, 
                                                           "proportion",
                                                           weight_method = "queen")
              
              PBP_df <- calculate_prevalence_gradient3D(proportion_grid_metrics,
                                                        "proportion",
                                                        show_AUC = F,
                                                        plot_image = F)
              
              
              metric_df_list[["PBSAC"]][pair_index, "PBSAC"] <- PBSAC
              metric_df_list[["PBP"]][pair_index, thresholds_colnames] <- PBP_df$prevalence
            } 
            
            entropy_grid_metrics <- calculate_entropy_grid_metrics3D(df, 
                                                                     n_splits,
                                                                     c(reference_cell_type, target_cell_type), 
                                                                     plot_image = F,
                                                                     feature_colname = "Cell.Type")
            
            if (is.null(entropy_grid_metrics)) {
              metric_df_list[["EBSAC"]][pair_index, "EBSAC"] <- NA
              metric_df_list[["EBP"]][pair_index, thresholds_colnames] <- NA
            }
            else {
              EBSAC <- calculate_spatial_autocorrelation3D(entropy_grid_metrics, 
                                                           "entropy",
                                                           weight_method = "queen")
              
              EBP_df <- calculate_prevalence_gradient3D(entropy_grid_metrics,
                                                        "entropy",
                                                        show_AUC = F,
                                                        plot_image = F)
              
              metric_df_list[["EBSAC"]][pair_index, "EBSAC"] <- EBSAC
              metric_df_list[["EBP"]][pair_index, thresholds_colnames] <- EBP_df$prevalence
            }    
          }
          else {
            metric_df_list[["PBSAC"]][pair_index, "PBSAC"] <- Inf
            metric_df_list[["PBP"]][pair_index, thresholds_colnames] <- Inf
            metric_df_list[["EBSAC"]][pair_index, "EBSAC"] <- Inf
            metric_df_list[["EBP"]][pair_index, thresholds_colnames] <- Inf
          }
          
          pair_index <- pair_index + 1
        }
      }
    }
    ### 2D analysis -----------------------------
    else {
      
      index <- n_cell_type_combinations * (i - 1) + 1 
      
      minimum_distance_data <- calculate_minimum_distances_between_cell_types2D(df,
                                                                                cell_types,
                                                                                show_summary = F,
                                                                                plot_image = F,
                                                                                feature_colname = "Cell.Type")
      
      minimum_distance_data_summary <- summarise_distances_between_cell_types2D(minimum_distance_data)
      
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "slice"] <- i
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "reference"] <- minimum_distance_data_summary$reference
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "target"] <- minimum_distance_data_summary$target
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "AMD"] <- minimum_distance_data_summary$mean
      
      # Need a new index for gradient-based data which increments after each target cell type
      pair_index <- n_cell_type_combinations * (i - 1) + 1 
      
      for (reference_cell_type in cell_types) {
        gradient_data <- calculate_all_gradient_cc_metrics2D(df,
                                                             reference_cell_type,
                                                             cell_types,
                                                             radii,
                                                             feature_colname = "Cell.Type")
        
        for (target_cell_type in cell_types) {
          print(paste(reference_cell_type, target_cell_type, sep = "/"))
          metric_df_list[["ANC"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["ACIN"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["COO"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["CK"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["CL"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["CG"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["MS"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["NMS"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["ANE"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, 
                                                                                         paste(reference_cell_type, target_cell_type, sep = ","))
          
          
          if (is.null(gradient_data)) {
            metric_df_list[["ANC"]][pair_index, radii_colnames] <- NA
            metric_df_list[["COO"]][pair_index, radii_colnames] <- NA
            metric_df_list[["CK"]][pair_index, radii_colnames] <- NA
            metric_df_list[["CL"]][pair_index, radii_colnames] <- NA
            metric_df_list[["CG"]][pair_index, radii_colnames] <- NA
            metric_df_list[["ACIN"]][pair_index, radii_colnames] <- NA
            metric_df_list[["MS"]][pair_index, radii_colnames] <- NA
            metric_df_list[["NMS"]][pair_index, radii_colnames] <- NA
            metric_df_list[["ANE"]][pair_index, radii_colnames] <- NA
          }
          else {
            metric_df_list[["ANC"]][pair_index, radii_colnames] <- gradient_data[["neighbourhood_counts"]][[target_cell_type]]
            metric_df_list[["COO"]][pair_index, radii_colnames] <- gradient_data[["co_occurrence"]][[target_cell_type]]
            metric_df_list[["CK"]][pair_index, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]] - gradient_data[["cross_K"]][["expected"]]
            metric_df_list[["CL"]][pair_index, radii_colnames] <- gradient_data[["cross_L"]][[target_cell_type]] - gradient_data[["cross_L"]][["expected"]]
            metric_df_list[["CG"]][pair_index, radii_colnames] <- gradient_data[["cross_G"]][[target_cell_type]][["observed_cross_G"]] - gradient_data[["cross_G"]][[target_cell_type]][["expected_cross_G"]]
            
            if (reference_cell_type != target_cell_type) {
              metric_df_list[["ACIN"]][pair_index, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
              metric_df_list[["MS"]][pair_index, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$mixing_score
              metric_df_list[["NMS"]][pair_index, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$normalised_mixing_score
              metric_df_list[["ANE"]][pair_index, radii_colnames] <- gradient_data[["neighbourhood_entropy"]][[target_cell_type]]
            }
          }
          if (reference_cell_type == target_cell_type) {
            metric_df_list[["ACIN"]][pair_index, radii_colnames] <- Inf
            metric_df_list[["MS"]][pair_index, radii_colnames] <- Inf
            metric_df_list[["NMS"]][pair_index, radii_colnames] <- Inf
            metric_df_list[["ANE"]][pair_index, radii_colnames] <- Inf
          }
          
          # Spatial heterogeneity metrics
          metric_df_list[["PBSAC"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          metric_df_list[["PBP"]][pair_index, c("slice", "reference", "target")] <- c(i, reference_cell_type, target_cell_type)
          
          metric_df_list[["EBSAC"]][pair_index, c("slice", "cell_types")] <- c(i, paste(reference_cell_type, target_cell_type, sep = ","))
          metric_df_list[["EBP"]][pair_index, c("slice", "cell_types")] <- c(i, paste(reference_cell_type, target_cell_type, sep = ","))
          
          if (reference_cell_type != target_cell_type) {
            proportion_grid_metrics <- calculate_cell_proportion_grid_metrics2D(df, 
                                                                                n_splits,
                                                                                reference_cell_type, 
                                                                                target_cell_type,
                                                                                plot_image = F,
                                                                                feature_colname = "Cell.Type")
            
            if (is.null(proportion_grid_metrics)) {
              metric_df_list[["PBSAC"]][pair_index, "PBSAC"] <- NA
              metric_df_list[["PBP"]][pair_index, thresholds_colnames] <- NA
            }
            else {
              PBSAC <- calculate_spatial_autocorrelation2D(proportion_grid_metrics, 
                                                           "proportion",
                                                           weight_method = "queen")
              
              PBP_df <- calculate_prevalence_gradient2D(proportion_grid_metrics,
                                                        "proportion",
                                                        show_AUC = F,
                                                        plot_image = F)
              
              
              metric_df_list[["PBSAC"]][pair_index, "PBSAC"] <- PBSAC
              metric_df_list[["PBP"]][pair_index, thresholds_colnames] <- PBP_df$prevalence
            } 
            
            entropy_grid_metrics <- calculate_entropy_grid_metrics2D(df, 
                                                                     n_splits,
                                                                     c(reference_cell_type, target_cell_type), 
                                                                     plot_image = F,
                                                                     feature_colname = "Cell.Type")
            
            if (is.null(entropy_grid_metrics)) {
              metric_df_list[["EBSAC"]][pair_index, "EBSAC"] <- NA
              metric_df_list[["EBP"]][pair_index, thresholds_colnames] <- NA
            }
            else {
              EBSAC <- calculate_spatial_autocorrelation2D(entropy_grid_metrics, 
                                                           "entropy",
                                                           weight_method = "queen")
              
              EBP_df <- calculate_prevalence_gradient2D(entropy_grid_metrics,
                                                        "entropy",
                                                        show_AUC = F,
                                                        plot_image = F)
              
              metric_df_list[["EBSAC"]][pair_index, "EBSAC"] <- EBSAC
              metric_df_list[["EBP"]][pair_index, thresholds_colnames] <- EBP_df$prevalence
            }    
          }
          else {
            metric_df_list[["PBSAC"]][pair_index, "PBSAC"] <- Inf
            metric_df_list[["PBP"]][pair_index, thresholds_colnames] <- Inf
            metric_df_list[["EBSAC"]][pair_index, "EBSAC"] <- Inf
            metric_df_list[["EBP"]][pair_index, thresholds_colnames] <- Inf
          }
          
          pair_index <- pair_index + 1
        }
      }
    }
  }
  return(metric_df_list)
}


