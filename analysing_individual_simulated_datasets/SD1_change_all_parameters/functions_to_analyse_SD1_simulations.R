# Code to generate SD1 and to perform SPIAT-3D analyis on SD1.
# Set seed for consistency
set.seed(678999821)


# Functions -----
# Generate data frame containing values of important parameters
generate_random_parameters <- function(
    n_simulations,
    bg_prop_A_range,
    bg_prop_B_range,
    E_radius_x_range,
    E_radius_y_range,
    E_radius_z_range,
    N_width_range,
    cluster_prop_A_range,
    ring_width_factor_range,
    cluster1_x_coord_range
) {

  # For each parameter, get 'n_simulations' values between specified min and max of each parameter  
  parameters_df <- data.frame(
    bg_prop_A = runif(n_simulations, bg_prop_A_range["min"], bg_prop_A_range["max"]),
    bg_prop_B = runif(n_simulations, bg_prop_B_range["min"], bg_prop_B_range["max"]),
    E_radius_x = runif(n_simulations, E_radius_x_range["min"], E_radius_x_range["max"]),
    E_radius_y = runif(n_simulations, E_radius_y_range["min"], E_radius_y_range["max"]),
    E_radius_z = runif(n_simulations, E_radius_z_range["min"], E_radius_z_range["max"]),
    N_width = runif(n_simulations, N_width_range["min"], N_width_range["max"]),
    cluster_prop_A = runif(n_simulations, cluster_prop_A_range["min"], cluster_prop_A_range["max"]),
    ring_width_factor = runif(n_simulations, ring_width_factor_range["min"], ring_width_factor_range["max"]),
    cluster1_x_coord = runif(n_simulations, cluster1_x_coord_range["min"], cluster1_x_coord_range["max"])
  )
  
  # Randomly choose which arrangement and shape to use
  arrangements <- c("mixed", "ringed", "separated")
  shapes <- c("ellipsoid", "network")
  
  parameters_df$arrangement <- sample(arrangements, size = nrow(parameters_df), replace = TRUE)
  parameters_df$shape <- sample(shapes, size = nrow(parameters_df), replace = TRUE)
  
  return(parameters_df)
}

# Analyse SD1 simulations
analyse_SD1_simulations <- function(parameters_df) {
  
  # Set defined parameters/values
  cell_types <- c('A', 'B')
  n_simulations <- nrow(parameters_df)
  
  radii <- seq(20, 100, 10)
  radii_colnames <- paste("r", radii, sep = "")
  
  n_splits <- 10
  thresholds <- seq(0.01, 1, 0.01)
  thresholds_colnames <- paste("t", thresholds, sep = "")
  
  # z-coords for slices
  bottom_z_coord_of_slices <- c(85, 95, 105, 115, 125, 135, 145, 155, 165, 175, 185, 195, 205)
  top_z_coord_of_slices <- bottom_z_coord_of_slices + 10
  n_slices <- length(bottom_z_coord_of_slices)
  
  # Generate simulation metadata using parameters data frame from 'generate_random_parameters' function
  generate_simulation_metadata <- function(parameters_df, simulation_index) {
    
    # Constant metadata for simulations
    simulation_metadata <- spe_metadata_background_template("random")
    simulation_metadata$background$n_cells <- 30000
    simulation_metadata$background$length <- 600
    simulation_metadata$background$width <- 600
    simulation_metadata$background$height <- 300
    simulation_metadata$background$minimum_distance_between_cells <- 10
    simulation_metadata$background$cell_types <- c("A", "B", "O") # Cell proportions will change later
    
    N_radius <- 125
    N_n_edges <- 20
    
    mixed_cluster_cell_types <- c("A", "B")
    mixed_cluster_centre_loc <- c(300, 300, 150)
    
    ringed_cluster_cell_type <- "A"
    ringed_cluster_cell_prop <- 1
    ringed_ring_cell_type <- "B"
    ringed_ring_cell_prop <- 1
    ringed_cluster_centre_loc <- c(300, 300, 150)
    
    separated_cluster1_cell_type <- "A"
    separated_cluster1_cell_prop <- 1
    separated_cluster1_y_z_centre_loc <- c(300, 150) # Append x-coord to it later
    separated_cluster2_shape <- "sphere"
    separated_cluster2_cell_type <- "B"
    separated_cluster2_cell_prop <- 1
    separated_cluster2_centre_loc <- c(450, 300, 150)
    separated_cluster2_radius <- 100
    
    # Alter background cell proportions
    simulation_metadata$background$cell_proportions <- c(parameters_df$bg_prop_A[simulation_index], 
                                                         parameters_df$bg_prop_B[simulation_index],
                                                         1 - parameters_df$bg_prop_A[simulation_index] - parameters_df$bg_prop_B[simulation_index]) # prop(O) = 1 - prop(A) - prop(B)
    
    shape <- parameters_df[simulation_index, "shape"]
    arrangement <- parameters_df[simulation_index, "arrangement"]
    
    # Determine what type of cluster to add to current metadata
    if (arrangement %in% c("mixed", "separated")) {
      simulation_metadata <- spe_metadata_cluster_template("regular", shape, simulation_metadata)
    }
    else if (arrangement == "ringed") {
      simulation_metadata <- spe_metadata_cluster_template("ring", shape, simulation_metadata)
    }
    
    # Modify shape parameters
    if (shape == "ellipsoid") {
      simulation_metadata$cluster_1$radii <- c(parameters_df$E_radius_x[simulation_index], parameters_df$E_radius_y[simulation_index], parameters_df$E_radius_z[simulation_index])
      simulation_metadata$cluster_1$axes_rotation <- c(runif(1, 0, 180), runif(1, 0, 180), runif(1, 0, 180)) # Choose random rotation angles
    }
    
    else if (shape == "network") {
      simulation_metadata$cluster_1$n_edges <- N_n_edges
      simulation_metadata$cluster_1$width <- parameters_df$N_width[simulation_index]
      simulation_metadata$cluster_1$radius <- N_radius
    }
    
    # Modify arrangemetn parameters
    if (arrangement == "mixed") {
      simulation_metadata$cluster_1$cluster_cell_types <- mixed_cluster_cell_types
      simulation_metadata$cluster_1$cluster_cell_proportions <- c(parameters_df$cluster_prop_A[simulation_index], 1 - parameters_df$cluster_prop_A[simulation_index])
      simulation_metadata$cluster_1$centre_loc <- mixed_cluster_centre_loc
    }
    else if (arrangement == "ringed") {
      simulation_metadata$cluster_1$cluster_cell_types <- ringed_cluster_cell_type
      simulation_metadata$cluster_1$cluster_cell_proportions <- ringed_cluster_cell_prop
      simulation_metadata$cluster_1$centre_loc <- ringed_cluster_centre_loc
      simulation_metadata$cluster_1$ring_cell_types <- ringed_ring_cell_type
      simulation_metadata$cluster_1$ring_cell_proportions <- ringed_ring_cell_prop
      simulation_metadata$cluster_1$ring_width <- parameters_df$ring_width_factor[simulation_index] * ifelse(shape == "network", 
                                                                                                             parameters_df$N_width[simulation_index],
                                                                                                             apply(parameters_df[, c("E_radius_x", "E_radius_y", "E_radius_z")], 1, mean))
    }
    else if (arrangement == "separated") {
      simulation_metadata$cluster_1$cluster_cell_types <- separated_cluster1_cell_type
      simulation_metadata$cluster_1$cluster_cell_proportions <- separated_cluster1_cell_prop
      simulation_metadata$cluster_1$centre_loc <- c(parameters_df$cluster1_x_coord[simulation_index], separated_cluster1_y_z_centre_loc)
      
      simulation_metadata <- spe_metadata_cluster_template("regular", "sphere", simulation_metadata)
      simulation_metadata$cluster_2$cluster_cell_types <- separated_cluster2_cell_type
      simulation_metadata$cluster_2$cluster_cell_proportions <- separated_cluster2_cell_prop
      simulation_metadata$cluster_2$centre_loc <- separated_cluster2_centre_loc
      simulation_metadata$cluster_2$radius <- separated_cluster2_radius
    }
    
    return(simulation_metadata)
  }
  
  # Function to get all slices from spe
  get_all_spe_slices <- function(spe, 
                                 bottom_z_coord_of_slices, 
                                 top_z_coord_of_slices) {
    
    spe_slices <- list()
    
    if (length(bottom_z_coord_of_slices) != length(top_z_coord_of_slices)) stop("Lengths of bottom_z_coords_of_slices and top_z_coords_of_slices should be equal.")
    number_of_slices <- length(bottom_z_coord_of_slices)
    
    for (i in seq(number_of_slices)) {
      bottom_z_coord_of_slice <- bottom_z_coord_of_slices[i]
      top_z_coord_of_slice <- top_z_coord_of_slices[i]
      z_coords_of_cells_in_spe <- spatialCoords(spe)[ , "Cell.Z.Position"]
      spe_for_slice <- spe[, bottom_z_coord_of_slice < z_coords_of_cells_in_spe & z_coords_of_cells_in_spe < top_z_coord_of_slice]
      spatialCoords(spe_for_slice) <- spatialCoords(spe_for_slice)[ , c("Cell.X.Position", "Cell.Y.Position")]
      
      spe_slices[[i]] <- spe_for_slice
    }
    
    return(spe_slices)
  }
  
  # Function to create empty metric df list
  create_empty_metric_df_list3D <- function(
    cell_types,
    n_simulations,
    radii_colnames,
    thresholds_colnames
  ) {
    n_cell_type_combinations <- length(cell_types)^2
    
    # Define AMD data frames as well as constants
    # n_slices + 1, to include the 3D value
    AMD_df_colnames <- c("simulation", "reference", "target", "AMD")
    AMD_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(AMD_df_colnames)))
    colnames(AMD_df) <- AMD_df_colnames
    
    # Define MS, NMS, ANC, ACIN, COO, ANE data frames as well as constants
    radii_colnames <- paste("r", radii, sep = "")
    
    MS_df_colnames <- c("simulation", "reference", "target", radii_colnames)
    MS_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(MS_df_colnames)))
    colnames(MS_df) <- MS_df_colnames
    
    NMS_df <- ANC_df <- ANE_df <- ACIN_df <- COO_df <- CK_df <- CL_df <- CG_df <- MS_df
    
    # Define SAC and prevalence data frames as well as constants
    thresholds_colnames <- paste("t", thresholds, sep = "")
    
    PBSAC_df_colnames <- c("simulation", "reference", "target", "PBSAC")
    PBSAC_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(PBSAC_df_colnames)))
    colnames(PBSAC_df) <- PBSAC_df_colnames
    
    PBP_df_colnames <- c("simulation", "reference", "target", thresholds_colnames)
    PBP_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(PBP_df_colnames)))
    colnames(PBP_df) <- PBP_df_colnames
    
    EBSAC_df_colnames <- c("simulation", "cell_types", "EBSAC")
    EBSAC_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(EBSAC_df_colnames)))
    colnames(EBSAC_df) <- EBSAC_df_colnames
    
    EBP_df_colnames <- c("simulation", "cell_types", thresholds_colnames)
    EBP_df <- data.frame(matrix(nrow = n_simulations * n_cell_type_combinations, ncol = length(EBP_df_colnames)))
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
  
  create_empty_metric_df_list2D <- function(
    cell_types,
    n_simulations,
    n_slices,
    radii_colnames,
    thresholds_colnames
  ) {
    n_cell_type_combinations <- length(cell_types)^2
    
    # Define AMD data frames as well as constants
    AMD_df_colnames <- c("simulation", "slice", "reference", "target", "AMD")
    AMD_df <- data.frame(matrix(nrow = n_simulations * n_slices * n_cell_type_combinations, ncol = length(AMD_df_colnames)))
    colnames(AMD_df) <- AMD_df_colnames
    
    # Define MS, NMS, ANC, ACIN, COO, ANE data frames as well as constants
    radii_colnames <- paste("r", radii, sep = "")
    
    MS_df_colnames <- c("simulation", "slice", "reference", "target", radii_colnames)
    MS_df <- data.frame(matrix(nrow = n_simulations * n_slices * n_cell_type_combinations, ncol = length(MS_df_colnames)))
    colnames(MS_df) <- MS_df_colnames
    
    NMS_df <- ANC_df <- ANE_df <- ACIN_df <- COO_df <- CK_df <- CL_df <- CG_df <- MS_df
    
    # Define SAC and prevalence data frames as well as constants
    thresholds_colnames <- paste("t", thresholds, sep = "")
    
    PBSAC_df_colnames <- c("simulation", "slice", "reference", "target", "PBSAC")
    PBSAC_df <- data.frame(matrix(nrow = n_simulations * n_slices * n_cell_type_combinations, ncol = length(PBSAC_df_colnames)))
    colnames(PBSAC_df) <- PBSAC_df_colnames
    
    PBP_df_colnames <- c("simulation", "slice", "reference", "target", thresholds_colnames)
    PBP_df <- data.frame(matrix(nrow = n_simulations * n_slices * n_cell_type_combinations, ncol = length(PBP_df_colnames)))
    colnames(PBP_df) <- PBP_df_colnames
    
    EBSAC_df_colnames <- c("simulation", "slice", "cell_types", "EBSAC")
    EBSAC_df <- data.frame(matrix(nrow = n_simulations * n_slices * n_cell_type_combinations, ncol = length(EBSAC_df_colnames)))
    colnames(EBSAC_df) <- EBSAC_df_colnames
    
    EBP_df_colnames <- c("simulation", "slice", "cell_types", thresholds_colnames)
    EBP_df <- data.frame(matrix(nrow = n_simulations * n_slices * n_cell_type_combinations, ncol = length(EBP_df_colnames)))
    colnames(EBP_df) <- EBP_df_colnames
    
    
    # Add all to list:
    metric_df_list <- list(AMD = AMD_df,
                           MS = MS_df,
                           NMS = NMS_df,
                           ACIN = ACIN_df,
                           ANE = ANE_df,
                           ANC = ANC_df,
                           CK = CK_df,
                           CL = CL_df,
                           CG = CG_df,
                           COO = COO_df,
                           PBSAC = PBSAC_df,
                           PBP = PBP_df,
                           EBSAC = EBSAC_df,
                           EBP = EBP_df)
    
    return(metric_df_list)
  }
  
  
  # Function analyse spe in 3D
  analyse_simulation3D <- function(spe, 
                                   cell_types, 
                                   radii, 
                                   thresholds, 
                                   n_splits, 
                                   simulation_index, 
                                   metric_df_list) {
    n_cell_type_combinations <- length(cell_types)^2
    radii_colnames <- paste("r", radii, sep = "")
    thresholds_colnames <- paste("t", thresholds, sep = "")
    
    index <- n_cell_type_combinations * (simulation_index - 1) + 1 
    
    minimum_distance_data <- calculate_minimum_distances_between_cell_types3D(spe,
                                                                              cell_types,
                                                                              show_summary = F,
                                                                              plot_image = F,
                                                                              feature_colname = "Cell.Type")
    
    minimum_distance_data_summary <- summarise_distances_between_cell_types3D(minimum_distance_data)
    
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "simulation"] <- simulation_index
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "reference"] <- minimum_distance_data_summary$reference
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "target"] <- minimum_distance_data_summary$target
    metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "AMD"] <- minimum_distance_data_summary$mean
    
    # Need a new index for 'pair' data which increments after each reference-target pair
    pair_index <- n_cell_type_combinations * (simulation_index - 1) + 1 
    
    for (reference_cell_type in cell_types) {
      gradient_data <- calculate_all_gradient_cc_metrics3D(spe,
                                                           reference_cell_type,
                                                           cell_types,
                                                           radii,
                                                           plot_image = F,
                                                           feature_colname = "Cell.Type")
      
      for (target_cell_type in cell_types) {
        print(paste(reference_cell_type, target_cell_type, sep = "/"))
        metric_df_list[["ANC"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["ACIN"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["COO"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["CK"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["CL"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["CG"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["MS"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["NMS"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["ANE"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, 
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
        metric_df_list[["PBSAC"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        metric_df_list[["PBP"]][pair_index, c("simulation", "reference", "target")] <- c(simulation_index, reference_cell_type, target_cell_type)
        
        metric_df_list[["EBSAC"]][pair_index, c("simulation", "cell_types")] <- c(simulation_index, paste(reference_cell_type, target_cell_type, sep = ","))
        metric_df_list[["EBP"]][pair_index, c("simulation", "cell_types")] <- c(simulation_index, paste(reference_cell_type, target_cell_type, sep = ","))
        
        if (reference_cell_type != target_cell_type) {
          proportion_grid_metrics <- calculate_cell_proportion_grid_metrics3D(spe, 
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
          
          entropy_grid_metrics <- calculate_entropy_grid_metrics3D(spe, 
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
    return(metric_df_list)
  }
  
  # Function analyse spe in 2D
  analyse_simulation2D <- function(spe_slices, 
                                   cell_types, 
                                   radii, 
                                   thresholds, 
                                   n_splits, 
                                   simulation_index, 
                                   metric_df_list) {
    n_cell_type_combinations <- length(cell_types)^2
    n_slices <- length(spe_slices)
    radii_colnames <- paste("r", radii, sep = "")
    thresholds_colnames <- paste("t", thresholds, sep = "")
    
    for (slice_index in seq_len(n_slices)) {
      spe <- spe_slices[[slice_index]]
      
      index <- n_cell_type_combinations * (simulation_index - 1) * n_slices +  n_cell_type_combinations * (slice_index - 1) + 1 
      
      minimum_distance_data <- calculate_minimum_distances_between_cell_types2D(spe,
                                                                                cell_types,
                                                                                show_summary = F,
                                                                                plot_image = F,
                                                                                feature_colname = "Cell.Type")
      
      minimum_distance_data_summary <- summarise_distances_between_cell_types2D(minimum_distance_data)
      
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "simulation"] <- simulation_index
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "slice"] <- slice_index
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "reference"] <- minimum_distance_data_summary$reference
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "target"] <- minimum_distance_data_summary$target
      metric_df_list[["AMD"]][index:(index + n_cell_type_combinations - 1), "AMD"] <- minimum_distance_data_summary$mean
      
      # Need a new index for 'pair' data which increments after each reference-target pair
      pair_index <- n_cell_type_combinations * (simulation_index - 1) * n_slices +  n_cell_type_combinations * (slice_index - 1) + 1 
      
      for (reference_cell_type in cell_types) {
        gradient_data <- calculate_all_gradient_cc_metrics2D(spe,
                                                             reference_cell_type,
                                                             cell_types,
                                                             radii,
                                                             plot_image = F,
                                                             feature_colname = "Cell.Type")
        
        for (target_cell_type in cell_types) {
          print(paste(reference_cell_type, target_cell_type, sep = "/"))
          metric_df_list[["ANC"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["ACIN"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["COO"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["CK"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["CL"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["CG"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["MS"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["NMS"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["ANE"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, 
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
          metric_df_list[["PBSAC"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          metric_df_list[["PBP"]][pair_index, c("simulation", "slice", "reference", "target")] <- c(simulation_index, slice_index, reference_cell_type, target_cell_type)
          
          metric_df_list[["EBSAC"]][pair_index, c("simulation", "slice", "cell_types")] <- c(simulation_index, slice_index, paste(reference_cell_type, target_cell_type, sep = ","))
          metric_df_list[["EBP"]][pair_index, c("simulation", "slice", "cell_types")] <- c(simulation_index, slice_index, paste(reference_cell_type, target_cell_type, sep = ","))
          
          if (reference_cell_type != target_cell_type) {
            proportion_grid_metrics <- calculate_cell_proportion_grid_metrics2D(spe, 
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
            
            entropy_grid_metrics <- calculate_entropy_grid_metrics2D(spe, 
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
    return(metric_df_list)
  }
  
  # Define metric df lists
  metric_df_list3D <- create_empty_metric_df_list3D(cell_types, n_simulations, radii_colnames, thresholds_colnames)
  metric_df_list2D <- create_empty_metric_df_list2D(cell_types, n_simulations, n_slices, radii_colnames, thresholds_colnames)
  
  for (simulation_index in seq_len(n_simulations)) {
    print(simulation_index)
    # Simulate spes in 3D, and get 2D slices
    simulation_metadata <- generate_simulation_metadata(parameters_df, simulation_index)
    
    spe3D <- simulate_spe_metadata3D(simulation_metadata, plot_image = F)
    
    spe_slices2D <- get_all_spe_slices(spe3D, bottom_z_coord_of_slices, top_z_coord_of_slices)
    
    # Analyse spes
    metric_df_list3D <- analyse_simulation3D(spe3D,
                                             cell_types,
                                             radii,
                                             thresholds,
                                             n_splits,
                                             simulation_index,
                                             metric_df_list3D)
    
    metric_df_list2D <- analyse_simulation2D(spe_slices2D,
                                             cell_types,
                                             radii,
                                             thresholds,
                                             n_splits,
                                             simulation_index,
                                             metric_df_list2D)
  }
  
  ## Turn gradient metrics into AUC and add to metric_df list
  get_AUC_for_radii_gradient_metrics <- function(y) {
    x <- radii
    h <- diff(x)[1]
    n <- length(x)
    
    AUC <- (h / 2) * (y[1] + 2 * sum(y[2:(n - 1)]) + y[n])
    
    return(AUC)
  }
  
  add_AUC_for_radii_gradient_metrics_to_metric_df_list3D <- function(metric_df_list) {
    gradient_radii_metrics <- c("MS", "NMS", "ACIN", "ANE", "ANC", "COO", "CK", "CL", "CG")
    
    for (metric in gradient_radii_metrics) {
      metric_AUC_name <- paste(metric, "AUC", sep = "_")
      
      if (metric %in% c("MS", "NMS", "ANC", "ACIN", "ANE", "COO", "CK", "CL", "CG")) {
        subset_colnames <- c("simulation", "reference", "target", metric_AUC_name)
      }
      else {
        stop("Unknown metric?")
      }
      
      df <- metric_df_list[[metric]]
      df[[metric_AUC_name]] <- apply(df[ , radii_colnames], 1, get_AUC_for_radii_gradient_metrics)
      df <- df[ , subset_colnames]
      metric_df_list[[metric_AUC_name]] <- df
    }
    
    # PBP_AUC 3D
    PBP_df <- metric_df_list[["PBP"]]
    PBP_df$PBP_AUC <- apply(PBP_df[ , thresholds_colnames], 1, function(y) {
      sum(diff(thresholds) * (head(y, -1) + tail(y, -1)) / 2)
    })
    PBP_AUC_df <- PBP_df[ , c("simulation", "reference", "target", "PBP_AUC")]
    metric_df_list[["PBP_AUC"]] <- PBP_AUC_df
    
    
    # EBP_AUC 3D
    EBP_df <- metric_df_list[["EBP"]]
    EBP_df$EBP_AUC <- apply(EBP_df[ , thresholds_colnames], 1, function(y) {
      sum(diff(thresholds) * (head(y, -1) + tail(y, -1)) / 2)
    })
    EBP_AUC_df <- EBP_df[ , c("simulation", "cell_types", "EBP_AUC")]
    metric_df_list[["EBP_AUC"]] <- EBP_AUC_df
    
    return(metric_df_list)
  }
  add_AUC_for_radii_gradient_metrics_to_metric_df_list2D <- function(metric_df_list) {
    gradient_radii_metrics <- c("MS", "NMS", "ACIN", "ANE", "ANC", "COO", "CK", "CL", "CG")
    
    for (metric in gradient_radii_metrics) {
      metric_AUC_name <- paste(metric, "AUC", sep = "_")
      
      if (metric %in% c("MS", "NMS", "ANC", "ACIN", "ANE", "COO", "CK", "CL", "CG")) {
        subset_colnames <- c("simulation", "slice", "reference", "target", metric_AUC_name)
      }
      else {
        stop("Unknown metric?")
      }
      
      df <- metric_df_list[[metric]]
      df[[metric_AUC_name]] <- apply(df[ , radii_colnames], 1, get_AUC_for_radii_gradient_metrics)
      df <- df[ , subset_colnames]
      metric_df_list[[metric_AUC_name]] <- df
    }
    
    # PBP_AUC 3D
    PBP_df <- metric_df_list[["PBP"]]
    PBP_df$PBP_AUC <- apply(PBP_df[ , thresholds_colnames], 1, function(y) {
      sum(diff(thresholds) * (head(y, -1) + tail(y, -1)) / 2)
    })
    PBP_AUC_df <- PBP_df[ , c("simulation", "slice", "reference", "target", "PBP_AUC")]
    metric_df_list[["PBP_AUC"]] <- PBP_AUC_df
    
    # EBP_AUC 3D
    EBP_df <- metric_df_list[["EBP"]]
    EBP_df$EBP_AUC <- apply(EBP_df[ , thresholds_colnames], 1, function(y) {
      sum(diff(thresholds) * (head(y, -1) + tail(y, -1)) / 2)
    })
    EBP_AUC_df <- EBP_df[ , c("simulation", "slice", "cell_types", "EBP_AUC")]
    metric_df_list[["EBP_AUC"]] <- EBP_AUC_df
    
    return(metric_df_list)
  }
  
  metric_df_list3D <- add_AUC_for_radii_gradient_metrics_to_metric_df_list3D(metric_df_list3D)
  metric_df_list2D <- add_AUC_for_radii_gradient_metrics_to_metric_df_list2D(metric_df_list2D)
  
  # Combine 3D and 2D metric df lists
  metric_df_list_combined <- list()
  
  metrics <- names(metric_df_list3D)
  for (metric in metrics) {
    metric_df_list3D[[metric]][["slice"]] <- 0 # Set slice index for 3D to be 0
    
    metric_df_combined <- rbind(metric_df_list3D[[metric]], metric_df_list2D[[metric]])
    
    metric_df_list_combined[[metric]] <- metric_df_combined
  }
  
  return(metric_df_list_combined)
}


# Running the functions ----
parameters_df <- generate_random_parameters(
  n_simulations = 10000,
  bg_prop_A_range = c("min" = 0, "max" = 0.10),
  bg_prop_B_range = c("min" = 0, "max" = 0.10),
  E_radius_x_range = c("min" = 75, "max" = 125),
  E_radius_y_range = c("min" = 75, "max" = 125),
  E_radius_z_range = c("min" = 75, "max" = 125),
  N_width_range = c("min" = 25, "max" = 35),
  cluster_prop_A_range = c("min" = 0.3, "max" = 0.7),
  ring_width_factor_range = c("min" = 0.1, "max" = 0.2) ,
  cluster1_x_coord_range = c("min" = 125, "max" = 175)
)

SD1_metric_df_list <- analyse_SD1_simulations(parameters_df)


