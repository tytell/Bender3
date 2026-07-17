# =============================================================================
# Muscle Mechanics Analysis Script for Bender CritterGripper Data
# Force-Velocity, Force-Length, and Power Analysis Pipeline
# =============================================================================
#
# This script analyzes muscle mechanics data from Bender HDF5 files to extract:
# - Force-Velocity (FV) properties and curves
# - Force-Length (FL) properties and curves  
# - Power output properties and curves
#
# Data Sources:
# - Input: HDF5 files from Bender CritterGripper system
# - Schema: Uses flat metadata/ + timeseries/ structure per context_jlab_cg_h5schema.md
#
# Output:
# - Individual trial plots (3-panel time series)
# - Summary plots for each protocol type
# - Extracted muscle properties (FV curves, FL curves, power curves)
#
# Dependencies: tidyverse, ggplot2, broom, fs, rhdf5, pracma, signal, patchwork

library(tidyverse)
library(ggplot2) 
library(broom)
library(fs)
library(rhdf5)
library(pracma)
library(patchwork)
library(cli)

# Source existing analysis functions
source("R/00_load_bender_flat.R")
source("R/03_analyze.R")

# =============================================================================
# Configuration and Constants
# =============================================================================

# Data directories
SOURCE_DATA_DIR <- "/Users/yjimenez/Library/CloudStorage/OneDrive-ProvidenceCollege/01_JimenezLab/01_PermanentArchive/bender_crittergripper/2026-07-14_bass16_bender"
OUTPUT_FIGURES_DIR <- "/Users/yjimenez/Library/CloudStorage/OneDrive-ProvidenceCollege/01_JimenezLab/02_ResearchHub/proj_crittergripper/figures/bass16_figures"

# Analysis constants
DEACTIVATION_WINDOW_SEC <- 0.5  # Continue passive subtraction for 0.5s after stimulus
STRAIN_CALCULATION_PARAMS <- list(
  half_thickness_default_mm = 5.0,  # Default if not in metadata
  muscle_depth_default_mm = 2.0     # Default if not in metadata  
)

# Create output directory
dir_create(OUTPUT_FIGURES_DIR, recurse = TRUE)

# =============================================================================
# Enhanced Data Loading Functions
# =============================================================================

#' Load and preprocess Bender data with muscle-specific enhancements
load_bender_enhanced <- function(filepath, protocol_type = NULL) {
  cli_alert_info("Loading {basename(filepath)}...")
  
  # Load using existing flat loader
  td <- load_bender_flat(filepath, loadtorques = c("x", "y", "z"))
  if (is.null(td)) {
    cli_alert_danger("Failed to load {basename(filepath)}")
    return(NULL)
  }
  
  # Extract additional metadata for muscle analysis
  h5 <- H5Fopen(filepath, "H5F_ACC_RDONLY")
  on.exit(try(H5Fclose(h5), silent = TRUE), add = TRUE)
  
  m_attrs <- tryCatch(h5readAttributes(h5, "/metadata"), error = function(e) list())
  m_a <- function(key) m_attrs[[key]]
  
  # Get muscle geometry parameters
  half_thickness_mm <- coalesce(
    m_a("measurement_specimen_local_body_height_millimeter"),
    STRAIN_CALCULATION_PARAMS$half_thickness_default_mm
  ) / 2  # Convert height to half-thickness
  
  muscle_depth_mm <- coalesce(
    m_a("measurement_target_muscle_depth_millimeter"),
    STRAIN_CALCULATION_PARAMS$muscle_depth_default_mm
  )
  
  # Add muscle-specific calculations
  td <- td %>%
    mutate(
      # Calculate predicted muscle strain based on curvature
      muscle_strain_predicted = ifelse(
        is.finite(curve.invm) & is.finite(muscle_depth_mm),
        curve.invm * muscle_depth_mm / 1000,  # Convert mm to m
        NA_real_
      ),
      
      # Add muscle geometry metadata
      half_thickness_mm = half_thickness_mm,
      muscle_depth_mm = muscle_depth_mm,
      
      # Protocol identification
      protocol_type = coalesce(protocol_type, unique(stimulus_type)[1])
    )
  
  # Add muscle mass estimates
  td <- add_muscle_mass(td)
  
  # Set cycle types for muscle analysis
  td <- set_cycle_types(td)
  
  # Estimate body torque for muscle analysis
  td <- estimate_body_torque(td, method = "xtorque")
  
  cli_alert_success("Loaded {nrow(td)} samples from {basename(filepath)}")
  return(td)
}

#' Batch load all files from data directory
load_all_trials <- function(data_dir = SOURCE_DATA_DIR, pattern = "\\.h5$") {
  cli_h1("Loading All Trial Data")
  
  files <- dir_ls(data_dir, regexp = pattern)
  if (length(files) == 0) {
    cli_alert_danger("No .h5 files found in {data_dir}")
    return(tibble())
  }
  
  cli_alert_info("Found {length(files)} HDF5 files")
  
  # Extract protocol types from filenames
  protocols <- str_extract(basename(files), "(?<=_)(isometric|isovelocity|dynamic|FL|FV)(?=\\.h5$)")
  
  all_data <- map2_dfr(files, protocols, ~{
    data <- load_bender_enhanced(.x, .y)
    if (!is.null(data)) {
      data$filepath <- .x
      data$protocol_filename = .y
    }
    data
  }, .id = "file_id")
  
  cli_alert_success("Loaded {nrow(all_data)} total samples from {length(files)} files")
  return(all_data)
}

# =============================================================================
# Muscle Force Calculation Functions
# =============================================================================

#' Calculate net muscle force with passive subtraction
calculate_muscle_force <- function(data, protocol_type = "dynamic") {
  cli_alert_info("Calculating muscle forces for {protocol_type} protocol...")
  
  # Ensure we have required columns
  if (!any(c("xtorque.Nm", "bodytorque.Nm") %in% names(data))) {
    cli_alert_danger("No torque columns found for muscle force calculation")
    return(data)
  }
  
  # Select primary torque column
  torque_col <- intersect(c("xtorque.Nm", "bodytorque.Nm"), names(data))[1]
  
  # Protocol-specific passive force calculation
  if (protocol_type == "dynamic") {
    # Dynamic: strictly passive loops for baseline
    data <- calculate_dynamic_muscle_force(data, torque_col)
  } else if (protocol_type == "isometric") {
    # Isometric: passive at-length hold prior to stimulus  
    data <- calculate_isometric_muscle_force(data, torque_col)
  } else if (protocol_type %in% c("isovelocity", "FV")) {
    # Isovelocity: passive blocks for baseline
    data <- calculate_isovelocity_muscle_force(data, torque_col)
  } else {
    # Generic approach
    data <- calculate_generic_muscle_force(data, torque_col)
  }
  
  return(data)
}

#' Dynamic protocol muscle force calculation
calculate_dynamic_muscle_force <- function(data, torque_col) {
  # For dynamic protocols, use passive cycles as baseline
  # Check if cycletype column exists
  if ("cycletype" %in% names(data)) {
    passive_baseline <- data %>%
      filter(cycletype == "pass" | (is.na(cycletype) & stim == "0")) %>%
      group_by(across(any_of(c("freq.Hz", "curvature.invm", "phase", "duty")))) %>%
      summarise(passive_torque_baseline = mean(.data[[torque_col]], na.rm = TRUE), .groups = "drop")
  } else {
    # Fallback: use stim == "0" periods
    passive_baseline <- data %>%
      filter(stim == "0") %>%
      group_by(across(any_of(c("freq.Hz", "curvature.invm", "phase", "duty")))) %>%
      summarise(passive_torque_baseline = mean(.data[[torque_col]], na.rm = TRUE), .groups = "drop")
  }
  
  # Join and calculate muscle force
  data %>%
    left_join(passive_baseline, by = intersect(names(.), names(passive_baseline))) %>%
    mutate(
      passive_torque_baseline = coalesce(passive_torque_baseline, 0),
      muscle_force_active = .data[[torque_col]] - passive_torque_baseline,
      
      # Continue passive subtraction for deactivation window
      time_since_stim_end = case_when(
        stim != "0" ~ 0,  # During stimulation
        TRUE ~ t.s - lag(ifelse(stim != "0", t.s, NA_real_), default = -Inf)
      ),
      
      muscle_force = case_when(
        stim != "0" ~ muscle_force_active,  # Active period
        !is.na(time_since_stim_end) & time_since_stim_end <= DEACTIVATION_WINDOW_SEC ~ muscle_force_active,  # Deactivation window
        TRUE ~ 0  # Pure passive period
      )
    )
}

#' Isometric protocol muscle force calculation  
calculate_isometric_muscle_force <- function(data, torque_col) {
  # Use pre-stimulus period as passive baseline
  baseline <- data %>%
    filter(t.s < -0.1) %>%  # Pre-stimulus period
    group_by(step_number) %>%
    summarise(passive_torque_baseline = mean(.data[[torque_col]], na.rm = TRUE), .groups = "drop")
  
  data %>%
    left_join(baseline, by = "step_number") %>%
    mutate(
      passive_torque_baseline = coalesce(passive_torque_baseline, 0),
      muscle_force = .data[[torque_col]] - passive_torque_baseline
    )
}

#' Isovelocity protocol muscle force calculation
calculate_isovelocity_muscle_force <- function(data, torque_col) {
  # Use designated passive blocks
  passive_baseline <- data %>%
    filter(stim == "0") %>%
    group_by(across(any_of(c("step_number", "velocity_deg_s")))) %>%
    summarise(passive_torque_baseline = mean(.data[[torque_col]], na.rm = TRUE), .groups = "drop")
  
  data %>%
    left_join(passive_baseline, by = intersect(names(.), names(passive_baseline))) %>%
    mutate(
      passive_torque_baseline = coalesce(passive_torque_baseline, 0),
      muscle_force = .data[[torque_col]] - passive_torque_baseline
    )
}

#' Generic muscle force calculation fallback
calculate_generic_muscle_force <- function(data, torque_col) {
  # Simple pre-stimulus baseline approach
  baseline_value <- data %>%
    filter(t.s < 0) %>%
    pull(.data[[torque_col]]) %>%
    mean(na.rm = TRUE)
  
  if (is.finite(baseline_value)) {
    data %>% mutate(muscle_force = .data[[torque_col]] - baseline_value)
  } else {
    data %>% mutate(muscle_force = .data[[torque_col]])
  }
}

# =============================================================================
# Individual Trial Visualization Functions
# =============================================================================

#' Create publication-quality compound plot for individual trials
plot_individual_trial <- function(data, output_dir = OUTPUT_FIGURES_DIR) {
  # Extract trial information
  trial_info <- data %>%
    slice_head(n = 1) %>%
    select(fishcode, trial, protocol_type, filename) %>%
    mutate(plot_title = glue::glue("{fishcode} - Trial {trial} ({protocol_type})"))
  
  # Filter to reasonable time window for plotting
  plot_data <- data %>%
    filter(t.s >= -1, t.s <= max(t.s, na.rm = TRUE)) %>%
    mutate(
      time_sec = t.s,
      # Convert torque to more readable units if needed
      torque_display = coalesce(bodytorque.Nm, xtorque.Nm, ytorque.Nm, ztorque.Nm, 0),
      muscle_force_display = coalesce(muscle_force, 0),
      strain_display = coalesce(muscle_strain_predicted, 0)
    )
  
  # Panel 1: Inertial-corrected torques vs time
  p1 <- plot_data %>%
    ggplot(aes(x = time_sec)) +
    geom_line(aes(y = torque_display), color = "blue", linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.7) +
    labs(
      x = "Time (s)",
      y = "Torque (N⋅m)",
      subtitle = "A. Inertial-corrected torque"
    ) +
    theme_minimal() +
    theme(
      plot.subtitle = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  # Panel 2: Calculated muscle forces vs time
  p2 <- plot_data %>%
    ggplot(aes(x = time_sec)) +
    geom_line(aes(y = muscle_force_display), color = "red", linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "solid", alpha = 0.5, color = "gray") +
    labs(
      x = "Time (s)", 
      y = "Muscle Force (N⋅m)",
      subtitle = "B. Calculated muscle force (active - passive)"
    ) +
    theme_minimal() +
    theme(
      plot.subtitle = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  # Panel 3: Predicted muscle strain vs time
  p3 <- plot_data %>%
    ggplot(aes(x = time_sec, y = strain_display)) +
    geom_line(color = "darkgreen", linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "solid", alpha = 0.5, color = "gray") +
    labs(
      x = "Time (s)",
      y = "Muscle Strain",
      subtitle = "C. Predicted muscle strain"
    ) +
    theme_minimal() +
    theme(
      plot.subtitle = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  # Combine panels
  combined_plot <- (p1 / p2 / p3) + 
    plot_annotation(
      title = trial_info$plot_title,
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  
  # Save plot
  output_filename <- glue::glue("{trial_info$fishcode}_trial_{trial_info$trial}_{trial_info$protocol_type}.png")
  output_path <- file.path(output_dir, output_filename)
  
  ggsave(
    output_path, 
    combined_plot, 
    width = 10, height = 8, 
    dpi = 300, 
    bg = "white"
  )
  
  cli_alert_success("Saved individual trial plot: {output_filename}")
  return(combined_plot)
}

# =============================================================================
# Force-Velocity Analysis Functions
# =============================================================================

#' Extract Force-Velocity properties from isovelocity data
analyze_force_velocity <- function(data) {
  cli_h2("Analyzing Force-Velocity Properties")
  
  # Filter to isovelocity protocols
  fv_data <- data %>%
    filter(protocol_type %in% c("isovelocity", "FV")) %>%
    filter(!is.na(muscle_force))
  
  if (nrow(fv_data) == 0) {
    cli_alert_warning("No isovelocity data found for FV analysis")
    return(tibble())
  }
  
  # Calculate FV properties per step
  fv_summary <- fv_data %>%
    group_by(fishcode, trial, step_number, protocol_type) %>%
    summarise(
      velocity_deg_s = first(anglevel.degps),  # Use actual angular velocity
      muscle_force_mean = mean(muscle_force, na.rm = TRUE),
      muscle_force_peak = max(abs(muscle_force), na.rm = TRUE),
      muscle_force_se = sd(muscle_force, na.rm = TRUE) / sqrt(n()),
      n_samples = n(),
      .groups = "drop"
    ) %>%
    filter(is.finite(velocity_deg_s), is.finite(muscle_force_mean))
  
  # Fit FV curve (Hill equation: F = F0 * (1 - v/vmax))
  if (nrow(fv_summary) > 3) {
    fv_fit <- tryCatch({
      # Simple linear approximation first
      lm_fit <- lm(muscle_force_mean ~ velocity_deg_s, data = fv_summary)
      
      tibble(
        F0_estimate = coef(lm_fit)[1],  # Force at zero velocity
        slope = coef(lm_fit)[2],        # Force-velocity slope
        r_squared = summary(lm_fit)$r.squared,
        p_value = summary(lm_fit)$coefficients[2, 4],
        n_points = nrow(fv_summary)
      )
    }, error = function(e) {
      cli_alert_warning("FV curve fitting failed: {e$message}")
      tibble(
        F0_estimate = NA_real_, slope = NA_real_, 
        r_squared = NA_real_, p_value = NA_real_, n_points = nrow(fv_summary)
      )
    })
  } else {
    fv_fit <- tibble(
      F0_estimate = NA_real_, slope = NA_real_, 
      r_squared = NA_real_, p_value = NA_real_, n_points = nrow(fv_summary)
    )
  }
  
  return(list(
    raw_data = fv_summary,
    curve_fit = fv_fit
  ))
}

# =============================================================================
# Force-Length Analysis Functions  
# =============================================================================

#' Extract Force-Length properties from isometric data
analyze_force_length <- function(data) {
  cli_h2("Analyzing Force-Length Properties")
  
  # Filter to isometric protocols
  fl_data <- data %>%
    filter(protocol_type %in% c("isometric", "FL")) %>%
    filter(!is.na(muscle_force))
  
  if (nrow(fl_data) == 0) {
    cli_alert_warning("No isometric data found for FL analysis")
    return(tibble())
  }
  
  # Calculate FL properties per step
  fl_summary <- fl_data %>%
    group_by(fishcode, trial, step_number, protocol_type) %>%
    summarise(
      target_angle_deg = first(angle.deg),  # Use actual angle as proxy
      muscle_length_proxy = first(angle.deg),  # Use actual angle as muscle length proxy
      muscle_force_mean = mean(muscle_force, na.rm = TRUE),
      muscle_force_peak = max(abs(muscle_force), na.rm = TRUE), 
      muscle_force_se = sd(muscle_force, na.rm = TRUE) / sqrt(n()),
      n_samples = n(),
      .groups = "drop"
    ) %>%
    filter(is.finite(muscle_length_proxy), is.finite(muscle_force_mean))
  
  # Fit FL curve (parabolic approximation)
  if (nrow(fl_summary) > 3) {
    fl_fit <- tryCatch({
      # Quadratic fit: F = a + b*L + c*L^2
      quad_fit <- lm(muscle_force_mean ~ muscle_length_proxy + I(muscle_length_proxy^2), 
                     data = fl_summary)
      
      # Find optimal length (maximum force)
      coeffs <- coef(quad_fit)
      optimal_length <- -coeffs[2] / (2 * coeffs[3])
      max_force <- predict(quad_fit, newdata = data.frame(muscle_length_proxy = optimal_length))
      
      tibble(
        F0_estimate = max_force,
        optimal_length = optimal_length,
        quadratic_coef = coeffs[3],
        r_squared = summary(quad_fit)$r.squared,
        n_points = nrow(fl_summary)
      )
    }, error = function(e) {
      cli_alert_warning("FL curve fitting failed: {e$message}")
      tibble(
        F0_estimate = NA_real_, optimal_length = NA_real_, 
        quadratic_coef = NA_real_, r_squared = NA_real_, n_points = nrow(fl_summary)
      )
    })
  } else {
    fl_fit <- tibble(
      F0_estimate = NA_real_, optimal_length = NA_real_, 
      quadratic_coef = NA_real_, r_squared = NA_real_, n_points = nrow(fl_summary)
    )
  }
  
  return(list(
    raw_data = fl_summary,
    curve_fit = fl_fit
  ))
}

# =============================================================================
# Power Analysis Functions
# =============================================================================

#' Calculate power output properties
analyze_power_output <- function(data) {
  cli_h2("Analyzing Power Output Properties")
  
  # Calculate instantaneous power
  power_data <- data %>%
    filter(!is.na(muscle_force), !is.na(angle.deg)) %>%
    arrange(fishcode, trial, t.s) %>%
    group_by(fishcode, trial) %>%
    mutate(
      # Calculate angular velocity (rad/s)
      angular_velocity_rad_s = c(NA, diff(angle.deg * pi/180) / diff(t.s)),
      
      # Calculate instantaneous power (W = torque * angular_velocity)
      instantaneous_power_W = muscle_force * angular_velocity_rad_s,
      
      # Mass-normalized power
      mass_specific_power_W_kg = if_else(
        is.finite(muscle_mass.kg) & muscle_mass.kg > 0,
        instantaneous_power_W / muscle_mass.kg,
        NA_real_
      )
    ) %>%
    ungroup()
  
  # Summarize power by cycles/steps
  if ("cycle" %in% names(power_data) && any(power_data$cycle > 0, na.rm = TRUE)) {
    # Dynamic protocols - by cycle
    power_summary <- power_data %>%
      filter(cycle > 0, !is.na(instantaneous_power_W)) %>%
      group_by(fishcode, trial, cycle, protocol_type) %>%
      summarise(
        mean_power_W = mean(instantaneous_power_W, na.rm = TRUE),
        peak_power_W = max(abs(instantaneous_power_W), na.rm = TRUE),
        mean_power_W_kg = mean(mass_specific_power_W_kg, na.rm = TRUE),
        peak_power_W_kg = max(abs(mass_specific_power_W_kg), na.rm = TRUE),
        frequency_Hz = first(freq.Hz),
        curvature_m_inv = first(curvature.invm),
        n_samples = n(),
        .groups = "drop"
      )
  } else if ("step_number" %in% names(power_data)) {
    # Step protocols
    power_summary <- power_data %>%
      filter(!is.na(step_number), !is.na(instantaneous_power_W)) %>%
      group_by(fishcode, trial, step_number, protocol_type) %>%
      summarise(
        mean_power_W = mean(instantaneous_power_W, na.rm = TRUE),
        peak_power_W = max(abs(instantaneous_power_W), na.rm = TRUE),
        mean_power_W_kg = mean(mass_specific_power_W_kg, na.rm = TRUE),
        peak_power_W_kg = max(abs(mass_specific_power_W_kg), na.rm = TRUE),
        n_samples = n(),
        .groups = "drop"
      )
  } else {
    # Fallback - by trial
    power_summary <- power_data %>%
      filter(!is.na(instantaneous_power_W)) %>%
      group_by(fishcode, trial, protocol_type) %>%
      summarise(
        mean_power_W = mean(instantaneous_power_W, na.rm = TRUE),
        peak_power_W = max(abs(instantaneous_power_W), na.rm = TRUE),
        mean_power_W_kg = mean(mass_specific_power_W_kg, na.rm = TRUE),
        peak_power_W_kg = max(abs(mass_specific_power_W_kg), na.rm = TRUE),
        n_samples = n(),
        .groups = "drop"
      )
  }
  
  return(list(
    raw_data = power_data,
    summary = power_summary
  ))
}

# =============================================================================
# Summary Visualization Functions
# =============================================================================

#' Create master summary plots for each protocol type
create_master_summary_plots <- function(data, output_dir = OUTPUT_FIGURES_DIR) {
  cli_h2("Creating Master Summary Plots")
  
  protocols <- unique(data$protocol_type)
  protocols <- protocols[!is.na(protocols)]
  
  for (protocol in protocols) {
    cli_alert_info("Creating summary plot for {protocol} protocol...")
    
    protocol_data <- data %>% filter(protocol_type == protocol)
    
    if (protocol %in% c("dynamic", "frequency_sweep")) {
      p <- create_dynamic_summary_plot(protocol_data, protocol)
    } else if (protocol %in% c("isometric", "FL")) {
      p <- create_isometric_summary_plot(protocol_data, protocol)
    } else if (protocol %in% c("isovelocity", "FV")) {
      p <- create_isovelocity_summary_plot(protocol_data, protocol)
    } else {
      cli_alert_warning("Unknown protocol type: {protocol}")
      next
    }
    
    # Save plot
    output_filename <- glue::glue("bass16_{tolower(protocol)}_summary.png")
    output_path <- file.path(output_dir, output_filename)
    
    ggsave(
      output_path, 
      p, 
      width = 12, height = 8, 
      dpi = 300, 
      bg = "white"
    )
    
    cli_alert_success("Saved summary plot: {output_filename}")
  }
}

#' Create summary plot for dynamic protocols
create_dynamic_summary_plot <- function(data, protocol) {
  # Calculate summary statistics
  summary_stats <- data %>%
    filter(!is.na(muscle_force)) %>%
    group_by(t.s) %>%
    summarise(
      mean_force = mean(muscle_force, na.rm = TRUE),
      se_force = sd(muscle_force, na.rm = TRUE) / sqrt(n()),
      lower_se = mean_force - se_force,
      upper_se = mean_force + se_force,
      n_trials = n(),
      .groups = "drop"
    ) %>%
    filter(between(t.s, -1, 5))  # Focus on relevant time window
  
  # Create plot
  ggplot() +
    # Individual data points (semi-transparent)
    geom_line(
      data = data %>% filter(between(t.s, -1, 5)),
      aes(x = t.s, y = muscle_force, group = interaction(trial, fishcode)),
      alpha = 0.3, color = "gray50", linewidth = 0.3
    ) +
    # Mean curve with SE ribbon
    geom_ribbon(
      data = summary_stats,
      aes(x = t.s, ymin = lower_se, ymax = upper_se),
      alpha = 0.3, fill = "blue"
    ) +
    geom_line(
      data = summary_stats,
      aes(x = t.s, y = mean_force),
      color = "blue", linewidth = 1.2
    ) +
    # Reference lines
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "solid", alpha = 0.5, color = "gray") +
    labs(
      title = glue::glue("Bass16 {str_to_title(protocol)} Protocol Summary"),
      subtitle = "Individual trials (gray) with mean ± SE (blue)",
      x = "Time (s)",
      y = "Muscle Force (N⋅m)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      panel.grid.minor = element_blank()
    )
}

#' Create summary plot for isometric protocols
create_isometric_summary_plot <- function(data, protocol) {
  # Analyze FL properties
  fl_results <- analyze_force_length(data)
  
  if (nrow(fl_results$raw_data) == 0) {
    return(ggplot() + 
           labs(title = glue::glue("No {protocol} data available")) + 
           theme_minimal())
  }
  
  # Create FL curve plot
  p1 <- ggplot(fl_results$raw_data, aes(x = muscle_length_proxy, y = muscle_force_mean)) +
    geom_point(size = 3, alpha = 0.7, color = "blue") +
    geom_errorbar(
      aes(ymin = muscle_force_mean - muscle_force_se, 
          ymax = muscle_force_mean + muscle_force_se),
      width = 0.1, alpha = 0.7, color = "blue"
    ) +
    labs(
      title = glue::glue("Bass16 {str_to_title(protocol)} - Force-Length Relationship"),
      x = "Muscle Length Proxy (deg)",
      y = "Mean Muscle Force (N⋅m)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  # Add fitted curve if available
  if (!is.na(fl_results$curve_fit$F0_estimate)) {
    # Create prediction data
    length_range <- range(fl_results$raw_data$muscle_length_proxy, na.rm = TRUE)
    pred_data <- tibble(
      muscle_length_proxy = seq(length_range[1], length_range[2], length.out = 100)
    )
    
    # Add fitted curve (this would need the actual model - simplified here)
    p1 <- p1 + 
      annotate(
        "text", 
        x = Inf, y = Inf, 
        label = glue::glue("F₀ = {round(fl_results$curve_fit$F0_estimate, 2)} N⋅m\nR² = {round(fl_results$curve_fit$r_squared, 3)}"),
        hjust = 1.1, vjust = 1.1,
        size = 3.5
      )
  }
  
  return(p1)
}

#' Create summary plot for isovelocity protocols  
create_isovelocity_summary_plot <- function(data, protocol) {
  # Analyze FV properties
  fv_results <- analyze_force_velocity(data)
  
  if (nrow(fv_results$raw_data) == 0) {
    return(ggplot() + 
           labs(title = glue::glue("No {protocol} data available")) + 
           theme_minimal())
  }
  
  # Create FV curve plot
  p1 <- ggplot(fv_results$raw_data, aes(x = velocity_deg_s, y = muscle_force_mean)) +
    geom_point(size = 3, alpha = 0.7, color = "red") +
    geom_errorbar(
      aes(ymin = muscle_force_mean - muscle_force_se, 
          ymax = muscle_force_mean + muscle_force_se),
      width = 0, alpha = 0.7, color = "red"
    ) +
    labs(
      title = glue::glue("Bass16 {str_to_title(protocol)} - Force-Velocity Relationship"),
      x = "Velocity (deg/s)",
      y = "Mean Muscle Force (N⋅m)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  # Add fitted line if available
  if (!is.na(fv_results$curve_fit$F0_estimate)) {
    p1 <- p1 + 
      geom_smooth(method = "lm", se = TRUE, color = "darkred", alpha = 0.3) +
      annotate(
        "text", 
        x = Inf, y = Inf, 
        label = glue::glue("F₀ = {round(fv_results$curve_fit$F0_estimate, 2)} N⋅m\nSlope = {round(fv_results$curve_fit$slope, 4)}\nR² = {round(fv_results$curve_fit$r_squared, 3)}"),
        hjust = 1.1, vjust = 1.1,
        size = 3.5
      )
  }
  
  return(p1)
}

# =============================================================================
# Main Analysis Pipeline
# =============================================================================

#' Main analysis pipeline function
run_muscle_mechanics_analysis <- function(data_dir = SOURCE_DATA_DIR, 
                                        output_dir = OUTPUT_FIGURES_DIR,
                                        create_individual_plots = TRUE,
                                        create_summary_plots = TRUE) {
  
  cli_h1("Muscle Mechanics Analysis Pipeline")
  cli_alert_info("Data source: {data_dir}")
  cli_alert_info("Output directory: {output_dir}")
  
  # Step 1: Load all data
  all_data <- load_all_trials(data_dir)
  
  if (nrow(all_data) == 0) {
    cli_alert_danger("No data loaded. Exiting.")
    return(invisible(NULL))
  }
  
  # Step 2: Process each protocol type
  protocols <- unique(all_data$protocol_type)
  protocols <- protocols[!is.na(protocols)]
  
  results <- list()
  processed_data <- list()  # Store processed data for summary plots
  
  for (protocol in protocols) {
    cli_h2("Processing {protocol} Protocol")
    
    protocol_data <- all_data %>% 
      filter(protocol_type == protocol) %>%
      calculate_muscle_force(protocol)
    
    # Store processed data
    processed_data[[protocol]] <- protocol_data
    
    # Step 3: Create individual trial plots
    if (create_individual_plots) {
      cli_alert_info("Creating individual trial plots...")
      
      trial_groups <- protocol_data %>%
        group_by(fishcode, trial) %>%
        group_split()
      
      walk(trial_groups, plot_individual_trial, output_dir)
    }
    
    # Step 4: Analyze muscle properties by protocol
    if (protocol %in% c("isometric", "FL")) {
      results[[protocol]] <- analyze_force_length(protocol_data)
    } else if (protocol %in% c("isovelocity", "FV")) {
      results[[protocol]] <- analyze_force_velocity(protocol_data)
    }
    
    # Step 5: Analyze power output
    power_results <- analyze_power_output(protocol_data)
    results[[paste0(protocol, "_power")]] <- power_results
  }
  
  # Step 6: Create summary plots
  if (create_summary_plots) {
    # Combine all processed data for summary plots
    combined_processed_data <- bind_rows(processed_data)
    create_master_summary_plots(combined_processed_data, output_dir)
  }
  
  # Step 7: Generate analysis summary
  cli_h2("Analysis Summary")
  cli_alert_success("Analysis completed!")
  cli_alert_info("Processed {length(protocols)} protocol types:")
  walk(protocols, ~cli_alert_info("  - {.x}"))
  cli_alert_info("Results saved to: {output_dir}")
  
  return(results)
}

# =============================================================================
# Execute Analysis (when script is run directly)
# =============================================================================

if (!interactive()) {
  # Run the full analysis pipeline
  results <- run_muscle_mechanics_analysis()
  
  # Save results to RDS file for later use
  saveRDS(results, file.path(OUTPUT_FIGURES_DIR, "muscle_mechanics_results.rds"))
  
  cli_alert_success("Analysis complete! Results saved to muscle_mechanics_results.rds")
}

# Example usage:
# source("muscle_mechanics_analysis.R")
# results <- run_muscle_mechanics_analysis()