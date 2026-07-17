# plot_trial_compound.R
# Required per-trial visualization (task spec): one publication-quality
# compound plot per individual trial file, three panels:
#   1. inertial-corrected torque vs. time
#   2. calculated muscle force vs. time
#   3. predicted muscle strain vs. time (passive vs. active differentiated)
#
# Protocol-agnostic by design: the pipeline (run_fv_fl_power_pipeline.R)
# is responsible for populating td$torque_inertia_corrected_Nm,
# td$muscle_force_Nm (NA where undefined, e.g. passive-only trials), and
# td$strain_active_pct / td$strain_passive_pct BEFORE calling this builder --
# this file only plots what it is given.

library(dplyr)
library(ggplot2)

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("plot_trial_compound.R requires the 'patchwork' package")
}

#' Stitch per-step local time axes into one continuous trial time axis using
#' the recorded inter-step rest gaps (index_step_rest_before_second) --
#' local copy of the equivalent helper in plot_torque_vs_time_batch.R,
#' duplicated here (not sourced) because that script runs unguarded
#' top-level file-listing code on source().
stitch_step_time_local <- function(td, filename) {
  if (!"step_number" %in% names(td)) return(td$t.s)
  rest <- tryCatch(
    as.numeric(rhdf5::h5read(filename, "/metadata/index_step_rest_before_second")),
    error = function(e) rep(0, length(unique(td$step_number)))
  )
  steps <- sort(unique(td$step_number))
  if (length(rest) < length(steps)) rest <- c(rest, rep(0, length(steps) - length(rest)))
  stitched <- numeric(nrow(td))
  t_off <- 0
  for (i in seq_along(steps)) {
    idx <- which(td$step_number == steps[[i]])
    t_local <- td$t.s[idx]
    t0 <- min(t_local, na.rm = TRUE)
    if (i > 1L) {
      gap <- rest[[i]]
      if (is.finite(gap)) t_off <- t_off + gap
    }
    stitched[idx] <- t_local - t0 + t_off
    t_off <- max(stitched[idx], na.rm = TRUE)
  }
  stitched
}


#' Downsample a data frame to a target display rate to keep ggplot/patchwork
#' rendering fast for multi-hundred-thousand-sample trials.
.thin_for_plot <- function(df, time_col, target_hz = 200, max_rows = 400000L) {
  if (nrow(df) <= max_rows) return(df)
  t_vals <- df[[time_col]]
  dt <- stats::median(diff(sort(t_vals)), na.rm = TRUE)
  if (!is.finite(dt) || dt <= 0) return(df)
  step <- max(1L, round((1 / target_hz) / dt))
  df[seq(1L, nrow(df), by = step), , drop = FALSE]
}


#' Build the required 3-panel compound plot for one trial.
#'
#' @param td Tibble with time_col, torque_inertia_corrected_Nm,
#'   muscle_force_Nm, strain_active_pct, strain_passive_pct.
#' @param time_col Name of the (already-stitched, if segmented) time column.
#' @param title,subtitle Plot title/subtitle strings.
#' @param muscle_force_note Optional caption shown under panel 2 (e.g. "passive-only trial -- no active muscle force").
#' @return A patchwork object (ggsave()-able).
build_trial_compound_plot <- function(td, time_col = "t.s", title = "", subtitle = "",
                                      muscle_force_note = NULL) {
  req <- c(time_col, "torque_inertia_corrected_Nm", "muscle_force_Nm",
           "strain_active_pct", "strain_passive_pct")
  missing <- setdiff(req, names(td))
  if (length(missing) > 0L) {
    cli::cli_abort("build_trial_compound_plot: missing column(s): {paste(missing, collapse=', ')}")
  }

  pd <- .thin_for_plot(td, time_col)
  pd <- pd[order(pd[[time_col]]), , drop = FALSE]

  p_torque <- ggplot(pd, aes(x = .data[[time_col]], y = .data$torque_inertia_corrected_Nm)) +
    geom_line(color = "#b91c1c", linewidth = 0.4) +
    labs(y = "Inertial-corrected\ntorque (N*m)", x = NULL) +
    theme_bw(base_size = 11)

  has_force <- any(is.finite(pd$muscle_force_Nm))
  p_force <- ggplot(pd, aes(x = .data[[time_col]], y = .data$muscle_force_Nm)) +
    geom_line(color = "#1d4ed8", linewidth = 0.6) +
    labs(y = "Muscle force\n(N*m)", x = NULL) +
    theme_bw(base_size = 11)
  if (!has_force) {
    note <- muscle_force_note %||% "no active stimulation in this trial -- Force_muscle undefined"
    p_force <- p_force +
      annotate("text", x = mean(range(pd[[time_col]], na.rm = TRUE)), y = 0,
               label = note, color = "gray40", size = 3.2, fontface = "italic")
  }

  strain_long <- dplyr::bind_rows(
    tibble::tibble(x = pd[[time_col]], y = pd$strain_active_pct,  state = "active"),
    tibble::tibble(x = pd[[time_col]], y = pd$strain_passive_pct, state = "passive")
  ) |> dplyr::filter(is.finite(.data$y))

  p_strain <- ggplot(strain_long, aes(x = .data$x, y = .data$y, color = .data$state, linetype = .data$state)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(name = NULL, values = c(active = "#b45309", passive = "#166534")) +
    scale_linetype_manual(name = NULL, values = c(active = "solid", passive = "dashed")) +
    labs(y = "Predicted muscle\nstrain (%)", x = "Time (s)") +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")

  patchwork::wrap_plots(p_torque, p_force, p_strain, ncol = 1) +
    patchwork::plot_annotation(title = title, subtitle = subtitle)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
