# plot_fatigue_check.R
# Diagnostic (not one of the originally-required deliverables): tests whether
# the monotonic (non-parabolic) isometric/isovelocity FL/FV data could be a
# fatigue artifact rather than a true force-length/velocity property, by
# plotting muscle_force_Nm against STIMULATION ORDER within each
# index_step_block_number block -- if force declines steadily with repeated
# stimulation regardless of the commanded strain/strain-rate at each step,
# that is consistent with fatigue/potentiation confounding the FL/FV shape,
# not a real length- or velocity-dependence.

library(dplyr)
library(ggplot2)

#' @param step_summary Combined step_summary across trials (must have
#'   trial_id, block_number, step_number, muscle_side, contraction_mode,
#'   shortening_strain_pct, muscle_force_Nm -- i.e. analyze_isometric()/
#'   analyze_isovelocity()$step_summary, row-bound across trials with a
#'   trial_id column, as built by run_fv_fl_power_pipeline.R).
build_fatigue_check_plot <- function(step_summary, title = "Fatigue check: muscle force vs. stimulation order") {
  df <- step_summary |>
    dplyr::filter(.data$muscle_side %in% c("left", "right"), is.finite(.data$muscle_force_Nm)) |>
    dplyr::mutate(block_number = dplyr::if_else(is.na(.data$block_number), 1L, .data$block_number)) |>
    dplyr::group_by(.data$trial_id, .data$block_number) |>
    dplyr::arrange(.data$step_number, .by_group = TRUE) |>
    dplyr::mutate(order_in_block = dplyr::row_number()) |>
    dplyr::ungroup()

  ggplot(df, aes(x = .data$order_in_block, y = .data$muscle_force_Nm,
                 color = .data$muscle_side, group = interaction(.data$trial_id, .data$muscle_side))) +
    geom_line(alpha = 0.7) +
    geom_point(aes(shape = .data$contraction_mode), size = 2.2) +
    facet_grid(trial_id ~ block_number, labeller = label_both) +
    scale_color_manual(values = c(left = "#1d4ed8", right = "#b91c1c"), name = "Muscle side") +
    labs(title = title,
         subtitle = "Flat/declining trend WITHIN a block regardless of strain(-rate) at each step suggests fatigue, not a true FL/FV shape",
         x = "Stimulation order within block", y = "Muscle force (N*m)", shape = "Contraction mode") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
}
