# plot_passive_subtraction_decision_tree.R
# PI directive (2026-07-27): save the passive-subtraction decision tree as
# an actual figure. Reflects the OPTIMIZED logic (PI directive same day:
# "abandon legacy and move forward with best practice") -- i.e. this is the
# NEW decision tree (`passive_force_Nm_optimizedbaseline`,
# build_segmented_step_summary()/03_analyze.R), not the old legacy-default
# one. See analysis_muscle_force_vector_log.md for the full investigation
# that justifies each box.
#
# Run with:  Rscript R/plot_passive_subtraction_decision_tree.R
# Output -> figs_diagnostic/passiveSubtraction_decisionTree_optimizedbaseline.png

suppressPackageStartupMessages({ library(ggplot2); library(tibble); library(cli) })

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
source(file.path(.pipeline_root, "paths_config.R"))
OUT_DIR <- FIGS_DIAGNOSTIC_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

box <- function(id, x, y, w, h, label, fill) tibble::tibble(id = id, x = x, y = y, w = w, h = h, label = label, fill = fill)

boxes <- dplyr::bind_rows(
  box("root",   6.0, 10.0, 3.4, 0.8, "Trial type?", "#e5e7eb"),

  box("dyn",     1.6, 8.4, 3.0, 0.8, "DYNAMIC\n(work loop)", "#dbeafe"),
  box("iso",     6.0, 8.4, 3.0, 0.8, "ISOMETRIC\n(held steps)", "#dbeafe"),
  box("isv",    10.4, 8.4, 3.0, 0.8, "ISOVELOCITY\n(ramp steps)", "#dbeafe"),

  box("dyn_q",   1.6, 6.8, 3.2, 0.9, "Passive cycles exist in\nthe SAME continuous trace?", "#f3f4f6"),
  box("iso_q",   6.0, 6.8, 3.4, 0.9, "No-stim control exists?\nNEVER -- excluded by\nprotocol design", "#f3f4f6"),
  box("isv_q",  10.4, 6.8, 3.4, 0.9, "Same-velocity no-stim\nstep exists in THIS trial?", "#f3f4f6"),

  box("dyn_yes", 0.3, 5.0, 2.6, 0.9, "YES", "#f3f4f6"),
  box("dyn_no",  3.1, 5.0, 2.4, 0.9, "NO (rare)", "#f3f4f6"),
  box("isv_yes", 9.2, 5.0, 2.4, 0.9, "YES", "#f3f4f6"),
  box("isv_no",  11.9, 5.0, 2.6, 0.9, "NO", "#f3f4f6"),

  box("dyn_method_y", 0.3, 3.1, 2.6, 1.1,
      "phase-matched\nhalfcycle subtraction\n(active vs. passive at\nsame %-of-cycle)", "#bbf7d0"),
  box("dyn_method_n", 3.1, 3.1, 2.4, 1.1,
      "stim_off fallback\n(active vs. stim-off\nsamples)", "#fef9c3"),

  box("iso_method", 6.0, 3.1, 3.4, 1.1,
      "baselineInterp\n(pre+post linear\ninterpolation at the\nactive window's own time)", "#bbf7d0"),

  box("isv_method_y", 9.2, 3.1, 2.4, 1.1,
      "velocity_matched\n(real no-stim step,\nsame commanded velocity)", "#bbf7d0"),
  box("isv_method_n", 11.9, 3.1, 2.6, 1.1,
      "baselineInterp\n(pre+post linear\ninterpolation)", "#bbf7d0"),

  box("legacy", 6.0, 0.9, 5.4, 0.9,
      "legacy (pre-only static baseline) -- ABANDONED 2026-07-27:\ndirection-dependent bias, confirmed in all 3 specimens (see log)",
      "#fecaca")
)

seg <- function(x0, y0, x1, y1) tibble::tibble(x0 = x0, y0 = y0, x1 = x1, y1 = y1)
arrows <- dplyr::bind_rows(
  seg(6.0, 9.6, 1.6, 8.8), seg(6.0, 9.6, 6.0, 8.8), seg(6.0, 9.6, 10.4, 8.8),
  seg(1.6, 8.0, 1.6, 7.25), seg(6.0, 8.0, 6.0, 7.25), seg(10.4, 8.0, 10.4, 7.25),
  seg(1.6, 6.35, 0.3, 5.45), seg(1.6, 6.35, 3.1, 5.45),
  seg(10.4, 6.35, 9.2, 5.45), seg(10.4, 6.35, 11.9, 5.45),
  seg(0.3, 4.55, 0.3, 3.65), seg(3.1, 4.55, 3.1, 3.65),
  seg(6.0, 6.35, 6.0, 3.65),
  seg(9.2, 4.55, 9.2, 3.65), seg(11.9, 4.55, 11.9, 3.65),
  seg(6.0, 2.55, 6.0, 1.35)
)

p <- ggplot() +
  geom_rect(data = boxes, aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2, fill = fill),
            color = "grey30", linewidth = 0.4) +
  geom_text(data = boxes, aes(x = x, y = y, label = label), size = 3.1, lineheight = 0.95) +
  geom_segment(data = arrows, aes(x = x0, y = y0, xend = x1, yend = y1),
               arrow = arrow(length = unit(0.12, "cm"), type = "closed"), color = "grey40", linewidth = 0.4) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(-0.5, 14), ylim = c(0.2, 10.6)) +
  labs(title = "Passive subtraction for active muscle force -- decision tree (OPTIMIZED, 2026-07-27)",
       subtitle = "Green = method used. Yellow = rare fallback. Red = ABANDONED (kept in figs for now, not deleted; see analysis_muscle_force_vector_log.md)") +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30"))

fout <- file.path(OUT_DIR, "passiveSubtraction_decisionTree_optimizedbaseline.png")
ggplot2::ggsave(fout, p, width = 13, height = 9, dpi = 150)
cli::cli_alert_success("Saved {fout}")
