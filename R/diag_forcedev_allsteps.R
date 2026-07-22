# diag_forcedev_allsteps.R
# READ-ONLY DIAGNOSTIC (2026-07-22, PI-requested). Sources the existing
# pipeline functions and only READS data + WRITES PNGs -- modifies no analysis.
#
# Purpose: sanity-check the active/passive/residual decomposition BY EYE. Unlike
# forcedevtiming (bass17, 3 hand-picked steps, legacy scalar force), this draws
# EVERY step of EVERY isometric / isovelocity trial for ALL fish (bass16/17/18,
# faceted), each step its own line colored by held strain (isometric) / strain
# rate (isovelocity). Y = the PRODUCTION vector muscle force
# (force_ts$muscle_force_vector_N, RAW sign -- follows bend direction), i.e. the
# exact quantity the FL/FV superplots sample. Stim window shaded.
#
# Why: the FL superplot's pooled concave-up is NOT universal across specimens.
# These traces make visible (a) that bass17's "muscle force" is a post-stim
# baseline DRIFT (climbs ~1 s after stim-off, not an activation twitch), and
# (b) that every isovelocity step has a large NON-ZERO pre-stim offset (the
# windowed-mean passive subtraction leaves the moving passive+inertia
# uncancelled). Both are baseline problems, not FL biology.
#
# Run: Rscript R/diag_forcedev_allsteps.R
# Canon token: forcedevtiming (the *_allsteps.png variant) -- see FIGURES_README.md.

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(stringr)
  library(ggplot2); library(cli)
})

.root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.root, f))
src("paths_config.R")
src("00_load_bender_flat.R")
src("01_calibrate.R")
src("02_deconvolve.R")
src("muscle_geometry.R")
src("fit_fv_fl.R")
src("03_analyze.R")
src("parse_trial_filename.R")
src("plot_force_vs_time.R")   # .smooth_trace_display_only
src("muscle_force_vector.R")  # attach_vector_muscle_force, compute_isovelocity_vector_batch

OUT_DIR <- FIGS_DIAGNOSTIC_DIR
DIRS <- list(bass16 = raw_source_dir(BASS16_RAW_SUBFOLDER),
             bass17 = raw_source_dir(BASS17_RAW_SUBFOLDER),
             bass18 = raw_source_dir(BASS18_RAW_SUBFOLDER))

load_one <- function(f) {
  td <- load_bender_flat(f, do_filter = TRUE, loadtorques = "x")
  tau <- deconvolve_bender(f, hub_path = NULL, verbose = FALSE)
  N <- min(nrow(td), length(tau)); td <- td[seq_len(N), , drop = FALSE]
  td$torque_inertia_corrected_Nm <- tau[seq_len(N)]
  attr(td, "Filename") <- f
  td
}

iso_all <- list(); isov_all <- list()
for (fish in names(DIRS)) {
  man <- parse_trial_directory(DIRS[[fish]])
  cli::cli_h2("{fish}")

  # ---- ISOMETRIC (own-step pre-stim baseline) ----
  for (fp in man$fullpath[man$protocol == "isometric"]) {
    res <- tryCatch(analyze_isometric(load_one(fp), filename = fp), error = function(e) NULL)
    if (is.null(res)) next
    res$trial_id <- basename(fp)
    v <- tryCatch(attach_vector_muscle_force(res, fp, "isometric"), error = function(e) NULL)
    if (is.null(v) || nrow(v$force_ts) == 0) next
    ts <- v$force_ts
    ts$step_number <- as.integer(str_extract(ts$unit_id, "\\d+"))
    ts <- left_join(ts, dplyr::select(v$step_summary, step_number, cond = shortening_strain_pct),
                    by = "step_number")
    ts$fish <- fish; ts$grp <- paste(ts$trial_id, ts$unit_id)
    iso_all[[length(iso_all) + 1L]] <- ts
  }

  # ---- ISOVELOCITY (cross-trial velocity-matched passive batch) ----
  isv <- man$fullpath[man$protocol == "isovelocity"]
  if (length(isv) > 0) {
    inputs <- lapply(isv, function(fp) {
      td <- load_one(fp)
      list(trial_id = basename(fp), filename = fp, res = analyze_isovelocity(td, filename = fp))
    })
    vb <- tryCatch(compute_isovelocity_vector_batch(inputs),
                   error = function(e) { cli::cli_warn(conditionMessage(e)); NULL })
    if (!is.null(vb) && nrow(vb$force_ts) > 0) {
      ts <- vb$force_ts
      ts$step_number <- as.integer(str_extract(ts$unit_id, "\\d+"))
      ss <- bind_rows(lapply(names(vb$step_summaries), function(tid) {
        s <- vb$step_summaries[[tid]]; s$trial_id <- tid; s
      }))
      key <- ss %>% dplyr::select(trial_id, step_number, cond = shortening_strain_pct) %>% distinct()
      ts <- left_join(ts, key, by = c("trial_id", "step_number"))
      ts$fish <- fish; ts$grp <- paste(ts$trial_id, ts$unit_id)
      isov_all[[length(isov_all) + 1L]] <- ts
    }
  }
}

ISO  <- bind_rows(iso_all)
ISOV <- bind_rows(isov_all)
cli::cli_alert_info("iso rows: {nrow(ISO)}, steps: {length(unique(ISO$grp))} | isov rows: {nrow(ISOV)}, steps: {length(unique(ISOV$grp))}")

prep <- function(D) {
  D %>% filter(is.finite(cond)) %>% group_by(grp) %>%
    arrange(t_rel, .by_group = TRUE) %>%
    mutate(y = .smooth_trace_display_only(muscle_force_vector_N)) %>% ungroup()
}
ISO <- prep(ISO); ISOV <- prep(ISOV)

build <- function(D, condlab, condunit, title) {
  dur <- median(D$stim_duration_s, na.rm = TRUE)
  ggplot(D, aes(t_rel, y, group = grp, color = cond)) +
    annotate("rect", xmin = 0, xmax = dur, ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.12) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = c(0, dur), linetype = "dashed", color = "grey40", linewidth = 0.3) +
    geom_line(linewidth = 0.5, alpha = 0.85) +
    scale_color_gradient2(low = "#2166ac", mid = "grey80", high = "#b2182b", midpoint = 0, name = condlab) +
    facet_wrap(~fish, ncol = 1, scales = "free_y") +
    labs(title = title,
         subtitle = sprintf("EVERY step = one line, colored by %s | orange = stim window (~%.2fs) | vector muscle force (production), RAW sign (follows bend direction)",
                            condunit, dur),
         x = "Time relative to stim onset (s)", y = "Vector muscle force (N)") +
    theme_bw(11) + theme(legend.position = "right", plot.subtitle = element_text(size = 8))
}

p_iso <- build(ISO, "strain (%)", "held strain (blue=shortened, red=lengthened)",
               "Within-step vector muscle force development -- ISOMETRIC, all steps")
ggsave(file.path(OUT_DIR, "forcedevtiming_isometric_allsteps.png"), p_iso, width = 9, height = 10, dpi = 140)

p_isov <- build(ISOV, "strain rate (%/s)", "strain rate",
                "Within-step vector muscle force development -- ISOVELOCITY, all steps")
ggsave(file.path(OUT_DIR, "forcedevtiming_isovelocity_allsteps.png"), p_isov, width = 9, height = 10, dpi = 140)

# -- per-fish shape verdict on the ISOMETRIC peak in-stim force (numbers for the log) --
cli::cli_h2("Isometric peak in-stim vector muscle force: per-fish shape vs held strain")
pk <- ISO %>% filter(t_rel >= 0, t_rel <= stim_duration_s) %>%
  group_by(fish, grp, cond) %>% summarise(peak = max(y, na.rm = TRUE), .groups = "drop")
for (fh in unique(pk$fish)) {
  d <- pk[pk$fish == fh, ]
  r_abs  <- suppressWarnings(cor(d$peak, abs(d$cond)))
  r_sign <- suppressWarnings(cor(d$peak, d$cond))
  verdict <- if (is.finite(r_abs) && is.finite(r_sign) && r_abs > abs(r_sign) + 0.1) "SYMMETRIC U (concave-up)"
             else if (is.finite(r_sign) && abs(r_sign) > r_abs + 0.1) "MONOTONIC with signed strain (NOT a U)"
             else "flat / mixed"
  cli::cli_alert("{fh}: cor(force,|strain|)={sprintf('%+.2f', r_abs)}  cor(force,signed)={sprintf('%+.2f', r_sign)} -> {verdict}")
}
cli::cli_alert_success("Saved 2 forcedevtiming *_allsteps diagnostic figures to {OUT_DIR}")
