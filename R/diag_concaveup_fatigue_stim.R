# diag_concaveup_fatigue_stim.R
# READ-ONLY DIAGNOSTIC (2026-07-22, PI-requested). Sources the pipeline and only
# READS data + WRITES a PNG -- modifies no analysis.
#
# Question: after the isometric baseline drift is removed, the in-stim PEAK force
# the FL superplot samples STILL rises with |bend|. Is that FATIGUE, or a
# per-step STIMULUS-strength effect (force-per-volt)? Both tested here:
#
#   STIM VOLTAGE -- RULED OUT. Every isometric step is a CONSTANT 5.00 V command
#   (recruitment blocks vary SIDE: left_unilateral / right_unilateral, NOT
#   amplitude). Force-per-volt is a no-op (nothing varies to divide out).
#
#   FATIGUE -- REAL, but does NOT create the concave-up. L0 force decays across
#   the session (bass17 L0 0.043 -> 0.011 N; bass18 similar). BUT within a SINGLE
#   block (fixed fatigue state, fixed 5 V), force rises monotonically with |bend|
#   from the block's own fresh L0 (bass17 within-block cor(F,|strain|) = +0.84,
#   +0.91). Fatigue would push the OPPOSITE way (each block STARTS at its fresh
#   L0), so it cannot produce the arms; it only drags the pooled L0 reps down
#   across the session, mildly deepening the trough (and destabilizing F0 for the
#   normalized superplot). The arms are the within-block force∝|bend| residual
#   (small difference of large passive) -- a passive-subtraction problem, not
#   fatigue or stim.
#
# This figure: |F| vs |strain|, ONE line per block, colored by session order, per
# fish -- each line rises L->R (within-block |bend| scaling at fixed fatigue),
# later blocks sit lower at |strain|=0 (fatigue).
#
# Run: Rscript R/diag_concaveup_fatigue_stim.R
# Canon token: concaveupfatiguestim -- see FIGURES_README.md.

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(stringr)
  library(ggplot2); library(rhdf5); library(cli)
})

.root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
for (f in c("paths_config.R","00_load_bender_flat.R","01_calibrate.R","02_deconvolve.R",
            "muscle_geometry.R","fit_fv_fl.R","03_analyze.R","parse_trial_filename.R",
            "plot_force_vs_time.R","muscle_force_vector.R")) source(file.path(.root, f))

OUT_DIR <- FIGS_DIAGNOSTIC_DIR
DIRS <- list(bass16 = raw_source_dir(BASS16_RAW_SUBFOLDER),
             bass17 = raw_source_dir(BASS17_RAW_SUBFOLDER),
             bass18 = raw_source_dir(BASS18_RAW_SUBFOLDER))

load_one <- function(f) {
  td <- load_bender_flat(f, do_filter = TRUE, loadtorques = "x")
  tau <- deconvolve_bender(f, hub_path = NULL, verbose = FALSE)
  N <- min(nrow(td), length(tau)); td <- td[seq_len(N), , drop = FALSE]
  td$torque_inertia_corrected_Nm <- tau[seq_len(N)]
  attr(td, "Filename") <- f; td
}

# per-step commanded stim amplitude = max|command| over the step (either channel)
read_step_stim_v <- function(filename) {
  h5 <- H5Fopen(filename, "H5F_ACC_RDONLY"); on.exit(try(H5Fclose(h5), silent = TRUE), add = TRUE)
  info <- tryCatch(h5ls(h5, recursive = TRUE), error = function(e) NULL); if (is.null(info)) return(tibble())
  sn <- sort(info$name[info$group == "/timeseries" & grepl("^step_\\d+$", info$name)])
  if (!length(sn)) return(tibble())
  map_dfr(sn, function(s) {
    b <- paste0("/timeseries/", s, "/")
    rd <- function(nm) { p <- paste0(b, nm); if (H5Lexists(h5, p)) tryCatch(as.numeric(h5read(h5, p)), error = function(e) NULL) else NULL }
    c1 <- rd("stim_channel1_command_volt"); c2 <- rd("stim_channel2_command_volt")
    v <- max(c(if (!is.null(c1)) max(abs(c1), na.rm = TRUE) else NA,
               if (!is.null(c2)) max(abs(c2), na.rm = TRUE) else NA), na.rm = TRUE)
    tibble(step_number = as.integer(sub("step_", "", s)), stim_v = v)
  })
}

rows <- list()
for (fish in names(DIRS)) {
  man <- parse_trial_directory(DIRS[[fish]])
  for (fp in man$fullpath[man$protocol == "isometric"]) {
    res <- tryCatch(analyze_isometric(load_one(fp), filename = fp), error = function(e) NULL); if (is.null(res)) next
    res$trial_id <- basename(fp)
    v <- tryCatch(attach_vector_muscle_force(res, fp, "isometric"), error = function(e) NULL); if (is.null(v)) next
    sv <- tryCatch(read_step_stim_v(fp), error = function(e) tibble())
    ss <- left_join(v$step_summary, sv, by = "step_number")
    ss <- ss[ss$muscle_side %in% c("left", "right"), ]
    if (!nrow(ss)) next
    ss$fish <- fish; ss$block <- paste(basename(fp), ss$recruitment)
    rows[[length(rows) + 1L]] <- ss[, c("fish","block","step_number","operating_point",
      "shortening_strain_pct","muscle_force_vector_N","activation_snr","stim_v")]
  }
}
D <- bind_rows(rows)
D$absF <- abs(D$muscle_force_vector_N); D$abs_strain <- abs(D$shortening_strain_pct)
D <- D %>% group_by(block) %>% mutate(block_order = min(step_number)) %>% ungroup()

stim_lab <- sprintf("stim = constant %.2f V on every step (recruitment = SIDE, not amplitude) -> force-per-volt is a no-op",
                    max(D$stim_v, na.rm = TRUE))
p <- ggplot(D, aes(abs_strain, absF, group = block, color = block_order)) +
  geom_line(linewidth = 0.6, alpha = 0.9) + geom_point(size = 1.3, alpha = 0.9) +
  scale_color_viridis_c(name = "session order\n(block start step)", option = "C") +
  facet_wrap(~fish, ncol = 1, scales = "free_y") +
  labs(title = "Concave-up cause: fatigue vs stim vs within-block |bend| scaling",
       subtitle = paste0("Each line = one block (fixed fatigue state + fixed 5 V). Lines RISE L->R = within-block force-|bend| scaling (the arms).\n",
                         "Later blocks (darker) sit lower at |strain|=0 = fatigue at L0. ", stim_lab),
       x = "|muscle strain| (%)", y = "Vector muscle force |F| (N)") +
  theme_bw(11) + theme(legend.position = "right", plot.subtitle = element_text(size = 8))
ggsave(file.path(OUT_DIR, "concaveupfatiguestim.png"), p, width = 9, height = 10, dpi = 140)

cli::cli_h2("Summary")
cli::cli_alert("stim_v range: {sprintf('%.2f-%.2f V', min(D$stim_v,na.rm=TRUE), max(D$stim_v,na.rm=TRUE))} (constant -> force-per-volt moot)")
for (fh in unique(D$fish)) {
  d <- D[D$fish == fh, ]
  wb <- sapply(split(d, d$block), function(b) if (nrow(b) >= 4 && sd(b$abs_strain) > 0) cor(b$absF, b$abs_strain) else NA)
  wb <- wb[is.finite(wb)]
  l0 <- d[abs(d$operating_point) < 0.5, ]
  fat <- if (nrow(l0) >= 3) cor(l0$absF, l0$step_number) else NA
  cli::cli_alert("{fh}: within-block cor(F,|strain|) mean={sprintf('%+.2f', mean(wb))} | L0 fatigue cor(F,order)={sprintf('%+.2f', fat)}")
}
cli::cli_alert_success("Saved concaveupfatiguestim.png to {OUT_DIR}")
