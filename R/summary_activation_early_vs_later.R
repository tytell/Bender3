# summary_activation_early_vs_later.R
# PI request 2026-07-24: compare ISOMETRIC L0 contraction kinetics (activation
# and relaxation times) between EARLY (preconditioning) and LATER (stable)
# trials, where both are available. Separate figure from
# isometric_L0_activation_kinetics.png.
#
# Reads the per-contraction times written by summary_isometric_l0_activation.R
# (data_processed/isometric_l0_activation_times.csv; unit_id embeds the
# chronological bender_NN), classifies each contraction early/later with the
# specimen-specific session-preconditioning cutoff (dynamic_trial_precondition.R:
# bass16=5, bass17=9, bass18=5; trial_num < cutoff = early), and plots
# early-vs-later boxplots (+ raw points) per specimen for activation and
# relaxation time, with the Coughlin & Carroll (2006) red-muscle band for
# reference.
#
# Run: Rscript R/summary_activation_early_vs_later.R  -> figs_summary/

suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork); library(cli) })

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
source(file.path(.pipeline_root, "paths_config.R"))
source(file.path(.pipeline_root, "dynamic_trial_precondition.R"))

OUT_DIR      <- FIGS_SUMMARY_DIR
DATA_OUT_DIR <- file.path(.crittergripper_root(), "02_processed", "data_processed")
CSV_IN       <- file.path(DATA_OUT_DIR, "isometric_l0_activation_times.csv")
PRECOND_COLORS <- c("early (preconditioning)" = "#d73027", "later (stable)" = "#4575b4")
COUGHLIN <- list(activation = c(mean = 78, sd = 22), relaxation = c(mean = 150, sd = 25))
COUGHLIN_WHITE <- list(activation = c(lo = 10, hi = 20, mid = 15),
                       relaxation = c(lo = 28, hi = 45, mid = 37))

read_retry <- function(p, tries = 5L) {
  for (i in seq_len(tries)) { d <- tryCatch(read.csv(p), error = function(e) NULL)
    if (!is.null(d) && nrow(d) > 0L) return(d); Sys.sleep(2) }
  stop("could not read ", p)
}

d <- read_retry(CSV_IN)
d$trial_num  <- extract_bender_trial_num(d$unit_id)
d$precondition <- classify_session_precondition(d$specimen, d$trial_num)
d <- dplyr::filter(d, !is.na(.data$precondition))
write.csv(d, file.path(DATA_OUT_DIR, "isometric_l0_activation_times_precondition.csv"), row.names = FALSE)

# TYPE CONTROL: restrict to dynamic L0 bookends. The isometric-protocol and
# isovelocity V=0 L0 contractions were run almost exclusively LATER in each
# session (see source x stage table below), so a raw all-source early-vs-later
# split would confound session stage with contraction TYPE (early = fast
# bookend twitches; later = slower 0.3-0.5 s tetani). The dynamic pre/post
# L0 bookends are the SAME twitch type and occur throughout the session
# (early AND later for all three specimens), giving a clean stage comparison.
cli::cli_h2("source x stage (why we restrict to bookends)")
print(table(d$specimen, d$source, d$precondition))
d <- dplyr::filter(d, .data$source == "dynamic L0 bookend")

cli::cli_h2("bookend counts by specimen x precondition")
print(dplyr::count(d, .data$specimen, .data$precondition, .drop = FALSE))
# Wilcoxon per specimen where BOTH groups have >= 4 finite values.
.wilcox_lab <- function(df, col) {
  s <- split(df[[col]][is.finite(df[[col]])], droplevels(df$precondition[is.finite(df[[col]])]))
  if (length(s) < 2L || any(lengths(s) < 4L)) return(NA_character_)
  p <- suppressWarnings(stats::wilcox.test(s[[1L]], s[[2L]])$p.value)
  sprintf("Wilcoxon p=%.3f", p)
}

.panel <- function(col, coug, white, title, ylab) {
  df <- dplyr::filter(d, is.finite(.data[[col]]))
  # only specimens with BOTH groups present (where the comparison is possible)
  present <- df |> dplyr::group_by(.data$specimen) |>
    dplyr::summarise(g = dplyr::n_distinct(.data$precondition), .groups = "drop")
  both <- present$specimen[present$g >= 2]
  df <- dplyr::filter(df, .data$specimen %in% both)
  band <- data.frame(lo = coug[["mean"]] - coug[["sd"]], hi = coug[["mean"]] + coug[["sd"]], mid = coug[["mean"]])
  ncnt <- df |> dplyr::count(.data$specimen, .data$precondition)
  ytop <- max(df[[col]], na.rm = TRUE)
  plab <- df |> dplyr::group_by(.data$specimen) |>
    dplyr::summarise(lab = .wilcox_lab(dplyr::pick(dplyr::everything()), col), .groups = "drop") |>
    dplyr::filter(!is.na(.data$lab))
  ggplot(df, aes(x = .data$specimen, y = .data[[col]], color = .data$precondition)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = white[["lo"]], ymax = white[["hi"]], fill = "grey45", alpha = 0.18) +
    geom_hline(yintercept = white[["mid"]], color = "grey35", linetype = "dashed", linewidth = 0.4) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = band$lo, ymax = band$hi, fill = "#7f0000", alpha = 0.10) +
    geom_hline(yintercept = band$mid, color = "#7f0000", linetype = "dashed", linewidth = 0.4) +
    geom_boxplot(fill = NA, outlier.shape = NA, width = 0.6, position = position_dodge(0.7)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7), size = 1.7, alpha = 0.75) +
    geom_text(data = ncnt, aes(label = paste0("n=", .data$n), y = -0.03 * ytop),
              position = position_dodge(0.7), size = 2.7, show.legend = FALSE, vjust = 1) +
    { if (nrow(plab)) geom_text(data = plab, aes(x = .data$specimen, y = 1.05 * ytop, label = .data$lab),
                                inherit.aes = FALSE, size = 3) } +
    scale_color_manual(values = PRECOND_COLORS, name = "session stage", drop = FALSE) +
    labs(title = title, x = NULL, y = ylab) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
}

pA <- .panel("activation_ms", COUGHLIN$activation, COUGHLIN_WHITE$activation, "A. Activation time (rise to 90% peak)", "Activation time (ms)")
pB <- .panel("relaxation_ms", COUGHLIN$relaxation, COUGHLIN_WHITE$relaxation, "B. Relaxation time (offset to 50% decay)", "Relaxation time (ms)")

fig <- (pA | pB) + patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = "L0 bookend-twitch kinetics: early (preconditioning) vs. later (stable) trials",
    subtitle = "Dynamic pre/post L0 BOOKEND twitches only (type-controlled: iso/isovel L0 ran almost only later, so all-source would confound stage with\ntwitch-vs-tetanus). Classified by the specimen session cutoff (bass16<5, bass17<9, bass18<5 = early). Red band = C&C 2006 red (slow) muscle; gray band = C&C 2006 white/fast (sternohyoideus/epaxial).",
    theme = theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 9))) &
  theme(legend.position = "bottom")

fout <- file.path(OUT_DIR, "isometric_L0_activation_earlyVsLater.png")
ggplot2::ggsave(fout, fig, width = 11, height = 6, dpi = 150)
cli::cli_alert_success("Saved {fout}")
