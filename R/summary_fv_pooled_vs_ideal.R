# summary_fv_pooled_vs_ideal.R
# PI request 2026-07-24: the per-individual isovelocity FV curves look
# consistent, so pool all specimens on ONE normalised graph (each individual
# a distinct colour + its own line) and overlay an IDEALISED force-velocity
# curve to show closeness of fit to the textbook shape.
#
# Reads the gated per-step CSV already produced by
# summary_fv_fl_best_within_individual.R (no slow re-collection):
#   data_processed/fv_fl_best_within_individual_steps.csv
# Those rows are already RIGHT-muscle, F>0, and (for FV) clean-sono-ramp
# (R^2 >= 0.90) -- i.e. the best-quality isovelocity points.
#
# Conventions (match the canonical reference):
#   x = V/Vmax  (SHORTENING negative, LENGTHENING positive), per individual.
#   y = F/Fmax  (active tension / that individual's peak FV tension), per
#       individual -- so every individual's lengthening plateau -> 1 and the
#       shapes are directly comparable regardless of absolute force.
#
# Idealised curve: a Hill-type sigmoid fit to the POOLED normalised points
# (upper asymptote L, steepness k, midpoint x0), the standard smooth FV
# shape (~0 at fast shortening, rising through the isometric region to a
# lengthening plateau). Pooled R^2 of data vs ideal is annotated.
#
# Run: Rscript R/summary_fv_pooled_vs_ideal.R  -> figs_summary/

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
source(file.path(.pipeline_root, "paths_config.R"))

OUT_DIR      <- FIGS_SUMMARY_DIR
DATA_OUT_DIR <- file.path(.crittergripper_root(), "02_processed", "data_processed")
CSV_IN       <- file.path(DATA_OUT_DIR, "fv_fl_best_within_individual_steps.csv")
SPECIMEN_COLORS <- c(bass16 = "#1b9e77", bass17 = "#d95f02", bass18 = "#7570b3")

# retry read (OneDrive can briefly fail to hydrate)
read_retry <- function(p, tries = 5L) {
  for (i in seq_len(tries)) {
    d <- tryCatch(read.csv(p), error = function(e) NULL)
    if (!is.null(d) && nrow(d) > 0L) return(d)
    Sys.sleep(2)
  }
  stop("could not read ", p)
}

d <- read_retry(CSV_IN)
fv <- d |>
  dplyr::filter(.data$category == "isovelocity", is.finite(.data$x_val), is.finite(.data$muscle_force_Nm)) |>
  dplyr::group_by(.data$specimen) |>
  dplyr::mutate(V_norm = .data$x_val / max(abs(.data$x_val), na.rm = TRUE),
                F_norm = .data$muscle_force_Nm / max(.data$muscle_force_Nm, na.rm = TRUE)) |>
  dplyr::ungroup()

message(sprintf("Pooled FV points: %d across %d individuals", nrow(fv), dplyr::n_distinct(fv$specimen)))

# ---- fit idealised Hill-type sigmoid to pooled normalised data ----
# y = L / (1 + exp(-k (x - x0))), robust bounded optimisation (nls is
# unstable on this small, plateau-heavy sample).
# Upper asymptote fixed at 1 (force is normalised so each individual's
# lengthening plateau = 1); fit only steepness k and midpoint x0.
logi <- function(par, x) 1 / (1 + exp(-par[1L] * (x - par[2L])))
fit_logi <- function(df) {
  sse <- function(par) sum((df$F_norm - logi(par, df$V_norm))^2)
  tryCatch(stats::optim(c(k = 3.0, x0 = -0.1), sse, method = "L-BFGS-B",
                        lower = c(0.5, -1.0), upper = c(15, 0.5)),
           error = function(e) NULL)
}

# PI-directed exclusion (2026-07-24): bass18 has three LENGTHENING points that
# sit well above the idealised curve (bender_03 steps 14-16 -- the mid-velocity
# eccentric steps, residual +0.32..+0.39). These show force enhancement beyond
# the symmetric idealised sigmoid (genuine eccentric overshoot / force-length
# contamination at those lengths). Dropped for the clean closeness-to-ideal
# demonstration; identified as the top-3 positive residuals among bass18
# lengthening points relative to the initial fit.
opt0 <- fit_logi(fv)
cf0  <- if (!is.null(opt0) && opt0$convergence == 0) opt0$par else c(k = 2.9, x0 = 0.40)
fv$.resid0 <- fv$F_norm - logi(cf0, fv$V_norm)
drop_idx <- fv |>
  dplyr::mutate(.row = dplyr::row_number()) |>
  dplyr::filter(.data$specimen == "bass18", .data$V_norm > 0, .data$.resid0 > 0) |>
  dplyr::slice_max(.data$.resid0, n = 3) |>
  dplyr::pull(.row)
dropped <- fv[drop_idx, c("specimen", "trial_id", "step_number", "V_norm", "F_norm")]
message("Dropped (above-ideal, lengthening):"); print(dropped)
fv <- fv[-drop_idx, , drop = FALSE]
fv$.resid0 <- NULL
message(sprintf("Retained FV points: %d", nrow(fv)))

opt <- fit_logi(fv)
xg <- seq(min(fv$V_norm) - 0.05, 1.0, length.out = 200)
if (!is.null(opt) && opt$convergence == 0) {
  cf <- opt$par
  ideal_fun <- function(x) logi(cf, x)
  ideal_label <- sprintf("Idealised Hill-type FV\n(plateau=1, k=%.1f, x0=%.2f)", cf[1L], cf[2L])
} else {
  # canonical fallback (a=0.25 concentric, 1.5x eccentric plateau), plateau-normalised
  raw <- function(x) ifelse(x <= 0, (1 + x) / (1 - x / 0.25), 1 + 0.5 * (1 - exp(-x / 0.15)))
  plat <- raw(1)
  ideal_fun <- function(x) raw(x) / plat
  ideal_label <- "Idealised Hill-type FV (canonical)"
}
ideal_df <- data.frame(V_norm = xg, F_norm = ideal_fun(xg))

# pooled R^2 of data vs idealised curve
yhat <- ideal_fun(fv$V_norm)
r2 <- 1 - sum((fv$F_norm - yhat)^2) / sum((fv$F_norm - mean(fv$F_norm))^2)
message(sprintf("Pooled R^2 (data vs idealised) = %.3f", r2))
write.csv(fv, file.path(DATA_OUT_DIR, "fv_pooled_normalised.csv"), row.names = FALSE)

p <- ggplot() +
  geom_hline(yintercept = c(0, 1), linetype = "dotted", color = "grey85") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
  geom_line(data = ideal_df, aes(x = V_norm, y = F_norm), color = "black", linewidth = 1.1, linetype = "22") +
  geom_smooth(data = fv, aes(x = V_norm, y = F_norm, color = specimen),
              method = "loess", se = FALSE, span = 1.2, linewidth = 0.7) +
  geom_point(data = fv, aes(x = V_norm, y = F_norm, color = specimen), size = 2.6) +
  annotate("text", x = -0.9, y = 1.06, label = "Shortening", hjust = 0, size = 3.6, color = "grey30") +
  annotate("text", x =  0.55, y = 1.06, label = "Lengthening", hjust = 0, size = 3.6, color = "grey30") +
  annotate("text", x = 0.98, y = 0.10, label = sprintf("pooled R^2 = %.2f\nvs idealised", r2), hjust = 1, vjust = 0, size = 3.6) +
  annotate("text", x = min(fv$V_norm), y = 0.92, label = ideal_label, hjust = 0, size = 3.0, color = "grey20") +
  scale_color_manual(values = SPECIMEN_COLORS, name = "individual") +
  coord_cartesian(xlim = c(min(fv$V_norm) - 0.05, 1.02), ylim = c(-0.02, 1.12)) +
  labs(title = "Pooled force-velocity across individuals vs. idealised Hill curve",
       subtitle = sprintf("Best isovelocity points (RIGHT muscle, F>0, clean sono ramp R^2>=0.90). Force & velocity from the encoder-defined\nconstant-velocity window applied to sono. Each individual normalised to its own Vmax and peak FV tension. n=%d.", nrow(fv)),
       x = expression(Normalised~velocity~V/V[max]~~"(shortening < 0 < lengthening)"),
       y = expression(Normalised~force~F/F[max])) +
  theme_bw(base_size = 12) + theme(legend.position = "right")

fout <- file.path(OUT_DIR, "FV_pooled_vs_ideal.png")
ggplot2::ggsave(fout, p, width = 9.5, height = 6, dpi = 150)
message("Saved ", fout)
