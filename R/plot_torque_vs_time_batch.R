# plot_torque_vs_time_batch.R
# Batch torque-vs-time PNGs using 00_load_bender_flat.R (flat-schema HDF5).

repo <- normalizePath("/Users/yjimenez/code/Bender3")
source(file.path(repo, "R/paths_config.R"))

# Defaults come from paths_config.R (single source of truth) -- see that
# file if the OneDrive folder layout ever moves again.
dropzone <- DROPZONE_DIR
outdir <- FIGS_TESTS_DIR
date_filter <- "2026-07-07"

source(file.path(repo, "R/00_load_bender_flat.R"))
library(rhdf5)
library(dplyr)
library(ggplot2)

stitch_step_time <- function(td, filename) {
  if (!"step_number" %in% names(td)) {
    return(td$t.s)
  }
  rest <- tryCatch(
    as.numeric(h5read(filename, "/metadata/index_step_rest_before_second")),
    error = function(e) rep(0, length(unique(td$step_number)))
  )
  steps <- sort(unique(td$step_number))
  if (length(rest) < length(steps)) {
    rest <- c(rest, rep(0, length(steps) - length(rest)))
  }
  stitched <- numeric(nrow(td))
  t_off <- 0
  for (i in seq_along(steps)) {
    s <- steps[[i]]
    idx <- which(td$step_number == s)
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

read_primary_axis <- function(filename) {
  attrs <- tryCatch(h5readAttributes(filename, "/metadata"), error = function(e) list())
  ax <- attrs[["daq_primary_bending_axis"]]
  if (is.null(ax) || !nzchar(as.character(ax[[1]]))) return("zTorque")
  as.character(ax[[1]])
}

read_concat_timeseries <- function(filename, ds_name, sampling_mode) {
  h5 <- H5Fopen(filename, "H5F_ACC_RDONLY")
  on.exit(try(H5Fclose(h5), silent = TRUE), add = TRUE)

  if (sampling_mode == "single_finite") {
    path <- paste0("/timeseries/", ds_name)
    if (H5Lexists(h5, path)) {
      return(as.numeric(h5read(h5, path)))
    }
    return(NULL)
  }

  ts_info <- tryCatch(h5ls(h5, recursive = TRUE), error = function(e) NULL)
  if (is.null(ts_info)) return(NULL)
  step_names <- sort(ts_info$name[
    ts_info$group == "/timeseries" & grepl("^step_\\d+$", ts_info$name)
  ])
  parts <- lapply(step_names, function(sname) {
    path <- paste0("/timeseries/", sname, "/", ds_name)
    if (H5Lexists(h5, path)) as.numeric(h5read(h5, path)) else NULL
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0L) return(NULL)
  unlist(parts, use.names = FALSE)
}

has_active_stimulus <- function(td, stim_s1 = NULL, stim_s2 = NULL) {
  vals <- c(
    if ("stim.V" %in% names(td)) td$stim.V else numeric(0),
    if (!is.null(stim_s1)) stim_s1 else numeric(0),
    if (!is.null(stim_s2)) stim_s2 else numeric(0)
  )
  any(vals > 0.1, na.rm = TRUE)
}

torque_plot_column <- function(td) {
  if ("ztorque.Nm" %in% names(td)) "ztorque.Nm" else "ztorque0.Nm"
}

prepare_ztorque_stim_overlay <- function(td, stim_s1 = NULL, stim_s2 = NULL) {
  td$stim_S1_V <- if ("stim.V" %in% names(td)) td$stim.V else if (!is.null(stim_s1)) stim_s1 else NA_real_
  if (!is.null(stim_s2)) {
    td$stim_S2_V <- stim_s2
  } else if (!"stim_S2_V" %in% names(td)) {
    td$stim_S2_V <- NA_real_
  }
  td$ztorque_plot <- td[[torque_plot_column(td)]]

  tq <- td$ztorque_plot
  tq_min <- min(tq, na.rm = TRUE)
  tq_max <- max(tq, na.rm = TRUE)
  tq_span <- tq_max - tq_min
  if (!is.finite(tq_span) || tq_span < 1e-9) {
    pad <- 0.05
    tq_min <- tq_min - pad
    tq_max <- tq_max + pad
    tq_span <- tq_max - tq_min
  }
  stim_max_V <- max(c(td$stim_S1_V, td$stim_S2_V), na.rm = TRUE)
  if (!is.finite(stim_max_V) || stim_max_V <= 0) stim_max_V <- 5.0

  stim_to_torque <- function(v) tq_min + (v / stim_max_V) * tq_span
  td$stim_S1_plot <- stim_to_torque(td$stim_S1_V)
  td$stim_S2_plot <- stim_to_torque(td$stim_S2_V)

  list(
    td = td,
    tq_min = tq_min,
    tq_max = tq_max,
    tq_span = tq_span,
    stim_max_V = stim_max_V
  )
}

thin_for_plot <- function(df, time_col = "t.s", target_hz = 100) {
  if (nrow(df) < 3000L) return(df)
  t_vals <- df[[time_col]]
  dt <- stats::median(diff(sort(t_vals)), na.rm = TRUE)
  if (!is.finite(dt) || dt <= 0) return(df)
  step <- max(1L, round((1 / target_hz) / dt))
  df[seq(1L, nrow(df), by = step), , drop = FALSE]
}

# Command angle is slow-varying; collapse to one sample per AI tick before display.
# Add bracketed phase annotations (segment + tick + label) to an angle panel.
# Bracket lines are drawn at y_bkt; labels at y_lbl. Both must be below the
# data area. Caller is responsible for expanding coord ylim to include them.
add_phase_brackets <- function(p, meta, protocol, y_bkt, y_lbl, tick_h) {
  is_iso <- grepl("isometric", protocol, ignore.case = TRUE)

  fin <- function(x) !is.null(x) && length(x) >= 1L && is.finite(x[[1L]])
  g   <- function(x) if (fin(x)) x[[1L]] else NA_real_

  t_pre_bs <- g(meta$t_pre_baseline_start_s)
  t_act_s  <- g(meta$t_active_start_s)
  t_act_e  <- g(meta$t_active_end_s)
  stim0    <- g(meta$stim_t0_s)
  stim1    <- g(meta$stim_t1_s)
  x_hi     <- g(meta$t_post_baseline_end_s)

  phases <- list()

  # bookend / approach ramp: start to t_pre_baseline_start
  if (is.finite(t_pre_bs) && t_pre_bs > 0.05) {
    phases <- c(phases, list(list(
      t0 = 0, t1 = t_pre_bs,
      lbl = if (is_iso) "bookend + ramp" else "bookend",
      col = "gray40"
    )))
  }

  # pre-stim baseline / approach ramp: t_pre_baseline_start to t_active_start
  if (is.finite(t_pre_bs) && is.finite(t_act_s) && t_act_s - t_pre_bs > 0.05) {
    phases <- c(phases, list(list(
      t0 = t_pre_bs, t1 = t_act_s,
      lbl = if (is_iso) "pre-stim hold" else "approach",
      col = "#1d4ed8"
    )))
  }

  # hold before stim (isometric only): t_active_start to stim_t0 if gap > 0.05 s
  if (is_iso && is.finite(t_act_s) && is.finite(stim0) && stim0 - t_act_s > 0.05) {
    phases <- c(phases, list(list(t0 = t_act_s, t1 = stim0, lbl = "hold", col = "#166534")))
  }

  # stim burst
  if (is.finite(stim0) && is.finite(stim1) && stim1 > stim0) {
    phases <- c(phases, list(list(t0 = stim0, t1 = stim1, lbl = "stim", col = "#b45309")))
  }

  # post-stim hold (isometric: stim_t1 to t_active_end)
  if (is_iso && is.finite(stim1) && is.finite(t_act_e) && t_act_e - stim1 > 0.05) {
    phases <- c(phases, list(list(
      t0 = stim1, t1 = t_act_e, lbl = "post-stim hold", col = "#166534"
    )))
  }

  # return ramp: t_active_end (or stim_t1) to x_hi
  ret_s <- if (is.finite(t_act_e)) t_act_e else (if (is.finite(stim1)) stim1 else NA_real_)
  if (is.finite(ret_s) && is.finite(x_hi) && x_hi - ret_s > 0.05) {
    phases <- c(phases, list(list(t0 = ret_s, t1 = x_hi, lbl = "return", col = "gray40")))
  }

  for (ph in phases) {
    t_mid <- (ph$t0 + ph$t1) / 2
    p <- p +
      annotate("segment",
               x = ph$t0, xend = ph$t1, y = y_bkt, yend = y_bkt,
               color = ph$col, linewidth = 0.85) +
      annotate("segment",
               x = ph$t0, xend = ph$t0,
               y = y_bkt - tick_h, yend = y_bkt + tick_h,
               color = ph$col, linewidth = 0.85) +
      annotate("segment",
               x = ph$t1, xend = ph$t1,
               y = y_bkt - tick_h, yend = y_bkt + tick_h,
               color = ph$col, linewidth = 0.85) +
      annotate("text",
               x = t_mid, y = y_lbl,
               label = ph$lbl, size = 2.4, color = ph$col,
               fontface = "bold", hjust = 0.5)
  }
  p
}

# Per-step LP for zoom display: removes stim electrical pickup spikes while
# preserving the slower bending-torque envelope during motion.
filter_step_torque_for_zoom <- function(
    df,
    raw_col = "ztorque0.Nm",
    out_col = "ztorque_plot",
    sampfreq = 1000,
    cutoff_hz = 20
) {
  if (!raw_col %in% names(df)) {
    if (out_col %in% names(df)) return(df)
    stop("missing torque column for zoom filter")
  }
  if (!requireNamespace("signal", quietly = TRUE)) {
    df[[out_col]] <- df[[raw_col]]
    return(df)
  }
  nyq <- 0.5 * sampfreq
  if (!is.finite(cutoff_hz) || cutoff_hz <= 0 || cutoff_hz >= nyq) {
    df[[out_col]] <- df[[raw_col]]
    return(df)
  }
  filt <- signal::butter(9L, cutoff_hz / nyq, type = "low")
  df[[out_col]] <- signal::filtfilt(filt, df[[raw_col]])
  df
}

stim_shade_fill <- function(stim_side) {
  side <- tolower(as.character(stim_side %||% "left"))
  if (side %in% c("right", "s2")) "#6b7280" else "#0f766e"
}

add_stim_shade_zoom <- function(p, stim_t0, stim_t1, stim_side, ymin = -Inf, ymax = Inf) {
  if (is.null(stim_t0) || is.null(stim_t1) ||
        !is.finite(stim_t0) || !is.finite(stim_t1) || stim_t1 <= stim_t0) {
    return(p)
  }
  side_lbl <- if (tolower(as.character(stim_side %||% "")) %in% c("right", "s2")) {
    "S2 right stim"
  } else {
    "S1 left stim"
  }
  p +
    annotate(
      "rect",
      xmin = stim_t0,
      xmax = stim_t1,
      ymin = ymin,
      ymax = ymax,
      fill = stim_shade_fill(stim_side),
      alpha = 0.18
    ) +
    annotate(
      "text",
      x = mean(c(stim_t0, stim_t1)),
      y = Inf,
      label = side_lbl,
      vjust = 1.4,
      size = 3.2,
      color = stim_shade_fill(stim_side)
    )
}

combine_zoom_panels <- function(p_angle, p_torque, title, subtitle) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("install.packages('patchwork') required for zoom plots with command angle")
  }
  patchwork::wrap_plots(p_angle, p_torque, ncol = 1, heights = c(1, 2)) +
    patchwork::plot_annotation(title = title, subtitle = subtitle)
}

build_ztorque_stim_overlay_plot <- function(
    td,
    time_col,
    title,
    subtitle,
    xlim = NULL,
    stim_t0 = NULL,
    stim_t1 = NULL,
    stim_side = NULL,
    zoom_style = FALSE,
    phase_meta = NULL,
    protocol = "unknown"
) {
  if (zoom_style) {
    if (!"ztorque_plot" %in% names(td)) {
      td$ztorque_plot <- td[[torque_plot_column(td)]]
    }
    td <- td |> dplyr::arrange(.data[[time_col]])
    tq     <- td$ztorque_plot
    tq_min <- min(tq, na.rm = TRUE)
    tq_max <- max(tq, na.rm = TRUE)
    tq_span <- tq_max - tq_min
    if (!is.finite(tq_span) || tq_span < 1e-9) {
      pad <- 0.05; tq_min <- tq_min - pad; tq_max <- tq_max + pad; tq_span <- tq_max - tq_min
    }

    has_angle  <- "angle.deg" %in% names(td) && any(is.finite(td$angle.deg))
    use_patchwork <- has_angle && requireNamespace("patchwork", quietly = TRUE)

    if (use_patchwork) {
      # --- angle data: direct subsample to ~50 Hz from the already-thinned td ---
      ang_raw  <- td |>
        dplyr::filter(is.finite(.data[[time_col]]), is.finite(.data$angle.deg)) |>
        dplyr::arrange(.data[[time_col]])
      n_ang    <- nrow(ang_raw)
      step_ang <- max(1L, as.integer(floor(n_ang / 175L)))
      ang_df   <- ang_raw[seq(1L, n_ang, by = step_ang), c(time_col, "angle.deg"), drop = FALSE]

      ang_min  <- min(ang_df$angle.deg, na.rm = TRUE)
      ang_max  <- max(ang_df$angle.deg, na.rm = TRUE)
      ang_span <- max(ang_max - ang_min, 1)
      # pad the data range; extra room below for bracket annotations
      ang_data_pad <- max(1.0, 0.04 * ang_span)
      has_brackets  <- !is.null(phase_meta) && is.list(phase_meta)
      bkt_space     <- if (has_brackets) (ang_span + ang_data_pad) * 0.45 else 0
      y_lo <- ang_min - ang_data_pad - bkt_space
      y_hi <- ang_max + ang_data_pad
      # bracket geometry relative to data range
      y_bkt  <- ang_min - ang_data_pad - bkt_space * 0.30
      y_lbl  <- ang_min - ang_data_pad - bkt_space * 0.72
      tick_h <- ang_span * 0.045
      # axis breaks only over data range (no brackets shown on y axis)
      bk <- pretty(c(ang_min, ang_max), n = 4L)
      bk <- bk[bk >= ang_min - ang_data_pad & bk <= ang_max + ang_data_pad]

      use_xlim <- !is.null(xlim) && length(xlim) == 2L && all(is.finite(xlim))

      p_angle <- ggplot(ang_df,
                        aes(x = .data[[time_col]], y = .data$angle.deg)) +
        geom_line(color = "#2563eb", linewidth = 0.7) +
        scale_y_continuous(name = "Command\nangle (deg)", breaks = bk) +
        labs(x = NULL) +
        coord_cartesian(
          xlim = if (use_xlim) xlim else NULL,
          ylim = c(y_lo, y_hi),
          clip = "off"
        ) +
        theme_bw(base_size = 11) +
        theme(axis.title.x = element_blank(),
              plot.margin = ggplot2::margin(t = 2, r = 4, b = 14, l = 4, unit = "mm"))

      p_angle <- add_stim_shade_zoom(p_angle, stim_t0, stim_t1, stim_side,
                                      ang_min - ang_data_pad, y_hi)

      if (has_brackets) {
        p_angle <- add_phase_brackets(p_angle, phase_meta, protocol,
                                       y_bkt, y_lbl, tick_h)
      }

      # --- torque panel ---
      p_torque <- ggplot(td, aes(x = .data[[time_col]], y = .data$ztorque_plot)) +
        geom_line(color = "#b91c1c", linewidth = 0.6) +
        scale_y_continuous(name = "zTorque (N*m)") +
        labs(x = "Time (s)") +
        coord_cartesian(
          xlim = if (use_xlim) xlim else NULL,
          ylim = c(tq_min, tq_max)
        ) +
        theme_bw(base_size = 11)
      p_torque <- add_stim_shade_zoom(p_torque, stim_t0, stim_t1, stim_side, tq_min, tq_max)

      return(combine_zoom_panels(p_angle, p_torque, title, subtitle))
    }

    p <- ggplot(td, aes(x = .data[[time_col]])) +
      geom_line(aes(y = ztorque_plot), color = "#b91c1c", linewidth = 0.6) +
      scale_y_continuous(name = "zTorque (N*m)", limits = c(tq_min, tq_max)) +
      labs(x = "Time (s)", title = title, subtitle = subtitle) +
      theme_bw(base_size = 11)
    p <- add_stim_shade_zoom(p, stim_t0, stim_t1, stim_side, tq_min, tq_max)
    if (!is.null(xlim) && length(xlim) == 2L && all(is.finite(xlim))) {
      p <- p + coord_cartesian(xlim = xlim)
    }
    return(p)
  }

  prep <- prepare_ztorque_stim_overlay(td)
  td <- prep$td |> dplyr::arrange(.data[[time_col]])
  tq_min <- prep$tq_min
  tq_max <- prep$tq_max
  tq_span <- prep$tq_span
  stim_max_V <- prep$stim_max_V

  p <- ggplot(td, aes(x = .data[[time_col]])) +
    geom_line(aes(y = ztorque_plot, color = "zTorque"), linewidth = 0.55) +
    geom_line(aes(y = stim_S1_plot, color = "S1 left command (V)"), linewidth = 0.45) +
    geom_line(aes(y = stim_S2_plot, color = "S2 right command (V)"), linewidth = 0.45) +
    scale_color_manual(
      name = NULL,
      values = c(
        "zTorque" = "#b91c1c",
        "S1 left command (V)" = "#0f766e",
        "S2 right command (V)" = "#6b7280"
      )
    ) +
    scale_y_continuous(
      name = "zTorque (N*m)",
      limits = c(tq_min, tq_max),
      sec.axis = sec_axis(
        transform = ~ (. - tq_min) * stim_max_V / tq_span,
        name = "Stimulus voltage (V)"
      )
    ) +
    labs(x = "Time (s)", title = title, subtitle = subtitle) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")

  if (!is.null(xlim) && length(xlim) == 2L && all(is.finite(xlim))) {
    p <- p + coord_cartesian(xlim = xlim)
  }
  p
}

read_step_metadata <- function(filename, step_number) {
  idx <- as.integer(step_number) - 1L
  out <- list(step_number = step_number)
  tryCatch({
    h5 <- H5Fopen(filename, "H5F_ACC_RDONLY")
    on.exit(try(H5Fclose(h5), silent = TRUE), add = TRUE)
    m_ds <- function(key) {
      path <- paste0("/metadata/", key)
      if (H5Lexists(h5, path)) tryCatch(h5read(h5, path), error = function(e) NULL) else NULL
    }
    scalar_at <- function(v) if (!is.null(v)) as.numeric(v)[[idx + 1L]] else NULL

    op    <- m_ds("index_step_operating_point")
    if (!is.null(op)) out$operating_point <- as.numeric(op)[[idx + 1L]]
    units <- m_ds("index_step_operating_point_units")
    if (!is.null(units)) out$operating_point_units <- as.character(units[[idx + 1L]])
    sides <- m_ds("index_step_block_stim_sides")
    if (!is.null(sides)) out$stim_side <- as.character(sides[[idx + 1L]])

    out$stim_t0_s              <- scalar_at(m_ds("index_step_stim_t0_second"))
    out$stim_t1_s              <- scalar_at(m_ds("index_step_stim_t1_second"))
    out$t_post_baseline_end_s  <- scalar_at(m_ds("index_step_t_post_baseline_end_second"))
    out$t_pre_baseline_start_s <- scalar_at(m_ds("index_step_t_pre_baseline_start_second"))
    out$t_pre_baseline_end_s   <- scalar_at(m_ds("index_step_t_pre_baseline_end_second"))
    out$t_active_start_s       <- scalar_at(m_ds("index_step_t_active_start_second"))
    out$t_active_end_s         <- scalar_at(m_ds("index_step_t_active_end_second"))
  }, error = function(e) NULL)
  out
}

select_extreme_stim_step <- function(filename) {
  td <- load_bender_flat(filename, do_filter = FALSE, loadtorques = "none")
  if (is.null(td) || !"step_number" %in% names(td)) return(NA_integer_)

  mode <- attr(td, "protocol_sampling_mode")
  stim_s2 <- read_concat_timeseries(filename, "stim_channel2_command_volt", mode)
  if (!is.null(stim_s2) && length(stim_s2) == nrow(td)) {
    td$stim_S2_V <- stim_s2
  }

  step_summ <- td |>
    dplyr::group_by(.data$step_number) |>
    dplyr::summarize(
      has_stim = any(.data$stim.V > 0.1, na.rm = TRUE) ||
        (if ("stim_S2_V" %in% names(td)) any(.data$stim_S2_V > 0.1, na.rm = TRUE) else FALSE),
      .groups = "drop"
    )

  op <- tryCatch(as.numeric(h5read(filename, "/metadata/index_step_operating_point")), error = function(e) NULL)
  steps <- tryCatch(as.integer(h5read(filename, "/metadata/index_step_number")), error = function(e) NULL)
  if (is.null(op) || is.null(steps) || length(op) != length(steps)) return(NA_integer_)

  cand <- tibble::tibble(step_number = steps, operating_point = op) |>
    dplyr::inner_join(step_summ, by = "step_number") |>
    dplyr::filter(.data$has_stim)

  if (nrow(cand) == 0L) return(NA_integer_)
  cand$step_number[[which.max(abs(cand$operating_point))]]
}

plot_ztorque_stim_extreme_step <- function(filename, outdir) {
  step_number <- select_extreme_stim_step(filename)
  if (is.na(step_number)) {
    message("[skip] no stimulated step found for extreme pick: ", basename(filename))
    return(invisible(NULL))
  }
  message("[info] extreme stimulated step for ", basename(filename), ": ", step_number)
  plot_ztorque_stim_one_step(filename, outdir, step_number = step_number)
}

plot_ztorque_stim_one_step <- function(filename, outdir, step_number = 2L) {
  td <- load_bender_flat(filename, do_filter = FALSE, loadtorques = "z")
  if (is.null(td) || nrow(td) == 0L || !"step_number" %in% names(td)) {
    message("[skip] not segmented or empty: ", basename(filename))
    return(invisible(NULL))
  }

  mode <- attr(td, "protocol_sampling_mode")
  if (!identical(mode, "segmented_finite")) {
    message("[skip] not segmented_finite: ", basename(filename))
    return(invisible(NULL))
  }

  stim_s2 <- read_concat_timeseries(filename, "stim_channel2_command_volt", mode)
  if (!is.null(stim_s2) && length(stim_s2) == nrow(td)) {
    td$stim_S2_V <- stim_s2
  }

  sd <- td |> dplyr::filter(.data$step_number == step_number)
  if (nrow(sd) == 0L) {
    message("[skip] step ", step_number, " missing in ", basename(filename))
    return(invisible(NULL))
  }

  if (!has_active_stimulus(sd)) {
    message("[skip] step ", step_number, " has no active stimulus in ", basename(filename))
    return(invisible(NULL))
  }

  meta <- read_step_metadata(filename, step_number)
  protocol <- unique(td$stimulus_type)
  protocol <- if (length(protocol)) protocol[[1]] else "unknown"
  specimen <- unique(td$fishcode)
  specimen <- if (length(specimen)) specimen[[1]] else "unknown"

  op_txt <- if (!is.null(meta$operating_point) && is.finite(meta$operating_point)) {
    paste0(meta$operating_point, " ", meta$operating_point_units %||% "")
  } else {
    ""
  }
  stim_txt <- if (!is.null(meta$stim_t0_s) && !is.null(meta$stim_t1_s)) {
    paste0(" | stim ", meta$stim_t0_s, "-", meta$stim_t1_s, " s")
  } else {
    ""
  }
  if (grepl("isovelocity", protocol, ignore.case = TRUE)) {
    stim_txt <- paste0(stim_txt, " | stim at CV segment onset (stim_onset_s=0)")
  }
  torque_col <- if ("ztorque0.Nm" %in% names(sd)) "ztorque0.Nm" else torque_plot_column(sd)
  filter_note <- " | 20 Hz LP per-step (zoom display)"
  x_hi <- if (!is.null(meta$t_post_baseline_end_s) && is.finite(meta$t_post_baseline_end_s)) {
    meta$t_post_baseline_end_s
  } else if (!is.null(meta$stim_t1_s) && is.finite(meta$stim_t1_s)) {
    meta$stim_t1_s + 2.0
  } else {
    max(sd$t.s, na.rm = TRUE)
  }
  sd_plot <- sd |> dplyr::filter(.data$t.s <= x_hi) |> dplyr::arrange(.data$t.s)
  sampfreq <- attr(td, "SampleFrequency.Hz")
  if (!is.finite(sampfreq) || sampfreq <= 0) sampfreq <- 1000
  sd_plot <- filter_step_torque_for_zoom(
    sd_plot,
    raw_col = torque_col,
    out_col = "ztorque_plot",
    sampfreq = sampfreq,
    cutoff_hz = 20
  )
  sd_plot <- thin_for_plot(sd_plot, time_col = "t.s")
  subtitle <- paste0(
    specimen, " | ", protocol, " | step ", step_number,
    if (nzchar(op_txt)) paste0(" | op ", op_txt) else "",
    if (!is.null(meta$stim_side)) paste0(" | stim side ", meta$stim_side) else "",
    stim_txt,
    " | local step time (trimmed before return-to-zero)",
    filter_note
  )

  # store x_hi into meta so add_phase_brackets can use it for the return bracket
  meta$t_post_baseline_end_s <- x_hi

  p <- build_ztorque_stim_overlay_plot(
    sd_plot,
    time_col = "t.s",
    title = paste0(basename(filename), " -- step ", step_number),
    subtitle = subtitle,
    xlim = c(0, x_hi),
    stim_t0 = meta$stim_t0_s,
    stim_t1 = meta$stim_t1_s,
    stim_side = meta$stim_side,
    zoom_style = TRUE,
    phase_meta = meta,
    protocol = protocol
  )

  out_base <- sub(
    "\\.h5$",
    sprintf("_step%03d_ztorque_stim_zoom.png", step_number),
    basename(filename)
  )
  out_path <- file.path(outdir, out_base)
  ggsave(out_path, p, width = 10, height = 7, dpi = 150)
  message("[saved] ", out_path)
  out_path
}

`%||%` <- function(x, y) if (is.null(x)) y else x

plot_ztorque_stim_vs_time <- function(filename, outdir) {
  td <- load_bender_flat(filename, do_filter = TRUE, loadtorques = "z")
  if (is.null(td) || nrow(td) == 0L) {
    message("[skip] empty: ", basename(filename))
    return(invisible(NULL))
  }

  mode <- attr(td, "protocol_sampling_mode")
  stim_s1 <- if ("stim.V" %in% names(td)) td$stim.V else NULL
  stim_s2 <- read_concat_timeseries(filename, "stim_channel2_command_volt", mode)
  if (!is.null(stim_s2) && length(stim_s2) != nrow(td)) {
    message("[warn] stim channel2 length mismatch in ", basename(filename))
    stim_s2 <- NULL
  }

  if (!has_active_stimulus(td, stim_s1 = stim_s1, stim_s2 = stim_s2)) {
    message("[skip] no active stimulus: ", basename(filename))
    return(invisible(NULL))
  }
  if (!"ztorque0.Nm" %in% names(td)) {
    message("[skip] no z torque: ", basename(filename))
    return(invisible(NULL))
  }

  td$t_stitch <- stitch_step_time(td, filename)
  protocol <- unique(td$stimulus_type)
  protocol <- if (length(protocol)) protocol[[1]] else "unknown"
  specimen <- unique(td$fishcode)
  specimen <- if (length(specimen)) specimen[[1]] else "unknown"
  trial <- unique(td$trial)
  trial <- if (length(trial) && is.finite(trial[[1]])) trial[[1]] else NA_real_

  subtitle <- paste0(
    specimen,
    if (is.finite(trial)) paste0(" | trial ", trial) else "",
    " | ", protocol,
    " | ", mode,
    if (mode == "segmented_finite") " (steps stitched with recorded rest gaps)" else "",
    " | stim overlay uses right-hand voltage axis | 50 Hz LP filtered"
  )

  p <- build_ztorque_stim_overlay_plot(
    td,
    time_col = "t_stitch",
    title = basename(filename),
    subtitle = subtitle
  )

  out_base <- sub("\\.h5$", "_ztorque_stim_vs_time.png", basename(filename))
  out_path <- file.path(outdir, out_base)
  ggsave(out_path, p, width = 12, height = 6, dpi = 150)
  message("[saved] ", out_path)
  out_path
}

plot_torque_vs_time <- function(filename, outdir) {
  td <- load_bender_flat(filename, do_filter = FALSE, loadtorques = c("x", "y", "z"))
  if (is.null(td) || nrow(td) == 0L) {
    message("[skip] empty: ", basename(filename))
    return(invisible(NULL))
  }

  primary <- read_primary_axis(filename)
  td$t_stitch <- stitch_step_time(td, filename)
  mode <- attr(td, "protocol_sampling_mode")
  protocol <- unique(td$stimulus_type)
  protocol <- if (length(protocol)) protocol[[1]] else "unknown"
  specimen <- unique(td$fishcode)
  specimen <- if (length(specimen)) specimen[[1]] else "unknown"
  trial <- unique(td$trial)
  trial <- if (length(trial) && is.finite(trial[[1]])) trial[[1]] else NA_real_

  torque_map <- c(
    xTorque = "xtorque0.Nm",
    yTorque = "ytorque0.Nm",
    zTorque = "ztorque0.Nm"
  )

  long <- bind_rows(lapply(names(torque_map), function(nm) {
    col <- torque_map[[nm]]
    if (!col %in% names(td)) return(NULL)
    tibble(
      t.s = td$t_stitch,
      torque = td[[col]],
      axis = nm,
      stim.V = if ("stim.V" %in% names(td)) td$stim.V else NA_real_
    )
  }))

  if (nrow(long) == 0L) {
    message("[skip] no torque columns: ", basename(filename))
    return(invisible(NULL))
  }

  long$axis <- factor(long$axis, levels = c("xTorque", "yTorque", "zTorque"))
  long$is_primary <- long$axis == primary

  stim_df <- long |>
    distinct(t.s, stim.V) |>
    filter(is.finite(stim.V))

  p <- ggplot(long, aes(x = t.s, y = torque, color = is_primary)) +
    geom_line(linewidth = 0.35, show.legend = FALSE) +
    scale_color_manual(values = c("TRUE" = "#b91c1c", "FALSE" = "#374151")) +
    facet_grid(axis ~ ., scales = "free_y", labeller = labeller(axis = function(x) paste0(x, " (N*m)"))) +
    labs(
      x = "Time (s)",
      y = NULL,
      title = basename(filename),
      subtitle = paste0(
        specimen,
        if (is.finite(trial)) paste0(" | trial ", trial) else "",
        " | ", protocol,
        " | ", mode,
        if (mode == "segmented_finite") " (steps stitched with recorded rest gaps)" else ""
      )
    ) +
    theme_bw(base_size = 11) +
    theme(strip.text.y = element_text(angle = 0))

  if (nrow(stim_df) > 0L && any(stim_df$stim.V > 0, na.rm = TRUE)) {
    p_stim <- ggplot(stim_df, aes(x = t.s, y = stim.V)) +
      geom_line(color = "#0f766e", linewidth = 0.35) +
      labs(x = "Time (s)", y = "Stim command (V)", title = NULL) +
      theme_bw(base_size = 10)
    if (requireNamespace("patchwork", quietly = TRUE)) {
      p <- patchwork::wrap_plots(p, p_stim, ncol = 1, heights = c(3, 1))
    }
  }

  out_base <- sub("\\.h5$", "_torque_vs_time.png", basename(filename))
  out_path <- file.path(outdir, out_base)
  ggsave(out_path, p, width = 12, height = 9, dpi = 150)
  message("[saved] ", out_path)
  out_path
}

pattern <- paste0("^", date_filter, "_.*\\.h5$")
files <- sort(list.files(dropzone, pattern = pattern, full.names = TRUE))
if (length(files) == 0L) {
  stop("No H5 files matching ", pattern, " in ", dropzone)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (identical(sys.nframe(), 0L)) {
if (isTRUE(getOption("bender.plot_torque_all", FALSE))) {
  saved <- lapply(files, plot_torque_vs_time, outdir = outdir)
  cat("Torque plots:", length(Filter(Negate(is.null), saved)), "saved to", outdir, "\n")
}

saved_stim <- lapply(files, plot_ztorque_stim_vs_time, outdir = outdir)
cat("Z-torque + stim plots:", length(Filter(Negate(is.null), saved_stim)), "saved to", outdir, "\n")

zoom_files <- c(
  "2026-07-07_rod40a_bender_02_isovelocity.h5",
  "2026-07-07_rod40a_bender_04_isometric.h5",
  "2026-07-07_rod40a_bender_05_isometric.h5"
)
saved_zoom <- lapply(zoom_files, function(fn) {
  path <- file.path(dropzone, fn)
  if (!file.exists(path)) {
    message("[skip] missing file for zoom: ", fn)
    return(invisible(NULL))
  }
  plot_ztorque_stim_extreme_step(path, outdir)
})
cat("Single-step zoom plots:", length(Filter(Negate(is.null), saved_zoom)), "saved to", outdir, "\n")
}
