library('hdf5r')
library('dplyr')
library('purrr')
# library('pracma')

get_stim_cycles <- function(data, filename) {
  # use just the first row, so that we get the single attribute values
  stim <-
    with(slice(data, 1), {
      tibble(cycle = seq(from=floor(min(data$cycle)), to=ceiling(max(data$cycle)))) %>%
        mutate(on_left = (cycle + stim_phase/100 + 0.25) / frequency.Hz,
               off_left = on_left + stim_duty/100/frequency.Hz,
               on_right = (cycle + stim_phase/100 + 0.75) / frequency.Hz,
               off_right = on_right + stim_duty/100/frequency.Hz,
               filename = first(filename))
    })
  
  stim <-
    stim %>% 
    group_by_all() %>% 
    mutate(V_left = mean(data$stim.V[(data$t.s >= on_left) & (data$t.s < off_left)]),
           V_right = mean(data$stim.V[(data$t.s >= on_right) & (data$t.s < off_right)]),
           torque_at_stim_left = first(data$torque.Nm[data$t.s >= on_left]),
           torque_at_stim_right = first(data$torque.Nm[data$t.s >= on_right])) %>% 
    ungroup()
  
#  stim$filename <- first(filename)
  
  stim <-
    filter(stim, (V_left > 0.1) | (V_right > 0.1))
  
  return(stim)
}

build_is_stim <- function(data) {
  stim1 <- get_stim_cycles(data, data$filename)
  data$isstim <- factor(0, levels = c("left", "right"))

  for (i in seq_len(nrow(cycledata))) {
    data$isstim[(data$t.s >= cycledata$on_left[i]) &
                  (data$t.s < cycledata$off_left[i])] <- "left"
    data$isstim[(data$t.s >= cycledata$on_right[i]) &
                  (data$t.s < cycledata$off_right[i])] <- "right"
  }
  return(data)
}

load_bender_data <- function(filename, cycle_offset=0.1) {
  tryCatch(
    {
      h5file <- H5File$new(filename, mode='r')
      
      # these have lots of rows
      data <- tibble(angle.deg = h5file[['/RawInput/Encoder']][],
                        torque.Nm = h5file[['/Calibrated/xTorque']][],
                        stim.V = h5file[['/RawInput/Left stim']][],
                        t.s = h5file[['/NominalStimulus/t']][],
                        t.norm = h5file[['/NominalStimulus/tnorm']][],
                        cycle = floor(t.norm-cycle_offset))

      # all of these things are single values that will get repeated for each row
      data <-
        cbind(data,
              tibble(filename = basename(filename),
                     fullpathname = file.path(normalizePath(dirname(filename)), filename),
                     stimulus_type = h5attr(h5file[['/ParameterTree/Stimulus']], 'Type'),
                     amplitude.deg = h5attr(h5file[['/NominalStimulus']], 'Amplitude'),
                     dclamp.m = h5attr(h5file[['/ParameterTree/Geometry']], 'Distance between clamps'),
                     width.m = h5attr(h5file[['/ParameterTree/Geometry/Cross-section']], 'width'),
                     height.m = h5attr(h5file[['/ParameterTree/Geometry/Cross-section']], 'height'),
                     sample_frequency.Hz = h5attr(h5file[['/RawInput']], 'SampleFrequency')))

      if (slice(data,1)$stimulus_type == "Sine") {
        data <- cbind(data,
                      tibble(frequency.Hz = h5attr(h5file[['/NominalStimulus']], 'Frequency'),
                             ncycles = h5attr(h5file[['/NominalStimulus']], 'Cycles'),
                             stim_phase = h5attr(h5file[['/ParameterTree/Stimulus/Parameters/Activation']], 'Phase'),
                             stim_start_cycle = h5attr(h5file[['/ParameterTree/Stimulus/Parameters/Activation']], 'Start cycle'),
                             stim_duty = h5attr(h5file[['/ParameterTree/Stimulus/Parameters/Activation']], 'Duty')))
        
        data <- mutate(data, duration.sec = ncycles / frequency.Hz)
      }
      else if (slice(data,1)$stimulus_type == "Frequency Sweep") {
        data <- cbind(data,
                      tibble(start_frequency.Hz = h5attr(h5file[['/NominalStimulus']], 'StartFrequency'),
                             end_frequency.Hz = h5attr(h5file[['/NominalStimulus']], 'EndFrequency'),
                             duration.sec = h5attr(h5file[['/NominalStimulus']], 'Duration'),
                      ))
      }
      
      data <- mutate(data, curve.invm = angle.deg * pi/180 / first(dclamp.m),
                     torque.Nm.orig = torque.Nm,
                     torque.Nm = torque.Nm.orig - mean(torque.Nm))
    },
    finally = {
      h5file$close_all()
    }
  )

  return(data)
}


get_avg_stiffness <- function(curve, torque) {
  s <- svd(cbind(curve, torque))
  if (all(s$v[,1] < 0)) {
    return(-s$v[2,1])
  }
  else {
    return(s$v[2,1])
  }
}

get_high_curve_stiffness <- function(curve, torque, sign=1) {
  c1 <- sign * curve
  t1 <- sign * torque
  return(median(t1[c1 >= quantile(c1, 0.99)], na.rm = TRUE) / quantile(c1, 0.99))
}

get_low_curve_stiffness <- function(curve, torque, sign=1,
                                    low_curve_vel=0.1, low_curve=0.1) {
  if (sum(!is.na(torque) & !is.na(curve)) == 0) {
    return(NA)
  }
  # quick central difference derivative
  dc <- sign * c(NA, curve[3:length(curve)] - curve[1:(length(curve)-2)], NA)
  
  low_dc <- quantile(abs(dc), low_curve_vel, na.rm = TRUE)
  low_curve <- quantile(abs(curve), low_curve, na.rm = TRUE)
  
  is_low <- (dc >= low_dc) & (abs(curve) <= low_curve)
  if (all(!is_low, na.rm = TRUE)) {
    return(NA)
  }
  else {
    trng <- diff(range(torque[is_low], na.rm = TRUE))
    crng <- diff(range(curve[is_low], na.rm = TRUE))
    
    return(trng / crng)
  }
}

get_cycle_data <- function(rawdata,
                           low_angvel = 10,
                           low_curve = 0.2) {
  cycledata <- rawdata %>%
    group_by(cycle) %>%
    summarize(
      work = pracma::trapz(-angle.deg * pi/180, torque.Nm),
      EI1 = get_avg_stiffness(curve.invm, torque.Nm.ctr),
      # stiffness when curvature is large = max torque / max curve
      # add some robustness by using the 99th percentile rather than max
      EIl1 = median(torque.Nm.ctr[curve.invm >= quantile(curve.invm, 0.99)], na.rm = TRUE) / quantile(curve.invm, 0.99),
      EIl2 = median(torque.Nm.ctr[curve.invm <= quantile(curve.invm, 0.01)], na.rm = TRUE) / quantile(curve.invm, 0.01))
      
  rawdata <- mutate(rawdata,
                    dangle.degs = (lead(angle.deg) - lag(angle.deg)) / (lead(t.s) - lag(t.s)))
  
  EIm1 <- rawdata %>%
    filter((dangle.degs > low_angvel) & (abs(curve.invm) < low_curve)) %>%    # gets high positive angular velocity and low curvature
    group_by(cycle) %>%
    summarize(EIm1 = diff(range(torque.Nm)) / diff(range(curve.invm))) # this is the slope of torque vs curvature when curvature is small
  
  EIm2 <- rawdata %>%
    filter((dangle.degs < -low_angvel) & (abs(curve.invm) < low_curve)) %>%    # gets high negative angular velocity and low curvature
    group_by(cycle) %>%
    summarize(EIm2 = diff(range(torque.Nm)) / diff(range(curve.invm))) # this is the slope of torque vs curvature when curvature is small
  
  cycledata <- cycledata %>%
    left_join(EIm1, by="cycle") %>%
    left_join(EIm2, by="cycle")

  stim <- tibble(cycle=seq())
  return(cycledata)
}

