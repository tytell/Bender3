library('hdf5r')
library('dplyr')
library('purrr')
# library('signal')
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

load_bender_data <- function(filename, cycle_offset=0, filterorder=7, filtercutoff=50.0) {
  tryCatch(
    {
      h5file <- H5File$new(filename, mode='r')
      
      # these have lots of rows
      data <- tibble(torque.Nm = h5file[['/Calibrated/xTorque']][],
                     t.s = h5file[['/NominalStimulus/t']][],
                     t.norm = h5file[['/NominalStimulus/tnorm']][],
                     cycle = floor(t.norm-cycle_offset))
      
      if ('Encoder' %in% names(h5file[['/RawInput']])) {
        data$angle.deg <- h5file[['/RawInput/Encoder']][]
      } else if ('Encoder' %in% names(h5file[['/Calibrated']])) {
        data$angle.deg <- h5file[['/Calibrated/Encoder']][]
      }
      
      if ('Left stim' %in% names(h5file[['/RawInput']])) {
        data$stim.V = h5file[['/RawInput/Left stim']][]
      } else if ('activation_monitor' %in% names(h5file[['/RawInput']])) {
        data$stim.V = h5file[['/RawInput/activation_monitor']][]
      }
      
      if ('ParameterTree' %in% names(h5file)) {
        # data from Bender2
        data <-
          data %>%
          bind_cols(
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
      } else {
        data <-
          data %>%
          bind_cols(
            tibble(filename = basename(filename),
                   fullpathname = file.path(normalizePath(dirname(filename)), filename),
                   stimulus_type = h5attr(h5file[['/NominalStimulus']], 'Type'),
                   dclamp.m = h5attr(h5file, 'ClampDistance_mm') / 1000,
                   width.m = h5attr(h5file, 'FishCrossSectionWidth_mm')/1000,
                   height.m = h5attr(h5file, 'FishCrossSectionHeight_mm')/1000,
                   sample_frequency.Hz = h5attr(h5file[['/RawInput']], 'SampleFrequency')))
        
        filetype <- h5attr(h5file[['/NominalStimulus']], 'Type')
        if (filetype == "Frequency Sweep") {
          # t.norm is in radians in frequency sweeps, just to be confusing
          dur <- h5attr(h5file[['/NominalStimulus']], 'Duration')
          data <- mutate(data,
                         t.norm = t.norm / (2*pi),
                         frequency.Hz = (lead(t.norm) - lag(t.norm)) / 
                           (lead(t.s) - lag(t.s)))
          ncyc <- data %>%
            filter(t.s < dur) %>%
            summarize(ncyc = last(t.norm)) %>%
            pull(ncyc)
          
          data <- mutate(data,
                         t.norm = case_when(t.norm < 0.25  ~  NA_real_,
                                            t.norm > ncyc - 0.25  ~  NA_real_,
                                            TRUE  ~  t.norm),
                         frequency.Hz = if_else(is.na(t.norm), NA_real_, frequency.Hz))
        }

        if ("Amplitude" %in% h5attr_names(h5file[['/NominalStimulus']])) {
          data$amplitude.deg = h5attr(h5file[['/NominalStimulus']], 'Amplitude')
        }
        data <- mutate(data, 
                       curve.invm = angle.deg * pi/180 / first(dclamp.m),
                       halfcycle = floor(t.norm*2)/2)
        
        sampfreq <- first(pull(data, sample_frequency.Hz))
        
        bf <- signal::butter(filterorder, filtercutoff / (0.5*sampfreq), type = 'low')
        data <- mutate(data,
                       torque.Nm.orig = torque.Nm,
                       torque.Nm = signal::filtfilt(bf, torque.Nm.orig),
                       curve.invm.orig = curve.invm,
                       curve.invm = signal::filtfilt(bf, curve.invm.orig))
      }
    },
    finally = {
      h5file$close_all()
    }
  )

  data <-
    data %>%
    mutate(torquebias.Nm = mean(torque.Nm[t.s < -0.5], na.rm = TRUE),
           torque.Nm.ctr = torque.Nm - torquebias.Nm,
           curverate.invms = (lead(curve.invm) - lag(curve.invm)) / (lead(t.s) - lag(t.s)))

  return(data)
}

load_cycle_data <- function(filename, tsdata) {
  tryCatch(
    {
      h5file <- H5File$new(filename, mode='r')

      filetype <- h5attr(h5file[['/NominalStimulus']], 'Type')
      if (filetype == 'Frequency Sweep') {
         cycledata <-
           tsdata %>%
           filter(t.s >= 0 & !is.na(halfcycle)) %>%
           group_by(halfcycle) %>%
           summarize(
             t_cycle.s = first(t.s),
             amplitude.deg = max(abs(angle.deg)),
             is_active = factor('passive', levels=c('passive', 'active')),
             activation_duty = 0,
             activation_phase = NA_real_,
             t_act_on.s = NA_real_,
             t_act_off.s = NA_real_,
             t_act_on.norm = NA_real_,
             t_act_off.norm = NA_real_,
             curvature.invm = amplitude.deg * pi/180 / first(tsdata$dclamp.m))
         
         # we use 0.5 / (stuff) here, because we're grouping things by half cycles
         cycledata <-
           cycledata %>%
           mutate(frequency.Hz = 0.5 / (lead(t_cycle.s) - t_cycle.s),
                  side = factor(case_when(halfcycle %% 1 < 0.1  ~  'left',
                                          halfcycle %% 1 >= 0.4  ~  'right'),
                                levels = c('left', 'right')))
         # the first value is usually messed up because of the ramp to position
         cycledata$frequency.Hz[1] = NA_real_
      }
      else if ('AmplitudeByCycle' %in% h5attr_names(h5file[['/NominalStimulus']])) {
        cycledata1 <- tibble(amplitude.deg = h5attr(h5file[['/NominalStimulus']], 'AmplitudeByCycle'),
                            frequency.Hz = h5attr(h5file[['/NominalStimulus']], 'FrequencyByCycle'),
                            is_active = h5attr(h5file[['/NominalStimulus']], 'IsActiveByCycle'),
                            activation_duty = h5attr(h5file[['/NominalStimulus']], 'ActualActivationDuty')) %>%
          mutate(is_active = factor(if_else(is_active, 'active', 'passive'), levels=c('passive', 'active')))
        
        cycledata1$activation_phase <- h5attr(h5file[['/NominalStimulus']], 'ActivationPhase')
        cycledata1 <-
          mutate(cycledata1, 
                 halfcycle = seq(0, length.out = nrow(cycledata1)),
                 side = factor('left', levels=c('left', 'right')))
        
        cycledata2 <- mutate(cycledata1,
                             halfcycle = halfcycle+0.5,
                             side = factor('right', levels=c('left', 'right')))
    
        cycledata <- bind_rows(cycledata1, cycledata2)
        
        tcycle <-
          tsdata %>%
          filter(t.s >= 0) %>%
          group_by(halfcycle) %>%
          summarize(t_cycle.s = min(t.s))
        
        cycledata <- left_join(cycledata, tcycle, by="halfcycle")
  
        cycledata <-
          cycledata %>%
          mutate(t_act_on.s = t_cycle.s + (0.25 + activation_phase)/frequency.Hz,
                 t_act_off.s = t_act_on.s + activation_duty/frequency.Hz,
                 t_act_on.norm = halfcycle + 0.25 + activation_phase,
                 t_act_off.norm = t_act_on.norm + activation_duty,
                 activation_phase = case_when(side == 'right'  ~  activation_phase+0.5,
                                              side == 'left'   ~  activation_phase),
                 curvature.invm = amplitude.deg * pi/180 / first(tsdata$dclamp.m))
        
        # Set the activation times to NA if there wasn't any activation in this cycle
        cycledata <-
          cycledata %>%
          mutate(across(c(t_act_on.s, t_act_off.s, activation_phase), 
                        ~ if_else(is_active == 'active', .x, NA_real_))) %>%
          arrange(halfcycle)
      }
      else {
        ncycles <- h5attr(h5file[['/NominalStimulus']], 'Cycles')
        
        cycledata1 <- tibble(
          cycle = seq(0, length.out = ncycles),
          halfcycle = seq(0, length.out = ncycles),
          amplitude.deg = h5attr(h5file[['/NominalStimulus']], 'Amplitude'),
          frequency.Hz = h5attr(h5file[['/NominalStimulus']], 'Frequency'),
          activation_duty = h5attr(h5file[['/NominalStimulus']], 'ActivationDuty'))

        actstartcycle = h5attr(h5file[['/NominalStimulus']], 'ActivationStartCycle')
        cycledata1 <- cycledata1 %>%
          mutate(is_active = case_when(cycle < actstartcycle  ~  'passive',
                                       TRUE  ~  'active'),
                 is_active = factor(is_active, levels=c('passive', 'active')),
                 activation_phase = h5attr(h5file[['/NominalStimulus']], 'ActivationStartPhase'),
                 side = factor('left', levels=c('left', 'right')))
        
        cycledata2 <- mutate(cycledata1,
                             halfcycle = halfcycle+0.5,
                             side = factor('right', levels=c('left', 'right')))
        
        cycledata <- bind_rows(cycledata1, cycledata2)
        
        tcycle <-
          tsdata %>%
          filter(t.s >= 0) %>%
          group_by(halfcycle) %>%
          summarize(t_cycle.s = min(t.s))
        
        cycledata <- left_join(cycledata, tcycle, by="halfcycle")
        
        cycledata <-
          cycledata %>%
          mutate(t_act_on.s = t_cycle.s + (0.25 + activation_phase)/frequency.Hz,
                 t_act_off.s = t_act_on.s + activation_duty/frequency.Hz,
                 t_act_on.norm = halfcycle + 0.25 + activation_phase,
                 t_act_off.norm = t_act_on.norm + activation_duty,
                 activation_phase = case_when(side == 'right'  ~  activation_phase+0.5,
                                              side == 'left'   ~  activation_phase),
                 curvature.invm = amplitude.deg * pi/180 / first(tsdata$dclamp.m))
        
        # Set the activation times to NA if there wasn't any activation in this cycle
        cycledata <-
          cycledata %>%
          mutate(across(c(t_act_on.s, t_act_off.s, activation_phase), 
                        ~ if_else(is_active == 'active', .x, NA_real_))) %>%
          arrange(halfcycle)
      }
    },
    finally = {
      h5file$close_all()
    }
  )
  
  return(cycledata)
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

get_avg_damping <- function(curverate, torque) {
  if (any(is.na(curverate)) | any(is.na(torque))) {
    return(NA_real_)
  }
    
  s <- svd(cbind(curverate, torque))
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

calculate_cycle_data <- function(rawdata,
                           low_curvevel = 1,
                           low_curve = 0.2) {
  rawdata <- rawdata %>%
    mutate(midcycle = floor((t.norm-0.25)*2)/2 + 0.25)
  
  cycledata <- rawdata %>%
    group_by(halfcycle) %>%
    summarize(
      work = pracma::trapz(-angle.deg * pi/180, torque.Nm),
      EI1 = get_avg_stiffness(curve.invm, torque.Nm.ctr),
      # stiffness when curvature is large = max torque / max curve
      # add some robustness by using the 99th percentile rather than max
      EIL = median(torque.Nm.ctr[abs(curve.invm) >= quantile(abs(curve.invm), 0.98)], na.rm = TRUE) / quantile(abs(curve.invm), 0.98),
      EIL = abs(EIL),
      etaI1 = get_avg_damping(curverate.invms, torque.Nm.ctr),
      #damping when the curvature rate is maximal
      curveratemax.invms = quantile(abs(curverate.invms), 0.98, na.rm=TRUE),
      etaIL = median(torque.Nm.ctr[abs(curverate.invms) >= curveratemax.invms], na.rm = TRUE) / 
        curveratemax.invms,
      etaIL = abs(etaIL))
      # EIl2 = median(torque.Nm.ctr[curve.invm <= quantile(curve.invm, 0.01)], na.rm = TRUE) / quantile(curve.invm, 0.01))
  
  # here we have to group differently, because we're looking for the middle of each cycle (from phase 0.25 to 0.75 and 0.75 back to 0.25),
  # which is when the velocities are high and the curvature is near zero
  EIM <- rawdata %>%
    filter((abs(curverate.invms) > low_curvevel) & (abs(curve.invm) < low_curve)) %>%    # gets high curvature velocity and low curvature
    group_by(midcycle) %>%
    summarize(
      # this is the slope of torque vs curvature when curvature is small 
      ## TODO: maybe think if this is the best algorithm
      EIM = diff(range(torque.Nm.ctr)) / diff(range(curve.invm)), 
      # this spans halfcycles, but we'll attribute the stiffness to the first half cycle
      halfcycle = first(midcycle)-0.25)
  
  # here we want small curvature rates, but large curvature
  etaIM <- rawdata %>%
    filter((abs(curverate.invms) < low_curvevel) & (abs(curve.invm) > low_curve)) %>%    # gets low curvature velocity and high curvature
    group_by(midcycle) %>%
    summarize(
      # this is the slope of torque vs curvature when curvature is small
      etaIM = diff(range(torque.Nm.ctr)) / diff(range(curverate.invms)), 
      # this spans halfcycles, but we'll attribute the stiffness to the first half cycle
      halfcycle = first(midcycle)-0.25)

  # EIm2 <- rawdata %>%
  #   filter((dangle.degs < -low_angvel) & (abs(curve.invm) < low_curve)) %>%    # gets high negative angular velocity and low curvature
  #   group_by(cycle) %>%
  #   summarize(EIm2 = diff(range(torque.Nm)) / diff(range(curve.invm))) # this is the slope of torque vs curvature when curvature is small
  
  cycledata <- cycledata %>%
    left_join(EIM, by="halfcycle") %>%
    left_join(etaIM, by="halfcycle")

  return(cycledata)
}

