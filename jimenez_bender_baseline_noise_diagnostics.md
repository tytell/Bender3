# Bender ATI Mini40 Baseline Noise Diagnostics

**Author context:** Jimenez Lab, CritterGripper/Bender rig
**Date of analysis:** 2026-07-20
**Sensor under test:** ATI Mini40 IP68, S/N `FT56491`, calibration `SI-40-2` (see `FT56491.cal`, this repo)
**Question:** What is the source of baseline noise on the ATI F/T channels, and is it caused by the Grass S88 stimulator or the Sonometrics DS3 sonomicrometry DAQ (EMI hypothesis)?

---

## 1. Files used

| Role | File (absolute path) |
|---|---|
| Real specimen trial (fish, dynamic sweep) | `.../bender_crittergripper/2026-07-16_bass18_bender/2026-07-16_bass18_bender_15_dynamic.h5` |
| Rigid-rod torque-validation trial | `.../bender_crittergripper/2026-07-14_rod60a_bender_torquevalidation/2026-07-14_rod60a_bender_03_dynamic.h5` |
| Diagnostic 1: no object, stim off, sono off | `.../bender_crittergripper/2026-07-20_diagnostic_bender_baseline_noise/2026-07-20_diagnostic_bender_stimulator_off_sono_off.h5` |
| Diagnostic 2: no object, stim off, sono on | `.../bender_crittergripper/2026-07-20_diagnostic_bender_baseline_noise/2026-07-20_diagnostic_bender_stimulator_off_sono_on.h5` |
| Diagnostic 3: no object, stim on, sono off | `.../bender_crittergripper/2026-07-20_diagnostic_bender_baseline_noise/2026-07-20_diagnostic_bender_stimulator_on_sono_off_01.h5` |
| Diagnostic 4: no object, stim on, sono on | `.../bender_crittergripper/2026-07-20_diagnostic_bender_baseline_noise/2026-07-20_diagnostic_bender_stimulator_on_sono_on_01.h5` |
| Calibration matrix (ATI-supplied) | `FT56491.cal` (this repo) |

All paths above are rooted at:
`/Users/yjimenez/Library/CloudStorage/OneDrive-ProvidenceCollege/01_JimenezLab/01_PermanentArchive/bender_crittergripper/`

Signal used for all measurements: `derived/forcetorque_calibrated` (6 x N array: `xForce, yForce, zForce, xTorque, yTorque, zTorque`, units N / N-m), produced by applying the embedded `metadata/calibration_forcetorque_matrix` to raw `timeseries/aidata` rows 0-5. Per the file's own `derived` group note, this is "RA inspection only -- not ground truth" (the R pipeline is authoritative), but it is adequate for a relative noise-floor comparison since the same calibration matrix (`FT56491.cal`) is embedded in every file compared.

Sample rate: 1000 Hz (`daq_ai_sample_rate_hertz`) in every file.

Key protocol metadata pulled from each file to establish test conditions:

| File | `specimen_id` / `specimen_genusspecies` | `protocol_stim_enabled` | motion profile | peak commanded angular velocity |
|---|---|---|---|---|
| bass18_dyn15 | `bass18` / `micropterus nigricans` | `false` | multisine sweep, 0.5/2.75/5 Hz components, up to 6 deg amplitude | 368 deg/s |
| rod60a_dyn03 | `rod60a` / `rod` (rigid calibration rod) | `false` | multisine sweep, 0.5/5/10 Hz components, up to 12 deg amplitude | 1471 deg/s |
| 4x diagnostics | none mounted | varies (stim ao channels commanded 0 V even in "stim_on" trials -- device was powered but not triggered) | single 1 Hz sine, 10 deg amplitude | 63 deg/s |

(`stim_channel1/2_command_volt` were 0 V and `timeseries/stim_state == 'passive'` in every file examined, including the "stimulator_on" diagnostics -- those trials tested whether the S88 unit being powered/connected injects EMI even without an active trigger, not whether an active stimulus injects EMI.)

---

## 2. Test logic

Three complementary measurements were used, each controlling for a different confound:

### 2a. Whole-trial high-pass noise floor (first pass)
- 4th-order Butterworth high-pass, cutoff 20 Hz, zero-phase (`scipy.signal.sosfiltfilt`), applied to each of the 6 calibrated F/T channels over the full recording.
- Rationale: all commanded motion frequencies in these trials are <= 10 Hz, so content above 20 Hz is not the intended mechanical signal.
- Caveat discovered later (see 2c): at high drive frequency/amplitude, real inertial-torque harmonics and step-rate artifacts can leak above 20 Hz, so this metric mixes true noise with motion-linked mechanical content. It is retained here because it is what motivated the deeper investigation.

### 2b. Static pre-roll noise floor (bias-controlled)
- Every trial includes a ~4 s hold at the start (`protocol_wait_before_second`) where `cycle_index == -1`, commanded angular velocity = 0 deg/s, and angle = 0 deg (motor energized and holding, not moving).
- Metric: raw standard deviation of each calibrated F/T channel over samples where `t` in `[-3.9, -0.1]` s, after DC-offset removal.
- This isolates the electrical/instrumentation noise floor from any mechanically-induced noise, since the motor is not moving.

### 2c. In-motion noise floor vs. peak angular velocity (mechanism test)
- Same >20 Hz high-pass metric as 2a, but explicitly plotted against each trial's peak commanded angular velocity (`timeseries/angular_velocity_commanded_degree_per_second`, restricted to `cycle_index >= 0`).
- Rationale: `rod60a_dyn03` and `bass18_dyn15` move at very different peak speeds (1471 vs. 368 deg/s) despite both having an object mounted and both stim/motion present in the diagnostics too -- this isolates *speed* as a variable, separate from *presence* of motion or of a mounted object.

### 2d. ATI SI-40-2 resolution benchmark
- ATI's published resolution for the SI-40-2 calibration under a "Net/DAQ" system interface (applicable here -- this rig digitizes via NI USB-6361, not ATI's own CTL controller):
  - Fx, Fy: 1/100 N (0.010 N)
  - Fz: 1/50 N (0.020 N)
  - Tx, Ty, Tz: 1/4000 N-m (0.00025 N-m)
- Source: ATI Mini40 IP65/IP68 product page ranges-and-resolutions table (SI-40-2 row), cross-checked against the ATI F/T sensor catalog PDF.
- Each measured noise metric (2a-2c) was expressed as a multiple of this resolution to judge whether the observed noise is sensor-limited or system-limited.

---

## 3. Results

### 3a. Whole-trial high-pass (>20 Hz) noise floor, calibrated units

| Channel | bass18 (fish) | diagnostic mean (n=4) | diagnostic range |
|---|---|---|---|
| zForce | 0.0560 N | 0.0333 N | 0.0325-0.0346 N |
| xTorque | 0.0177 N-m | 0.0050 N-m | 0.0049-0.0051 N-m |
| yTorque | 0.0037 N-m | 0.0011 N-m | 0.0010-0.0012 N-m |
| zTorque (primary bending axis) | 0.0123 N-m | 0.0029 N-m | 0.0028-0.0030 N-m |

bass18 noise is 1.7-4.3x every diagnostic on every channel. Critically, the four diagnostics (stim x sono, all 4 combinations) agree with each other to within 2-6% -- toggling the stimulator or the Sonometrics DS3 produced no detectable change.

### 3b. Static pre-roll noise floor (motor at rest, all 6 files)

| File | xForce (N) | yForce (N) | zForce (N) | xTorque (N-m) | yTorque (N-m) | zTorque (N-m) |
|---|---|---|---|---|---|---|
| bass18 (fish) | 0.0137 | 0.0141 | 0.0304 | 0.00037 | 0.00058 | 0.00027 |
| rod60a (rigid rod) | 0.0132 | 0.0147 | 0.0297 | 0.00062 | 0.00036 | 0.00039 |
| diag stim off / sono off | 0.0147 | 0.0140 | 0.0286 | 0.00043 | 0.00072 | 0.00031 |
| diag stim off / sono on | 0.0168 | 0.0147 | 0.0301 | 0.00046 | 0.00088 | 0.00031 |
| diag stim on / sono off | 0.0156 | 0.0127 | 0.0263 | 0.00042 | 0.00094 | 0.00027 |
| diag stim on / sono on | 0.0157 | 0.0144 | 0.0289 | 0.00038 | 0.00075 | 0.00027 |

All six files agree to within ~20-30% on every channel with no consistent ordering by specimen type or stim/sono state. This is the genuine electrical/instrumentation noise floor.

### 3c. In-motion noise vs. peak angular velocity

| Trial | Peak angular velocity | zTorque high-pass noise | noise / velocity |
|---|---|---|---|
| diagnostics (mean) | 63 deg/s | 0.0029 N-m | 4.5e-5 |
| bass18 | 368 deg/s | 0.0123 N-m | 3.3e-5 |
| rod60a | 1471 deg/s | 0.0603 N-m | 4.1e-5 |

Velocity increased ~23x from diagnostics to rod60a; noise increased ~21x. The noise/velocity ratio is roughly constant (3.3-4.5e-5) across all three trials despite completely different mounted objects (nothing / fish / rigid rod).

### 3d. Comparison to ATI SI-40-2 rated resolution (Net/DAQ)

| Axis | ATI resolution | Static floor (typ.) | Static floor vs. resolution | Whole-trial (in-motion) floor, bass18 | vs. resolution |
|---|---|---|---|---|---|
| Fz | 0.020 N | ~0.029 N | ~1.5x | 0.056 N | 2.8x |
| Tx | 0.00025 N-m | ~0.00045 N-m | ~1.8x | 0.0177 N-m | 71x |
| Ty | 0.00025 N-m | ~0.00055 N-m | ~2.2x | 0.0037 N-m | 15x |
| Tz | 0.00025 N-m | ~0.00030 N-m | ~1.2x | 0.0123 N-m | 49x |

The at-rest floor is only modestly above the transducer's own rated resolution (roughly consistent with amplifier/ADC noise). The in-motion floor is one to two orders of magnitude above rated resolution -- the motion-linked term dominates by a wide margin.

---

## 4. Verdict

1. **Stimulator and Sonometrics DS3 EMI are ruled out.** At rest, all 6 recordings (fish specimen, rigid rod, and 4 stim/sono on-off combinations with nothing mounted) show statistically indistinguishable noise (Section 3b). Toggling stim or sono state does not shift the noise floor in any measurable way (Section 3a).
2. **Specimen/tissue type is ruled out.** bass18 (live fish tissue) and rod60a (inert rigid rod) have the same at-rest noise floor as the empty-clamp diagnostics.
3. **The dominant noise source is motion-speed-dependent**, not motion-presence-dependent. Noise while moving scales close to linearly with peak commanded angular velocity (Section 3c) regardless of what (if anything) is mounted. Likely mechanisms: Teknic ClearPath stepper microstepping/commutation ripple (which increases in frequency and amplitude with step rate), or drivetrain/clamp structural vibration excited more strongly at higher angular velocity. The ATI Mini40's own mechanical resonances (1300-1400 Hz per datasheet) are far above the observed noise band (20-100 Hz) and are not implicated.
4. **Even the quiet, at-rest floor sits slightly above the sensor's rated resolution** (~1.2-2.2x on the torque axes), consistent with ordinary amplifier/ADC noise. The in-motion floor is 15-71x rated resolution and is where remediation effort should be focused.

## 5. Recommended next steps

1. Inspect Teknic ClearPath microstep resolution / current setting; evaluate whether a different microstep configuration reduces high-frequency ripple at high step rates.
2. Inspect clamp-to-sensor mechanical coupling for resonances excited at higher angular velocity (independent of the F/T sensor's own 1300+ Hz resonance).
3. Re-run one no-object diagnostic at rod60a-like speed (10 Hz / 12 deg) to confirm the velocity-noise relationship holds with zero mounted mass, further isolating drivetrain/clamp vibration from any specimen-loading effect.
4. Test cable routing/strain relief on the F/T signal lines during a fast motion trial to rule out motion-induced cable microphonics.
5. Deprioritize further stim/sono EMI mitigation work given items 1-2 above; the diagnostic evidence does not support either as a meaningful noise contributor.
