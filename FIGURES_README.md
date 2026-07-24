# Figure naming convention (adopted 2026-07-19)

Figures live outside the repo, on the shared OneDrive hub:
`proj_crittergripper/02_processed/`. This file documents the naming scheme so
a filename tells you what it is without opening it. It replaces ad hoc
suffixes (`_interpbaseline`, `_vector`, `_L0_sono`, `_snr_passing`) that had
grown inconsistently across scripts.

All folder paths below are single-sourced in `R/paths_config.R` -- if the
OneDrive layout moves again, update that file only (this doc + the
`R/*.R` scripts read from it, not from hardcoded strings).

## Folder layout
Renamed 2026-07-19/2026-07-20/2026-07-21 to match the `.cursorrules` "File
placement" convention (`figures/` -> `02_processed/figs_*`;
`bass_summary_figures/` -> `figs_summary/`, matching the "cross-specimen or
pooled summary figures" rule). Flattened 2026-07-21 (PI-directed): every
folder below is FLAT, no subfolders anywhere -- `trial_plots/`/
`summary_plots/` per-fish subfolders were removed; trial-plot and
summary-plot filenames already have distinct naming shapes
(`{bassID}_bender_{NN}_{protocol}...` vs. `{signal}_{protocol}...`), so
merging them into one folder creates no collisions.
- `figs_{bassID}/` — flat: one compound trial plot per trial file
  (`{bassID}_bender_{NN}_{protocol}[_{method}].png`) PLUS that fish's
  pooled summary plots, side by side in the same folder.
- `figs_diagnostic/` — flat: cross-individual pooled superplots, AND
  (2026-07-21) diagnostic/decision-making plots that compare filters or
  calculations rather than reporting a final result (see "Diagnostic vs.
  individual vs. summary" below) -- topic carried in the filename, not a
  subfolder, since these plots are per-decision, not per-fish.
- `figs_summary/` — flattened, curated copies (only SNR-passing
  material) pulled from all three fish for side-by-side browsing; every
  filename here is prefixed with the bass ID since the folder mixes fish.

## Diagnostic vs. individual vs. summary (adopted 2026-07-21)
Three tiers, not just two -- see `analysis_muscle_force_vector_log.md` for
the point-selection design this feeds:
- **Diagnostic** (`figs_diagnostic/`, flat, filename-tagged by topic, not
  per-fish): plots that compare filters/calculations to make a decision
  (e.g. empirical vs. geometric u_hat, legacy vs. vector force, own-step
  vs. interpolated baseline, smoothing cutoffs). Diagnostic is the FIRST
  step of the analysis (PI directive, 2026-07-21) -- everything else is
  downstream of these decisions, so this tier is meant to be comprehensive:
  practically any subjective choice that determines the shape of all future
  data belongs here, not just the choices already made. Once a decision is
  locked, it's recorded in a manifest (`00_records/`, not yet built) rather
  than re-litigated per plot. Current tokens: `musclepullmethod`
  (empirical vs. geometric u_hat line-of-action), `muscleforceestimate`
  (legacy zTorque-only vs. 6-axis vector force), `forceextractionmethod`
  (BUILT 2026-07-21, `bass16_forceextractionmethod.png`,
  `R/diag_force_extraction_baseline.R` -- one representative step per trial
  type. REBUILT 2026-07-22 (was the old full-window MEAN as-of 2026-07-21):
  now shows what CURRENT production does = Method D (narrow-window mean of
  RAW samples centered on the active window's own smoothed-trace peak,
  duration-guarded), SHADING the actual averaging window in green (matching
  `muscleforcemethodcompare.png`) so the averaged samples are a visible
  REGION, not just a horizontal reference value; the old full-window MEAN is
  kept as a secondary blue-dashed comparison line specifically so the
  dynamic-bookend panel's DURATION-GUARD fallback (Method D -> plain mean for
  the ~0.05s bursts) is legible. Distinct from `muscleforceestimate`, which
  is about WHICH axes/method, not which time-window statistic.
  BUG FIXED same day: `.raw_window()` filtered `td` (the WHOLE trial file,
  16 steps) by time-range alone, with no `step_number` filter -- but a
  segmented_finite file's `t.s` RESETS TO 0 AT THE START OF EVERY STEP
  (confirmed empirically, NOT a global file-wide clock), so the "one
  representative step" traces were actually all 16 steps' data concatenated
  together onto the same relative time axis (n_samples was ~801 real vs.
  ~12816-19224 reported -- ~16x inflated). This is what the earlier
  "isometric ripple, plausible unfused tetanus" finding below was based on
  and is now RETRACTED -- the real single-step isometric trace has a much
  smaller ripple than that (PI correctly doubted the original
  overplotting-density explanation and asked for it to be re-investigated).
  Dynamic-bookend panels were never affected (single_finite files have one
  continuous monotonic `t.s` and no `step_number` column at all).),
  `passivebaselinemethod` (BUILT 2026-07-21, `bass16_passivebaselinemethod.png`,
  same script -- STATIC (own-step pre-stim window) vs. INTERPOLATED
  (pre-/post-stim linear interpolation) side by side for isometric and
  isovelocity, PLUS the dynamic L0 bookend's fundamentally different
  runtime-shrunk-window mechanism with no post-baseline counterpart, all
  three trial types in one figure. ENHANCED 2026-07-22 so the flat baseline
  line's FIT to the real passive data is judgeable by eye: the raw samples
  INSIDE each baseline window are now overlaid as bold distinct-colored points
  (pre = red, post = purple) on top of the smoothed trace, plus a zoomed inset
  per window -- the pale raw trace was always plotted but a ~0.2-0.4s window on
  a multi-second axis made its in-window scatter invisible. The dynamic panel
  also carries an explicit in-panel annotation that NO post-baseline window
  exists for single_finite files -- a real limitation, not a rendering gap, and
  deliberately NOT fabricated for symmetry.), `torquesmoothingmethod`
  (this planned item is now covered by `torquesmoothing`, BUILT -- see the
  Diagnostic tier entry below; the two names differ because this list
  predates that entry and was never reconciled),
  `sonosmoothingmethod`, `sonotiminglag`, `forcedevtiming` (dual-tagged,
  see below), `fatiguetimeline` (dual-tagged), `snrfiltereffect`,
  `signconventioncheck` (added 2026-07-21, NOT YET BUILT -- confirms
  left/right stim produces the anatomically-expected force sign; per the
  2026-07-22 sign-convention resolution below, only `force_zTorque_N` is
  actually side-mirrored/anatomically-fixed, so this check applies to that
  axis specifically, not to all four),
  `forcetorquecalcheck` (added 2026-07-21 -- ATI Mini40 `.cal` file
  raw-voltage-to-N/N*m conversion checked for saturation/clipping/drift
  across the dataset), `lengthsignalsource` (added 2026-07-21 -- pooled
  cross-fish version of the per-fish strainValid* commanded/encoder/sono
  cross-checks, decided once instead of per-fish), `geometrysensitivity`
  (added 2026-07-21 -- how much the muscle-force estimate shifts under
  plausible clamp-geometry/moment-arm assumption ranges), `motorstepartifact`
  (added 2026-07-21 -- whether stepper motor step/dir pulses leak
  synchronous noise into F/T or sono channels). The 5 tokens added
  2026-07-21 are DESIGNED/NAMED but not yet built -- no script/figure exists
  for them yet; do not treat their presence here as confirmation the plot
  has been generated.
  `ytorquesignexamples` (BUILT 2026-07-22, `ytorquesignexamples.png`,
  `R/diag_ytorque_sign_examples.R` -- 3 real positive + 3 real negative
  `force_yTorque_N` examples, all SNR >= 4, raw Ty trace with the actual
  passive/active means the pipeline used marked; built as evidence for the
  sign-convention resolution below) and `axissnrcomparison` (BUILT
  2026-07-22, same script -- SNR compared across all 4 vector-force axes on
  the SAME baseline noise floor, pooled and split by trial category) are
  BUILT, distinct from the 5 designed-but-not-built tokens just above.
  CORRECTED same day (PI follow-up on isovelocity's y-torque timing, see
  `ytorqueinertialtiming` below): v1 of `ytorquesignexamples.png` captured
  isovelocity's non-V0 rows via a STATIC pre-stim baseline (motionless vs.
  moving -- guaranteed motion-linked distortion), not the REAL angle-matched
  path (`compute_isovelocity_vector_batch()`, what actually feeds
  `FV_isovelocity_uhatBoth.png`). Regenerated with the correct angle-matched
  computation -- reported percentages dropped from 73-83% negative (wrong,
  motion-contaminated) to 54% negative (correct, angle-matched) for
  isovelocity's moving steps, much closer to a coin flip than static holds'
  97% positive -- see `ytorqueinertialtiming` for why.
  `ytorqueinertialtiming` (BUILT 2026-07-22, `ytorqueinertialtiming.png` +
  `ytorqueinertialtiming_stats.png`, `R/diag_ytorque_inertial_timing.R` --
  PI hypothesis test: "negative y torques in isovelocity happen well after
  the stimulus... velocity ramp suggests inertial noise?" CONFIRMED: the
  SAME Ty deflection, at the SAME time relative to the SAME velocity-profile
  feature, appears in BOTH a stimulated ramp AND a completely unstimulated
  no-stim ramp of the same commanded speed (bass16 step 12 vs. step 20,
  within-trial) -- a kinematic feature of the ramp's motion, not a
  stimulus-locked muscle event. Aggregate check across 123 ramps, 3 fish:
  median lag between the Ty extremum and the angular-velocity extremum is
  ~0.07s WHETHER OR NOT stim is present -- angle-matched passive subtraction
  reduces but does not fully cancel this, leaving a residual in
  `force_yTorque_N` for moving isovelocity steps.
  `apparatusinertiafit` (BUILT 2026-07-22, `apparatusinertiafit.png`,
  `R/diag_apparatus_inertia.R` -- PRE-WIRING evidence for generalizing the
  inertial correction (uniaxial -> multi-axial) and for fixing Gap 1 (apparatus
  MOI is a silent no-op: `calibration_inertia_apparatus_moi_gram_millimeter_
  squared` is NaN in every real trial, so `02_deconvolve.R` drops the apparatus
  term -- which is ~10x the specimen term at real geometries). On the 11-trial
  empty-apparatus corpus: per-channel (all 6 F/T axes) empirical inertia I vs
  angular acceleration (which axes actually carry an I*alpha term), how well
  geometry (aor, width) explains each channel's I (F4 vs F5, PI decision 1), and
  validation of the stored zTorque F4 fit. NOTE on PI decision 3: this corpus
  excites only ONE rotational DOF (bending), so it can identify a per-channel
  I-vs-alpha VECTOR, not a full cross-axis inertia tensor.) and
  `multiaxialinertiacompare` (BUILT 2026-07-22, `multiaxialinertiacompare.png`,
  same script -- on real specimen trials, per-axis raw torque vs
  (raw - sign*I*alpha), with the correction sign auto-chosen on the most
  inertia-dominated window (PI decision 2, since the JSON's stored sign was never
  rig-verified), per-axis variance-reduction quantified, and the yTorque
  residual re-examined against `ytorqueinertialtiming` (PI decision 4). Both are
  PRE-WIRING diagnostics: they change nothing in the production correction path.)
  `specimeninertiacompare` (BUILT 2026-07-22, `specimeninertiacompare.png`,
  `R/diag_specimen_inertia.R` -- PI follow-up to `apparatusinertiafit`: "run a
  similar analysis on specimen inertia using lxwxh frustum dimensions in metadata
  and provided tissue density." Two estimates: (1) ANALYTIC -- solid elliptical
  cylinder from `measurement_specimen_local_body_height/width_millimeter`,
  `measurement_clamp_separation_millimeter`, and
  `measurement_specimen_density_gram_per_cubic_millimeter`, full inertia tensor
  about the sensor origin (bending term validated vs stored
  `calibration_inertia_specimen_moi`; yTorque needs the product of inertia
  I_yz = -m*y_cm*d, ZERO unless the specimen CoM is mediolaterally offset); (2)
  EMPIRICAL -- on passive no-stim ramps, per-channel `channel ~ angle + alpha`
  separates elastic (angle) from inertial (alpha). RESULT: yTorque's alpha
  partial-R^2 on passive ramps is ~0.0003 (three fish) -- yTorque carries NO
  meaningful I*alpha term from apparatus OR specimen, and the corrected traces sit
  exactly on the raw (panel b). So the `ytorqueinertialtiming` residual is NOT an
  inertial (I*alpha) effect; consistent with that diagnostic's own finding that
  the Ty extremum tracks the VELOCITY feature, it is velocity-correlated
  (viscous/damping or centrifugal ~omega^2), which a multi-axial I*alpha
  correction will not fix. yForce/xTorque/zTorque DO carry real inertial terms
  (alpha partial-R^2 0.15-0.25). PRE-WIRING; changes nothing in production.)
  `muscleforcemethodcompare` (BUILT 2026-07-22,
  `muscleforcemethodcompare.png`, `R/diag_force_extraction_methods_compare.R`
  -- PI follow-up on `forceextractionmethod`: "max in window is probably NOT
  the best way... it should be calculated from the smoothed black line."
  4 candidate scalar-extraction methods on the same 3 representative traces
  (isometric/isovelocity/dynamic bookend): A = MEAN full window (current),
  B = MAX raw sample, C = peak of the smoothed trace, D = narrow-window
  (0.15s) mean of RAW samples centered on the smoothed peak's timing. The
  dynamic bookend panel shows A (0.09) substantially underestimating the
  visually obvious peak (~0.22-0.3) because the window's post-peak decay
  dilutes the mean -- C (0.225) and D (0.184) both track it much better,
  D being the more robust/less noise-sensitive of the two.) and
  `muscleforcemethodsensitivity` (BUILT 2026-07-22, same script -- aggregate
  across 92 real SNR-eligible isometric steps, 3 fish: method D tracks A
  almost exactly (near-perfect correlation) for isometric's sustained FLAT
  holds, while B and C both show large, noisy departures -- i.e. the
  extraction-method choice barely matters for isometric specifically, and
  the dynamic-bookend problem above is about TRANSIENT (rise-then-decay)
  windows, not a universal flaw in method A). Both are prototypes (PI has
  not yet picked a replacement method); `forceextractionmethod`'s own
  blue "MEAN" line was also re-labeled the same day to clarify it is the
  ACTIVE window's mean, NOT a baseline (a separate quantity, shown in
  `passivebaselinemethod`) -- this was a real point of PI confusion, not
  just a plot cosmetic.

  **PI decision 2026-07-22: adopted Method D** (narrow-window mean of RAW
  samples, centered on the ACTIVE window's own smoothed-trace peak, peak
  SEARCH restricted to the true stim duration) as the production
  active-force extraction method, REPLACING the plain full-window mean --
  in BOTH the 6-axis vector path (`.mfv_window_peak_means()`,
  `muscle_force_vector.R`, feeding all `FL_isometric_uhatBoth.png`/
  `FV_isovelocity_uhatBoth.png`/`FLsuperplot`/`FVsuperplot` outputs) and the
  legacy zTorque-only path (`.legacy_peak_window_mean()`, `03_analyze.R`'s
  `active_force_Nm`, feeding `muscle_force_Nm`/the frequency-sweep power
  outputs). The PASSIVE/baseline window is unchanged (still a plain mean --
  a steady reference has no "peak" to chase). Full pipeline rerun clean on
  all 3 fish (bass16/17/18) after the switch.

  **Bug found + fixed same day, before promoting Method D**: the dynamic
  trials' L0 bookend contractions (`extract_dynamic_l0_bookends.R`) have a
  genuinely different regime from isometric/isovelocity's stim windows --
  their stim bursts are only ~0.05s (vs ~0.5-1s+ elsewhere), so the
  smoothed trace is often STILL RISING at the search window's own edge,
  pinning the found "peak" to that edge; the fixed-width 0.15s narrow
  averaging window then balloons ~3x past the search boundary into the
  deactivation tail, which can land on a large transient unrelated to the
  tiny bookend pulse (observed >2N spurious spikes on bass18, vs a real
  ~0.1-0.8N range elsewhere, first surfaced as an unexplained spike at 0%
  strain in a regenerated `FLsuperplot`). Fixed by a DURATION guard (not a
  sample-count guard): both `.mfv_window_peak_means()` and
  `.legacy_peak_window_mean()` now fall back to the plain full-window mean
  whenever the search window itself is narrower than the 0.15s narrow
  window, so short/fast bursts (dynamic bookends) get Method A while
  long/sustained or ramping windows (isometric, isovelocity, dynamic
  bookends' own longer siblings if any) get Method D. Re-verified against
  the exact bass18 bookend that produced the spurious spike (now matches
  the plain mean exactly) and reran the full pipeline + both superplots
  again to confirm the spike is gone.
  `strainValidCmd`/`angleValid`/`strainValidSonoEnc`/`strainValidSonoCmd`
  RELOCATED here from the individual tier 2026-07-21 (PI-directed): these
  are rig/system-behavior checks (does the motor track its command? does
  curvature-derived strain track real motion?), not biological results --
  they belong in diagnostic, not `figs_{bassID}/`. Filenames are
  `{bassID}_strainValidCmd.png` etc. (prefixed, unlike most diagnostic
  tokens, since these are still per-fish content living in a flat
  cross-individual folder). `lengthsignalsource` (above) is the FUTURE
  pooled-across-all-fish version of this same check, decided once instead
  of per-fish -- not yet built, distinct from these 4 per-fish files.
- **Individual** (`figs_{bassID}/`): trial-level plots, PLUS per-fish
  aggregated plots (all trials for that fish, trial identity kept visible
  via color/shape) using whichever method won its diagnostic decision.
  `forcedevtiming` and `fatiguetimeline` are dual-tagged -- diagnostic
  output that is ALSO a legitimate individual-tier visualization (within-
  step force development timing; near-L0 force vs. real elapsed session
  time), so they may appear in both places. `forcedevtiming` is
  OVERLAY-ONLY as of 2026-07-21 (PI feedback: the faceted per-step variant
  wasn't useful and was dropped, code and files both). ADDED 2026-07-22
  (PI-requested, to sanity-check the active/passive/residual decomposition
  by eye):   `forcedevtiming_isometric_allsteps.png` /
  `forcedevtiming_isovelocity_allsteps.png` (`R/diag_forcedev_allsteps.R`,
  read-only) -- cross-fish (bass16/17/18 faceted), EVERY step drawn as its
  own line colored by held strain (isometric) / strain rate (isovelocity),
  y = PRODUCTION vector muscle force (`force_ts$muscle_force_vector_N`, RAW
  sign), stim window shaded. These plot the exact quantity the FL/FV
  superplots sample, so the concave-up and the baseline-drift / non-zero
  pre-stim offset that drive it are directly visible per step.
  `isometricbaselinedrift` (BUILT 2026-07-22, `isometricbaselinedrift.png`,
  `R/diag_isometric_baseline_drift.R`, read-only, cross-fish) -- isometric
  baseline-DRIFT test: per step, muscle force under the CURRENT static
  pre-stim baseline (left) vs a pre->post INTERPOLATED baseline (right; same
  linear scheme as `passive_force_Nm_interp`). Shows the ~6 s viscoelastic
  creep the static baseline leaves in (force drifts up for seconds after
  stim-off, scaling with |bend|) and that interpolation returns it to ~0
  after the stim transient. Worked in the single-axis inertia-corrected
  torque domain (carries pre/post windows), NOT the 6-axis vector the FL
  superplot samples -- the drift MECHANISM is torque-level; the interp does
  NOT flatten the in-stim PEAK the superplot samples (peak is early, where
  static ~= interp).
  `concaveupfatiguestim` (BUILT 2026-07-22, `concaveupfatiguestim.png`,
  `R/diag_concaveup_fatigue_stim.R`, read-only, cross-fish) -- tests whether
  the residual in-stim concave-up is FATIGUE or a per-step STIMULUS effect.
  Both ruled out as the cause: stim is a CONSTANT 5.00 V on every step
  (recruitment blocks vary SIDE, not amplitude -> force-per-volt is a no-op);
  fatigue is REAL (L0 force decays across the session, bass17 cor(F,order)
  =-0.66) but cannot produce the arms because WITHIN a single block (fixed
  fatigue + fixed stim) force still rises with |bend| from that block's own
  fresh L0 (bass17 within-block cor(F,|strain|)=+0.88). Each line = one block,
  colored by session order. The arms are the within-block force-|bend|
  residual (passive-subtraction problem), not fatigue/stim.
  `isopassivemodels` (BUILT 2026-07-22, PROTOTYPE, `R/diag_isometric_passive_models.R`,
  read-only, cross-fish, 4 files `isopassivemodels_{1_relaxshape,2_leverage,
  3_zeromusclecontrol,4_flshape}.png`) -- prototypes improved isometric passive
  subtraction and visualizes the reasoning. Compares M0 (static pre-stim mean,
  current) vs M1 (pre->post linear interp, pointwise) vs M2 (viscoelastic
  relaxation loess-fit over quiescent pre+post samples, subtracted POINTWISE),
  all projected on the geometric u_hat (= muscle_force_vector_geom_N). Plot 3 is
  a MODEL-FREE proof (a no-activation quiescent window minus M0 reproduces
  121-169% of the bass16/17 |bend| slope -> the concave-up is mostly a stale-
  baseline artifact). Pointwise M2 removes it without a spurious bell
  (cor(F,|strain|) M0->M2: bass16 +0.19->+0.03, bass17 +0.93->+0.25, bass18
  +0.57->+0.42; bass18's monotonic rise is genuine). CAVEAT: M2 interpolates the
  passive across the unobservable ~0.3 s stim gap; low-force fish are below the
  passive-drift floor and must be magnitude/SNR-gated, not read as a shape. NOT
  a production change -- prototype pending PI decision.
  `isovpassivemodels` (BUILT 2026-07-22, DIAGNOSTIC, `R/diag_isovelocity_passive_models.R`,
  read-only, cross-fish, 3 files `isovpassivemodels_{1_rampshape,2_fvpayoff,
  3_rampstruct}.png`) -- the ISOVELOCITY counterpart to `isopassivemodels`,
  answering how the moving-trial passive logic compares to the isometric fix.
  The PRE-FIX production (compute_isovelocity_vector_batch) subtracted an
  ANGLE+signed-VELOCITY-matched no-stim ramp -- the right raw material -- but
  COLLAPSED it to a scalar window-MEAN, then subtracted that from the Method-D
  (peak) active. Because the active window SWEEPS through angle, the passive
  varies a LOT across it (range median 1.8/2.2/2.4 N bass16/17/18, up to ~5.6 N),
  so the flat mean was a poor stand-in for the passive at the peak's own angle:
  pointwise-minus-mean muscle force differs by up to +-5.6 N. Panel 2 (FV
  payoff): the window-MEAN passive manufactured a CONCAVE-UP FV in low-force
  bass16/17 (same artifact as the FL concave-up), which POINTWISE angle-matched
  subtraction FLATTENS to ~0; bass18 goes from a flat FV to a genuine, cleanly
  POSITIVE curve across the FULL velocity range -- NO negative overshoot
  anywhere (corrected 2026-07-22: the earlier "overshoots negative at high |v|"
  finding was an artifact of (a) `velocity matching` originally only covering an
  EXACT commanded velocity -- bass18's OTHER trials moved at velocities with NO
  stim-off ramp anywhere in the corpus and fell to a meaningless static baseline,
  since fixed with a nearest-same-sign-velocity fallback (`passive_source ==
  "angle_matched_nearest_v"`, see analysis log), and (b) this diagnostic's own
  reimplementation misusing `uhat_geometric(velocity)` -- now uses the fixed
  longitudinal u_hat production actually uses for isovelocity, matching
  `muscle_force_vector_geom_N` exactly). LOGIC vs isometric: isometric
  passive varies only in TIME (relaxation, 1 d.o.f., bracketed -> M2 time-fit);
  isovelocity varies in ANGLE (elastic, large) + velocity + direction, so the
  time-relaxation fit does NOT transfer -- the analog fix is POINTWISE
  angle-matched subtraction (subtract the ramp sample-by-sample by angle, then
  Method D on the delta). IMPLEMENTED IN PRODUCTION 2026-07-22 (PI-approved):
  .mfv_ramp_passive_pointwise() replaced the mean-collapse; velocity matching
  falls back exact -> within-fish nearest-same-sign-velocity -> static (last
  resort, added 2026-07-22 per PI direction to use the stim-off ramps as the
  key ingredient "cooked into the experimental design"). This diagnostic calls
  the SAME production helper (not a reimplementation) and builds an equivalent
  fish-wide no-stim-ramp library, so it remains the canonical BEFORE/AFTER
  record, bit-for-bit consistent with production. Rebuilt FVsuperplot
  geometric-u_hat FV: concave-up artifact gone, all bass18 points positive.
  TARGET SHAPE (PI-clarified 2026-07-22): FV should be a Hill hyperbola
  (monotonic-decreasing with shortening velocity, eccentric > isometric >
  concentric), NOT a bell/peak-at-V=0 -- that's the FL target, not FV's. A
  pairwise eccentric-vs-concentric check (SNR-passing points, 127/255 %/s) is
  directionally suggestive for bass18 -- SEE `isovhillcheck` BELOW: once the
  V=0 isometric anchor is actually plotted, the full curve is NOT a clean
  monotonic Hill curve (confounded by cross-trial fatigue) -- "bass18
  reproduces the Hill relationship" was an overclaim, corrected in the
  analysis log the same day. empirical-u_hat FV stays U-shaped (empirical
  direction is unstable for these low-force moving steps -- unrelated to the
  passive-baseline fix, still pending a PI decision, see analysis log "flag 2").
  `isovhillcheck` (BUILT 2026-07-22, DIAGNOSTIC, `R/diag_isovelocity_hillcheck.R`,
  read-only, cross-fish, 1 file `isovhillcheck.png`) -- PI-requested plot after
  pushing back on a verbal Hill-consistency claim above ("I'm not seeing how
  bass18 resembles the Hill-type relationship"). Plots isovelocity's own V=0
  hold + moving ramps (pointwise angle-matched passive), all points shown,
  alpha-flagged by 4-tier confidence (never dropped), summary line = mean of
  confident + confidently_small points per velocity. ORIGINAL FINDING
  (ratio-only gate): bass18 right is a "W", not a Hill hyperbola -- eccentric
  plateau (~1.65-1.76 N) drops sharply to 0.56 N at V=0, then oscillates
  (1.37/1.03/1.83 N at 127/255/382). REVISED 2026-07-22 (SNR-magnitude
  conflation audit IMPLEMENTED, see `snrmagnitudeaudit` below): the V=0 mean
  is now 0.94 N (n=7, up from n=3) once the stronger trial's confidently-small
  V=0 reps are correctly included instead of dropped as noise -- CONFIRMING
  those reps are real, but only PARTIALLY resolving the notch (0.94 N is
  still below the eccentric plateau and most concentric points). bass16/17
  show a PEAK AT V=0 with decline on both sides (a tent) -- neither Hill's
  plateau-then-decline nor flat noise; unexplained. Honest state: the
  eccentric>concentric pairwise ordering remains directionally suggestive for
  bass18 but is NOT a settled demonstration of the Hill relationship, even
  under the magnitude-aware gate -- a proper test needs each trial's OWN V=0
  vs its OWN moving steps (within-trial only), not pooled across trials --
  not yet built.
  `snrmagnitudeaudit` (BUILT 2026-07-22, DIAGNOSTIC, `R/diag_snr_magnitude_audit.R`,
  read-only, cross-fish, 3 files `snrmagnitudeaudit_{1_forceVsNoiseFloor,
  2_quadrantCounts,3_isovelocityV0notch}.png`) -- PI-requested audit of whether
  the ratio-only `activation_snr >= MFV_UHAT_SNR_MIN` gate (used throughout the
  pipeline -- see analysis_muscle_force_vector_log.md's 2026-07-22 "SNR-based
  confidence gating audit" entry for the full site-by-site writeup) conflates
  "elevated noise floor" with "genuinely small real force." Classifies every
  finalized step/hold/ramp (full 3-fish corpus, isometric + isovelocity V=0 +
  isovelocity moving + dynamic bookend) into 4 quadrants by crossing SNR-pass
  (`activation_snr >= snr_min`) with magnitude-pass (`|force_geom_N| >=
  baseline_force_noise_N`, the SAME two-condition test `mfv_gate_f0()` already
  uses for the F0 denominator) -- "confident" (both pass), "confidently small"
  (SNR fails, magnitude passes -- the conflation case), "unstable magnitude"
  (SNR passes, magnitude fails -- what `mfv_gate_f0` was built to catch), and
  "unconfirmable" (both fail). Plot 1: force vs noise floor directly, log-log,
  quadrant-colored, by category. Plot 2: quadrant counts by category (42
  "confidently small" points and 13 "unstable magnitude" points exist across
  the corpus outside the F0-denominator context `mfv_gate_f0` already
  protects). Plot 3: re-examines the `isovhillcheck` V=0 notch directly --
  bass18's isovelocity V=0 holds are 8/9 "confidently small" + 1/9 "confident",
  ZERO "unconfirmable" (including the strong trial's 1.99/1.31 N reps the
  ratio-only gate used to exclude, which are consistent with its eccentric
  plateau) -- under a magnitude-aware scheme none of bass18's V=0 data is
  noise-indistinguishable, unlike bass16/17 which each have some genuinely
  unconfirmable V=0 reps. This diagnostic itself remains READ-ONLY and
  unchanged; the 4-tier scheme it motivated was IMPLEMENTED in production
  2026-07-22 (see `analysis_muscle_force_vector_log.md`'s "implementation"
  entry the same day for the full site list -- `mfv_confidence_tier()` /
  `MFV_CONFIDENCE_ALPHA` in `muscle_force_vector.R`, consumed by
  `plot_muscle_force_vector.R`, `superplot_fl_pooled_snr_passing.R`,
  `export_snr_summary_figures.R`, `diag_isovelocity_hillcheck.R`,
  `plot_fatigue_timeline.R`, and `superplot_fl_fv_tiers.R`).
  `torquesmoothing` (BUILT earlier, EXTENDED 2026-07-22 (twice, same day),
  DIAGNOSTIC, `R/diag_torque_smoothing.R`, read-only, now runs per-specimen
  via `BENDER3_BASS_ID` env var (bass16/17/18 all run -- was bass16-only,
  "a few examples," before the second extension), 3 files per specimen
  `diag_torque_smoothing_{timedomain,spectrum,examples}.png`) -- PI question,
  prompted by dissatisfaction with the FL/FV pattern strength: "we might need
  to adjust [the smoothing cutoff] per experiment type or angular
  acceleration." Compares raw calibrated torque vs. inertia-corrected vs.
  30/15/8 Hz low-pass options, per step, plus a pooled power spectrum of the
  unsmoothed corrected signal, BY PROTOCOL (isometric = static hold, no
  commanded motion, so any high-frequency content there is unambiguously
  noise; isovelocity = actively moving the whole window, so the same content
  could ALSO be genuine motion-coupled signal).
  FINDING, now confirmed on ALL 3 fish (was bass16-only): isovelocity's
  active-window power spectrum is elevated across the ENTIRE 0-100 Hz band --
  **QUANTIFIED median isovelocity/isometric power ratio: bass16 9783x,
  bass17 3909x, bass18 3947x** (~3.5-4 orders of magnitude, broadband, not
  one narrow spike) -- NOT the same noise floor scaled up, and consistent
  across specimens. Peak frequency also consistent: isometric 18-20 Hz,
  isovelocity 12 Hz, all 3 fish. New `_examples.png` (2 isometric + 2
  isovelocity median-power steps per specimen, legible large panels) shows
  WHY directly: for isometric, all filter options (raw through 8 Hz) are
  visually indistinguishable except for shaving a single-sample outlier
  spike -- the smoothing choice barely matters. For isovelocity, the filter
  choice changes the trace SHAPE dramatically -- large (order ~0.3-0.8 N*m),
  several-cycle oscillations right at motion onset that 8/15 Hz collapse
  into a smooth curve while 30 Hz/unsmoothed preserve the ringing -- whether
  that ringing is real mechanical response or an artifact to remove is
  exactly the open question a single pipeline-wide cutoff can't resolve for
  both protocols at once. Does NOT yet test the "angular acceleration" half
  of the question directly (would need the commanded/encoder bend-angle 2nd
  derivative correlated against the torque residual, per step) -- flagged as
  a follow-up, not done here. PRE-DECISION diagnostic: does not change
  `NOISE_FILTER_HZ` / `DISPLAY_SMOOTH_HZ` or any other production filter
  setting; a per-protocol (or per-motion-state) cutoff is now evidenced as
  plausibly necessary, confirmed across all 3 fish, but not yet decided or
  implemented.
  `fltiers` / `fvtiers` (BUILT 2026-07-22, `R/superplot_fl_fv_tiers.R`, read-only
  re-aggregation of the two pooled builders, 3 files
  `fltiers_1_within_trial.png`, `fvtiers_1_within_trial.png`,
  `fltiers_2_within_protocol_isometriconly.png`) -- FL/FV at the three POOLING
  TIERS so every claim's pooling level is explicit. within-trial (each line =
  one trial's own step series, RAW geometric force, NO cross-trial
  normalization) is the cleanest shape view and exposes that bass17's isometric
  forces are ~0.01 N (at the noise floor). within-protocol isometric-only pooled
  FL is the un-mixed companion to the across-protocol FLsuperplot (grand mean is
  concave-down, peak near L0). FV's within-protocol pool already IS the
  FVsuperplot (isovelocity moving only); the across-protocol pool is the existing
  FLsuperplot/FVsuperplot. Geometric u_hat throughout. Points alpha-flagged by
  4-tier confidence (added 2026-07-22, SNR-magnitude conflation audit -- this
  file previously had NO confidence flagging at all; lines stay at fixed
  alpha, only points are tiered).
  `precondition*`/`isovelocity_*power*`/`isometric_*tension*` (BUILT
  2026-07-22/23/24, DIAGNOSTIC, cross-fish, read-only, several scripts --
  the "dynamic sono validation is tight for isovelocity/isometric but weak
  for dynamic" investigation, see `analysis_muscle_force_vector_log.md`
  for the full multi-day writeup) --
  `dynamic_trial_precondition.R` (shared, not a plotter): hard specimen-
  specific "early (preconditioning)" vs. "later (stable)" trial-order
  cutoff (bass16=5, bass17=9, bass18=5), chosen so trial-mean sono-strain
  offset falls and stays below 1.5 pct-points -- root cause traced to
  early-trial tissue preconditioning, not signal processing (calibration/
  crosstalk, 247 Hz staircase, mechanical vibration, and true phase lag
  were all ruled out first). `diag_precondition_power_check.R`
  (`dynamic_precondition_{offset,power}_vs_trialorder.png`,
  `dynamic_precondition_power_vs_curvature.png`) -- shows the SAME
  early-trial pattern in muscle power output (active dynamic cycles: mean
  power/offset both fall from trial ~1-4 onward). `diag_precondition_
  power_vs_offset.R` (dynamic, `dynamic_precondition_{mean,max}power_
  vs_offset.png` + 2 by-specimen facet versions) -- power DIRECTLY
  correlated against offset (not just both against trial order): pooled
  r=0.73 (mean)/0.63 (max), n=30, holds WITHIN every specimen individually
  (r=0.55-0.98), i.e. not a between-fish confound. `diag_precondition_
  power_vs_offset_isovelocity.R` (`isovelocity_{mean,max}power_
  vs_offset.png`) -- independent confirmation in a 2nd protocol, per-step
  torque x angular-velocity power (same conversion muscle_geometry.R
  already used for a Hill fit's Vmax, applied per-step instead): r=0.92
  (mean)/0.82 (max), n=11. `diag_precondition_tension_vs_offset_
  isometric.R` (`isometric_{mean,max}tension_vs_offset.png`) -- isometric
  has no external power by design (motor doesn't move), so uses SPECIFIC
  TENSION instead (same lever-arm torque->force->N/cm^2 conversion,
  per-step); only 5 isometric trials exist corpus-wide, correlation not
  significant (n=5), and isometric offsets are all small (<1.3 pct-points)
  unlike dynamic's early-trial 2-10 pct-point offsets.
  `diag_precondition_passive_sono_fidelity.R` (`dynamic_precondition_
  passive_sonoValidation_earlyVsLater.png`) -- PI question: does PASSIVE
  fidelity also improve early->later? PASSIVE here means STRICT cycle-level
  (PI-clarified 2026-07-24: no stim on EITHER side, anywhere in the whole
  cycle -- `step_activity == "purely passive"`, not the looser per-sample
  "no/left stim" window used in an earlier version of this script, which
  still mixed in left-stim cycles that have real muscle-driven bending).
  Yes for pointwise r (0.915->0.975) but the MEAN-LEVEL OFFSET stays small
  and FLIPS SIGN in both groups (-0.45 early, +0.16 later) vs. active's
  early offset of +4.92 pct-points, one direction, ~30x larger -- the
  dramatic bias is specific to active (right-side) stimulation, not
  generic early-session mechanics.
  Sanity check (not a separate figure): commanded p2p strain amplitude for
  these dynamic trials is ~5% (median 5.04%); early-trial OFFSETS alone
  are 92-200% of that commanded amplitude (later trials: 10-25%) -- the
  early bias is the same order of magnitude as the entire commanded
  motion, not a subtle measurement artifact.
  `diag_precondition_calibration_gain_check.R` (`dynamic_precondition_
  calibration_{gain,offset}_vs_trialorder.png`) -- PI question: can a
  fixed calibration factor recover true strain from the offset? Fits
  sono~predicted PER TRIAL (dynamic, active/right-stim) to get a gain
  (slope) and offset (intercept). NO: early-trial gain is not even a
  stable relationship (ranges -0.24 to 1.66 trial-to-trial within the
  same specimen, some trials near-zero/negative r) -- it's the trial-order
  decay, not a fixed error. Later/stable gain is tighter WITHIN a specimen
  (bass17 1.24-1.28, bass18 0.99-1.08, bass16 0.93-1.24) but DIFFERS
  BETWEEN specimens (median 1.03/1.27/1.08) -- no universal factor; a
  per-specimen gain on later/stable trials only is plausible for
  bass17/bass18, not for early trials.
  `diag_sono_vs_geometric_dynamic_power.R` (`dynamic_sonoVsGeometric_
  power_vs_trialorder.png`, `..._power_scatter.png`, `..._power_boxplot.png`)
  -- PI reframing: early dynamic trials show REAL muscle behavior at lengths
  NOT prescribed by the motor/geometry; can sono recover valid power from
  them instead of excluding them? Computes dynamic cycle power two ways from
  the SAME measured torque -- GEOMETRIC (commanded-angle kinematics, the
  pipeline's existing method) vs. SONO (directly measured muscle length,
  zero-phase Butterworth LOW-PASS FILTERED at 40 Hz on the full-trial
  continuous series -- "condition the sono signal" per PI direction,
  2026-07-24, using the best-justified candidate from `diag_sono_
  smoothing.R`'s earlier comparison -- THEN decimated to the true
  ~241-247 Hz DS3 rate before differentiating), RIGHT-STIM cycles only, ALL
  trials (early + later, no exclusion applied here). RESULT: the early-trial
  power "inflation" (geometric median 38.3 W/kg early vs. 0.22 W/kg later,
  ~170x) nearly VANISHES under the sono method (-2.49 vs. -0.02 W/kg) --
  strong evidence the inflation is a geometric-model artifact, not real
  muscle output; unchanged by the filtering follow-up. The 40 Hz filter DID
  fix the peak-power noise flagged in the first pass (bass17 trial 4:
  sono max_peak 1,059 -> 151.8 W/kg, now comparable to geometric's 167 W/kg;
  per-cycle peak-power r rises to 0.965 in later trials). REMAINING CAVEAT:
  per-cycle AVG-power agreement between the two methods stays weak-to-
  negative even in later trials (r=-0.27 early, r=-0.68 later) -- filtering
  fixed peak-power noise but not cycle-by-cycle mean-power tracking; both
  methods independently sit near a shared near-zero floor in later trials.
  The boxplot (`_power_boxplot.png`, 2x2: mean/peak rows x geometric/sono
  columns, n=34 early / 89 later cycles) visualizes the trial-median table
  directly -- still a first-pass comparison, not a finished replacement
  pipeline. UPDATE 2026-07-24 (PI: "some early trials have much higher cycle
  power, much closer to Coughlin -- look into that"): added a red Coughlin
  (2000) bass power reference band (14.4+/-1.9 W/kg, derived -- see
  `summary_coughlin2000_bass_comparison.R`) to all 4 facets. CONFIRMS the
  effect is the SAME geometric artifact documented above: in the
  mean/geometric panel, early-trial power sits ABOVE the Coughlin band
  (looks "closer to/exceeding" a real literature benchmark); in the
  mean/sono panel (real, measured muscle length), the SAME early trials'
  median goes slightly NEGATIVE and nowhere near the Coughlin band -- the
  apparent match disappears (and reverses) once corrected with real muscle
  kinematics. Peak power still shows some earlyvs.later elevation even
  after sono correction (median ~40 vs. ~5 W/kg), but well below the
  geometric method's early-trial peak inflation.
  `diag_precondition_sono_length_activeVsPassive.R`
  (`dynamic_precondition_sonoLengthExcess_{vs_trialorder,boxplot}.png`) --
  ground-truths the early-trial "excess shortening" using STRICTLY the
  active-vs-passive difference (PI request), reusing `calc_muscle_torque()`'s
  own phase-matched act-minus-(averaged)-passive machinery UNCHANGED, but
  pointed at the 40 Hz-filtered sono STRAIN instead of a torque channel.
  Output ("sono strain excess") does not depend on L0, the commanded-angle
  assumption, the encoder, or any cross-trial calibration -- it is a
  within-trial, phase-matched comparison against each trial's own passive
  baseline. RESULT: reproduces the same early > later pattern (early
  median 4.81%, mean 5.29%, n=36 cycles; later median 0.20%, mean 0.03%,
  n=89 cycles), decaying toward zero at each specimen's own precondition
  cutoff -- the strongest available confirmation that the excess shortening
  is a real, trial-order-decaying biomechanical property of the muscle, not
  an artifact of any one length/velocity model. Also clarified in the same
  session: the torque source for this ENTIRE investigation (`deconvolve_
  bender()`) is the single primary-bending-axis ("z-torque") channel, NOT
  the multi-axis empirical uHat vector from `muscle_force_vector.R` (a
  separate pipeline feeding the `_uhatBoth` FL/FV plots) -- corrects an
  earlier mis-statement in the decision log, no code changed.
- **Summary** (`figs_summary/`): GENUINELY cross-fish content only --
  either one pooled-across-all-fish panel (`FLsuperplot_*`) or side-by-side
  per-individual panels in one figure (`specimen_comparison_specific_properties.png`),
  never a single-fish plot with just a `{bassID}_` filename prefix.
  Cleaned out 2026-07-21 (PI-directed, "completely unnecessary and
  unhelpful ... should appear in the bass## folders"): `export_snr_summary_figures.R`
  used to copy per-fish trial plots and per-fish `_snrPass` variants in here
  prefixed with the specimen ID -- those were pure duplicates (or, for
  `_snrPass`, single-fish-only content that never belonged here) and are now
  generated/kept ONLY in `figs_{bassID}/`, unprefixed. Raw / individual-mean
  / group-mean cross-fish tiers are still open (point-selection method
  design, see the decision log) -- do not assume a fitted curve vs. binned
  mean without checking that file first.
  `pooled_strainValidSonoEnc_*.png` (BUILT 2026-07-22/23, `R/plot_sono_
  strain_validation_pooled.R`, part of the sono-validation/preconditioning
  investigation above) -- pools measured (sono) vs. predicted (encoder)
  RIGHT-muscle strain across all 3 specimens: one file per protocol family
  (`_isometric`/`_isovelocity`/`_dynamic`/`_frequency_sweep`, active vs.
  passive panels), `_dynamic_later.png` (dynamic restricted to "later
  (stable)" trials only, dynamic_trial_precondition.R cutoff), and two
  ALL-protocol combined figures (protocol_family x column facet):
  `_allProtocols_later.png` (columns = per-sample windowed active/passive
  state) and `_allProtocols_later_stepActivity.png` (columns = per-STEP/
  CYCLE any-stim-anywhere, isolating residual post-stim force/strain from
  truly-never-stimulated units -- revealed isometric has ZERO purely-
  passive steps corpus-wide, and dynamic's stimulated-cycle fit improves,
  not degrades, once the full cycle is included: r 0.905->0.932). Sono is
  DECIMATED (1-in-~4) to the true ~241-247 Hz DS3 update rate before
  plotting/pooling -- the 1 kHz AI clock oversamples the DS3 by ~4x.
  `precondition_power_tension_earlyVsLater.png` (BUILT 2026-07-24,
  `R/summary_precondition_power_tension_earlyVsLater.R`) -- PI question:
  does the early/later exclusion leave too little power/tension data?
  6-panel grid (protocol x mean/max, independent y-axes, n labeled per
  x-tick): dynamic loses close to half its trials to the cutoff (14 early
  / 16 later) and the excluded early trials are the highest-power ones
  in the dataset (power-side signature of the same slip effect); but
  isovelocity (1 early / 10 later) and isometric (1 early / 4 later) are
  barely affected -- almost none of those trials were run early in a
  session to begin with, so there's little to lose either way. CAVEAT:
  the early/later trial-number cutoff was fit and validated on DYNAMIC
  trials only; applying it to isovelocity/isometric trial numbers here is
  a reasonable but not independently re-validated extrapolation (same
  specimen/session, same chronological trial counter).
  `sonoLengthExcess_activeVsPassive_byTrialOrder.png` (BUILT 2026-07-24,
  `R/diag_precondition_sono_length_activeVsPassive.R` -- same script also
  writes the `figs_diagnostic/dynamic_precondition_sonoLengthExcess_*`
  copies) -- cross-fish (3-panel) trial-order view of the ground-truthed
  active-minus-passive sono strain excess. Excess = active-cycle sono
  strain minus PHASE-MATCHED average passive-cycle sono strain, same trial
  (reuses `calc_muscle_torque()`'s act-minus-passive logic, pointed at sono
  strain instead of torque). Independent of L0/commanded-angle/encoder/
  cross-trial calibration. Shows the early > later excess decaying to zero
  at each specimen's own precondition cutoff (early median 4.81% vs. later
  0.20%) -- the headline confirmation that early-trial excess shortening is
  a real, trial-order-decaying muscle property, promoted here as the
  cross-fish summary of that finding.
  `sono_conditioning_not_fatigue.png` (BUILT 2026-07-24,
  `R/summary_sono_conditioning_not_fatigue.R`) -- PI request: show, from the
  data, that good sono requires CONDITIONING (a modest dose), not a fatigued
  muscle. 2 panels. (A) |active-vs-passive sono strain excess| vs. cumulative
  dynamic conditioning cycles delivered before each trial -- reads off the
  dose (excess drops into a <1 pct-pt "good" band after ~10-20 cycles;
  lower bound, other protocols also condition). (B) DECOUPLING scatter:
  |excess| vs. peak active force capacity (max |muscle_torque.Nm| per trial,
  the standard fatigue readout, as %% of each specimen's own session max) --
  low excess occurs across the FULL force range. Decisive case: bass17 ends
  its session at ~99%% of peak strength with the excess already resolved
  (dose to good ~12 cycles at ~89%% force), so stabilization is a
  conditioning-dose effect, not force loss. HONEST CAVEATS: n=3; bass16 and
  bass18 are noisier/sparser and their sessions ended at lower force
  (possibly real fatigue for those two) -- bass17, the specimen with the
  most sustained dynamic protocol, carries the clean decoupling.
  `FL_isometric_sonoVsCommanded.png` / `FV_isovelocity_sonoVsCommanded.png`
  (BUILT 2026-07-24, `R/prototype_fv_fl_sono_length.R`) -- PI request 2A:
  prototype FL/FV curves with the x-axis (length/velocity) recomputed from
  the REAL muscle length (40 Hz-filtered sono) instead of the commanded
  operating_point; force (y) unchanged (active-minus-passive primary-axis
  torque). RIGHT muscle only (sono's only channel). RESULT: for these
  protocols the commanded and sono x-axes AGREE tightly (isometric r=0.992
  in value AND scale; isovelocity r=0.972 in rank) -- i.e. the existing
  commanded-based FL/FV curves are already faithful here, because
  isometric/isovelocity were run under good conditioning; the dramatic
  sono-vs-commanded divergence was specific to EARLY DYNAMIC trials, not
  these. FV PROTOTYPE CAVEAT: sono strain-RATE comes out ~10x smaller in
  magnitude than the commanded rate (partly real muscle-tendon-grip series
  compliance, partly a still-imprecise steady-ramp window) -- do NOT fit a
  sono-based Vmax until the constant-velocity ramp segment is isolated; FL
  (a position, not a rate) has no such issue.
  `FL_FV_best_within_individual.png` (BUILT 2026-07-24,
  `R/summary_fv_fl_best_within_individual.R`) -- PI request: reproduce the
  canonical textbook FL/FV shapes (see the Hill reference) from the
  BEST-quality points, WITHIN A SINGLE INDIVIDUAL, in standard sign
  conventions. Three fixes vs. the messy zTorqueVsUhat diagnostic:
  (1) FORCE = the sign-FOLDED active-minus-(velocity-matched)passive
  single-axis tension (`muscle_force_Nm`), always POSITIVE -- the earlier
  "negative forces" were the UNfolded `force_zTorque_N`, a display artifact.
  (2) AXES in canonical convention: FL x = length strain with shorter muscle
  LEFT (negative) / longer RIGHT (positive) = -(shortening-positive sono
  strain); FV x = velocity with SHORTENING negative / LENGTHENING positive,
  normalised to V/Vmax. (3) THE UNLOCK for FV -- the constant-velocity window
  is detected from the COMMANDED angle (|omega| >= 60% of the step peak,
  longest contiguous run) and that SAME time window is applied to sono for
  BOTH the strain-rate (velocity) AND the active force (phase-matched). Using
  the mean force over the whole ramp had made the concentric limb look
  ANTI-Hill (force rising with shortening speed) because it mixed force-length
  effects; the CV-window force removes that and the concentric limb recovers
  the expected Hill shape. BEST-POINT gates (principled, then rationalised):
  RIGHT muscle only; genuine activation (F>0); FV only -- the sono ramp must
  be a clean straight line over the CV window (R^2 >= 0.90). Individual chosen
  by best combined FL+FV coverage = bass18 (confirmed by the all-individuals
  overview below: bass18 is the ONLY specimen whose isovelocity set spans BOTH
  a shortening and a lengthening side -- bass16/bass17 are eccentric-only, so
  they cannot form a full FV curve). FL additionally uses the SINGLE
  best-activated isometric trial within the individual (bass18 = bender_04);
  pooling both isometric trials had smeared the curve vertically because the
  second trial (bender_11) was only weakly activated. RESULT: Panel B (FV) is textbook
  -- force near zero at fast shortening (V ~ -0.6), ~0.45 at isometric (V=0),
  plateauing ~0.95 during lengthening (concentric-rise -> eccentric plateau).
  Panel A (FL) is an inverted-U peaking near -8% strain, declining on both
  limbs. CAVEAT/why per-panel normalisation: in bass18 the isometric holds
  produced anomalously little force (F0 ~= 1/17 of the eccentric peak, a
  weak-isometric artifact, NOT a physiological feature), so a shared F0
  normalisation would crush FL to a sliver; each panel is therefore
  normalised to its own within-individual peak to show the SHAPES faithfully.
  Per-step values (all 3 specimens, incl. whole-window force for comparison)
  in `data_processed/fv_fl_best_within_individual_steps.csv`.
  `FV_pooled_vs_ideal.png` (BUILT 2026-07-24,
  `R/summary_fv_pooled_vs_ideal.R`) -- PI request: pool the (surprisingly
  consistent) isovelocity FV points across ALL individuals on one graph,
  each specimen its own colour + loess line, and overlay an IDEALISED FV
  curve to show closeness of fit. Reads the gated per-step CSV from
  `summary_fv_fl_best_within_individual.R` (RIGHT muscle, F>0, clean sono
  ramp R^2>=0.90; force & velocity from the encoder-defined constant-velocity
  window applied to sono). Each individual normalised to its OWN Vmax and
  peak FV tension so lengthening plateaus -> 1 and shapes are comparable
  regardless of absolute force. Idealised curve = a Hill-type logistic
  (upper asymptote fixed at 1, steepness k and midpoint x0 fit to the pooled
  points by bounded optimisation). PI-directed exclusion: three bass18
  LENGTHENING points (bender_03 steps 14-16, residual +0.32..+0.39 above the
  ideal -- genuine eccentric overshoot / force-length contamination at those
  mid-velocity lengths) are dropped for the clean demonstration (identified
  as the top-3 positive residuals among bass18 lengthening points). RESULT:
  all three specimens trace the same concentric-rise -> eccentric-plateau
  shape; pooled R^2 vs the idealised curve = 0.83 (n=24, k=3.8, x0=0.50).
  NOTE the shortening (V<0) side is bass18-only
  -- bass16/bass17 isovelocity points that pass the gates are lengthening-only
  -- so the concentric limb is sparsely sampled; the midpoint sits right of
  isometric because bass18's near-isometric shortening force is elevated
  relative to a symmetric sigmoid.   Pooled normalised points saved to
  `data_processed/fv_pooled_normalised.csv`.
  `isometric_L0_activation_kinetics.png` (BUILT 2026-07-24,
  `R/summary_isometric_l0_activation.R`) -- PI request: normalised force-vs-time
  of ISOMETRIC NEAR-ZERO (L0) contractions pooled across every source they
  occur in, plus activation/relaxation-time boxplots vs. Coughlin & Carroll
  (2006, CBP-A 145:533, Fig. 2) red (slow) muscle. Sources of L0 (~0 deg)
  stimulated contractions: (1) isometric-protocol near-zero step(s), (2)
  isovelocity V=0 step, (3) dynamic pre/post-cycling L0 bookends (left+right,
  via `detect_dynamic_l0_bookends()`). Force is peak-normalised,
  baseline-subtracted, sign-folded contraction-positive. TOP panel: thin
  lines = individual contractions coloured by specimen, thick = per-specimen
  mean drawn only while >=50% of that specimen's contractions still
  contribute (the ~54 ms bookend twitches -- the majority, 120/159 -- end
  near 0.4 s, so this stops the mean cleanly instead of letting it jump onto
  the 0.3-0.5 s tetani-only tail). BOTTOM: boxplots (+ raw points) of
  activation time (rise to 90% of peak) and relaxation time (offset to 50%
  decay) per specimen, with a red shaded band = C&C red (slow) muscle
  mean+/-SD read off Fig. 2 (TA ~78 ms, TR ~150 ms) AND a gray shaded band =
  C&C white/fast feeding muscle (sternohyoideus + epaxial, TA ~10-20 ms,
  TR ~28-45 ms) -- so ours sit between the two references. RESULT: bass feeding-muscle L0
  kinetics are FASTER than red muscle -- activation medians ~45-80 ms, relaxation
  ~50-90 ms (mostly below the C&C red band) -- consistent with fast-twitch
  feeding muscle. CAVEAT (on the figure): ours are TETANIC rise/half-relaxation,
  C&C are single-twitch TA/TR, so the comparison is qualitative (fast-vs-slow).
  Bookends capped at 0.35 s post-stim AND before the next stim event of any
  kind, to avoid pre-cycling motor-motion contaminating the tiny twitch tail.
  Per-contraction times + durations in
  `data_processed/isometric_l0_activation_times.csv`.
  `isometric_L0_activation_kinetics_bookendsOnly.png` (BUILT 2026-07-24,
  same script) -- PI request: a bookend-ONLY twin of the figure above so the
  force-vs-time traces and the kinetics boxplots share ONE comparable
  contraction type. Restricted to the dynamic pre/post L0 bookend twitches
  (drops the isometric-protocol steps + isovelocity V=0 tetani). Same two
  reference bands (C&C red slow + white/fast). The traces are now all ~54 ms
  twitches ending near 0.4 s, and the boxplots are on a matched (compressed)
  y-scale instead of being squashed by the long tetani tails.
  `isometric_L0_activation_earlyVsLater.png` (BUILT 2026-07-24,
  `R/summary_activation_early_vs_later.R`) -- PI request: compare L0 activation
  (and relaxation) time between EARLY (preconditioning) and LATER (stable)
  session stages. Reads the per-contraction times above, classifies each by
  the specimen session cutoff (`dynamic_trial_precondition.R`:
  classify_session_precondition; bass16<5, bass17<9, bass18<5 = early).
  TYPE-CONTROLLED: restricted to dynamic pre/post L0 BOOKEND twitches only,
  because the isometric-protocol and isovelocity V=0 L0 contractions ran
  almost exclusively LATER (source x stage table printed by the script), so an
  all-source split would confound session stage with contraction type (fast
  bookend twitches vs slower 0.3-0.5 s tetani). Bookends are the one L0 source
  present in BOTH stages for all three specimens. Boxplots (+ points) per
  specimen x stage, Wilcoxon p per specimen (both groups n>=4), Coughlin red
  (slow) + gray white/fast bands for reference. RESULT: within the same twitch type, LATER twitches
  relax slower / more variably (bass17 p=0.006, bass18 p<0.001) and bass18
  also activates slightly slower (p=0.002), while early twitches are fast and
  tight -- a within-session slowing consistent with mild fatigue, NOT the
  early-fast artifact the confounded all-source version would have implied.
  Classified per-contraction table (all sources) in
  `data_processed/isometric_l0_activation_times_precondition.csv`.
  [DIAGNOSTIC tier] `isometric_meantension_vs_offset.png` /
  `isometric_maxtension_vs_offset.png` (`R/diag_precondition_tension_vs_
  offset_isometric.R`, `figs_diagnostic/`) -- specific tension (N/cm^2) vs.
  sono-strain offset, per isometric trial (right-muscle steps only; PI
  direction, 2026-07-24: use tension not power for isometric, since the
  motor doesn't move during activation so power is ~zero by design). Only
  5 isometric trials exist across all 3 specimens -- correlation is
  descriptive only (n=5, not significant). UPDATE 2026-07-24 (PI: "I agree
  about the CSA [geometric estimate]. Use the CSA in this table as an
  estimate for now" -- `01_inputs/bass_csa_measurements/
  bass_csa_measurements.xlsx`): tension now divides by
  `MEASURED_RED_MUSCLE_CSA_CM2` (0.55 cm^2, `muscle_geometry.R`) instead of
  the 3%-of-body-oval geometric guess. That constant is an image-analysis
  CSA measurement from a REFERENCE specimen ("bass07", NOT bass16/17/18),
  averaged over the middle 3 of 5 body "chunks" ASSUMED to bracket ~50%L
  (no chunk-length/position metadata exists to verify this). Effect was
  modest (~1.2-1.5x higher tension per trial) because the prior geometric
  CSA (~0.65-0.83 cm^2, specimen-specific) was already in a similar range
  to the new fixed 0.55 cm^2 -- our tension is STILL ~11-20x below
  Coughlin (2000)'s 18.64 N/cm^2, so CSA method was NOT the dominant driver
  of that gap; some other factor (activation level, series compliance,
  measurement method) likely also contributes.
  `coughlin2000_bass_power_work_tension_comparison.png` (BUILT 2026-07-24,
  `R/summary_coughlin2000_bass_comparison.R`) -- PI request: compare OUR
  measured mass-specific power (W/kg), work per cycle (J/kg) and specific
  tension (kN/m^2) against published largemouth bass RED MUSCLE values from
  Coughlin (2000, J. Exp. Biol. 203:617-629) -- a DIFFERENT reference paper
  than the "Coughlin & Carroll 2006" activation/relaxation band used above.
  User's 10-question clarification form was skipped, so the RECOMMENDED
  default from each question was used (documented in the script header).
  Panel A: the SAME dynamic L0 bookend twitches as
  `isometric_L0_activation_kinetics_bookendsOnly.png` (reused via the new
  `data_processed/isometric_l0_activation_traces.csv` export -- see that
  script), x-axis tightened to the twitch itself (-0.2 to 0.45 s, vs. the
  0-1.2 s span before), with shaded bands for stim-ON, the activation window
  (onset->90% peak) and the relaxation window (offset->50% decay) plus their
  median durations. Panels B/C/D: boxplots (+ points) of, respectively,
  isometric specific tension (trial-max, N/cm^2 -> kN/m^2), dynamic mean
  cycle power (W/kg) and dynamic mean work per cycle (J/kg) -- LATER
  (stable) trials only, dynamic power/work restricted to right-stim cycles
  (this repo's "active" convention) -- each with a red Coughlin (2000)
  reference band. PROVENANCE CAVEATS (also on the figure): Coughlin's
  tension (186.4+/-33.6 kN/m^2, N=18) is POSITION-POOLED (no significant
  position effect, so no single 50%L value exists); the work value
  (3.56+/-0.47 J/kg, N=6) is a SINGLE data point at 0.572L ("MID",
  closest to 50%L) for the FASTEST bass swimming speed tested (2.4 L/s,
  4.05 Hz) -- not a range, and embedded on-figure text in the PDF (all
  other position x speed combinations exist only as bar/scatter graphics
  and were not digitized); power (14.4+/-1.9 W/kg) is DERIVED as
  work x frequency (Coughlin's own definition), not independently reported.
  Coughlin's numbers reflect SUBMAXIMAL in-vivo swimming stimulation, not
  the animal's maximum -- an apples-to-oranges risk against our closer-to-
  maximal characterization protocols. RESULT: our isometric specific
  tension (later trials: ~0.15-0.81 N/cm^2 = 1.5-8.1 kN/m^2) is roughly
  15-50x BELOW Coughlin's band, most likely because our tension uses a
  GEOMETRIC whole-body-oval CSA estimate (`compute_muscle_mass_and_csa()`)
  rather than Coughlin's HISTOLOGICAL live-fiber-area measurement on an
  isolated single-myotome bundle. Our dynamic power/work land much closer
  to (and for some trials, above) the single Coughlin reference point.
  No bath/room temperature is logged anywhere in this rig's HDF5 metadata;
  assumed ambient ~20-22 C (uncontrolled) vs. Coughlin's controlled 20 C
  bath. Comparison data in
  `data_processed/coughlin2000_bass_comparison_data.csv`.
  UPDATE 2026-07-24: tension (panel B) now uses `MEASURED_RED_MUSCLE_CSA_CM2`
  (0.55 cm^2, `muscle_geometry.R`) instead of the geometric body-oval
  estimate -- see that entry's log addendum below. Tension improved only
  modestly (still ~11-20x below Coughlin, not fully explained by the CSA
  guess alone) -- do not describe this figure as "resolving" the tension gap.
  [DIAGNOSTIC tier] `FL_FV_candidates_all_individuals.png`
  (`R/summary_fv_fl_best_within_individual.R`, `figs_diagnostic/`) -- the
  all-individuals companion used to justify the
  bass18 choice: FL (top) and FV (bottom) faceted by specimen, points
  coloured by trial, per-specimen loess overlaid. Shows at a glance that (a)
  bass18's FL peak comes entirely from the well-activated bender_04 while
  bender_11 is a weak second trial, (b) bass17's FL is a clean descending
  limb only, (c) bass16's FL is trial-mixed/messy, and (d) only bass18 has
  isovelocity points on BOTH the shortening and lengthening sides.
  [DIAGNOSTIC tier] `FL_isometric_zTorqueVsUhat.png` /
  `FV_isovelocity_zTorqueVsUhat.png` (BUILT 2026-07-24,
  `R/diag_fv_fl_ztorque_vs_uhat.R`, saved to `figs_diagnostic/`) -- PI
  request 2B: FL/FV built with the two FORCE-reconstruction methods as
  separate panels. Both forces (N) come from `attach_vector_muscle_force()`:
  `force_zTorque_N` (simple single-axis) vs. `muscle_force_vector_N` (uHat
  empirical vector); commanded x-axis held constant so only the force method
  varies. FINDINGS: (1) the two columns carry OPPOSITE sign conventions
  (compare magnitude/trend, not sign); (2) uHat COLLAPSES toward zero for
  bass17 (low-SNR empirical direction) and is generally noisier/compressed,
  while the single-axis zTorque robustly recovers clean monotonic FL and
  structured FV trends for all 3 specimens -- concrete support for the
  earlier conclusion that the simple single-axis (zTorque/primary-bending-
  axis) force is the more robust default for these data. RIGHT muscle only.
  [DIAGNOSTIC tier] `tension_zTorqueVsUhat_scatter.png` (BUILT 2026-07-24,
  same script) -- PI follow-up: "how does zTorque muscle TENSION compare
  to uHat tension? double check the math." Converts the same 114 right-
  muscle steps to specific tension (`|force_N| / MEASURED_RED_MUSCLE_
  CSA_CM2`, 0.55 cm^2, `muscle_geometry.R` -- one shared CSA for both
  methods, since CSA is a muscle property, not a force-method property);
  compared on MAGNITUDE (`|force|`) since the two raw force columns use
  opposite sign conventions by construction (confirmed: 85% sign-disagree
  for isometric steps, 62% for isovelocity). Faceted isometric/isovelocity,
  1:1 dashed line, Coughlin (2000) 18.64+/-3.36 N/cm^2 band on the x-axis.
  RESULT: agreement is weak AND specimen-dependent, not a stable
  relationship -- magnitude correlation ranges from r=0.89 (bass16
  isometric) to r=-0.58 (bass18 isometric, inverted) to r=-0.08 (bass18
  isovelocity, none); the uHat/zTorque ratio swings from ~0.10 (bass17,
  matching its known uHat SNR collapse) to ~5.4 (bass18 isovelocity) --
  too unstable to treat as a fixed calibration factor between methods.
  Pooled, zTorque tension is systematically ~2x (mean) to ~3x (max) larger
  than uHat's. CAVEAT (on the figure and in the log): the isovelocity
  tension values here (up to several x the Coughlin band) are NOT a
  trustworthy Coughlin comparison -- this script's isovelocity force is a
  pre-fix, local PEAK-WINDOW MEAN (`MFV_PEAK_WINDOW_S=0.15s`,
  `muscle_force_vector.R`), the same windowing `summary_fv_fl_best_
  within_individual.R` later found produced an "anti-Hill" inflation
  artifact before its constant-velocity-window fix; that fix was not
  retrofitted here (out of scope -- this is a force-method comparison, not
  a length/velocity-model one). The ISOMETRIC-only panel (no motion
  confound) is the trustworthy read, and it shows both methods sitting
  far below Coughlin regardless of force method -- the CSA-driven tension
  gap documented elsewhere is not resolved by switching force methods.
  Per-step data: `data_processed/fv_fl_ztorque_vs_uhat_tension.csv`.

## Filename tokens
`{signal}_{protocol}[_{method}][_{filter}].png` inside per-fish folders.
`figs_summary/` has no fixed prefix convention -- it only holds genuinely
cross-fish content (see "Summary" above), so names describe the comparison
itself (`FLsuperplot_...`, `specimen_comparison_...`), not a
`{bassID}_`-prefixed single-fish plot. Per-fish tokens always appear in the
order above; omit a bracketed token when it doesn't apply (e.g. `legacy` method is implicit and
dropped in a few original-workflow names — see table below for the exact
names in use).

- **protocol**: `isometric` | `isovelocity` | `dynamic` | `freqsweep`
- **signal**: what's plotted — `FL`, `FV`, `FVl0` (FV sampled at sono-confirmed
  L0 crossing), `forceTime`, `uhatCompare`, `fatigueCheck`,
  `powerDynamic`, `stiffnessDamping`, `FLsuperplot`.
  `FLsuperplot_*` REVISED 2026-07-21 (PI-directed, "the rule is that FL
  superplot only contains moments or steps where V = 0"): the isovelocity
  continuous-ramp sweep (angle-matched pointwise passive subtraction) was
  REMOVED outright -- pooling it folded real force-VELOCITY behavior into a
  force-LENGTH plot regardless of how carefully the passive baseline was
  subtracted. Three V=0-only sources are pooled instead: isometric holds,
  isovelocity's OWN embedded V=0 (`isometric_zero`) holds (previously used
  only as an internal F0 reference, never plotted), and the pre-/post-
  cycling L0 stim bookends every dynamic trial brackets its sinusoidal sweep
  with (motor stationary at 0 deg -- `R/extract_dynamic_l0_bookends.R`).
  Fixing this ALSO resolved the earlier "raw peak DOWN, normalized peak UP"
  confusion at strain=0% (see `analysis_muscle_force_vector_log.md`) -- the
  trough/spike artifact was intrinsic to the sweep's mid-ramp zero-angle
  crossings, not a sign-convention bug, and disappeared once the sweep was
  dropped.
  `FLsuperplot_..._normalized[_snrPass].png` (added 2026-07-21, Gate A
  exploration -- see `analysis_muscle_force_vector_log.md`): same pooling as
  the corresponding non-`_normalized` file, each point expressed as F/F0
  (that trial+side's own SNR-passing L0/V0 force, or for a dynamic trial its
  own pre+post bookend mean) instead of raw absolute N, to correct
  cross-trial force-scale differences before pooling. Kept SIDE BY SIDE with
  the raw file, NOT a replacement. bass17 now contributes points to the
  normalized version too (previously zero) because dynamic bookends give it
  an SNR-trustworthy L0 reference its isometric/isovelocity trials lacked.
  `FVsuperplot_isovelocity_pooled[_normalized].png` (BUILT 2026-07-22,
  `R/superplot_fv_pooled.R`, PI-requested "similar to FL superplots, raw and
  normalized") -- pooled Force-Velocity superplot, same 3-tier connect-the-
  dots convention as `FLsuperplot` (thin=trial, medium=individual mean,
  thick black=group mean), x-axis swapped to `shortening_strain_pct` treated
  as a STRAIN RATE (%/s, muscle-centric: + = concentric/shortening, - =
  eccentric/lengthening -- NOT raw commanded `operating_point` sign, which
  can differ by side/motor-direction at the same contraction_mode).
  V=0 anchor pools isometric L0 reps + isovelocity's own V=0 holds + dynamic
  L0 bookends (same 3 sources as `FLsuperplot`'s V=0-only rule, just serving
  as the anchor here instead of the whole dataset); the velocity-dependent
  part is isovelocity's actual moving (concentric/eccentric) ramps, via the
  REAL angle-matched batch (`compute_isovelocity_vector_batch()`,
  angle_matched(/_cross_trial/_nearest_v) sources ONLY -- see
  `ytorquesignexamples`'s 2026-07-22 correction for why static-baseline
  isovelocity rows are excluded. `_nearest_v` (added 2026-07-22, see
  `isovpassivemodels`) covers steps whose EXACT commanded velocity has no
  stim-off ramp anywhere in the corpus -- confirmed to remove 100% of a
  prior negative-overshoot artifact in bass18, see analysis log "flag 1
  ROOT-CAUSED + FIXED"). F0 (normalized companion) falls back from trial+side to
  fish+side when a trial's own V=0 doesn't clear SNR_MIN (isovelocity's
  embedded V=0 holds are usually too low-SNR to self-normalize -- observed
  SNR ~0.1-1.9 vs. isometric's dedicated holds). LIMITATION (documented in
  the figure's own subtitle): force is MEAN-over-active-window, sampled over
  whatever strain range the ramp swept through -- NOT sampled at a fixed L0
  crossing the way a textbook FV curve is defined; a stricter alternative
  (`.mfv_fv_l0_crossing()` / `fv_l0`, sono-confirmed, right-side-only)
  already exists in the codebase but would drop most of the data for a
  first pass, so it is a candidate refinement, not silently substituted in.
  `angleValid`/`strainValidCmd`/`strainValidSonoEnc`/`strainValidSonoCmd`
  MOVED to the diagnostic tier 2026-07-21 (PI-directed) -- see "Diagnostic
  vs. individual vs. summary" above; no longer in `figs_{bassID}/`.
  `strainValidMeas` REMOVED 2026-07-21 (PI-directed) as redundant with
  `strainValidCmd`'s isometric panel. `powerDynamicMassSpec` REMOVED
  2026-07-21 (PI-directed) -- merged INTO `powerDynamic` as a second
  (patchwork-stacked) panel in the same file rather than a separate one.
- **method**: which force/baseline calculation —
  - `legacy` = original single-axis (zTorque only) calculation, no û.
    MOSTLY REMOVED 2026-07-21 (PI-directed: "legacy figures muddy the
    waters") once a vector (`uhatBoth`) or `baselineInterp` alternative
    existed for that signal/protocol combo. TWO exceptions remain
    (flagged for PI confirmation, not yet resolved, because no
    alternative exists for them): `fatigueCheck_isovelocity_legacy.png`
    and `forceTime_dynamic_legacy.png`.
  - `baselineInterp` = isometric per-step baseline computed via
    cross-step interpolation instead of that step's own pre-stim window
    -- NOT affected by the `legacy` removal above (this is the separate,
    still-open `passivebaselinemethod` diagnostic decision, not the
    force-calculation-method decision)
  - `uhatEmp` / `uhatGeom` = 6-axis vector force via empirical /
    geometric line-of-action
  - `uhatBoth` = both û methods faceted side by side in one file
- **filter**: `snrPass` = regenerated with low-confidence points/steps
  dropped entirely. REVISED 2026-07-22 (SNR-magnitude conflation audit,
  `analysis_muscle_force_vector_log.md`) -- the drop criterion is now the
  4-tier confidence factor (`mfv_confidence_tier()`, `muscle_force_vector.R`),
  NOT a bare `activation_snr < 3` ratio: a row is kept if its confidence tier
  is `confident` OR `confidently_small` (i.e. its own force clears its own
  `baseline_force_noise_N` even when the SNR ratio doesn't), and dropped only
  if `unstable_magnitude` or `unconfirmable`. The one exception is
  `uhatCompare_empVsGeom_snrPass` (a u_hat DIRECTION/angle comparison, which
  has no force magnitude to test) -- that one file is still ratio-only by
  design. Omitted = unfiltered (all 4 tiers are still alpha-flagged in-plot
  where applicable, never dropped).

## Old name -> new name
Historical migration record (2026-07-19 renaming pass) -- the `trial_plots/`/
`summary_plots/` subfolder references below describe the layout AS IT WAS
THEN, not the current one (flattened 2026-07-21, see "Folder layout" above).
The filename tokens themselves are still current.

### Trial plots (was in `trial_plots/`, filename root = `{bassID}_bender_{NN}_{protocol}`, unchanged)
| Old | New |
|---|---|
| `{tid}.png` | unchanged |
| `{tid}_interpbaseline.png` | `{tid}_baselineInterp.png` |

### Summary plots (`summary_plots/`, per-fish, no bassID prefix needed — folder scopes it)
| Old | New |
|---|---|
| `summary_isometric_FL.png` | `FL_isometric_legacy.png` |
| `summary_isometric_FL_interpbaseline.png` | `FL_isometric_baselineInterp.png` |
| `summary_isometric_FL_vector.png` | `FL_isometric_uhatBoth.png` |
| `summary_isovelocity_FV.png` | `FV_isovelocity_legacy.png` |
| `summary_isovelocity_FV_vector.png` | `FV_isovelocity_uhatBoth.png` |
| `summary_isovelocity_FV_vector_L0_sono.png` | `FVl0_isovelocity_uhatBoth.png` |
| `uhat_empirical_vs_geometric.png` | `uhatCompare_empVsGeom.png` |
| `force_vs_time_isometric.png` | `forceTime_isometric_legacy.png` |
| `force_vs_time_isovelocity.png` | `forceTime_isovelocity_legacy.png` |
| `force_vs_time_isometric_vector.png` | `forceTime_isometric_uhatBoth.png` |
| `force_vs_time_isovelocity_vector.png` | `forceTime_isovelocity_uhatBoth.png` |
| `force_vs_time_dynamic_by_{bname}.png` | `forceTime_dynamic_by_{bname}.png` (bname = frequency/amplitude/phase/duty) |
| `fatigue_check_isometric.png` | `fatigueCheck_isometric_legacy.png` |
| `fatigue_check_isovelocity.png` | `fatigueCheck_isovelocity_legacy.png` |
| `fatigue_check_isometric_interpbaseline.png` | `fatigueCheck_isometric_baselineInterp.png` |
| `strain_validation_commanded.png` | `strainValidCmd.png` |
| `strain_validation_measured.png` | `strainValidMeas.png` |
| `strain_validation_sono_vs_encoder.png` | `strainValidSonoEnc.png` |
| `strain_validation_sono_vs_commanded.png` | `strainValidSonoCmd.png` |
| `angle_validation.png` | `angleValid.png` |
| `summary_dynamic_power.png` | `powerDynamic.png` |
| `summary_dynamic_power_massspecific.png` | `powerDynamicMassSpec.png` |
| `summary_frequency_sweep_stiffness_damping.png` | `stiffnessDamping_freqsweep.png` |

### Cross-individual superplots (`figs_diagnostic/`, flat)
| Old | New |
|---|---|
| `superplot_FL_pooled_bass16_17_18.png` | `FLsuperplot_isometric_isovelocity_pooled.png` |
| `superplot_FL_pooled_bass16_17_18_snr_passing.png` | `FLsuperplot_isometric_isovelocity_pooled_snrPass.png` |

### `figs_summary/` (flat, curated, bassID-prefixed; only SNR-passing material)
Trial plots copied as-is (already `bassID_bender_NN_protocol[_baselineInterp].png`).
Legacy summary copies (no filtering possible) and regenerated SNR-passing
summaries both follow `{bassID}_{new-summary-name-above}[_snrPass].png`, e.g.
`bass18_FL_isometric_uhatBoth_snrPass.png`,
`bass17_uhatCompare_empVsGeom_snrPass.png`.

## Adding a new figure
Pick tokens from the lists above in order; add a new token to this file in
the same change if you introduce a new `signal`/`method`/`filter` value.
Don't invent a new suffix style — extend this scheme.
