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
  deliberately NOT fabricated for symmetry.), `torquesmoothingmethod`,
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
  hold + moving ramps (pointwise angle-matched passive), all points shown
  (SNR-failing alpha-flagged, never dropped), summary line = mean of
  SNR-passing points per velocity. FINDING: bass18 right is a "W", not a Hill
  hyperbola -- eccentric plateau (~1.65-1.76 N) drops sharply to 0.56 N at V=0,
  then oscillates (1.37/1.03/1.83 N at 127/255/382). The V=0 notch traces to a
  CROSS-TRIAL FATIGUE CONFOUND: its 3 SNR-passing contributors are all from one
  (weaker) trial-set while the stronger trial's V=0 values (1.69/2.41 N, fully
  consistent with the eccentric plateau) fail the SNR gate (brief/embedded
  hold). bass16/17 show a PEAK AT V=0 with decline on both sides (a tent) --
  neither Hill's plateau-then-decline nor flat noise; unexplained. Honest
  state: the eccentric>concentric pairwise ordering remains directionally
  suggestive for bass18 but is NOT a settled demonstration of the Hill
  relationship -- a proper test needs each trial's OWN V=0 vs its OWN moving
  steps (within-trial only), not pooled across trials -- not yet built.
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
  FLsuperplot/FVsuperplot. Geometric u_hat throughout.
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
- **filter**: `snrPass` = regenerated with sub-threshold
  (`activation_snr < 3`) points/steps dropped entirely. Omitted = unfiltered
  (low-confidence points are still alpha-flagged in-plot where applicable,
  not dropped).

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
