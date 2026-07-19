# Muscle force vector analysis — decision log

Purpose: preserve the *why* behind the 6-axis vector-force work (bass16/17/18)
that isn't obvious from the code alone. `muscle_force_vector_physics.md` is
the model/physics reference (coordinate frame, moment arms, formulas,
troubleshooting data). This file is the decision trail: what was asked, what
was tried, what changed and why, and what's still open. Update it (don't
replace it) the next time this analysis line moves.

## Why this started
QC plots showed a stim-locked DC shift in **yTorque** correlated with muscle
shortening (sono length decreasing), even though FL/FV force had always been
computed from **zTorque alone** with an assumed muscle depth. The muscle runs
longitudinally; a longitudinal (X-axis) force projects onto yTorque via the
dorsoventral offset, not onto zTorque. That meant real muscle force was very
likely leaking out of the force calculation, while work-loop power (computed
from zTorque/bending only, not X-axis force) still matched published values —
i.e. the power estimate could be right while the FL/FV force estimate was
structurally blind to the muscle's actual line of action.

## Foundational decisions (from the Q1-Q6 exchange)
- **Sensor frame:** X = anterior-positive, Y = mediolateral (sign flips with
  facing direction — see `daq_positive_motor_direction`), Z = ventral-positive.
- **Target quantity:** scalar force in Newtons along the muscle's line of
  action û, i.e. projected-moment / moment-arm — not just a shape-only
  moment proxy.
- **Moment arms:** z-torque arm = half-body-width minus muscle depth (muscle
  sits medial to the body surface); y-torque arm = the fixed dorsoventral
  clamp offset from metadata. Both live in `muscle_geometry.R`.
- **û(θ) strategy — both, not either/or:** empirical (active-minus-passive
  wrench direction at each isometric step) as primary, geometric (beam
  tangent from commanded curvature) as an independent cross-check. Agreement
  validates û without leaning on the "active-minus-passive is uncontaminated"
  assumption; divergence is itself diagnostic.
- **Step selection for building û:** rank by activation magnitude, human
  reviews/vetoes the top candidates rather than a silent auto-pick.

## What changed as the build progressed (and why)
1. **Empirical û was side-dependent-contaminated** (diverged from the
   geometric/beam-tangent prediction because of a lateral reaction force not
   related to the muscle). Resolved by driving the actual force *projection*
   off geometric û (physically motivated, side-independent) while still
   computing empirical û for comparison/diagnostic plots — both are kept as
   parallel outputs (`_uhat_emp` vs `_uhat_geom` / `_vector` facets), never
   silently dropped.
2. **Sign convention, twice.** First pass: standardized left/right muscle
   sign as "tension-positive," which flipped bass17's FL curve concave-up.
   That was wrong: concentric contraction is not "more extreme" than
   eccentric in an absolute sense — it's only more extreme **in the
   direction of the passive baseline torque**. Fixed by re-deriving the sign
   from whether the active wrench reinforces or opposes the passive-baseline
   torque direction (`tension_relative_to_passive` in
   `.mfv_finalize_step`, `muscle_force_vector.R`), not from an assumed
   absolute tension convention. This restored the expected concave-down
   hill-type shape.
3. **Isovelocity passive-baseline matching widened deliberately**, in this
   order: same-trial angle-matched passive ramp -> cross-trial angle-matched
   passive **within the same individual only** (never across fish) -> flag
   if no angle-overlapping passive exists anywhere. "Same individual only"
   is a hard rule, not a fallback of convenience — different fish have
   different passive stiffness.
4. **Sono added as an independent confirmation, not a primary signal.**
   40 Hz low-pass (`mfv_load_sono_lp40`), used only to confirm the FV L0
   crossing sampled from the force trace — sono is right-muscle-only, so it
   can corroborate but never substitute for the force-based measurement.
5. **bass17 isometric step 15 looked like a beautiful hill-type FL curve,
   but the ascending/descending legs exceeded the L0 magnitude** — a
   red flag. No-stim drift check (early vs. late pre-stim baseline) showed
   this trial is dominated by passive creep, not true active force. By
   contrast bass18's steps hold clean, force-like shapes with the *same*
   baseline logic applied — so the problem is trial-specific signal quality,
   not a flaw in the baseline method itself. This is why **SNR-based
   filtering**, not per-trial detrending, was chosen as the fix: detrending
   would risk quietly manufacturing a shape.
6. **`activation_snr` became the universal confidence gate** — used to (a)
   decide which steps are trustworthy enough to build empirical û from, (b)
   alpha-flag "low confidence" points in every vector summary plot rather
   than silently drop them, and (c) gate which trials/steps get copied or
   regenerated into `bass_summary_figures/` (majority-of-steps-pass rule for
   trial plots; point/step-level drop for regenerated summary plots).

## Validated / cross-checked
- Empirical vs. geometric û agreement checked directly
  (`uhat_empirical_vs_geometric.png` / `uhatCompare_empVsGeom.png`).
- One clearly-active bass18 step run through the moment arm and compared to
  specific-tension x CSA (Coughlin range) as a sanity check — reported as a
  check, not tuned to hit the range.
- SNR>=3 threshold cross-checked against physiologically expected muscle
  tension (see `muscle_force_vector_physics.md` troubleshooting log):
  real activation should clear it by 50-350x; failing steps are genuinely
  weak/contaminated, not victims of an overly strict threshold.
- Ruled out (not the cause of per-fish noise variability): calibration file
  choice (FT56491 vs FT56492), ADC gain/quantization, stim_monitor
  crosstalk. Session-specific mounting/setup remains the leading unproven
  hypothesis for bass17's ~2x higher noise floor.

## Where things live (code map)
- `muscle_force_vector.R` — core: baseline subtraction, û construction
  (empirical + geometric), wrench->force solve, sign standardization,
  sono LP/L0, per-step SNR.
- `plot_muscle_force_vector.R` — all vector-force plots, confidence-alpha
  flagging, emp-vs-geom faceting.
- `superplot_fl_pooled.R` / `superplot_fl_pooled_snr_passing.R` — pooled
  cross-individual FL superplot (isometric + isovelocity instantaneous
  sweep), unfiltered and SNR-filtered.
- `export_snr_summary_figures.R` — per-fish SNR-gated export pass into
  `bass_summary_figures/`.
- `muscle_force_vector_physics.md` — physics model + troubleshooting data
  (SNR/calibration/gain/noise investigation).
- This file — decision trail only; no formulas, no data.

## Open / not yet built
**Dynamic-trial pre/post isometric-hold extraction (planned, not yet
built).** Dynamic (`single_finite`, sweep) trials contain flat/motionless
windows bracketing the sweep (confirmed directly: e.g.
`bass16_bender_01_dynamic.h5` is flat from t=-4s to -0.17s, sweeps to 5.5s,
flat again to 9.33s) and `stim_side` can be active in those flat windows.
Plan: (1) detect flat pre/post windows per dynamic file (angle derivative
~0 for a sustained span) and any stim event inside them; (2) treat each such
window as an isometric-equivalent step (reuse the existing per-step
baseline+û+SNR pipeline, `.mfv_finalize_step`, unchanged) at whatever
angle/length it sits at; (3) feed the resulting points into the pooled FL
(and FV, if a hold has a brief ramp) plots as additional data,
`activation_snr`-gated like everything else. Explicitly **not** in scope:
extracting a continuous instantaneous force-length trajectory from the
moving/sweeping middle of a dynamic trial — only the flat bracketing holds
are usable as isometric-equivalent points.

**Persistent good/bad step registry.** `activation_snr` already exists per
step; what's missing is a durable per-fish CSV (step_id, protocol, SNR,
pass/fail) so "give me the good steps" doesn't require re-deriving it from
a plot each time. Deferred by request — revisit once the dynamic pre/post
feature above is settled, since it will add new step rows to register too.

**Figure naming.** See `FIGURES_README.md` (new) for the naming convention
adopted 2026-07-19 and what replaced the old `_interpbaseline`/`_vector`/
`_L0_sono` suffixes.
