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
   regenerated into `figs_summary/` (majority-of-steps-pass rule for
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

## Gate A exploration: trial-to-trial force-scale normalization (2026-07-21)
**Motivation.** `FLsuperplot_isometric_isovelocity_pooled.png` pools RAW
absolute force across trials (e.g. bass16 has 2 isometric + 2 isovelocity
trials at IDENTICAL nominal operating points, recorded at different session
times) -- PI feedback: "currently very unhelpful." Hypothesis: cross-trial
force-SCALE differences (fatigue, clamp re-seating, etc. between repeat
trials) are indistinguishable from real FL shape once pooled raw, and are a
major contributor.

**Two things tried, two different outcomes:**
1. *Fatigue-timeline (time-based) normalization* -- force / that trial's own
   MAX force, tested on `plot_fatigue_timeline.R`. Mixed/net-negative
   result (fixed bass16's mid-session rebound, but hurt bass17/bass18's
   already-clean signals) -- ABANDONED for that plot, reverted. Wrong
   reference point (arbitrary max) for a TIME diagnostic.
2. *FL/FV pooling (length-based) normalization* -- force / that
   TRIAL+SIDE's own L0 (isometric) or V=0 (isovelocity) force, i.e. F/F0,
   the classic muscle-physiology normalization, built into
   `superplot_fl_pooled.R`/`superplot_fl_pooled_snr_passing.R`. First
   attempt (F0 computed from ALL L0/V0 reps regardless of confidence) blew
   up into +-1000 spikes -- dividing by a noisy near-zero F0 amplifies
   noise. FIX: gate F0 itself by the SAME `activation_snr >= MFV_UHAT_SNR_MIN`
   threshold used everywhere else (not a new invented number) -- only use an
   L0/V0 rep as F0 if IT is independently trustworthy. Result: clean,
   symmetric, physiologically-plausible FL peak shape, bass16/bass18
   tracking consistently.

**Cost, and the decision:** SNR-gating F0 means a trial/fish with NO
trustworthy L0/V0 rep contributes NO normalized points at all -- bass17
drops out ENTIRELY (traced: its isometric L0 activation_snr is 1.18-2.36
across all 4 reps, its isovelocity V=0 holds also never clear 3.0), which
tracks with the already-documented "~2x higher noise floor" for bass17
above. PI decision (2026-07-21): ACCEPT this exclusion (a curve built from
an untrustworthy reference point shouldn't count), and KEEP raw + normalized
as two side-by-side files for now (`FLsuperplot_..._pooled.png` /
`FLsuperplot_..._pooled_normalized.png`, same for the `_snrPass` pair) --
normalized is NOT yet promoted to replace raw as canonical.
(UPDATE 2026-07-21, later same day: dynamic L0 bookends were added as a
THIRD F0 source alongside isometric/isovelocity -- bass17 now clears the
SNR gate via its bookends even though its isometric/isovelocity L0/V0 reps
still don't, so this specific exclusion no longer applies. See the
"RESOLVED" item further down for the fuller V=0-only rework.)

**New Gate A framing this opened up:** F/F0 normalization, done safely,
acts as an IMPLICIT quality gate distinct from a fatigue-based one -- a
trial is only usable if ITS OWN reference point is independently
trustworthy, not just whether individual points are SNR-passing. Whether
this becomes (part of) Gate A's actual architecture, vs. staying a parallel
diagnostic lens, is still open.

**Raw-peak-DOWN vs. normalized-peak-UP at strain=0% (2026-07-21, PI-flagged
"very confusing," STILL OPEN).** Root-caused as far as code hygiene goes,
NOT yet resolved as a modeling question:
1. *Drift fixed, no visual change (as predicted):* `superplot_fl_pooled.R`'s
   `extract_isovelocity_sweep()` had its OWN hand-copied implementation of
   the tension-sign decision, never updated with the `is_l0_step` guard
   added to `.mfv_finalize_step()` (muscle_force_vector.R) on 2026-07-21.
   Unified both into one shared helper, `.mfv_tension_sign()`
   (muscle_force_vector.R). Regenerated both superplots after the fix --
   the raw plot's trough at strain=0% (group mean ~-0.25, bass18 spike to
   ~-1.9) is UNCHANGED. Expected: sweep steps already skip
   `operating_point == 0` holds upstream, so `is_l0_step` never fires for
   them regardless -- this refactor closes a real drift risk but was never
   going to touch this specific artifact.
2. *Actual mechanism (diagnosed, not yet fixed):* an isovelocity sweep step
   picks ONE fixed `tsign` for its ENTIRE ramp, from that step's own
   AGGREGATE passive-torque projection (usually unambiguous, since a real
   sweep's mean is dominated by its large-angle ends). That same fixed sign
   then gets applied to every sample along the ramp, INCLUDING near the
   geometric zero-angle crossing, where the projection lever `rxu_e` itself
   can swing through zero/sign-change purely from geometry -- independent
   of whether the muscle's real tension is changing. Different steps/
   trials/fish land on different dominant signs, so pooled+binned points
   right at the crossing (where the true signal is smallest, hence most
   sensitive to whichever fixed sign won) scatter/trough. Meanwhile F0
   (isometric_zero HOLD, routed through `.mfv_finalize_step()`'s
   `is_l0_step`-forced-ambiguous rule) is mathematically forced
   non-negative (`tension_sign * F_moment_N == abs(F_moment_N)` there) --
   so F0 is a stable, always-positive denominator, and isometric's OWN
   strain=0% points normalize to ~1.0 near-TRIVIALLY (F/F0 where F IS
   (close to) F0, same window/mean/sign method) -- this is very likely the
   real source of the normalized plot's "peak UP at exactly 0," not a
   biological signal.
3. *RESOLVED 2026-07-21, by removing the mechanism instead of patching it:*
   PI reframed the question entirely -- "the rule is that FL superplot only
   contains moments or steps where V = 0" -- rather than deciding how to
   patch the sweep's per-step-fixed `tsign` near its zero-crossing (option 3
   above), the isovelocity sweep extraction itself
   (`extract_isovelocity_sweep()`) was DELETED from `superplot_fl_pooled.R`.
   It is superseded by three V=0-only sources: isometric holds (unchanged),
   isovelocity's own embedded `isometric_zero` holds
   (`extract_isovelocity_zero_points()` -- previously used only as F0,
   never plotted as data), and dynamic trials' pre-/post-cycling L0 stim
   bookends (`extract_dynamic_l0_points()`, using
   `R/extract_dynamic_l0_bookends.R` + the same
   `attach_vector_muscle_force(..., td6_override=)` pattern already proven
   by the fatigue-timeline's dynamic integration). Since none of these three
   sources ever has a mid-ramp geometric zero-crossing (isometric holds a
   fixed nonzero-or-zero angle; V=0 holds/bookends are AT 0 deg, not
   crossing through it), the artifact's root cause (item 2 above) cannot
   occur anymore -- confirmed by re-running both superplots: the raw plot's
   strain=0% trough/spike is GONE (group mean now a clean local peak, no
   sign-driven scatter), and the normalized plot's near-1.0 peak at 0% is
   still present but no longer confusable with the trough -- it's the
   expected, honest "F/F0 at the reference point is ~1 by construction"
   behavior, now sitting alone without a mismatched raw-plot trough
   contradicting it. Bonus: dynamic bookends give bass17 an SNR-passing L0
   reference it never had from isometric/isovelocity alone, so bass17 now
   contributes to the normalized plot too (previously zero points, see the
   `figs_summary/` note above).

**L0 "spike" fully RESOLVED 2026-07-22, by DELETING the passive-relative
`tension_sign`/`.mfv_tension_sign()` mechanism entirely** (not patching its
`is_l0_step` special case further -- see item 3 above, which left a residual
near-1.0 spike/discontinuity at strain=0% because `is_l0_step` forced a
DIFFERENT sign rule at L0 than everywhere else). PI: "resolve now." Empirical
check across all 328 real finalized steps (3 fish x isometric + isovelocity-
V0 + dynamic bookends x muscle_force_vector_N/muscle_force_vector_geom_N/
force_force_channel_N/force_yTorque_N/force_zTorque_N, captured by
instrumenting the REAL `.mfv_finalize_step()` calls, not reimplemented) found:
- `force_zTorque_N` cleanly MIRRORS by anatomical side (left one sign, right
  the other, consistent across all 3 fish) -- for the SAME reason the legacy
  zTorque-only path already corrects for via `force_sign = rec_lidx *
  lidx_pos_motor` (`resolve_step_contraction()`, `muscle_geometry.R`).
  Unrelated to L0 or concentric/eccentric.
- `muscle_force_vector_N`/`muscle_force_vector_geom_N`/`force_force_channel_N`
  are mostly positive (91.5%/94.5%) with a REAL negative minority concentrated
  in isovelocity (16-21%) -- not a side effect, not an L0 artifact: a
  passive-relative correction would have been ERASING genuine eccentric-vs-
  concentric signal to force these positive.
- `force_yTorque_N` splits almost entirely on STATIC (isometric hold, or
  isovelocity's own V=0 embedded holds -- 96.4%/89.7% positive) vs. ACTIVELY
  MOVING (isovelocity's concentric/eccentric ramp, genuinely rotating during
  the stim window -- 73.3%/83.3% negative) -- see `ytorquesignexamples.png`
  (`figs_diagnostic/`) for 3 real positive + 3 real negative examples (all
  SNR >= 6), and `axissnrcomparison.png` for SNR compared across all 4 axes
  on the SAME noise floor (medians 3.8-5.4, roughly comparable, zTorque
  slightly highest).
**New rule (replaces the deleted mechanism):** two DIFFERENT treatments for
two DIFFERENT behaviors, not one `tension_sign` forced onto all five outputs.
`force_zTorque_N` gets the FIXED per-side `force_sign` correction (computed
from file geometry, works identically for dynamic bookends' synthetic step
rows). The other four are reported in their RAW, as-computed sign -- no
per-step decision, no L0 special case, because there is no more "ambiguous
vs. reinforcing/opposing" question left to answer for them. Regenerated
`FLsuperplot_isometric_isovelocity_pooled(_normalized).png`: the raw plot's
values are now all-positive (previously alternated sign, described as
"confusing" by PI) and the normalized plot's L0 region is now a smooth
trough (classic FL-curve shape), not a spike/discontinuity.

**Isovelocity's y-torque follow-up, same day: CONFIRMED as (partly)
uncorrected inertial coupling, not purely a stimulus-driven signal.** PI
observation: "y torque is positive for isometric stims and the rise is
coincident with muscle stimulation. On the other hand, negative y torques
in isovelocity happen well after the stimulus and most likely in relation
to the velocity ramp. Velocity ramp suggests inertial noise?"

First, a CORRECTION to the numbers above: `ytorquesignexamples.png`'s
isovelocity examples (and the "73.3%/83.3% negative" figure) were captured
via `extract_isovelocity_zero_points()`'s internal per-step loop, which
computes EVERY step (including moving ones) with a STATIC pre-stim baseline
before discarding everything except V=0 -- comparing a fast-moving active
window against a motionless pre-stim window is guaranteed to show a huge
motion-linked deviation, unrelated to the REAL isovelocity path
(`compute_isovelocity_vector_batch()`, angle-matched passive, what actually
feeds `FV_isovelocity_uhatBoth.png`). Re-run with the CORRECT angle-matched
computation (313 real finalized steps, angle-matched-only for isovelocity's
moving rows): isovelocity's moving steps are 54% negative (was 73-83%,
wrong) vs. static holds' 97% positive (was 96-97%, unaffected by the bug --
static holds were always computed correctly). 54% is much closer to a coin
flip than a real directional signal.

Second, the MECHANISM test (`ytorqueinertialtiming.png`,
`R/diag_ytorque_inertial_timing.R`): a representative active step
(bass16_bender_16_isovelocity, step 12, left, eccentric, -335.6 deg/s) and
its matched WITHIN-TRIAL no-stim ramp of the same commanded speed (step 20)
show the SAME Ty dip at the SAME time (2.129s vs. 2.128s) relative to the
SAME angular-velocity-profile feature (both peak at 2.057s) -- one of these
traces has ZERO stimulation. This is a KINEMATIC feature of the ramp's own
motion, not a stimulus-locked muscle event. Aggregate check across all 123
isovelocity ramps (moving steps only, both stimulated and no-stim), 3 fish:
median |t(Ty extremum) - t(|angular velocity| extremum)| = 0.071s (stim) vs.
0.072s (no-stim) -- essentially IDENTICAL lag with or without stimulation
(`ytorqueinertialtiming_stats.png`).

Conclusion: `force_yTorque_N`'s isovelocity sign is NOT purely a real
muscle-force signal being correctly captured -- angle-matched passive
subtraction reduces but does not fully cancel a genuine inertial-coupling
artifact from the ramp's own acceleration profile (the active ramp's dip is
~3.7x larger than the matched no-stim ramp's, so subtraction leaves a real
residual). This does NOT change the 2026-07-22 sign-convention decision
above (raw sign, no per-step correction is still correct -- there is no
"true" sign to correct TOWARD when part of the signal is motion artifact,
not muscle force) but DOES mean `force_yTorque_N` specifically should be
treated with LOWER confidence during isovelocity's actively-moving phase
than during any static hold -- a caveat for whoever next uses this axis for
anything quantitative (not just sign-checking).

**Muscle-force extraction method (mean-over-window vs. alternatives),
2026-07-22, PI follow-up on `bass16_forceextractionmethod.png`:** "the max
in window is probably NOT the best way... it should be calculated from the
smoothed black line." Built 4 candidates (`muscleforcemethodcompare.png`):
A = current (MEAN over the full active window), B = MAX raw sample,
C = peak of the smoothed trace, D = narrow-window (0.15s) mean of RAW
samples centered on the smoothed peak's timing. `bass16_forceextractionmethod.png`'s own dynamic-bookend
panel is the clearest case: the smoothed trace rises to a real peak
(~0.22-0.3 N*m) then fades WITHIN the same active window (twitch/tetanic
decline) -- A dilutes this to 0.09, while C (0.225) and D (0.184) both
track the visible peak much better, D being the more noise-robust of the
two (a single smoothed-line value vs. an averaged plateau).

Aggregate check (`muscleforcemethodsensitivity.png`, 92 real SNR-eligible
isometric steps, 3 fish) adds an important qualifier: for ISOMETRIC's
sustained, FLAT holds specifically, method D tracks A almost exactly
(near-1:1 correlation) -- the extraction-method choice barely matters
there, because there is no post-peak decay to dilute. B and C both show
large, noisy departures from A even for isometric (single-sample noise for
B; onset-transient sensitivity for C -- the smoothed trace's peak search
can catch the trace still settling into its plateau right at stim onset,
found and fixed for the compare figure by restricting C/D's peak SEARCH
window to the stim duration itself, excluding the deactivation tail, but
the underlying onset-transient sensitivity is a real property of C, not
fully eliminated). **Net finding: the dynamic-bookend problem is about
TRANSIENT (rise-then-decay) windows specifically, not a universal flaw in
method A** -- isovelocity's moving-ramp windows likely share this same
transient character (visibly non-flat in the compare figure's isovelocity
panel) and are a second candidate for the same fix; a firm method choice
is still open, this is a comparison, not yet a decision.

Also clarified the same day: `bass16_forceextractionmethod.png`'s blue
"MEAN" line is the mean of the ACTIVE (stim) window -- i.e. method A itself,
what currently gets read off as the "active" reading before baseline
subtraction. It is NOT a baseline (a source of PI confusion, since
`bass16_passivebaselinemethod.png` also has a blue "MEAN" line, but of the
PRE-STIM window instead) -- the extraction figure's subtitle now says so
explicitly.

**FV superplot built 2026-07-22** (`FVsuperplot_isovelocity_pooled(_normalized).png`,
`R/superplot_fv_pooled.R`, PI-requested "similar to FL superplots, raw and
normalized"). Mirrors the FL superplot's V=0-only rule, axis-swapped: pools
a V=0 anchor (isometric L0 reps + isovelocity's own V=0 holds + dynamic L0
bookends, same 3 sources FL already uses) plus isovelocity's actual moving
ramps via the REAL angle-matched batch (angle_matched(/_cross_trial) only,
same exclusion as the y-torque correction above). X-axis is
`shortening_strain_pct` (a strain RATE, %/s, for isovelocity) rather than
raw commanded `operating_point` -- verified empirically that two steps at
the same |op| but opposite raw sign can share the SAME contraction_mode
(mirrored motor directions across sides), which raw op would have pooled
as opposite directions; `shortening_strain_pct` already resolves this
muscle-centrically (concentric always +, eccentric always -), matching
`03_analyze.R`'s own documented intent ("FL x-axis / FV x-axis use
predicted MUSCLE strain / strain-rate"). Normalized F0 falls back from
trial+side to fish+side when a trial's own V=0 doesn't clear SNR_MIN
(isovelocity's embedded V=0 holds are usually too low-SNR to self-
normalize, observed SNR ~0.1-1.9 -- found while building this). Both fish
branches (concentric and eccentric) land ABOVE the V=0 (isometric) force in
the pooled result -- NOT the classic Hill-curve shape (concentric should
fall BELOW isometric, only eccentric should rise above it) -- plausibly
because the moving-ramp force estimate carries some of the same
uncorrected inertial-coupling contamination just diagnosed for
`force_yTorque_N` above (which would elevate both ramp directions
similarly, roughly independent of sign). Flagged in the figure's own
subtitle as a limitation, not silently smoothed over; PROTOTYPE, not yet
canonical.

**Method D adopted as the production active-force extraction method,
2026-07-22 (PI decision after reviewing `muscleforcemethodcompare.png`/
`muscleforcemethodsensitivity.png` above): "Agreed. Let's go with D."**
Replaces the plain full-window MEAN in BOTH real pipeline paths that
previously shared it (`diag_force_extraction_baseline.R`'s header already
noted they were "identical MEAN-over-window logic, just a different window
variable name" -- so both needed the same fix, not just the diagnostic):
- `muscle_force_vector.R`'s `.mfv_window_peak_means()` -- the 6-axis vector
  path's ACTIVE window (isometric's `attach_vector_muscle_force()` and
  isovelocity's `compute_isovelocity_vector_batch()`); feeds
  `FL_isometric_uhatBoth.png`, `FV_isovelocity_uhatBoth.png`,
  `FLsuperplot`, `FVsuperplot`, and every SNR-summary export.
- `03_analyze.R`'s `.legacy_peak_window_mean()` -- the legacy zTorque-only
  path's `active_force_Nm` (`build_segmented_step_summary()`); feeds
  `muscle_force_Nm`/`muscle_force_Nm_interp` and everything downstream of
  those (frequency-sweep power, the pre-vector-force FL/FV outputs).
The PASSIVE/baseline window in both files is UNCHANGED (still the plain
mean via `.mfv_window_means()` -- a steady reference has no "peak" to
chase; only the ACTIVE window's extraction method changed).

**Bug found + fixed the SAME day, before trusting Method D in production**:
regenerating `FLsuperplot_isometric_isovelocity_pooled.png` with Method D
surfaced a large, unexplained spike at 0% strain (some individual points
>2N, vs a real ~0.1-0.8N range everywhere else) that had NOT been present
before the switch. Root cause traced to the dynamic trials' L0 bookend
contractions (`extract_dynamic_l0_bookends.R`) specifically: their stim
bursts are only ~0.05s long (verified directly on bass18's raw files --
`stim_t1_s - stim_t0_s` = 0.053-0.054s for all 4 bookends checked), a
fundamentally different regime from isometric/isovelocity's ~0.5-1s+
windows that `muscleforcemethodsensitivity.png` actually validated Method
D against. On a burst this short, the smoothed trace is often STILL
RISING at the search window's own trailing edge (verified: for a bass18
bookend with a spurious 2.3N output, the smoothed search-window max was
only 0.48 -- the peak was pinned to the search edge, not a genuine
interior maximum), so the found "peak time" sits right at the boundary;
Method D's fixed 0.15s narrow averaging window (bounded by the FULL
active+deactivation window, not the tiny search window) then extends ~3x
past that boundary into the deactivation tail, where it can land on a
large transient (verified directly: raw xforce on that same bookend
reached 5.16 N a few hundred ms into the deactivation window, versus <0.5N
everywhere within the actual 55-sample stim burst) unrelated to the pulse
itself -- most likely a mechanical/motor-related artifact (candidate
`motorstepartifact` territory), not a real muscle force.

Fixed with a DURATION guard, not a sample-count guard (the existing "too
few finite samples" fallback wasn't the problem -- there were plenty of
samples, just not enough elapsed TIME for the trace to reach a genuine
interior peak): both `.mfv_window_peak_means()` and
`.legacy_peak_window_mean()` now fall back to the plain full-window mean
(Method A) whenever the search window's own duration is shorter than the
0.15s narrow-averaging width, on the reasoning that a narrow window wider
than the entire search range can't possibly stay confined to it. Re-verified
directly against the exact bass18 bookend that produced the 2.3N/3.27N
spurious values -- with the guard, `.mfv_window_peak_means()` now returns
values identical to the plain mean (0.21-0.79N, matching the physiologically
plausible range) for all 4 of that trial's bookends. Confirmed the guard
does NOT change isometric/isovelocity's own Method D behavior (their
search windows are comfortably >0.15s). Full pipeline + SNR-summary export
+ both superplots rerun clean on all 3 fish after the fix; the spurious
spike is gone from `FLsuperplot`, and only Method A ever fires for dynamic
bookends specifically -- isometric and isovelocity's genuinely long/
sustained or ramping windows still get Method D as intended.

**`FLsuperplot_..._normalized.png` is CONCAVE-UP (U-shaped) -- ROOT-CAUSED
2026-07-22, NOT yet patched (prototype, PI decision pending).** PI question:
why is the NORMALIZED FL superplot U-shaped (F/F0 high at both strain
extremes, low near L0) -- inverted from the textbook bell curve? Investigated
all four candidate causes with instrumented per-point data (dumped `pooled`
to CSV, 241 points, 3 fish; also re-ran the whole superplot with a temporary
`BENDER3_MFV_FORCE_PLAIN_MEAN` env toggle forcing Method A, since removed):
- **RAW superplot is NOT U-shaped** (group-mean empirical force 0.08-0.44 N,
  a gentle peak near L0), so the force estimates are not "wrong everywhere" --
  the U is introduced by F/F0 NORMALIZATION only. This immediately rules out a
  global extraction/baseline/sign problem.
- **(d) Sign: RULED OUT.** Raw empirical force is 94-100% positive in EVERY
  strain region (near-L0 <2.5%: 94.2% pos; mid: 100%; extreme >=15%: 100%) --
  no concentric/eccentric sign cancellation near L0.
- **(b) Passive baseline: RULED OUT** as the driver -- isometric's own-step
  baselines yield all-positive, reasonable raw forces; the U is not in the raw
  curve, so a strain-dependent baseline bias is not producing it.
- **(c) F0 normalization + (a) Method D, INTERACTING: CONFIRMED cause.** The
  extreme strain bins (|strain| >= 15%) are populated ONLY by isometric steps,
  and at |strain| >= 20% almost entirely by ONE trial, `bass17_bender_15_
  isometric` (bass17 = the documented ~2x-noise-floor fish). Its off-L0 forces
  (~0.07-0.17 N) get divided by its own near-noise-floor L0 F0, giving F/F0 =
  6-17x -> the U's arms. Near-L0 bins are dominated by dynamic bookend +
  isovelocity-V0 points that self-normalize to ~1 by construction, pinning the
  trough. The TRIGGER is the Method D switch: `bass17_bender_15`'s L0 reps have
  `activation_snr` = 1.18/1.39/2.27/2.36 under Method A (all BELOW the 3.0 F0
  gate -> F0 = NA -> bass17 dropped from the normalized plot) but 3.25/3.26/
  5.09/5.36 under Method D (all ABOVE 3.0 -> admitted as F0, with force down to
  0.0006 N geometric). Method D roughly DOUBLES the SNR of marginal low-SNR
  flat holds because the peak-centered narrow window sits on the local noise
  peak instead of averaging the whole window -- pushing bass17's L0 across the
  gate. Extreme/trough ratio of the normalized empirical group mean: 4.8x
  under Method D vs 1.6x under Method A (geometric: 11.5x vs 1.7x) -- i.e. the
  U is overwhelmingly a Method-D-admitted-F0 artifact, not biology.
- **Confirmed by direct comparison**, not inference: the raw plot (Method D)
  is fine; the U appears only in the F/F0 companion; and forcing Method A
  removes it by (correctly) excluding bass17's noise-floor L0 as F0.
NOT patched: `FLsuperplot`/`FVsuperplot` are explicitly prototype. The clean
fix is a threshold decision for the PI -- the `activation_snr >= 3` gate on F0
measures activation-above-noise, not absolute force magnitude, so a near-zero-
force L0 hold can still pass it and become a pathological denominator. Candidate
fixes (propose one, PI to choose): (1) also gate F0 by absolute magnitude (F0
must exceed the baseline force-noise floor, not just SNR); (2) require n_fish>1
per strain bin before drawing a group-mean point (the extreme arms are
single-trial bass17); (3) treat the RAW superplot as the trustworthy one and
keep normalized as a flagged prototype. The RAW plot is unaffected and remains
the reliable read for now.

**REFINED / PARTLY REVISED 2026-07-22 (PI-requested decomposition + all-steps
overlay tests).** The PI was skeptical of the "(b) passive baseline: RULED OUT"
call above (baselines not settled). Two read-only tests were run:
- **eff_arm-vs-angle: RULED OUT as a cause.** `eff_arm^2 = d^2 + r_m^2*cos^2(theta/2)`
  is dominated by the dorsoventral offset `d` (~91-100 mm), so it barely moves
  with bend (0.998-1.000x out to 37 deg; denominator inflation <=0.3%). The
  geometric moment-arm does not create the concave-up.
- **active/passive/residual decomposition (geometric u_hat, isometric, 92 steps,
  `/tmp/diag_decomp.R`): the earlier "(b) ruled out" is TOO STRONG.** The
  reported muscle force is a SMALL residual (3-13% of |passive|) sitting on a
  large, bend-scaling, passive-dominated projection (active_proj ~= passive_proj,
  both up to -4 to -5.7 N at extremes). For bass17 the residual is a NEAR-CONSTANT
  ~3% of |passive| at every bend -> it tracks passive magnitude (subtraction
  leakage), not independent activation. So a strain-dependent passive-subtraction
  error DOES contribute to the raw concave-up; it is not purely an F/F0 artifact.
- **all-steps force-vs-time overlays (`R/diag_forcedev_allsteps.R`,
  `forcedevtiming_*_allsteps.png`) show WHY the residual is untrustworthy:**
  bass17 isometric traces climb MONOTONICALLY for ~1 s AFTER stim-off (baseline
  DRIFT / viscoelastic creep, not a 0.3 s twitch); every isovelocity step has a
  large NON-ZERO pre-stim muscle force (+-3 N) because the windowed-MEAN passive
  subtraction leaves the moving passive+inertia uncancelled pointwise.
- **The concave-up is NOT universal across specimens** (per-fish quadratic +
  magnitude/signed correlation on the isometric residual): **bass17** = clean
  SYMMETRIC U (cor(F,|strain|)=+0.93, cor(F,signed)=+0.01; and it is drift-shaped
  -> baseline artifact); **bass18** = MONOTONIC, force rises toward LENGTHENING
  (cor(F,signed)=+0.78, L0 0.27 N -> +25% 0.56 N; passive-tension-on-stretch
  signature, NOT a U); **bass16** = essentially FLAT (cor(F,|strain|)=+0.08). The
  pooled FL U therefore emerges from heterogeneous per-fish shapes with the
  extreme bins bass17-dominated, not a shared curve. Interpretation is blocked
  until the isometric baseline drift and isovelocity pointwise subtraction are
  fixed (item-2 work).

**Isometric drift + FATIGUE + STIM tested 2026-07-22 (PI-requested).**
- **Baseline drift (`R/diag_isometric_baseline_drift.R`,
  `isometricbaselinedrift.png`):** the ~6 s isometric hold's raw torque creeps
  monotonically (viscoelastic relaxation, amplitude ~|bend|); the static
  pre-stim baseline leaves it in, so force drifts up for seconds after stim-off.
  A pre->post INTERPOLATED baseline (existing `passive_force_Nm_interp` scheme)
  removes it cleanly (post-baseline residual 0.005-0.010 -> 0.000). BUT it does
  NOT flatten the in-stim PEAK the FL superplot samples (peak is early, where
  static ~= interp), so the concave-up survives drift correction.
- **Stim voltage: RULED OUT.** Every isometric step commands a CONSTANT 5.00 V
  (recruitment blocks = left_unilateral/right_unilateral SIDE selection, NOT
  amplitude recruitment). Force-per-volt normalization is a no-op (cor(F,|strain|)
  unchanged to 2 dp for all fish). No per-step stim variation exists to correct.
- **Fatigue: REAL but NOT the cause of the concave-up.** L0 force decays across
  the session (bass17 L0 0.043->0.011 N, cor(F,step_order)=-0.66; bass18 -0.73).
  DECISIVE test (`R/diag_concaveup_fatigue_stim.R`, `concaveupfatiguestim.png`):
  WITHIN a single block (same fatigue state, same 5 V) force rises monotonically
  with |bend| from that block's own fresh L0 (bass17 within-block cor(F,|strain|)
  =+0.84/+0.91). Fatigue pushes the OPPOSITE way (each block starts at its fresh
  L0), so it cannot create the arms; it only drags the pooled L0 reps down across
  the session, mildly deepening the trough AND destabilizing F0 for the
  normalized superplot. A fatigue correction (e.g. per-block fresh-L0 as F0)
  addresses F0 stability, a SEPARATE problem, not the arms.
- **Conclusion:** the residual concave-up = within-block force-proportional-to-
  |bend| = the passive-subtraction residual (small difference of large passive),
  confirmed independent of drift, fatigue, and stim. The real fix is accurate
  (pointwise/interp) passive subtraction (item-2).

**Improved isometric passive subtraction PROTOTYPED 2026-07-22 (PI-requested,
`R/diag_isometric_passive_models.R`, `isopassivemodels_1..4.png`).** Compared
three passive models on the geometric-u_hat projected force (= muscle_force_
vector_geom_N): M0 static pre-stim mean (current); M1 pre->post linear interp
(pointwise); M2 viscoelastic relaxation loess-fit over quiescent pre+post
samples, subtracted POINTWISE then Method D.
- **MODEL-FREE proof the concave-up is a stale-baseline artifact:** a NO-
  activation quiescent window minus M0 (zero-muscle control) reproduces 121% /
  169% / 61% of the bass16 / bass17 / bass18 |force|-vs-|strain| slope. The M0
  pre-stim baseline is sampled ~0.5-1 s before the active window on a relaxing
  passive, so it is STALE; the leftover drift scales with |bend| and masquerades
  as the FL arms.
- **Payoff cor(F,|strain|) M0->M1->M2:** bass16 +0.19->+0.11->+0.03; bass17
  +0.93->+0.91->+0.25; bass18 +0.57->+0.52->+0.42. M1 (linear interp) barely
  helps; pointwise M2 removes most of the artifactual concave-up WITHOUT
  creating a spurious bell. bass18's residual monotonic rise toward lengthening
  survives all baselines -> genuine (not the artifact).
- **CAVEATS / open decision:** (1) M2 interpolates the passive across the ~0.3 s
  stim gap (unobservable); it assumes the contraction does not discontinuously
  perturb the passive. (2) For the low-force fish (bass17) the true force is
  below the passive-drift floor -> its FL SHAPE is unresolvable and must be
  magnitude/SNR-gated, not reported as flat-or-bell. Production change (porting a
  pointwise relaxation-aware passive into the vector path + a magnitude gate) is
  NOT yet made -- PI decision pending on M2 vs M1 and the gating rule.

## Where things live (code map)
- `muscle_force_vector.R` — core: baseline subtraction, û construction
  (empirical + geometric), wrench->force solve, sign standardization,
  sono LP/L0, per-step SNR. ACTIVE-window extraction is
  `.mfv_window_peak_means()` (Method D, added 2026-07-22, duration-guarded
  fallback to the plain mean); PASSIVE/baseline window stays
  `.mfv_window_means()` (plain mean, unchanged).
- `03_analyze.R` — legacy zTorque-only path. ACTIVE-window extraction is
  `.legacy_peak_window_mean()` (Method D, added 2026-07-22, same
  duration-guarded fallback as above, independently defined since this
  file isn't always co-sourced with `muscle_force_vector.R`).
- `plot_muscle_force_vector.R` — all vector-force plots, confidence-alpha
  flagging, emp-vs-geom faceting.
- `superplot_fl_pooled.R` / `superplot_fl_pooled_snr_passing.R` — pooled
  cross-individual FL superplot (isometric + isovelocity instantaneous
  sweep), unfiltered and SNR-filtered.
- `superplot_fv_pooled.R` — pooled cross-individual FV superplot (added
  2026-07-22), reuses `superplot_fl_pooled.R`'s dynamic-bookend helper
  rather than re-deriving it.
- `diag_ytorque_sign_examples.R` — `ytorquesignexamples`/`axissnrcomparison`
  diagnostics (2026-07-22 correction: real angle-matched isovelocity data,
  not the static-baseline artifact v1 used).
- `diag_ytorque_inertial_timing.R` — `ytorqueinertialtiming`(`_stats`)
  diagnostics (2026-07-22, added): active-vs-matched-passive-ramp timing
  evidence for isovelocity's y-torque inertial-coupling artifact.
- `diag_force_extraction_methods_compare.R` — `muscleforcemethodcompare`/
  `muscleforcemethodsensitivity` diagnostics (2026-07-22, added):
  MEAN-vs-MAX-vs-smoothed-peak-vs-narrow-window extraction method
  comparison.
- `diag_forcedev_allsteps.R` — `forcedevtiming_*_allsteps` diagnostics
  (2026-07-22, added, PI-requested): cross-fish per-step vector-muscle-force
  vs time overlays colored by strain / strain rate, to sanity-check the
  active/passive/residual decomposition by eye (reveals bass17 isometric
  baseline drift + isovelocity non-zero pre-stim offset).
- `diag_isometric_baseline_drift.R` — `isometricbaselinedrift` diagnostic
  (2026-07-22, added): static vs pre->post interpolated baseline per isometric
  step, showing the ~6 s viscoelastic creep and that interp removes it but not
  the sampled in-stim peak.
- `diag_concaveup_fatigue_stim.R` — `concaveupfatiguestim` diagnostic
  (2026-07-22, added): rules out stim-voltage (constant 5 V) and fatigue as the
  concave-up cause via the within-block force-|bend| test.
- `export_snr_summary_figures.R` — per-fish SNR-gated export pass into
  `figs_summary/`.
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
