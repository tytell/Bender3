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

**M2 IMPLEMENTED IN PRODUCTION 2026-07-22 (PI chose M2; gating deferred).**
`.mfv_isometric_relaxation_passive()` (R/muscle_force_vector.R) replaces the
static pre-stim window mean for the ISOMETRIC vector passive: per channel it fits
the viscoelastic relaxation to the step's quiescent samples (pre-stim window ..
stim onset, plus stim_t1+relaxation_s .. post-baseline end) and evaluates it
inside the stim gap at EACH channel's own active-peak time (found the same way
.mfv_window_peak_means() finds it), so the passive is subtracted at the instant
the active value is sampled -- the scalar-per-channel equivalent of pointwise
subtraction. Per-channel fallback to the old static mean when the fit is
unreliable (<8 quiescent samples / search shorter than the peak window / loess
failure). Only `attach_vector_muscle_force()`'s isometric loop + the new helper
changed; `.mfv_finalize_step()` untouched, so isovelocity (angle-matched passive)
is byte-for-byte unchanged. Verified on the production path (all 3 fish,
isometric-only, every step used relaxation_fit, no fallbacks): cor(Fgeom,
|strain|) bass16 +0.02, bass17 -0.01, bass18 +0.43 -- the strong isometric
concave-up (bass17 was +0.93) is gone; bass18's genuine monotonic rise survives.
(bass17's -0.01 vs the prototype's projected-g +0.25 is the expected per-channel-
vs-projected difference for a near-noise-floor fish whose shape is unresolvable.)
FL/FV superplots regenerated: the RAW isometric arms are removed and the
NORMALIZED plot is flat ~1 -- EXCEPT one pathological bass17 strain-0 point
(near-zero F0 -> F/F0 ~ -64) that drags the y-axis. That single outlier is the
DEFERRED low-force gating decision, not a passive-subtraction failure. ts NOTE:
the isometric force_ts display trace still subtracts the peak-time passive as a
CONSTANT pointwise (finalize unchanged); the sampled SCALAR is the corrected
quantity, the trace is display-only.

**Negative-force question answered + ts trace ported to pointwise 2026-07-22
(UNCOMMITTED, R/muscle_force_vector.R).** PI asked whether negative force
values represent muscle "pushing away" physiologically. Traced the specific
pathological point (bass17's normalized-plot outlier): it is NOT from the
isometric protocol (all isometric L0 reps positive, 0.003-0.04 N) -- it is from
`bass17_bender_16_isovelocity`'s own embedded V=0 holds (last isovelocity trial,
heavily fatigued), where force_geom_N = -0.0005/+0.0007 N (left) and
+0.0012/-0.0033 N (right), all at the fish's documented noise floor, held at
V=0/L0 (not at a lengthening/shortening extreme where an eccentric-resistance
story could apply). Skeletal muscle cannot physically push (only generate
tension); a real physiological "push" is not possible here. Conclusion: the
sign is measurement noise on a near-zero, fatigued-out contraction, NOT a
biomechanical push -- consistent with the module header's existing finding that
the negative minority is "not explained by ... concentric/eccentric position."
The fix is the already-deferred MAGNITUDE-based F0 floor (reject F0 built from
near-noise-floor reps), not a sign-based rule -- still deferred per PI ("decide
later").

Separately, PORTED the M2 pointwise passive into the `ts` (force_ts) DISPLAY
trace, closing the gap flagged above: `.mfv_isometric_relaxation_passive()` now
also returns `fits_T` (the per-channel loess fit objects for xtorque/ytorque/
ztorque, reused not refit), which `.mfv_finalize_step()` accepts as an optional
`passive_curve_fits` argument and evaluates POINTWISE at every ts sample's own
time (falling back to the constant `pass$T[i]` per channel when no fit is
available). Isovelocity/dynamic callers omit the argument (defaults to NULL),
so their behavior is byte-for-byte unchanged. Verified: (1) scalar cor(Fgeom,
|strain|) unchanged from the M2 scalar commit (bass16 +0.02, bass17 -0.01,
bass18 +0.43) -- the ts port touches ONLY the display trace, not the sampled
scalar; (2) the fitted curve is confirmed genuinely time-varying, not frozen at
the constant (one bass17 step: range [-0.051793,-0.051541] vs constant
-0.051730, a real ~2.3e-4 N*m drift across the window); (3) regenerated
`forcedevtiming_isometric_allsteps.png` -- bass17's traces now decay back to
~0 within ~0.5-0.7 s of stim-off instead of climbing for ~1 s (the drift
artifact documented above is now visibly gone from the DISPLAY trace too, not
just the scalar). NOT YET COMMITTED per PI instruction.

**Magnitude-based F0 floor IMPLEMENTED 2026-07-22 (PI-approved; the deferred
low-force gate, now decided).** Root-caused the surviving normalized-FL outlier
(F/F0 ~ -64): NOT the isometric protocol (all its L0 reps 0.003-0.04 N) but
`bass17_bender_16_isovelocity`'s own embedded V=0 holds (last, fatigued
isovelocity trial), force_geom_N ~ +-0.0005-0.003 N at the noise floor, held at
V=0/L0 (NOT a lengthening extreme, so not an eccentric "push"). Skeletal muscle
cannot push; the negative sign there is noise on a near-zero fatigued
contraction, not physiology -- so the fix is a MAGNITUDE floor, not a sign rule.
New helper `mfv_gate_f0()` (R/muscle_force_vector.R) keeps an L0/V=0 rep as an
F0 contributor only if it passes BOTH activation_snr >= MFV_UHAT_SNR_MIN AND
|force| >= its own baseline_force_noise_N (the pre-stim ||force-channel|| SD,
now surfaced as a step_summary column). Rationale: the SNR gate qualifies the
force VECTOR's amplitude, but the F/F0 DENOMINATOR is the PROJECTED
muscle_force_vector_N, which can be near-zero even at high SNR -> a magnitude
floor on the denominator itself is required. Wired into all F0 blocks: 3 in
superplot_fl_pooled.R (isometric L0, isovelocity V=0, dynamic bookend), 2 in
superplot_fv_pooled.R (trial- and fish-level); the SNR-passing FL companion
inherits it by sourcing. Result: FL normalized F0 available for 42/46 trials (4
near-noise-floor trials correctly dropped); the -64 outlier is GONE, the plot
now reads as a sane peak-near-L0 curve (~0.3-1.5 F/F0). RAW plots unaffected
(F0 gate only touches normalized). Plumbing touched: .mfv_finalize_step row,
.mfv_empty_out_cols, .mfv_assign_row (add baseline_force_noise_N); FV extractors
carry the column through their transmute. NOT YET COMMITTED (continuing the
prior leave-uncommitted instruction).

**Isovelocity passive baseline DIAGNOSED 2026-07-22 (PI-requested "how does the
logic compare to isometric").** New READ-ONLY diagnostic
`R/diag_isovelocity_passive_models.R` (canon token `isovpassivemodels`, 3 PNGs).
Current production (compute_isovelocity_vector_batch -> .mfv_interp_ramp_onto)
subtracts an ANGLE-matched no-stim ramp at the same signed velocity -- the
correct raw material (it carries angle-elastic + velocity-viscous + inertial in
one) -- but COLLAPSES it to a scalar window-MEAN and subtracts that from the
Method-D (peak) active. Because the active window SWEEPS through angle, the
passive varies enormously across it (range median 2.1/3.0/3.1 N bass16/17/18,
up to ~6 N), so a single mean is a poor stand-in for the passive at the peak's
own angle: pointwise-vs-production muscle force differs by up to +-4.6 N. FV
payoff: the window-MEAN passive manufactures a CONCAVE-UP FV in low-force
bass16/17 (the SAME artifact class as the FL concave-up); pointwise angle-matched
subtraction flattens it to ~0 (those fish sit at the noise floor); bass18 (real
force) goes from a flat production FV to a plausible bell but overshoots negative
at high |v| -- residual inertial-transient / angle-alignment error (a flagged
2nd-order limit). COMPARISON TO ISOMETRIC: isometric passive varies only in TIME
(viscoelastic relaxation, 1 d.o.f., bracketed by quiescent pre/post samples ->
the M2 time-fit solved it cleanly); isovelocity passive varies in ANGLE (elastic,
large) + velocity (viscous) + direction (hysteresis/transients), so the
time-relaxation fit does NOT transfer. The direct analog of the isometric
pointwise fix is to subtract the angle-matched ramp POINTWISE (sample-by-sample
by angle), then take Method D on the delta -- same family of fix (kill a
stale/averaged baseline that fakes a velocity/strain-correlated shape), one extra
complication (motion transients) that isometric doesn't have. NO production change
yet -- diagnostic + comparison pending PI decision on whether to port the
pointwise angle-matched subtraction into compute_isovelocity_vector_batch.

**Isovelocity POINTWISE passive IMPLEMENTED + 3-tier FL/FV plots BUILT
2026-07-22 (PI-approved "proceed" on the proposed order).** (1) Production change
in compute_isovelocity_vector_batch: the mean-collapsing .mfv_interp_ramp_onto
was replaced by .mfv_ramp_passive_pointwise, which interpolates the
velocity+signed-direction-matched no-stim ramp onto EACH of the step's own
samples by ANGLE (full-length td6 vectors). .mfv_window_peak_means gained an
optional passive_pw arg (subtract the pointwise passive per-channel BEFORE peak
finding -> Method D runs on the ACTIVE-minus-passive delta; pass is then zero).
.mfv_finalize_step gained passive_pw_T so the force_ts display trace subtracts
the SAME pointwise passive per sample (scalar + display on one baseline).
Velocity matching, the within->cross-trial->static fallback, and the isometric
path are all UNCHANGED (attach_vector_muscle_force still uses the M2 relaxation
fit -- correct: it only handles isovelocity V=0 holds, which are non-moving).
Verified on bass18: runs clean, angle_matched steps now pointwise; the
static_baseline_fallback rate is inherited from the pre-existing matching logic
(unchanged), not introduced. (2) Rebuilt FVsuperplot: the GEOMETRIC-u_hat FV is
now a bell (concave-up GONE, peak near V=0), matching the diagnostic's
prediction; the empirical-u_hat FV stays U-shaped (empirical direction is
unstable for these low-force moving steps -- a separate known limitation), and a
negative overshoot persists at high |eccentric v| (the flagged inertial-transient
residual). FL superplot unchanged (its isovelocity points are V=0 holds, off the
moving path). (3) New canonical tiers script R/superplot_fl_fv_tiers.R (tokens
fltiers/fvtiers): within-trial FL+FV (each line one trial, RAW geom force, no
normalization -- exposes bass17 isometric forces ~0.01 N at the noise floor) and
within-protocol isometric-only pooled FL (un-mixed companion to the
across-protocol FLsuperplot; grand mean concave-down, peak near L0). Per-fish
run_fv_fl_power_pipeline rebuilt for all 3 fish. NOT YET COMMITTED (continuing
the leave-uncommitted instruction).

**Post-commit flags investigated 2026-07-22 (PI-directed "commit and push, then
address two flags").** Two residual concerns from the pointwise-passive rebuild:
(1) a negative-overshoot at extreme eccentric velocity in the FV superplots,
(2) the empirical-u_hat FV staying U-shaped. Initial pass (before the corpus-gap
fix below) ruled out turnaround-transient and low-SNR-noise as the overshoot's
cause (all strongly-negative points had HIGH SNR, median 11.8) and characterized
the empirical-u_hat bias (`dF/||dF||` is always-positive-biased at extreme
velocity, masking the true eccentric sign) -- concluding geometric/longitudinal
u_hat should be the reported-primary FV method. This left flag 1 open pending a
PI decision, since no simple fix (const-velocity restriction, SNR gate) resolved
it.

**Flag 1 ROOT-CAUSED + FIXED 2026-07-22 (PI-directed: "use the stim-off ramps to
calculate position AND velocity-specific isovelocity baseline -- that's the KEY
ingredient, cooked into the experimental design").** Verified the pointwise
angle-matched subtraction (above) DOES do exactly that -- but only for an EXACT
velocity match (tol 1.0 deg/s). Root cause: bass18's stim-off ("no-stim") ramps
only cover a SUBSET of its commanded velocities (107/214/320 deg/s from
`bender_03`); its OTHER 5 isovelocity trials move at a DIFFERENT, non-overlapping
set (142/285/427 deg/s) with NO stim-off ramp anywhere in the corpus at those
speeds. Those steps fell through to the STATIC single-angle pre-stim baseline --
meaningless for a swept ramp -- which was confirmed to be the ENTIRE overshoot:
`angle_matched` bass18 rows were 100% positive (median geom force +1.51 N, min
+0.77 N) while `static_baseline_fallback` rows were the overshoot itself (median
-0.36 N, min -5.86 N). FIX (`R/muscle_force_vector.R`,
`compute_isovelocity_vector_batch`): added a fallback tier between the exact
match and the static baseline -- the NEAREST same-sign-velocity stim-off ramp
(still angle-matched via `.mfv_ramp_passive_pointwise`, so its overlap guard
still applies; loops through candidates in order of increasing |Δv| until one
clears the guard). New `passive_source` value: `"angle_matched_nearest_v"`.
Verified via the real `compute_isovelocity_vector_batch()` output: bass18's
previously-fallback steps (142/285/427 deg/s) are now 100% positive too (median
geom force +1.64 N, min +0.72 N, i.e. IN LINE with the 107/214/320 series, as
physically expected) -- confirmed across the FULL corpus with `passive_source`
breakdown, not spot-checked. Also fixed two bugs found while reconciling the
`isovpassivemodels` diagnostic against this corrected production behavior: (a)
the diagnostic built its no-stim library PER-TRIAL, not per-fish, so it couldn't
even ATTEMPT the nearest-velocity fallback for trials with zero no-stim steps of
their own -- rebuilt as a fish-wide library (mirrors production's
`passive_library`); (b) the diagnostic's own force-projection reimplementation
called `uhat_geometric(op)` with `op` = a VELOCITY (deg/s), not a bend angle --
`uhat_geometric()` expects a bend ANGLE, so at high |v| it silently computed
cos/sin of a nonsense "angle" and flipped sign, manufacturing a SEPARATE spurious
negative dip that was never in production; fixed by switching to the SAME fixed
longitudinal u_hat = (1,0,0) `.mfv_finalize_step()` actually uses for
`muscle_force_vector_geom_N` on non-isometric categories. After both fixes, all
3 fish show clean pointwise FV curves with NO negative excursion anywhere in the
corpus -- the earlier "overshoots negative at high |v| -- residual
inertial-transient" conclusion was an artifact of the corpus-gap bug + the
diagnostic's own uhat-misuse bug, not a real limitation of pointwise angle-matched
subtraction. Rebuilt: `run_fv_fl_power_pipeline.R` (all 3 fish, refreshes
`FV_isovelocity_uhatBoth.png` etc.), `superplot_fv_pooled.R`,
`superplot_fl_pooled.R`, `superplot_fl_fv_tiers.R`, `diag_isovelocity_passive_
models.R` (canonical `isovpassivemodels` PNGs + header FINDING rewritten).
Flag 2 (empirical-u_hat U-shape) is UNCHANGED by this fix (it's a property of
the empirical direction estimate at low force, not the passive baseline) --
still pending a PI decision on primary-vs-caveat reporting.

**TARGET SHAPE CORRECTION for FV, 2026-07-22 (PI-clarified).** Several entries
above (and `isovpassivemodels`'s header) described the fixed bass18 pointwise FV
curve as "a plausible bell" / "bell-shaped" -- WRONG terminology, left as-is
above for the historical record but corrected going forward. FV is not supposed
to look like FL (a bell/symmetric peak at some optimum) -- the physiological
target is a Hill hyperbola: monotonic-DECREASING force with increasing
shortening (concentric) velocity, with eccentric (lengthening) force >=
isometric >= concentric, NOT a symmetric hump centered at V=0. Checked the
post-fix pointwise data against that target (SNR-passing right-side points,
grouped by |velocity|, `activation_snr >= 3`):

| \|v\| (%/s) | bass18 concentric F (N) | bass18 eccentric F (N) | Hill-consistent (ecc > con)? |
|---|---|---|---|
| 127 | 1.37 | 1.68 | yes |
| 255 | 1.03 | 1.65 | yes, gap widens |
| 382 | 1.83 | 1.76 | no -- concentric slightly exceeds eccentric |

bass16/bass17 show NO consistent eccentric-vs-concentric ordering at any
velocity (both fish's values are small, sign-unstable, and drift toward/past
zero with |v| -- noise floor, not signal). So bass18 pointwise isn't just the
best-LOOKING curve among the three fish -- it is the ONLY one that reproduces
the correct Hill-type SIGN relationship (eccentric > concentric), which is a
much stronger validation of the pointwise angle-matched fix than shape-matching
by eye. The one exception (382 %/s, where concentric slightly exceeds eccentric)
is the single highest velocity tested and also the nearest-velocity-fallback
point from the flag-1 fix above -- flagged for a closer look, not yet explained
(candidates: genuine small-N noise (n=3 concentric vs n=2 eccentric trials at
that speed) or a residual velocity-dependent passive/inertial mismatch at the
extreme end of the nearest-velocity fallback's angle-overlap tolerance).
Corrected the "bell" language in `FIGURES_README.md`'s `isovpassivemodels`
entry and in `R/diag_isovelocity_passive_models.R`'s own header FINDING to
state the Hill-hyperbola target explicitly, per the "keep comments/names
truthful" rule -- code/docs describing CURRENT understanding must not use
FL's bell language for FV's target shape.

**CORRECTION to the above, same day (PI pushback: "I'm not seeing how bass18
resembles the Hill-type relationship. Did you produce the plot yet?").** The
table above compared concentric-vs-eccentric PAIRS at each |v| and never
plotted/included the V=0 isometric anchor -- an incomplete test, and the PI
was right to push back before accepting "Hill-consistent" on that basis alone.
Built `R/diag_isovelocity_hillcheck.R` (canon token `isovhillcheck`, 1 PNG) to
actually plot the full curve (isovelocity's own V=0 hold + moving ramps,
SNR-passing summary line, all points shown). Result: bass18 right is NOT a
clean monotonic Hill curve -- it's a "W": eccentric plateau (~1.65-1.76 N,
-127 to -382 %/s) -> drops SHARPLY to 0.56 N at V=0 -> back up to 1.37 N (127)
-> down to 1.03 N (255) -> back up to 1.83 N (382). Root cause of the V=0
notch: its 3 SNR-passing contributors are all from the SAME weaker trial-set
(bender_06/07/08/09/10) while the STRONGER trial (bender_03) has V=0 values
(1.69/2.41 N) fully consistent with the eccentric plateau but FAILS the SNR
gate (its V=0 hold is brief/embedded, not sustained) -- so pooling V=0 and
moving points ACROSS trials mixes different fatigue states and manufactures
the notch; it is not itself evidence against Hill, but it also means the
pairwise eccentric>concentric comparison above was not a full/clean test
either. bass16/17 show a PEAK AT V=0 with decline on both sides (a tent) --
neither Hill's plateau-then-decline nor flat noise; unexplained, flagged.
REVISED CONCLUSION: the eccentric>concentric ordering at 127/255 %/s remains
directionally suggestive for bass18, but "bass18 pointwise reproduces the
Hill relationship" was an OVERCLAIM -- the honest state is "directionally
suggestive, confounded by cross-trial fatigue when V=0 is included, not yet
settled." A proper test needs each trial's OWN V=0 vs its OWN moving steps
(within-trial only), not pooled across trials -- not yet built.

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

## SNR-based confidence gating audit: does the ratio conflate noise with
## small-but-real force? (2026-07-22)

**Trigger.** PI audit request, prompted by a live example already logged
above (`isovhillcheck` V=0 notch): `activation_snr` is a RATIO
(`||dF_force|| / baseline_force_noise_N`); a low ratio can mean EITHER (a)
elevated noise, OR (b) genuinely small true force -- the ratio alone cannot
tell these apart, and a ratio-only gate risks silently erasing/de-weighting a
real physiological finding (weak/fatigued muscle, a brief V=0 hold, a point
near L0/V0 where near-zero force is the CORRECT answer). One fix already
exists for exactly this (`mfv_gate_f0()`, 2026-07-22, requires BOTH
`snr >= snr_min` AND `abs(force) >= noise`) but is NOT applied everywhere
else SNR gating happens. This entry is audit + diagnosis + proposal only --
**no production gating changed** (per `.cursorrules`: propose one fix at a
time, ask before changes touching >2 functions).

### 1-2. Every SNR gating site, what it does, and the conflation risk

| Site (file::function) | What it does with a failing point | Conflation risk (genuinely-small-but-real treated as noise)? |
|---|---|---|
| `muscle_force_vector.R::.mfv_finalize_step()` (`use_emp <- ... snr >= snr_min`, ~line 656) | SOURCE SELECTION, not a drop: below `snr_min`, the DIRECTION estimate falls back from empirical u_hat to geometric/longitudinal u_hat. The force itself is still computed and reported. | Low risk by design intent -- direction stability genuinely needs a ratio test (a small force's direction estimate really is noisier than a large one's, independent of whether the magnitude itself is real). Flagged so this ratio-only pattern isn't blindly copied elsewhere as if it were a magnitude judgment. |
| `muscle_force_vector.R::mfv_gate_f0()` (lines 96-119) | HARD DROP from the F0 (denominator) pool only: requires BOTH conditions. | **This is the reference fix**, not a risk site -- included as the baseline the other sites should match. |
| `muscle_force_vector.R::rank_isometric_steps_by_activation()` (`clearly_active` column, ~line 1071) | Nothing dropped -- a ratio-only boolean surfaced in a printed console table for human eyeball/veto (PI A6). | Frames the human reviewer's mental model on ratio alone; a genuinely-small-but-real step reads as "not clearly active" even if its magnitude is trustworthy. |
| `plot_muscle_force_vector.R::.mfv_confidence_factor()` / `MFV_CONFIDENCE_ALPHA` (lines 47-55) | ALPHA-FLAG only (never drops): 2-level ratio-only factor feeding `.build_vector_xy_summary()` (`FL_isometric_uhatBoth.png`, `FV_isovelocity_uhatBoth.png`) and `build_summary_plot_FV_L0_vector()` (`FVl0_isovelocity_uhatBoth.png`). | A confidently-small point is rendered visually IDENTICAL (same faint alpha) to a point that is pure noise -- no way to tell them apart by eye even though the plot never drops either. |
| `superplot_fl_pooled.R` / `superplot_fv_pooled.R` (F0 blocks) | F0 denominator: already uses `mfv_gate_f0()` (fixed). The PLOTTED NUMERATOR points (per-trial/per-fish/per-group means) have **no SNR or magnitude gating/flagging at all**. | Inconsistent in the opposite direction: F0 is doubly-gated but the force values folded into group means get NO confidence weighting whatsoever -- a pure-noise point and a confident point contribute equally to the mean. |
| `superplot_fl_pooled_snr_passing.R` | HARD DROP, ratio-only (`activation_snr >= SNR_MIN`), regenerates the entire `_snrPass` FL superplot with everything below threshold excluded. | **Direct hit**: a confidently-small real point is dropped identically to pure noise -- exactly the risk this audit is checking for. |
| `export_snr_summary_figures.R::snr_filter()` | HARD DROP, ratio-only, for the whole `_snrPass` per-fish figure family (`FL_isometric_uhatBoth_snrPass.png`, `FV_isovelocity_uhatBoth_snrPass.png`, `FVl0_isovelocity_uhatBoth_snrPass.png`, `uhatCompare_empVsGeom_snrPass.png`, `forceTime_*_uhatBoth_snrPass.png`). | Same direct-hit risk as above, across 5 figure families at once. |
| `superplot_fl_fv_tiers.R` | **No gating or flagging of any kind** -- every finite geometric-u_hat point plotted at equal opacity. | Not a drop risk, but an inconsistency with the rest of the pipeline's "never silently drop, always flag" policy: a low-SNR/noise point looks visually IDENTICAL to a high-confidence one here, worse than the 2-level alpha sites. |
| `diag_isovelocity_hillcheck.R` (`confident <- activation_snr >= snr_min`) | Point-level: ALPHA-FLAG only (never drops raw points). Summary-line level: the per-(fish,side,velocity) MEAN -- what actually carries the "Hill or not" scientific claim -- is computed from `confident`-only (ratio-only) points. | **This is the site that produced the documented V=0 notch** -- an effective drop at the level that matters for the conclusion, even though no raw point is deleted. See item 4 below. |
| `plot_fatigue_timeline.R` "Gate A" (`snr_pass <- activation_snr >= snr_min`) | Point-level: ALPHA-FLAG only (never drops). Trend line: the `loess` smooth is NOT restricted to SNR-passing points -- unlike `hillcheck`, it pools ALL points (pass+fail) with EQUAL weight. | The opposite inconsistency from `hillcheck`: here a genuinely unconfirmable (pure-noise) point gets full, unweighted influence on the trend, while there a confidently-small real point gets zero influence. Neither site treats "confidently small" as its own category. |

Also checked and found NOT to be gating sites (informational only): `diag_ytorque_sign_examples.R`'s `axissnrcomparison` plot draws the `MFV_UHAT_SNR_MIN` line as a reference annotation only; `diag_concaveup_fatigue_stim.R` carries `activation_snr` through as a column without filtering on it. The log's own "Gate A exploration" entry above mentions a "majority-of-steps-pass rule for trial plots" gating which trial plots get copied into `figs_summary/` -- this could NOT be found in current code; it most likely describes `export_snr_summary_figures.R`'s pre-2026-07-21 behavior, removed when the `figs_summary/` copy mechanism was cleaned out (see "Summary" in `FIGURES_README.md`) -- flagged as a stale reference, not a live site.

### 3. Proposed consistent scheme (NOT implemented)

Generalize `mfv_gate_f0()`'s two-condition test into a shared 4-level tier,
one function, sourced everywhere (mirrors how `mfv_gate_f0` is already the
one shared F0 gate -- "one source of truth," not four ad hoc re-derivations):

```r
mfv_confidence_tier <- function(force, snr, noise, snr_min = MFV_UHAT_SNR_MIN, k = 1.0) {
  snr_pass <- is.finite(snr) & snr >= snr_min
  mag_pass <- is.finite(force) & is.finite(noise) & abs(force) >= k * noise
  dplyr::case_when(
    snr_pass & mag_pass   ~ "confident",
    !snr_pass & mag_pass  ~ "confidently_small",   # the conflation case
    snr_pass & !mag_pass  ~ "unstable_magnitude",   # what mfv_gate_f0 already catches
    TRUE                  ~ "unconfirmable")
}
```

- **confident** (both pass): treat as now -- full weight/opacity, eligible as
  F0, eligible to seed empirical u_hat.
- **confidently_small** (SNR fails, magnitude passes): a REAL finding, not
  noise -- INCLUDE in summary/mean lines and F0 pools (already true for F0 via
  `mfv_gate_f0`), render at an intermediate alpha distinct from both
  "confident" and "unconfirmable" (currently lumped into "low confidence"
  alongside pure noise), but keep EXCLUDED from seeding empirical u_hat
  direction (a low-ratio point's DIRECTION is still unreliable even when its
  magnitude is real -- direction and magnitude are different questions, see
  item 1's u_hat-source-selection row).
- **unstable_magnitude** (SNR passes, magnitude fails): the case
  `mfv_gate_f0` was built to catch (high activation on the vector, near-zero
  on the projected scalar) -- exists OUTSIDE the F0 context too (13 points in
  the corpus below, "isovelocity_moving" and "isometric" categories). Should
  get the SAME caution `mfv_gate_f0` already applies (exclude from any
  denominator-like or direction-seeding role) wherever it currently isn't.
- **unconfirmable** (both fail): keep the CURRENT faint/never-drop treatment.

Concrete per-site changes this implies (proposal only, not yet made): (a) add
a 4-level `MFV_CONFIDENCE_ALPHA`-equivalent for the point-level flag sites
(`plot_muscle_force_vector.R`, `superplot_fl_fv_tiers.R` which currently has
none at all); (b) change `superplot_fl_pooled_snr_passing.R` /
`export_snr_summary_figures.R::snr_filter()`'s drop criterion from
`activation_snr >= snr_min` to `tier %in% c("confident","confidently_small")`
(drop ONLY `unstable_magnitude`/`unconfirmable`) -- this is the change most
directly responsive to the audit's opening question; (c) change
`diag_isovelocity_hillcheck.R`'s summary-line filter the same way; (d) decide
whether `plot_fatigue_timeline.R`'s trend line should match (currently pools
everything unweighted -- the opposite problem); (e) swap
`rank_isometric_steps_by_activation()`'s `clearly_active` boolean for the
4-level tier so the PI A6 human review sees the magnitude distinction too.
A magnitude multiplier other than `k=1` (matching `mfv_gate_f0`'s existing
default) was not explored -- `k=1` is proposed for consistency with the one
fix that already exists, not because alternatives were tested.

### 4. isovhillcheck V=0 notch, re-examined under the proposed scheme

Built `R/diag_snr_magnitude_audit.R` (see item 5) and ran it against the real
3-fish corpus. Direct answer: **under a magnitude-aware scheme, NONE of
bass18's isovelocity V=0 holds are unconfirmable.** All 9 finalized bass18 V=0
reps: 1 "confident" (bender_10, right, force 0.4284 N, SNR 3.37, noise
0.2342 N) + **8 "confidently_small"** (SNR 1.63-2.97, ALL below `snr_min`, but
EVERY ONE independently clears its own noise floor) -- including
`bass18_bender_03_isovelocity` right side's 1.9944 N (SNR 2.09, noise 1.0663)
and 1.3127 N (SNR 1.63, noise 0.8030) reps: exactly the strong trial the log
already identified as "consistent with the eccentric plateau but fails the
SNR gate." Zero bass18 V=0 rows are "unconfirmable" or "unstable_magnitude."
This confirms the log's diagnosis directly: the V=0 notch is a ratio-gating
artifact, not evidence of a real dip at V=0 -- the excluded points are large,
real, and consistent with the eccentric plateau; they fail SNR purely because
`bender_03`'s own noise floor (0.80-1.07 N) happens to be higher than the
weaker trials' (0.16-0.41 N), which a ratio test conflates with "the force
itself is unreliable." Including `confidently_small` points in
`hillcheck`'s summary line (proposal 3c above) would pull the V=0 anchor UP
toward the eccentric plateau, at least partially resolving the "W" shape --
not yet done (diagnostic only).

bass16/17 are a genuine MIX, unlike bass18: bass16's 8 V=0 rows are 5
confident / 1 confidently_small / 2 unconfirmable; bass17's 12 rows are 6
confident / 2 confidently_small / 1 unstable_magnitude / 3 unconfirmable. So
their already-documented "tent" shape (peak at V=0, decline both sides,
unexplained in the `isovhillcheck` entry above) is NOT purely a ratio-gating
artifact the way bass18's notch is -- some of their V=0 data genuinely cannot
be distinguished from noise on EITHER test, magnitude-aware or not.

Corpus-wide (all categories, all 3 fish, `snrmagnitudeaudit_2`): 42
"confidently_small" points and 13 "unstable_magnitude" points exist outside
the F0-denominator context `mfv_gate_f0` already protects -- isometric holds
also confidently_small mostly the noise-floor-limited category
("unconfirmable" = 31/92 isometric points, consistent with the already-
documented bass17 ~2x-noise-floor finding), while isovelocity_V0 (11/29),
isovelocity_moving (15/87), and dynamic_bookend (15/120) each carry a
non-trivial confidently_small fraction the current ratio-only gates at those
sites (items 6, 7, 9 in the table above) would treat identically to noise.

### 5. Diagnostic plots built (read-only, no production change)

`R/diag_snr_magnitude_audit.R` (canon token `snrmagnitudeaudit`, 3 files in
`figs_diagnostic/`, added to `FIGURES_README.md` in this same change):
1. `snrmagnitudeaudit_1_forceVsNoiseFloor.png` -- `|muscle_force_vector_geom_N|`
   vs. `baseline_force_noise_N` directly (log-log), quadrant-colored, faceted
   by category (isometric / isovelocity_V0 / isovelocity_moving /
   dynamic_bookend), reference lines at ratio=1 (magnitude floor) and
   ratio=`snr_min` (illustrative -- see the script's own caveat: `activation_snr`'s
   real numerator is `||dF_force||`, the 3-axis force-vector magnitude, NOT
   the projected scalar `force_geom_N` plotted here -- the same numerator
   mismatch `mfv_gate_f0`'s header already documents for the F0 case).
2. `snrmagnitudeaudit_2_quadrantCounts.png` -- quadrant counts by category,
   full corpus, to size the issue beyond one anecdote (numbers in item 4).
3. `snrmagnitudeaudit_3_isovelocityV0notch.png` -- the isovhillcheck V=0
   re-examination in item 4, log-log per fish.

All three are read-only re-aggregations of `attach_vector_muscle_force()` /
`compute_isovelocity_vector_batch()` / `.dynamic_l0_trial_rows()` output
(same production helpers `superplot_fl_pooled.R`/`superplot_fv_pooled.R`
already call) -- no new force computation, no production file touched.

### Status: STOPPING here per instruction

This entry is audit + diagnosis + proposal only. **No change was made to any
production gate** (`mfv_gate_f0`, `.mfv_confidence_factor`,
`superplot_fl_pooled_snr_passing.R`, `export_snr_summary_figures.R`,
`diag_isovelocity_hillcheck.R`, `plot_fatigue_timeline.R`, or
`rank_isometric_steps_by_activation()`). Per `.cursorrules` ("propose one fix
at a time," "if a fix requires changes in more than two functions, stop and
ask"), the proposal in item 3 touches more than two functions/files and needs
a PI decision on: (i) whether to adopt the 4-tier scheme at all, (ii) which
specific sites should change behavior first (the two hard-drop sites,
`superplot_fl_pooled_snr_passing.R` and `export_snr_summary_figures.R`, are
the highest-leverage/lowest-risk starting point since they are pure filters,
not part of the force-computation path), and (iii) the magnitude multiplier
`k` (proposed `k=1`, matching the existing `mfv_gate_f0`, but not
independently validated).

---

## SNR-magnitude conflation fix: IMPLEMENTED (2026-07-22)

PI direction: "Implement proposal. Would that make 4-tier scheme redundant or
still useful? If useful, keep that also." Answering the question first, then
the implementation.

### Is the 4-tier scheme redundant once implemented?

**No -- it IS the implementation, not a separable extra.** The proposal in
the entry above was never "add a 2-way real-vs-noise split" -- it was
specifically that a 2-way split (pass/fail) can't distinguish 3 genuinely
different situations a ratio-only gate collapses together: a point can (1)
pass both tests (`confident`), (2) fail the ratio but still be real
(`confidently_small` -- the conflation this whole audit exists to fix), or
(3) pass the ratio but have an unstable/near-zero PROJECTED force
(`unstable_magnitude` -- the case `mfv_gate_f0()` was already built to catch
for the F0 denominator, generalized here to every other site). Collapsing
(2) and (3) back into one "was low-SNR, now overridden" bucket would
re-introduce a NEW conflation of its own: a confidently_small point is safe
to average into a summary/F0/mean (it's a real, if quiet, force), but an
unstable_magnitude point is specifically NOT safe to use as a divisor or
direction seed (a near-zero denominator with high SNR on the ||dF|| vector
is a numerically fragile input) -- see `mfv_confidence_tier()`'s doc comment
in `muscle_force_vector.R`. Corpus-wide these are NOT rare edge cases: 42
`confidently_small` and 13 `unstable_magnitude` points exist outside the
F0-denominator context alone (item 4 above). So the 4 tiers are the minimum
needed to implement the proposal correctly, not a redundant elaboration of
a simpler 2-tier fix -- **kept, and is now the single shared implementation
(`mfv_confidence_tier()` / `MFV_CONFIDENCE_TIERS` / `MFV_CONFIDENCE_ALPHA`,
all in `muscle_force_vector.R`)** every other file below consumes.

### What changed (implements items 3(a)-3(e) from the proposal above)

`muscle_force_vector.R` (new shared machinery + 2 additive columns):
- Added `MFV_CONFIDENCE_TIERS` (ordered level vector), `mfv_confidence_tier()`
  (the function proposed in item 3, unchanged from the draft), and
  `MFV_CONFIDENCE_ALPHA` (4-level: confident=0.85, confidently_small=0.55,
  unstable_magnitude=0.35, unconfirmable=0.15) -- moved here from
  `plot_muscle_force_vector.R` (REPLACING that file's old 2-level version) so
  every consumer shares ONE definition; `muscle_force_vector.R` is sourced
  before every file that needs it in every entry point checked.
- `.mfv_finalize_step()`: added `baseline_force_noise_N` to the `ts` (force-
  vs-time) and `fv_l0` (sono L0-crossing) output tibbles -- both were missing
  the noise column entirely, which is why the summary/plot-level fixes below
  needed this additive change first (the tier function needs `force`, `snr`,
  AND `noise` together; `activation_snr` alone was already there, `noise`
  was not).
- `rank_isometric_steps_by_activation()`: `clearly_active` (ratio-only
  boolean) -> `confidence_tier` (4-level factor) in the printed PI A6 table.
- `superplot_fl_pooled.R`'s 3 extraction functions (isometric, isovelocity
  V=0, dynamic L0 bookend): added `baseline_force_noise_N` passthrough (was
  computed upstream but not carried into `pooled`).

`plot_muscle_force_vector.R` (item 3a, point-level alpha):
- `.mfv_long_methods()`: carries `baseline_force_noise_N` through the
  emp/geom pivot.
- `.mfv_confidence_factor()`: now a thin wrapper on `mfv_confidence_tier()`
  (was its own ad hoc 2-level `if_else`).
- `.build_vector_xy_summary()` (feeds `FL_isometric_uhatBoth.png` /
  `FV_isovelocity_uhatBoth.png`) and `build_summary_plot_FV_L0_vector()`
  (feeds `FVl0_isovelocity_uhatBoth.png`): confidence now computed from
  force+SNR+noise together (was SNR alone); subtitles report all 4 tier
  counts instead of one "n low-confidence" number.

`superplot_fl_pooled_snr_passing.R` (item 3b, hard-drop #1):
- Drop criterion: `activation_snr >= SNR_MIN` -> `tier_emp %in% c("confident",
  "confidently_small") | tier_geom %in% ...` (row survives if EITHER force
  estimate for that point qualifies, mirroring the existing design's single
  shared `activation_snr` across both method facets). Re-run against the
  real 3-fish corpus: **194/241 points now pass (was fewer under the pure
  ratio test) -- 21 of the 194 are kept ONLY via the magnitude floor** (would
  have been dropped as "noise" by the old ratio-only filter).

`export_snr_summary_figures.R` (item 3b, hard-drop #2):
- Added `tier_filter(df, force_cols)` alongside the existing `snr_filter()`.
  `tier_filter()` now feeds `iso_vec`/`isv_vec` (both force columns, OR
  logic), `fv_l0_all` (`force_at_L0_N`), and every `force_ts_vec_all` family
  (`muscle_force_vector_N`) -- i.e. every `_snrPass` figure that has an
  actual force magnitude to test. `snr_filter()` (unchanged, ratio-only) is
  KEPT deliberately for `uhat_tbl_all` (`uhatCompare_empVsGeom_snrPass.png`)
  ONLY -- that table compares u_hat DIRECTIONS (angles), which have no
  magnitude counterpart, and direction confidence is intentionally
  ratio-only by design (item 1's `.mfv_finalize_step()` finding). This is the
  one deliberate exception to "everything becomes tier-based," not an
  oversight.

`diag_isovelocity_hillcheck.R` (item 3c, the V=0-notch site):
- Summary-line filter: `confident <- activation_snr >= snr_min` ->
  `confidence_tier %in% c("confident", "confidently_small")`. Point alpha
  upgraded from 2-level TRUE/FALSE to the full 4-tier `MFV_CONFIDENCE_ALPHA`.
  Re-run against real data: **bass18 right's V=0 mean rose from 0.56 N (n=3)
  to 0.94 N (n=7)** once the stronger trial's confidently_small V=0 reps
  (SNR 2.08/2.55, but real magnitude) were correctly included -- confirming
  the audit's diagnosis that those reps are real, not noise. The notch is
  PARTIALLY resolved, not eliminated: 0.94 N is still below the eccentric
  plateau (1.65-2.33 N) and most concentric points (1.37-3.63 N) -- see the
  file's own updated FINDING comment and the regenerated `isovhillcheck.png`
  for the honest, non-overclaiming conclusion. This directly answers item 4's
  open question: **a magnitude-aware scheme changes bass18's V=0 points from
  "excluded as noise" to "included as real, but the curve is still not a
  clean monotonic Hill decline"** -- resolving the ratio-gating ARTIFACT
  without resolving the underlying scientific question (which needs the
  within-trial-only test the file already flags as not yet built).

`plot_fatigue_timeline.R` (item 3d, decided): restricted the loess trend
line to `confident + confidently_small` points (was: every point, unweighted,
regardless of tier) -- brought into line with `hillcheck`'s summary-line
policy so the two sites no longer make OPPOSITE choices (one used to
under-weight-by-omission, the other over-weight-by-inclusion). Point-level
alpha upgraded to the full 4-tier factor.

`superplot_fl_fv_tiers.R` (item 3a, second point-flag site): added 4-tier
point alpha to all 3 plots (`fltiers_1`, `fvtiers_1`, `fltiers_2`) -- this
file previously had NO confidence flagging of any kind, the inconsistency
item 2's table called out. Lines are left at their existing fixed alpha
(a per-segment tier isn't meaningful for a connect-the-dots line); only the
point layer is tiered.

### Not changed (explicitly, not by omission)

- `superplot_fl_pooled.R`/`superplot_fv_pooled.R`'s plotted per-trial/per-
  fish/per-group MEAN points (the "numerator" issue flagged in item 2's
  table) are left as-is: those are already-averaged bin means, not raw
  points, and the proposal's lettered items (3a-3e) did not include a
  decided fix for weighting a bin mean by tier -- flagged here as still open,
  not silently dropped from consideration.
- `.mfv_finalize_step()`'s `use_emp` (u_hat source selection) stays
  ratio-only, by design (item 1's finding, reaffirmed above for
  `uhat_tbl_all`).
- The magnitude multiplier is `k=1` everywhere (matching `mfv_gate_f0`'s
  existing default), not independently re-validated.

### Verification

All edited files parse clean (`Rscript -e 'parse(file=...)'`, zero errors).
Ran end-to-end against the real 3-fish corpus:
`superplot_fl_pooled_snr_passing.R`, `superplot_fl_fv_tiers.R` (which also
re-runs `superplot_fl_pooled.R`/`superplot_fv_pooled.R`), and
`diag_isovelocity_hillcheck.R` -- all completed without error and produced
the point counts and V=0 mean quoted above. Also ran
`export_snr_summary_figures.R` (bass16, default env vars, which re-runs the
FULL `run_fv_fl_power_pipeline.R` for that specimen) end-to-end without
error: `FL_isometric_uhatBoth_snrPass.png` 27/32, `FV_isovelocity_uhatBoth_
snrPass.png` 27/48, `FVl0_isovelocity_uhatBoth_snrPass.png` 3/13,
`uhatCompare_empVsGeom_snrPass.png` (ratio-only, unchanged) 54/64,
`forceTime_isometric_uhatBoth_snrPass.png` 18041/48000 rows,
`forceTime_isovelocity_uhatBoth_snrPass.png` 20313/48032 rows -- confirming
`tier_filter()`'s per-family force-column wiring (`iso_vec`/`isv_vec` OR
logic, `fv_l0_all`'s `force_at_L0_N`, and per-family `muscle_force_vector_N`
in `force_ts_vec_all`) all resolve correctly against live pipeline output;
the same run also exercised `plot_fatigue_timeline.R`'s revised trend-line
gate (`bass16_fatiguetimeline.png`, 60 L0 contractions) without error.

**UPDATE 2026-07-22 (same day): bass17/bass18 caught up.** Re-ran
`export_snr_summary_figures.R` for bass17 and bass18 (env-var pattern, full
per-specimen pipeline + `_snrPass` regeneration) -- both completed without
error. bass17: `FL..._snrPass` 5/32, `FV..._snrPass` 43/72, `FVl0..._snrPass`
4/11, `uhatCompare..._snrPass` 42/80, force-vs-time rows 2190/48000
(isometric) + 35992/72048 (isovelocity). bass18: `FL..._snrPass` 28/28 (ALL
pass -- consistent with item 4's finding that none of bass18's data is
unconfirmable), `FV..._snrPass` 35/44, `FVl0..._snrPass` 17/28,
`uhatCompare..._snrPass` 45/64, force-vs-time rows 24371/47628 (isometric) +
28899/54036 (isovelocity). All 3 specimens' per-specimen figures are now
current with the 4-tier scheme.

---

## Torque smoothing cutoff: does it need to vary by protocol/motion state? (2026-07-22)

PI, dissatisfied with FL/FV pattern strength, asked whether the raw 6-axis
(3 force + 3 torque) timeseries are filtered appropriately, and specifically
whether the smoothing cutoff needs to differ "per experiment type or angular
acceleration." Context: `load_bender_flat()` already applies a 9th-order
Butterworth `filtfilt` low-pass at LOAD time to all 6 raw channels (cutoff =
`cutoffmult` (8x, default) times the protocol's driving frequency, or 50 Hz
if that's not finite) -- so the signals are NOT literally unfiltered. The
open question is whether that (or the separate post-inertial-correction
smoothing choice `diag_torque_smoothing.R` exists to pick) is RIGHT for every
protocol, or whether isometric (static hold) and isovelocity (continuously
moving) need different treatment.

Extended `R/diag_torque_smoothing.R` (previously ISOMETRIC ONLY) to also
process isovelocity steps, tagged `protocol` through both output plots, and
re-ran it (bass16 only -- "a few examples," per PI's "let's run it on a few
examples"). See `FIGURES_README.md`'s `torquesmoothing` entry for the full
writeup; summary:

- **Isovelocity's active-window power spectrum is ~2-3 orders of magnitude
  ABOVE isometric's, across the ENTIRE 0-100 Hz band** (broadband elevation,
  not one narrow spike) -- isometric n=32 steps peaks near 18 Hz at a low
  absolute level; isovelocity n=48 steps peaks broadly near 10-18 Hz and
  stays elevated out past 75 Hz.
- The time-domain panel shows why: isovelocity torque oscillates for the
  WHOLE step (motion-coupled, not a brief transient), while isometric is
  quiet except right at stim onset.
- **This directly supports the PI's hypothesis for the "per experiment
  type" half of the question**: a single cutoff hand-picked from
  isometric-only data (this diagnostic's original scope, before today) is
  not obviously right for isovelocity -- either it's too tight (blurs real
  motion-coupled torque) or too loose (leaves isovelocity's much larger
  broadband content essentially unfiltered), and there was previously no
  isovelocity data in this diagnostic to even ask the question.
- The **"angular acceleration" half is NOT yet tested** -- that needs the
  commanded/encoder bend-angle signal's 2nd derivative correlated against
  the torque residual, per step (to check whether the elevated isovelocity
  power tracks acceleration events specifically, e.g. ramp start/stop, vs.
  being uniform throughout the moving window). Not built in this pass --
  flagged as the natural next diagnostic.

**No production filter changed** (`NOISE_FILTER_HZ`, `DISPLAY_SMOOTH_HZ`,
`load_bender_flat()`'s `cutoffmult`) -- this remains a read-only diagnostic
extension. A protocol-dependent (or motion-state-dependent) cutoff is now
evidenced as plausibly necessary, not decided; per `.cursorrules`, changing
a setting consumed pipeline-wide needs an explicit PI decision before
implementation.

### UPDATE (same day): confirmed on all 3 fish + quantified + example traces

PI: "run 1 [of the 5-step plan]. BUT make sure to include broadband-power gap
and examples of different filters for different samples." Three changes to
`R/diag_torque_smoothing.R`:

1. Parameterized `BASS_ID`/`SOURCE_DIR`/`OUTPUT_DIR` via `BENDER3_BASS_ID`
   env var (same pattern as `export_snr_summary_figures.R`) -- was
   hardcoded to bass16.
2. The broadband power gap is now a NUMBER, not an eyeballed log-scale plot:
   median isovelocity/isometric power ratio across every shared 0-100 Hz
   frequency bin. **bass16 = 9783x, bass17 = 3909x, bass18 = 3947x** --
   consistent order of magnitude (~3.5-4 decades) on all 3 fish, so this is
   NOT a bass16 idiosyncrasy. (Revises the earlier "2-3 orders of magnitude"
   visual estimate upward slightly, now that it's measured, not eyeballed.)
   Peak frequency also consistent: isometric 18-20 Hz, isovelocity 12 Hz,
   all 3 fish.
3. New Plot 3 (`diag_torque_smoothing_examples.png`): 2 isometric + 2
   isovelocity steps per specimen, chosen as the MEDIAN-total-variance step
   within their own protocol (representative, not cherry-picked), plotted
   large enough to actually see what each cutoff does. Result, same pattern
   on bass16 AND bass18 (visually spot-checked both): **isometric traces are
   visually IDENTICAL across every filter option** (raw through 8 Hz --the
   only thing any cutoff removes is a single-sample outlier spike, unrelated
   to the smoothing question). **Isovelocity traces are NOT** -- a large
   (~0.3-0.8 N*m), multi-cycle ringing right at motion onset that 8/15 Hz
   collapses into a smooth curve while 30 Hz/unsmoothed preserve in full.
   This is the direct visual answer to "why would the cutoff need to
   differ": for isometric the choice is inconsequential; for isovelocity it
   changes the trace's qualitative shape, and whether that ringing is a real
   mechanical response (worth preserving) or an artifact (worth removing) is
   exactly the question a single shared cutoff can't answer for both
   protocols simultaneously.

All 3 specimens re-ran end-to-end without error (`diag_torque_smoothing_{timedomain,spectrum,examples}.png`
regenerated in each `figs_bass{16,17,18}/`). Still READ-ONLY / PRE-DECISION
-- no production filter changed. Angular-acceleration correlation (plan step
2) not yet built.

---

## Dynamic-trial sono-strain gain traced to early-trial tissue preconditioning (2026-07-22)

PI flagged that `bass17_strainValidSonoEnc.png` shows a tight predicted-vs-
measured strain relationship for isometric/isovelocity but a "tall, narrow"
(excess-shortening) scatter for dynamic trials, and asked whether this is
real biology/slip or a signal-conditioning artifact.

**Ruled out** (each with its own diagnostic script, `R/diag_sono_*`):
calibration/crosstalk (isometric/isovelocity from the same sensors are
clean), the DS3's ~241-242 Hz staircase artifact (de-staircasing via
rolling-mean/linear-interp did not change the amplitude-ratio pattern),
mechanical vibration (no correlation of residual vs. peak angular
velocity/acceleration), true phase lag (~8 ms once a sign bug in a diagnostic
script -- NOT in the production `attach_sono_strain()` -- was fixed; the
production `strainValidSonoEnc.png` was never wrong), and muscle-depth
mis-estimation alone (sensitivity analysis: cannot produce the observed
~1.15x+ gain from a plausible 0-3 mm depth range).

**Root cause identified**: pooling all specimens/protocols
(`R/diag_sono_strain_validation_pooled.R`) showed the degradation is
concentrated almost entirely in dynamic, ACTIVE (right-stim) samples.
Disaggregating dynamic-active samples by trial order
(`R/diag_sono_progressive_slip.R`) revealed a clear decay: large excess
shortening (offset up to ~9.6 percentage points, r as low as -0.11) in each
specimen's earliest dynamic trials, decaying to near-zero offset and
r > 0.9 within that specimen's own session. Cross-validated: isometric/
isovelocity trials run EARLY in a session (bass18, trial 3) show the same
degradation; the same protocols run LATER (bass16, bass17, after prior
dynamic loading) do not. This is tissue preconditioning (a real, one-time-
per-specimen biomechanical settling), not a protocol-specific artifact.

**Cutoff decision (PROPOSED, NOT YET IMPLEMENTED -- awaiting PI go-ahead)**:
a single trial-count cutoff does not transfer across specimens -- bass16 and
bass18 stabilize by trial 5, bass17 not until trial 9 (using "|mean offset|
stays < 1.5 percentage points for every subsequent trial" as the stability
test; an r > 0.85 test agrees for bass16/bass18 but is more lenient for
bass17, stabilizing at trial 6 -- offset was preferred as the stricter/safer
criterion since downstream power-output work depends on absolute strain
accuracy, not just tracking shape). Proposed per-specimen lookup: bass16
exclude trials 1-4 (keep 5+), bass17 exclude 1-8 (keep 9+), bass18 exclude
1-2 (keep 5+). Effect: labels dynamic+active trials only as
early(preconditioning)/later(stable); does not touch isometric, isovelocity,
or passive samples; does not delete or overwrite any existing data --
intended as an opt-in join key (specimen + trial_num) for scripts that need
to exclude preconditioning trials (e.g., sono validation, and an upcoming
early-vs-later power-output comparison for comparable bending/stim).

### IMPLEMENTED (same day): hard cutoff + power check + cleaned pooled figure

PI approved the cutoff after confirming what "1.5 percentage-point offset"
physically means (a signed bias in %L0, not a correlation). Three new files:

- **`R/dynamic_trial_precondition.R`** (new, shared): the canonical
  `DYNAMIC_PRECONDITION_CUTOFF_TRIALNUM` lookup (bass16=5, bass17=9,
  bass18=5), `extract_bender_trial_num()`, and `classify_dynamic_precondition()`.
  Any script needing this exclusion should source this file rather than
  re-deriving the cutoff.
- **`R/diag_precondition_power_check.R`** (new): (a) re-plots the offset-vs-
  trial-order data directly from the already-vetted pooled sono CSV (not
  recomputed) with the hard cutoff drawn as a vline --
  `figs_diagnostic/dynamic_precondition_offset_vs_trialorder.png`; (b) computes
  per-cycle muscle power (`avg_power.W`, via the same `calc_muscle_torque()`/
  `summarize_muscle_cycles()` recipe `run_fv_fl_power_pipeline.R` uses) for
  every dynamic active cycle, tags each with trial_num/precondition/
  curvature.invm/freq.Hz, and saves `figs_diagnostic/dynamic_precondition_
  power_vs_{trialorder,curvature}.png` plus per-cycle/per-trial CSVs in
  `data_processed/`.
  **Finding (descriptive, not yet significant per-specimen -- small n):**
  the SAME early-trial elevation pattern seen in sono offset also appears in
  mechanical power (computed from encoder/torque, NOT sono, so this is an
  independent signal): all 3 specimens show numerically higher mean cycle
  power in early/preconditioning trials than in later/stable trials
  (bass16 -0.0034 W, bass17 -0.012 W, bass18 -0.051 W, all "later" minus
  "early", none individually significant at n=5-13 trials/specimen --
  `lm(avg_power.W ~ precondition + curvature + freq)` per specimen). Direction
  is consistent across all 3 specimens, suggesting preconditioning may affect
  more than just the sono-strain reading, but this needs a properly powered
  follow-up (more trials or a mixed-effects pooled model) before treating it
  as confirmed -- flagged as open, not claimed as resolved.
- **`R/plot_sono_strain_validation_pooled.R`** (extended, not overwritten):
  `.build_family_plot()` now accepts an optional `df`/`title_suffix` so the
  dynamic family can be re-plotted restricted to "later (stable)" trials only,
  alongside (not replacing) the existing all-trials dynamic panel. New output:
  `figs_summary/pooled_strainValidSonoEnc_dynamic_later.png` +
  `data_processed/sono_strain_validation_pooled_dynamic_later_samples.csv`.
  **Result: dynamic-active pointwise r rises from 0.447 (all trials) to 0.903
  (later-only)**, essentially recovering isometric/isovelocity-level
  tightness -- direct confirmation that the original dynamic-trial scatter
  (`bass17_strainValidSonoEnc.png`, the figure that started this whole
  investigation) is explained by early-trial preconditioning, not a
  protocol-wide sono/encoder artifact.

All 3 new/edited scripts ran end-to-end without error against the full
bass16/17/18 corpus.

### FOLLOW-UP (same day): combined all-protocol figure + sono oversampling fix

PI-requested: one multi-panel figure in the style of the original per-
specimen `<specimen>_strainValidSonoEnc.png` (facet by protocol_family),
extended to pool all 3 specimens, split by active/passive, and restricted to
"later (stable)" trials -- plus a request to confirm sono duplicate/
oversampled points aren't being double-counted (~241 Hz true DS3 rate vs.
1 kHz AI digitization).

**Sono oversampling check -- important correction to the initial approach.**
First attempt de-duplicated by exact-value repeats (assuming a zero-order-
hold staircase: N identical consecutive mm values per true DS3 update).
Checked empirically before trusting it (bass16 trial 1): only ~2% of
consecutive-sample diffs were exactly zero -- the raw calibrated sono signal
changes by small, continuously-varying amounts on nearly EVERY 1 kHz sample,
not in flat held steps. Exact-duplicate dedup was therefore the wrong tool
(it would have kept ~98% of the oversampled rows, doing almost nothing).
Most likely explanation: the DS3's analog output stage smooths/reconstructs
its own ~241-247 Hz updates before the 1 kHz ADC digitizes it, so there is
no literal staircase to detect by value equality even though the true
independent-information rate is still only ~241-247 Hz.
**Fix: decimate, don't dedupe** -- `plot_sono_strain_validation_pooled.R`'s
`.process_trial()` now keeps 1 sample in every `round(1000 / daq_sono_
internal_sample_rate_hertz)` (typically 1-in-4), applied to the OUTPUT rows
only (after stim-state labeling on the full-rate td, so `.detect_stim_
events()` cannot miss a single-sample stim pulse to a decimated-away row).
This reduced pooled n from ~3.6M to ~0.9M rows (all-trials) -- roughly the
expected 4x, confirming the fix is doing real work this time. Applies to ALL
4 existing pooled figures (`pooled_strainValidSonoEnc_{isometric,isovelocity,
dynamic,frequency_sweep}.png` and `_dynamic_later.png`), not just the new one
-- regenerated all of them; r/RMSE values did not materially change (as
expected -- decimating oversampled, non-independent points shouldn't bias a
correlation, just its inflated n), confirming this is a data-hygiene fix, not
a result-changing one.

**Protocol-agnostic precondition classification (new):** added
`classify_session_precondition()` (`dynamic_trial_precondition.R`) -- the
SAME specimen+trial_num cutoff as `classify_dynamic_precondition()`, applied
to every protocol family, not just dynamic. Justified because trial_num is a
per-specimen SESSION-chronology counter spanning every protocol type, and the
cutoff represents a one-time tissue-settling effect over the session, not a
dynamic-specific artifact (supported by the earlier cross-protocol
observation that early-session isometric/isovelocity trials show the same
degradation). `pooled` now carries `trial_num`/`precondition` columns
unconditionally.

**New figure:** `figs_summary/pooled_strainValidSonoEnc_allProtocols_later.png`
-- `facet_grid(protocol_family ~ active_passive)`, colored by specimen,
"later (stable)" trials only across all families. All 8 cells:
r = 0.905-0.989 (isometric passive/active 0.988/0.989, isovelocity
0.972/0.967, dynamic 0.973/0.905, frequency_sweep passive-only 0.965) --
i.e. once preconditioning trials are excluded and sono oversampling is
corrected, EVERY protocol family and stim state now sits in the same tight
0.90-0.99 r range, closing out the original investigation
(`bass17_strainValidSonoEnc.png`'s dynamic-panel scatter).
Also saved: `data_processed/sono_strain_validation_pooled_allprotocols_
later_samples.csv`.

### FOLLOW-UP (2026-07-23): step/cycle-level activity split (residual-force check)

PI hypothesis: the per-SAMPLE windowed `active_passive` split (burst window
+ relaxation tail = "active", everything else = "passive") could still leave
residual post-stimulation force/strain sitting in the "passive" bucket --
i.e. samples that occur later in the SAME designed step/cycle as a real stim
burst, after the window closes but before full mechanical recovery, get
counted as if that unit had never been stimulated at all.

Added a coarser, unit-level split as an alternate figure: for each designed
experimental unit (one isometric/isovelocity STEP via `step_number`, or one
dynamic/frequency_sweep bending CYCLE via `cycle_index` -- `step_number`
doesn't exist as a column for single_finite trials, `cycle_index` collapses
to 1 value for segmented ones, so the two are genuinely disjoint, not
redundant), flag whether stim fired ANYWHERE in that whole unit (any side),
computed from the raw per-sample `stim` column on the FULL-rate signal
(`.process_trial()`'s new `step_has_stim`, via `stats::ave(..., FUN = max)`)
-- NOT the windowed label, and computed BEFORE the sono decimation step so a
short burst can't be lost to either. Label: "purely passive (no stim in
step/cycle)" vs. "contains stimulation".

`.build_allprotocols_plot()` generalized to take `col_var`/`col_caption` so
both the sample-level and unit-level figures share one implementation.

New file: `figs_summary/pooled_strainValidSonoEnc_allProtocols_later_
stepActivity.png`.

**Result:** isometric's "purely passive" cell is now EMPTY (n=0) -- every
isometric step in the corpus has SOME stim scheduled somewhere in it (no
truly never-activated FL steps), which the old per-sample split couldn't
show since it only ever excluded/included individual windowed samples, never
revealed that the whole-step population itself has no null case. Dynamic's
"contains stimulation" cell grew from n=8,577 (old, sample-windowed
"active (right stim)") to n=52,851 (new, whole-cycle) AND got tighter
(r=0.905 -> 0.932, RMSE=0.96 -> 0.82) -- i.e. once every sample belonging to
a stimulated cycle is included (not just the narrow burst+relaxation
window), the fit doesn't degrade, it improves, meaning the earlier windowed
"active" bucket was actually the noisier/smaller slice, not the residual-
contaminated one. Isovelocity and frequency_sweep both show real,
comparably-sized populations on both sides of this split (r=0.97-0.98
throughout), i.e. no similar all-or-nothing pattern outside isometric.

FOLLOW-UP (2026-07-24): normalized power vs. sono-strain offset (direct
check on "is high power just the slipping artifact?")

PI concern: `dynamic_precondition_offset_vs_trialorder.png` and
`dynamic_precondition_power_vs_trialorder.png` (from
`diag_precondition_power_check.R`) both show the SAME early-trial pattern
side by side, implying (but not directly showing) that the trials with the
highest measured power are exactly the trials where the muscle appears to
be slipping (highest sono-strain offset). New script
`R/diag_precondition_power_vs_offset.R` tests this directly by correlating
per-trial normalized power against per-trial offset (joins the two CSVs
`diag_precondition_power_check.R` already writes, at the TRIAL level --
not per-cycle, to avoid reconciling two different cycle-index conventions
between the sono and power pipelines).

**Result: confirmed, and it is not a between-specimen confound.**
- Pooled (n=30 dynamic trials, all 3 specimens): trial-mean cycle-averaged
  power (W/kg) vs. trial-mean sono offset (%L0), r=0.727, p=5.3e-6.
  Trial-max instantaneous peak power (W/kg) vs. same offset, r=0.630,
  p=1.9e-4.
- WITHIN each specimen separately (controls for one fish simply producing
  more power AND happening to have bigger early-trial offsets by
  coincidence): bass16 r=0.68 (n=13), bass17 r=0.55 (n=12), bass18 r=0.98
  (n=5, most extreme). All positive, all substantial -- the relationship
  holds fish-by-fish, not just in the pooled cloud.
- Visually (`dynamic_precondition_meanpower_vs_offset.png`): every "later
  (stable)" trial (circle) sits at offset~0, power~0. Every trial with
  power >5 W/kg is an "early (preconditioning)" trial (triangle), and the
  two highest-power trials in the whole dataset (bass17 ~27 W/kg, bass18
  ~27 W/kg) are also among the highest-offset trials (7-8 pct-points).

**Interpretation:** this reinforces (does not merely repeat) the existing
early/later cutoff -- it shows the artifact is not confined to the strain
channel. The apparently high power numbers from early trials are not an
independent positive finding; they are very likely inflated by the same
excess/spurious shortening that inflates the sono-strain offset (more
measured length change per cycle -> more computed work/power for the same
commanded motion and stim). This is an ADDITIONAL reason (beyond sono
validation) to exclude "early (preconditioning)" trials from any power
analysis, not just from strain-validation analysis.

New outputs:
- `figs_diagnostic/dynamic_precondition_meanpower_vs_offset.png`
- `figs_diagnostic/dynamic_precondition_maxpower_vs_offset.png`
- `figs_diagnostic/dynamic_precondition_power_vs_offset_byspecimen.png`
- `data_processed/dynamic_precondition_power_vs_offset_by_trial.csv`
- `data_processed/dynamic_precondition_power_vs_offset_within_specimen_cor.csv`

FOLLOW-UP (2026-07-24, same day): peak-power version + isovelocity extension

PI asks: (1) "perform a similar analysis using peak power instead of cycle
average" and (2) "perform a peak/avg power analysis on the other protocol
types as well."

(1) Added `dynamic_precondition_maxpower_vs_offset_byspecimen.png` to
`diag_precondition_power_vs_offset.R` -- the peak-power counterpart of the
existing mean-power by-specimen facet plot (Plot 3). Trial-max
`peak_power.Wkg` vs. offset, faceted by specimen, one within-specimen r per
panel. Pooled: r=0.630 (already known from Plot 2); already-reported
within-specimen r_max values (bass16=0.74, bass17=0.63, bass18=0.99) are now
visualized directly per specimen, not just tabulated.

(2) Which OTHER protocol types can actually support a "peak/avg power"
analysis -- decided by physical definition, not code convenience:
- isometric: the motor does NOT move during the activation window (that is
  the literal definition of "isometric"). Mechanical power there is ~zero
  BY DESIGN, not by measurement -- there is no meaningful power signal to
  correlate against offset. FORCE (muscle_force_Nm, already computed by
  analyze_isometric()) is the correct isometric analog of "how hard is the
  muscle working" -- flagged to PI as a separate possible follow-up, not
  computed here (different physical quantity, not power).
- frequency_sweep: passive-only in every category grouping seen throughout
  this whole investigation (active_mask always FALSE) -- no stimulation is
  ever delivered, so calc_muscle_torque()'s cycletype=="act" filter (the
  existing per-cycle power framework) returns zero rows. No muscle-driven
  power exists to analyze.
- isovelocity: DOES have genuine per-step mechanical power, and the
  codebase already has a VETTED torque -> Watts conversion (muscle_geometry
  .R::add_specific_properties_to_fit(): "real muscle power at that
  operating point is peak_power_composite * dclamp_m/(100*r_m)", where the
  composite is F0(muscle_force_Nm, a TORQUE despite the "_force" name) *
  Vmax(shortening_strain_pct, %/s for isovelocity per 03_analyze.R's
  build_segmented_step_summary() comment) -- torque x angular_velocity =
  power). New script `R/diag_precondition_power_vs_offset_isovelocity.R`
  applies this IDENTICAL formula PER STEP (each step's own measured
  muscle_force_Nm x shortening_strain_pct, not just the Hill fit's single
  extrapolated Vmax) -- no new physics, a finer-grained application of an
  already-vetted conversion. Restricted to muscle_side=="right" (the
  sono-instrumented muscle), matching the dynamic analysis's side
  convention. Required two additional sources not already loaded by this
  script's pattern: `plot_force_vs_time.R` (.smooth_trace_display_only(),
  needed by build_segmented_step_summary()'s peak-window force extraction)
  and `fit_fv_fl.R` (fit_force_velocity_curve(), needed by
  analyze_isovelocity() even though this script never uses the resulting
  Hill fit) -- both were silent dependency gaps until the script was
  actually run, now documented as required sources in the file header.

**Result:** isovelocity offsets are much smaller in magnitude than dynamic's
(mostly -0.2% to +0.5%, i.e. near the sono noise floor) except one bass18
trial at +2.4%. Trial-mean power (W/kg) vs. offset: pooled r=0.919,
p=6.3e-5, n=11. Trial-max |power| vs. offset: r=0.820, p=0.002, n=11. The
single high-offset bass18 trial (+2.4%, ~90 W/kg mean / 603 W/kg peak) is a
high-leverage point, but the remaining 10 trials still trend the same
direction (e.g. bass17's three trials, all near-zero/slightly-negative
offset, all show negative mean power -- net energy absorption, physically
consistent with near-zero-offset "no slip" trials). Same qualitative
conclusion as dynamic: higher offset tracks higher measured power, now
confirmed in a SECOND, independently-computed protocol family using a
different power-extraction method (per-step torque x commanded angular
velocity, vs. dynamic's calc_muscle_torque()/summarize_muscle_cycles()
work-loop integration) -- this is not an artifact specific to one power
pipeline.

New outputs:
- `figs_diagnostic/dynamic_precondition_maxpower_vs_offset_byspecimen.png`
- `figs_diagnostic/isovelocity_meanpower_vs_offset.png`
- `figs_diagnostic/isovelocity_maxpower_vs_offset.png`
- `data_processed/isovelocity_power_percycle.csv` (per-step, all sides)
- `data_processed/isovelocity_power_vs_offset_by_trial.csv`

FOLLOW-UP (2026-07-24, same day): isometric using TENSION (not power) + offset
magnitude sanity check against commanded amplitude

PI: "Please run analysis for isometric using tension" -- following directly
from the isometric-has-no-power-by-design reasoning above. New script
`R/diag_precondition_tension_vs_offset_isometric.R` computes per-step
SPECIFIC TENSION (N/cm^2) via the exact same lever-arm conversion
muscle_geometry.R::add_specific_properties_to_fit() already applies to a
fit's single F0 anchor point (Force0_N = muscle_force_Nm / r_m;
specific_tension_Ncm2 = Force0_N / csa_muscle_cm2), applied PER STEP.
Right-muscle steps only, same convention as the power scripts. Required the
same `plot_force_vs_time.R` source-dependency fix found in the isovelocity
script.

**Result:** only 5 isometric trials exist in the WHOLE corpus (bass16=2,
bass17=1, bass18=2), so this is a very small-n check. Mean tension vs.
offset: r=0.323, p=0.596 (not significant). Max tension vs. offset:
r=-0.408, p=0.495 (not significant, opposite sign). Neither correlation is
meaningful with n=5 -- flagged in the plot subtitle itself ("interpret
cautiously"), and the script explicitly suppresses the fitted trend line
whenever n<4. Notably, isometric offsets themselves are ALL small (0.57-
1.24 pct-points, one trial slightly negative) -- nothing like dynamic's
2-10 pct-point early-trial offsets -- consistent with isometric never
showing the same large-magnitude preconditioning effect dynamic does (the
motor isn't moving, so there's much less opportunity for progressive
mechanical settling/slip to manifest as a strain-offset signal).

New outputs:
- `figs_diagnostic/isometric_meantension_vs_offset.png`
- `figs_diagnostic/isometric_maxtension_vs_offset.png`
- `data_processed/isometric_tension_percycle.csv` (per-step, all sides)
- `data_processed/isometric_tension_vs_offset_by_trial.csv`

PI: "2-10% offset seems high, especially if commanded strain was something
like 5%." Checked directly against the vetted pooled sono CSV (not the
older/incomplete amp_ratio_all_dynamic_trials.csv, which only had right-stim
coverage for 5/30 dynamic trials -- likely predates full bass17/18
collection): commanded (encoder-predicted) peak-to-peak strain amplitude
for the SAME active-right dynamic trials used in the offset analysis is
remarkably uniform, median 5.04% (range 4.8-6.6%) -- confirming the PI's
5% estimate almost exactly. Expressing each trial's OFFSET as a percentage
OF ITS OWN commanded p2p amplitude:
  - "early (preconditioning)" trials (n=14): median 92%, mean |ratio| 104%,
    max 200% of the commanded p2p amplitude -- i.e. in the worst early
    trials, the MEAN-level bias alone is comparable to or LARGER than the
    entire commanded oscillation's peak-to-peak swing.
  - "later (stable)" trials (n=16): median 9.9%, mean |ratio| 12.1%, max
    25.4% -- an order of magnitude smaller, but still non-trivial (roughly
    a tenth of the commanded swing).
This confirms the PI's intuition was correct and quantifies it: the early-
trial offsets are not a subtle measurement bias, they are on the same
scale as the entire commanded motion, consistent with real gross mechanical
settling/slip rather than a small calibration artifact. The "later
(stable)" residual offset, while much smaller, is not zero and should be
kept in mind as an ~10%-of-amplitude noise floor even after excluding early
trials.

FOLLOW-UP (2026-07-24, same day): does passive fidelity ALSO improve
early->later, or is the effect active-stim-specific?

PI: "Is the sono fidelity in passive trials ALSO better than earlier ones?"
New script `R/diag_precondition_passive_sono_fidelity.R` -- a pure re-slice
of the already-vetted pooled sono CSV by precondition status, no new
processing, giving the PASSIVE (no/left stim) counterpart of the ACTIVE
(right-stim) numbers `diag_precondition_power_check.R` already reported.

**Result: yes, pointwise fidelity improves, but the MEAN-LEVEL BIAS (the
thing the whole cutoff was built around) does NOT appear in passive samples
at all, even early.**

| samples | precondition | n | r | RMSE (%) | offset (pct-pts) |
|---|---|---|---|---|---|
| passive (no/left stim) | early | 47,773 | 0.610 | 2.25 | **+0.017** |
| passive (no/left stim) | later | 160,462 | 0.973 | 1.62 | +0.086 |
| active (right stim) | early | 4,235 | 0.293 | 6.34 | **+4.92** |
| active (right stim) | later | 8,577 | 0.905 | 0.957 | +0.426 |

Passive r rises 0.61->0.97 (less scatter/noise in later trials, same
direction as active) but the passive OFFSET stays ~0 in BOTH early and
later trials (+0.02 vs +0.09 pct-points -- both negligible). Active's early
offset is +4.92 pct-points -- two orders of magnitude larger than passive's
early offset. Per-specimen passive breakdown: bass18's early passive r is
the worst of all (0.246), consistent with bass18 also showing the most
extreme active-offset/power anomalies elsewhere in this investigation.

**Interpretation:** the dramatic MEAN-LEVEL bias (the "slip") is
specifically an ACTIVE-STIMULATION phenomenon -- it only appears when the
muscle is actually pulling on the grips, not just from general early-
session mechanical unsettledness. Passive samples DO get noisier early
(more scatter, lower r) -- consistent with a generally less-settled
early-session mechanical state -- but with no directional bias, because
there's no muscle tension to drive slip in one direction. This is
independent confirmation that the offset cutoff is capturing a real,
stimulation-linked mechanical effect, not just generic early-session sensor
noise that would show up equally regardless of stimulation state.

New outputs:
- `figs_diagnostic/dynamic_precondition_passive_sonoValidation_earlyVsLater.png`
- `data_processed/dynamic_precondition_sono_fidelity_active_vs_passive.csv`
- `data_processed/dynamic_precondition_sono_fidelity_passive_by_specimen.csv`

CORRECTION (same day): "passive" must be no stim on EITHER side, for the
WHOLE cycle -- not the per-sample "no/left stim" window above.

PI clarification: "For comparing early versus late trials, passive means
no left OR right stim at all during the cycle." The `active_passive`
bucket used just above is a per-SAMPLE windowed split where only right-stim
counts as "active" (since sono only instruments the right muscle) -- so its
"passive" bucket still includes LEFT-stim cycles, which have real
muscle-driven bending (just not from the muscle sono is watching), not a
clean no-activation baseline. Re-ran with the correct CYCLE-LEVEL
definition: `step_activity == "purely passive (no stim in step/cycle)"`
(stim never fired on either side, anywhere in that whole bending cycle --
already computed in the pooled CSV, `plot_sono_strain_validation_pooled.R`).

| samples | precondition | n | r | RMSE (%) | offset (pct-pts) |
|---|---|---|---|---|---|
| STRICT passive (no stim either side, whole cycle) | early | 15,988 | 0.915 | 1.63 | -0.447 |
| STRICT passive (no stim either side, whole cycle) | later | 116,188 | 0.975 | 1.84 | +0.158 |
| active (right stim) | early | 4,235 | 0.293 | 6.34 | +4.92 |
| active (right stim) | later | 8,577 | 0.905 | 0.957 | +0.426 |

With the correct definition, STRICT-passive fidelity is already quite good
even in EARLY trials (r=0.915, close to later's 0.975) -- much better than
the old loose passive number (r=0.610) suggested, because that number was
partly measuring left-stim-driven bending dynamics, not sensor noise. The
offset stays small in both groups AND FLIPS SIGN (-0.447 early ->
+0.158 later), consistent with residual scatter around zero rather than a
directional bias. This makes the contrast with active's early offset
(+4.92, one direction, ~30x larger) sharper, not weaker: the massive,
one-directional early bias is specific to right-side stimulation, not a
generic property of early-session mechanics. Conclusion from the earlier
(superseded) loose-passive comparison is unchanged, just now on firmer
footing. `diag_precondition_passive_sono_fidelity.R` now uses the strict
definition as primary and reports the old loose number only as a labeled
comparison row.

---

FOLLOW-UP (2026-07-24, same day): can a fixed calibration factor recover
the "true" strain from the active-trial offset?

PI: "Based on the offset, can we just add a calibration factor to truly
get the desired strain? If so, does the active data suggest such a factor
would be consistent across conditions or variable (and therefore kinda
useless)?"

New script `R/diag_precondition_calibration_gain_check.R` -- fits
`strain_sono_pct ~ strain_pred_encoder_right_pct` PER TRIAL (dynamic,
active/right-stim only), extracting slope (a hypothetical multiplicative
GAIN correction) and intercept (a hypothetical additive OFFSET correction)
separately, then checks whether either is stable enough across trials to
serve as a single fixed calibration factor.

**Result: no single factor works. Early trials aren't a fixed error at
all (the fit itself is barely a real relationship), and even the
"later/stable" residual gain differs by specimen.**

| trial group | n trials | median gain (slope) | sd | range | median offset (intercept) | sd |
|---|---|---|---|---|---|---|
| early (preconditioning) | 14 | 0.66 | 0.61 | -0.24 to 1.66 | 6.09 | ~4-5 |
| later (stable) | 16 | 1.08 | 0.13 | 0.93 to 1.28 | 0.40 | ~0.6 |

Early trials: gain swings from -0.24 to 1.66 trial-to-trial within the SAME
specimen/session (bass16 trial 1: slope=0.07, r=0.04 -- essentially no
linear relationship at all; bass17 trial 4: slope=1.66; bass18 trials 1-2:
NEGATIVE slope, r<0). A "calibration factor" presupposes a stable linear
relationship to correct -- there isn't one in early trials. This is the
regression-level signature of the same trial-order decay documented above
(early trials trend from near-zero/negative gain toward ~1 as trial number
increases), not a fixed sensor/geometry error.

Later (stable) trials: gain is far more self-consistent WITHIN a specimen
but DIFFERS BETWEEN specimens:

| specimen | n trials | median gain | sd gain | median offset | sd offset |
|---|---|---|---|---|---|
| bass16 | 9 | 1.03 | 0.124 | 0.38 | 0.76 |
| bass17 | 4 | 1.27 | 0.020 | 0.34 | 0.27 |
| bass18 | 3 | 1.08 | 0.050 | 0.49 | 0.16 |

bass17's later-trial gain is remarkably tight (1.24-1.28, sd=0.02) --
consistent with a small, real, specimen-specific geometric mis-estimate
(e.g. muscle depth/body width measured slightly off for that individual)
that a PER-SPECIMEN multiplicative correction could plausibly remove.
bass18 is also fairly tight (sd=0.05). bass16 is looser (sd=0.12, includes
trials 10-12 at 1.20-1.24 vs trials 5-8 at 0.93-0.99) -- a specimen-wide
factor would still leave meaningful trial-to-trial residual there.

**Bottom line for the PI's question:**
1. NO universal calibration factor (one number for the whole dataset) is
   defensible -- gain isn't even a fixed number for the SAME specimen in
   early trials, and later-trial gain differs specimen-to-specimen
   (1.03-1.27).
2. A PER-SPECIMEN multiplicative gain applied ONLY to later/stable trials
   is plausible for bass17 (and reasonably so for bass18), but would be a
   band-aid on top of an underlying geometric uncertainty already explored
   (muscle-depth/body-width sensitivity, ruled insufficient ALONE to
   explain the full early-trial gain but plausible for this smaller
   residual). It should not be applied to early trials, where slippage
   dominates and no fixed factor exists.
3. Simplest, already-implemented alternative that requires no assumed
   correction model: keep using the existing hard exclusion cutoff
   (`dynamic_trial_precondition.R`) for anything requiring accurate
   absolute strain, rather than trying to correct the excluded trials.

New outputs:
- `figs_diagnostic/dynamic_precondition_calibration_gain_vs_trialorder.png`
- `figs_diagnostic/dynamic_precondition_calibration_offset_vs_trialorder.png`
- `data_processed/dynamic_precondition_calibration_gain_by_trial.csv`
- `data_processed/dynamic_precondition_calibration_gain_consistency_by_precondition.csv`
- `data_processed/dynamic_precondition_calibration_gain_consistency_later_by_specimen.csv`
