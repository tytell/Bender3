# Muscle force vector physics (bending + LOA projection)

## Coordinate frame (ATI Mini40, sensor frame)
X = longitudinal, + = anterior. Y = mediolateral, + = right (flips if the fish
faces the opposite way — check `daq_specimen_lateral_index_on_positive_motor_side`
/ `daq_positive_motor_direction`). Z = dorsoventral, + = ventral.

## Bending geometry (two different moment arms on the same axis)
Bending is lateral (in the XY plane), so joint angle is read out as **zTorque**.
Two different radii convert bending-related quantities on this one axis:
- **Passive body bending arm:** half the local body width,
  `r_body = measurement_specimen_local_body_width_millimeter / 2`.
- **Muscle's z-torque arm:** `r_m = r_body - muscle_depth`, using
  `measurement_target_muscle_depth_millimeter` (falls back to an assumed 1 mm
  when unset/zero).

Muscle tension also acts at a **dorsoventral offset** from the sensor origin,
`d = measurement_clamp_offset_vertical_millimeter` (89 mm in the current
files), which converts longitudinal muscle force into **yTorque**.

## Muscle force vector model
Muscle tension is one scalar force *F* along a unit line-of-action
û = (ux, uy, uz) — mostly longitudinal, but rotating in the XY plane with bend
angle θ (the posterior half moves with the motor while the anterior half stays
fixed), with a possible small Z tilt if the fibers aren't perfectly
horizontal. So F·û gives force components Fx = F·ux, Fy = F·uy, Fz = F·uz —
all three of which the sensor can register directly as forces — and, via
torque = r × F at the muscle's attachment point r = (0, r_m, d):

```
tau_y =  d   * Fx      (dorsoventral arm x longitudinal force)
tau_z = -r_m * Fx      (mediolateral arm x longitudinal force)
tau_x =  r_m * Fz - d * Fy      (cross-axis leakage)
```

tau_y and tau_z each give an **independent estimate of F*ux** (F, when
û ~= X): `F*ux = tau_y / d = -tau_z / r_m`. Agreement between these two is a
built-in consistency check on the model. tau_x, Fy, and Fz capture off-axis
leakage from û's rotation/tilt and any mounting imperfection — diagnostic,
not part of the primary force estimate.

## Baseline logic
All of the above is applied to **active-minus-passive** wrench components
(componentwise baseline subtraction: pre-stim window for isometric,
angle-matched passive for isovelocity) before solving for F — a raw or
uncorrected wrench mixes passive bending stiffness into every channel.

## Troubleshooting log: resolution, SNR, calibration/gain, per-trial variability (2026-07-18)

**Activation SNR and what tension it actually requires.** `activation_snr =
||dF_force|| / noise`, where `dF_force` is the active-minus-passive
force-channel (Fx,Fy,Fz) vector for a step and `noise` is that step's
pre-stim force-channel SD. Because the muscle's force vector IS `F*u_hat`
(unit vector), `||dF_force||` is already directly comparable to muscle
tension in Newtons — no moment-arm conversion needed to answer "how much
tension clears SNR >= 3". Measured per-fish baseline noise (isometric
steps): bass16 31 mN, bass17 69 mN, bass18 24 mN (medians). The
corresponding tension floor (3x noise) is 0.07-0.21 N, or 0.09-0.32 N/cm^2
once normalized by each fish's estimated muscle CSA — 50-350x below the
typical vertebrate specific-tension range (~15-30 N/cm^2). Conclusion: real
muscle activation should clear SNR >= 3 by a wide margin; a step failing to
clear it is evidence of weak/absent recruitment or contamination (e.g.
drift), not an unreasonably strict threshold.

**Calibration file swap (FT56491 vs FT56492) does not change SNR.**
Re-applied FT56492.cal's (SI-20-1, "high resolution") matrix to the same raw
voltage in place of the embedded FT56491.cal (SI-40-2) matrix, across all 92
isometric steps (bass16/17/18): median SNR shifted -0.19%, per-fish pass
rates (SNR >= 3.0) were identical to the decimal. Expected: this is a
post-hoc linear rescale of already-digitized voltage (coefficients differ by
only ~1-9%); SNR is scale-invariant under a uniform per-channel rescale.
FT56492 is also a **different sensor serial number**, not a second
calibration mode of the unit mounted on the rig, so the swap isn't even a
valid physical substitute regardless. Do not chase calibration-file changes
as an SNR fix.

**ADC gain (±10V range) is also not the bottleneck.** All 8 AI channels
share one NI USB-6361 ADC range, `min_val=-10.0, max_val=10.0`
(`bender_functions.py` ~line 2607) — the real user-adjustable "gain" knob
(the ATI amplifier's own gain is fixed by hardware, unrelated to the `.cal`
file). The theoretical quantization-only noise floor from the 16-bit ADC and
the calibration matrix's channel sensitivities is ~1.9 mN combined
(Fx,Fy,Fz), only 3-8% of the observed baseline noise (24-69 mN) across
bass16/17/18. Halving the range to +/-5V (doubling resolution) would only
cut that 3-8% component in half, predicting <0.3% change in total noise. The
observed noise floor is dominated by something other than ADC bit depth.

**Stim_monitor channel: tested and not implicated.** Checked whether
`stim_monitor` (ai7) wiring leaks noise into the F/T channels (DAQ mux
crosstalk or shared electrical path) versus general EMI from the Grass S88's
proximity. Across 92k pre-stim baseline samples: `stim_monitor` itself is
quiescent (SD ~8 mV, tightly bounded, no signal to leak), its correlation
with Fx/Fy/Fz noise at the same instants is near-zero (r = 0.07, 0.01,
-0.08), and its own noise is nearly identical across all three fish
(7.8-8.2 mV) even though force-channel noise varies ~2x between fish
(bass17 highest). An aggressor with no signal, no correlation, and no
cross-fish variability of its own is not a plausible driver.

**Per-trial/per-fish variability is the real open question.** bass17's ~2x
higher baseline noise varies **between recording sessions**, not from any
fixed cause (calibration, gain, stim-monitor wiring all ruled out or shown
negligible above). That points toward session-specific factors —
mounting/clamping, cable routing, grounding differences at setup — rather
than anything fixable in software. Combined with the earlier finding that
bass17's isometric FL curve is dominated by passive creep (drift), the
working hypothesis is a less stable physical mount that session, not a
sensor/DAQ configuration problem. Not yet tested: whether noise correlates
with mounting hardware, cable reconnection, or nearby equipment power state
per session (needs session-log cross-referencing, not just h5 data).
