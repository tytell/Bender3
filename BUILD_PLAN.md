# BUILD_PLAN — Faithful scup-mechanics port + dynamic QC layer for Bender3

## Authority & scope
- **Deliverable of this plan:** completing a self-contained, native, flat-schema port of the
  validated scup muscle-mechanics / QC core into Bender3, a trial-type dispatcher, an
  auditable dynamic-trial QC layer, and a scup-vs-Bender3 numerical equivalence check.
- **Reuse strategy:** PORT, not `source()`. scup is a validated SPEC. Bender3 re-implements
  each function as its own code against the flat canonical schema
  (`context_jlab_cg_h5schema.md`). No scup runtime dependency, no fixed cross-project paths.
  The single exception is the **validation harness**, which `source()`s scup deliberately as a
  test oracle (never on the Bender3 run path).
- **First milestone:** solid QC for **dynamic** (`single_finite`) trials only. Segmented
  (`segmented_finite`) trials are detected + routed but analysis is DEFERRED (stubs exist).
- **Faithfulness:** identical math to scup, only the schema changes. Any place a faithful port
  is impossible under the flat schema is called out as an OPEN QUESTION (§7).

---

## 1. Repo inventory & gaps

### 1.1 What Bender3 already has
- [`R/00_load_bender_flat.R`](R/00_load_bender_flat.R): flat-schema loader. Reads `metadata/`
  (attrs + index datasets) and `timeseries/` (flat for `single_finite`, `step_NNN/` for
  `segmented_finite`); attaches per-cycle design (`index_cycle_*`), calibrated torque from
  `derived/forcetorque_calibrated`, bias subtraction, Butterworth LP. Produces the
  scup-compatible column names (`t.s`, `t.norm`, `angle.deg`, `enc.deg`, `stim`, `cycle`,
  `halfcycle`, `freq.Hz`, `amp.deg`, `duty`, `phase`, `curve.invm`, `curvature.invm`,
  `curverate.invms`, `dclamp.m`, `fishlength.m`, `fishmass.g`, `fishcode`, `stimulus_type`,
  `xtorque0.Nm`/`xtorque.Nm`, `is_active_by_cycle`).
- [`R/01_calibrate.R`](R/01_calibrate.R): applies embedded `calibration_forcetorque_matrix`
  to raw `aidata` → `derived/forcetorque_calibrated`; sono interp. (Bender3's replacement for
  the legacy pre-computed `/Calibrated/*`.)
- [`R/02_deconvolve.R`](R/02_deconvolve.R): inertial deconvolution
  `tau_corrected = tau - (I_total*alpha + bias)`.
- [`R/03_analyze.R`](R/03_analyze.R): already ports a large slice of scup —
  `set_cycle_types`, `detect_trial_type`/`detect_muscle_group_keys`, `estimate_muscle_mass`/
  `add_muscle_mass`, `add_nominal_curvature_if_missing`, `filter_cycles_if_set`/
  `filter_transition_halfcycles`/`filter_poster_kinematics`/`filter_poster_activation`,
  `calc_work`, stiffness/damping (`.avg_mech`/`.high_mech`/`.low_mech` + wrappers),
  `get_mechanics_by_half_cycle`, `interpolate_even_phase`, `estimate_body_torque` (transform
  only), `detect_passive_mode`, `calc_muscle_torque`, `add_muscle_instantaneous`,
  `summarize_muscle_cycles`, PARTIAL Coughlin + PARTIAL sinusoid QC, and
  `analyze_isometric`/`analyze_isovelocity` stubs.
- [`R/plot_torque_vs_time_batch.R`](R/plot_torque_vs_time_batch.R): batch torque/stim figures
  (dynamic + segmented) — reusable for the subjective QC figures.

### 1.2 Gaps vs. the validated scup stack (what must be ported/completed)
| Area | scup source | Bender3 status | Action |
|------|-------------|----------------|--------|
| Muscle work-time / phase-power | `calc_muscle_work_time.R` | MISSING | PORT full |
| Coughlin exclusion + audit | `filter_coughlin_power.R` | PARTIAL (`flag_coughlin_power_violations`; no trial-drop, no audit/CSV) | COMPLETE |
| Sinusoid step quality | `filter_cycle_quality.R` | PARTIAL (`score_step_sinusoid_quality`; no harmonic fit / roughness / 4 fail-rules / exclusion / CSV) | COMPLETE |
| QC orchestrator | `apply_analysis_filters.R` | MISSING | PORT |
| Power-smoothness audit | `power_smoothness_audit.R` | MISSING | PORT (objective audit + subjective panels) |
| Torque QC checks | `estimate_body_torque.R` (`check_good_torque`, `check_symmetric_torque`, `check_torque_consistency`, `exclude_bad_files`) | MISSING (only transform ported) | PORT |
| Single-trial orchestrator | `process_muscle_trial.R` | MISSING | PORT as `run_dynamic_trial()` |
| Dispatcher | (scup infers from `Type` string) | MISSING | NEW `route_trial()` on `protocol_sampling_mode` |
| Geometry test region | `geometry_test_region.R` | MISSING | PORT (flag missing `bendlocation`, §7) |
| Notebook-grid filters | `filter_activation.R` (`filter_muscle_to_trial_grid`, cohort gates, `load_trial_sheet`) | MISSING | OUT OF SCOPE (§8) — depends on the scup notebook CSV, not the H5 |
| Per-cycle nominal parse | `nominal_stimulus.R` | SUPERSEDED by flat `index_cycle_*` | Port-as-equivalence-note only |

---

## 2. Port map (scup → Bender3)

Proposed module layout keeps the numbered-stage convention and isolates QC. New files under
`R/`. "Identical" = same math; "Schema change" = only the field access differs.

| scup function (file) | Bender3 target (file) | What stays identical | What changes for flat schema |
|----------------------|-----------------------|----------------------|------------------------------|
| `set_cycle_types` (`set_cycle_types.R`) | `set_cycle_types` (`03_analyze.R`) ✓done | act/pass + half-cycle order logic | `is_active_by_cycle` now from `index_cycle_active`; `stim` from `stim_side` |
| `detect_trial_type`/`detect_muscle_group_keys` (`detect_trial_type.R`) | same (`03_analyze.R`) ✓done | layout thresholds (`n_freq/duty/phase<=1`), group keys | protocol strings snake_case (`constant_frequency`, `frequency_amplitude_combo`, `frequency_sweep`) — verify vs actual `protocol_type` (§7) |
| `estimate_muscle_mass`/`add_muscle_mass` (`estimate_muscle_mass.R`) | same (`03_analyze.R`) ✓done | CSA regression (`0.001023`, `0.0001249`, `9.043e-7`), density `1060` | inputs `fishlength.m`/`dclamp.m` from `measurement_*` |
| `calc_work` (`calc_work.R`) | same (`03_analyze.R`) ✓done | `pracma::trapz`, radians guard `2*pi` | none |
| `stiffness_damping.*` (`stiffness_damping.R`) | `.avg/.high/.low_mech` + wrappers (`03_analyze.R`) ✓done | PCA slope, `qhigh=0.99`, `qlow=0.08`, 25% time exclusion | none |
| `detect_passive_mode`/`calc_muscle_torque` (`calc_muscle_torque.R`) | same (`03_analyze.R`) mostly done | active−passive template, `phase_bin_width=0.02`, torque priority | ADD scup `sideward`/`contraction`/`velo.rads` columns to reach full parity (§7) |
| `add_muscle_instantaneous`/`summarize_muscle_cycles` (`calc_muscle_metrics.R`) | same (`03_analyze.R`) ✓done | insta power = torque·dθ/dt; trapz work; mass-normalize | mass col `muscle_mass.kg` |
| `estimate_body_torque` transform (`estimate_body_torque.R`) | same (`03_analyze.R`) ✓done | `-xtorque`, `+ztorque`, `yforce*dclamp` | torque cols from `derived/` |
| `halfcycle_side`, `apply_muscle_work_sign`, `prepare_muscle_work_time`, `summarize_muscle_work_time_levels`, `prepare_muscle_cycle_phase_power`, `summarize_muscle_cycle_phase_power_levels` (`calc_muscle_work_time.R`) | `R/muscle_work_time.R` (NEW) | side thresh `0.01`, `cycles_keep=3:6`, `freq_keep=c(1,3,5,7)`, `time_bin_width=0.02`, `phase_shift=0.25`, `phase_bin_width=0.05` | none (operates on ported columns) |
| `coughlin_steady_state_power_limits`, `summarize_step_peak_power`, `audit_coughlin_power_trials`, `filter_trials_coughlin_power` (`filter_coughlin_power.R`) | `R/qc_coughlin.R` (NEW) | mean `133`, sd `19`, `hi=152` one-sided, trial-level drop, `lo=114` computed-unused | audit CSV lands at hub `qualitycontrol/` |
| `sinusoid_quality_defaults`, `score_step_sinusoid_quality`, `summarize_step_sinusoid_quality`, `audit_sinusoid_quality_steps`, `filter_low_quality_steps` (`filter_cycle_quality.R`) | `R/qc_cycle_quality.R` (NEW) | harmonic `r_sin/r_cos/r_harmonic`, roughness `sd(diff)/sd`, 4 fail rules + all thresholds, `min_n=20`, step-key rounding (freq 2, curve 1) | measured curvature = `curve.invm`; step key = `curvature.invm` |
| `apply_figure_quality_filters` (`apply_analysis_filters.R`) | `R/qc_apply.R` (NEW) | Coughlin-then-sinusoid order | writes combined audit to hub |
| `audit_trial_power_roughness`, `rank_trials_by_power_smoothness`, `prepare_*`, `plot_trial_power_*` (`power_smoothness_audit.R`) | `R/qc_power_smoothness.R` (NEW) | roughness = `median(abs(diff2))/peak`, `min_peak_Wkg=20`, `smoothness_score`, panel params | none |
| `check_good_torque`, `check_symmetric_torque`, `check_torque_consistency`, `exclude_bad_files` (`estimate_body_torque.R`) | `R/qc_torque.R` (NEW) | badness scoring (`20/10/5/1`, `warn_threshold=2`, `0.8`, `0.5`) | torque cols from flat |
| `process_muscle_trial` (`process_muscle_trial.R`) | `run_dynamic_trial()` (`R/run_dynamic_trial.R`, NEW) | load→cycle types→(body torque)→muscle torque→mass→summarize | uses `load_bender_flat()`; calibrate/deconvolve upstream |
| (none — new) | `route_trial()` (`R/dispatch.R`, NEW) | — | routes by `protocol_sampling_mode` + `protocol_type` |
| `compute_test_region` (`geometry_test_region.R`) | `R/geometry.R` (NEW) | `L = fishlen − bendloc − clamp/2` | `bendlocation` MISSING in flat schema (§7) |
| `nominal_stimulus.R` helpers | (no port) | activation/curvature-from-amplitude equivalence | SUPERSEDED by `index_cycle_*`; keep as validation cross-check only |

---

## 3. Schema mapping (flat field → quantity each ported function needs)

| Analysis quantity (scup col) | Bender3 flat source | Path/kind | Notes |
|------------------------------|---------------------|-----------|-------|
| Time `t.s` | `time_second` | `timeseries/` dataset | per-step local for segmented |
| Norm time `t.norm` | `time_normalized` | `timeseries/` dataset | drives `cycle`/`halfcycle` |
| Commanded angle `angle.deg` | `angle_commanded_degree` | `timeseries/` | kinematics |
| Encoder angle `enc.deg` | `angle_measured_degree` | `timeseries/` | specimen-frame |
| Calibrated torque `xtorque0.Nm`→`xtorque.Nm` | `derived/forcetorque_calibrated` col 4 (xTorque) | `derived/` (from `01_calibrate.R`) | bias(`t.s<-0.1`)+Butterworth in loader |
| Stim state `stim` (0/L/R/B) | `stim_side` | `timeseries/` categorical | remapped in loader |
| Stim voltage `stim.V` | `stim_channel1_command_volt` | `timeseries/` | COMMAND, not monitor (§7) |
| Freq `freq.Hz` | `index_cycle_frequency_hertz` | `metadata/` dataset | per-cycle join via `cycle_index` |
| Amplitude `amp.deg` | `index_cycle_motor_amplitude_degree` | `metadata/` | → curvature fallback |
| Duty / phase `duty`/`phase` | `index_cycle_activation_duty`/`_phase` | `metadata/` | |
| Active-by-cycle `is_active_by_cycle` | `index_cycle_active` | `metadata/` | drives `cycletype` |
| Nominal curvature `curvature.invm` | `index_cycle_operating_point`(+`_units`) or amp→`amp*pi/180/dclamp` | `metadata/` | matches scup fallback |
| Clamp span `dclamp.m` | `measurement_clamp_separation_millimeter`/1000 | `metadata/` attr | |
| Fish length `fishlength.m` | `measurement_specimen_bodylength_millimeter`/1000 | `metadata/` attr | SL fallback |
| Fish mass `fishmass.g` | `measurement_specimen_body_mass_gram` | `metadata/` attr | |
| Sample rate | `daq_ai_sample_rate_hertz` | `metadata/` attr | tibble attr `SampleFrequency.Hz` |
| Protocol id `stimulus_type` | `protocol_type` | `metadata/` attr | dispatcher key (§7) |
| Sampling mode | `protocol_sampling_mode` | `metadata/` attr | `single_finite`/`segmented_finite` |
| MOI / bias (deconvolve) | `calibration_inertia_*` | `metadata/` attrs | 02_deconvolve |

**OPEN-QUESTION fields (missing / ambiguous — see §7):** `bendlocation` (no flat equivalent →
`geometry_test_region` blocked); measured stim monitor vs command voltage; exact `protocol_type`
string values for dynamic sub-types; species-appropriateness of scup-specific muscle-mass and
Coughlin constants.

---

## 4. Validation plan (scup ⇄ Bender3 equivalence)

### 4.1 Overlap trials — legacy→flat converter (DECIDED)
- Build a one-off, validation-only converter `tools/convert_legacy_to_flat.R` (NOT on the
  Bender3 run path). Input: a handful of validated **scup** legacy `.h5` trials (dynamic:
  combo / constant-frequency / sweep). Output: flat-schema `.h5` twins under
  `.../qualitycontrol/validation/flat_twins/`.
- Converter mapping (legacy → flat):
  - `/NominalStimulus/{t,tnorm,Position,Velocity}` → `timeseries/{time_second,time_normalized,angle_commanded_degree,angular_velocity_commanded_degree_per_second}`; `/Calibrated/Encoder` → `angle_measured_degree`.
  - Per-cycle nominal (`Frequencies`,`Curvatures`,`CyclesPerStep`,`IsActiveByCycle`, or `AmplitudeByCycle`/`FrequencyByCycle`, or constant-freq `Amplitude`/`Frequency`/`Cycles`/`ActivationStartCycle`) → `metadata/index_cycle_*` using the scup `nominal_stimulus.R` expansion rules verbatim.
  - `stim` windows (`Lonoff`/`Ronoff` or reconstructed) → `timeseries/stim_side` + `stim_channel1_command_volt`.
  - Root geometry attrs → `metadata/measurement_*`; `SampleFrequency` → `daq_ai_sample_rate_hertz`.
  - **Calibrated torque carry-over:** write the legacy `/Calibrated/xTorque…zForce` directly into `derived/forcetorque_calibrated` (correct column order via a synthetic `daq_ai_channel_map`). This isolates ANALYSIS parity from calibration differences.
  - `protocol_sampling_mode = single_finite`; `protocol_type` = mapped dynamic string.

### 4.2 Two-tier comparison (scup executed for real — DECIDED)
- **Tier 0 (calibration parity, informational):** Bender3 `01_calibrate.R` re-derives calibrated
  F/T from raw+matrix (only if a raw twin is also written) vs the legacy `/Calibrated`. Reported,
  not gating for M1 (M1 carries calibrated torque over, §4.1).
- **Tier 1 (analysis parity, GATING):** run scup (`process_muscle_trial`) on the ORIGINAL legacy
  file and Bender3 (`run_dynamic_trial`) on the flat twin. Align pre-processing exactly: same LP
  `cutoff` (use `15` Hz), same bias window (`t.s < -0.1`), same `cycles_keep`.

### 4.3 Compared quantities & tolerances
| Quantity | Column(s) | Tolerance (pass) |
|----------|-----------|------------------|
| Filtered torque | `xtorque.Nm` | max abs Δ ≤ `1e-6` N·m OR rel ≤ `1e-4` |
| Instantaneous curvature / rate | `curve.invm`, `curverate.invms` | rel ≤ `1e-6` |
| Cycle / half-cycle indexing | `cycle`, `halfcycle`, `cycletype`, `halfcycletype`, `halfcycleorder` | EXACT |
| Muscle torque | `muscle_torque.Nm` | rel ≤ `1e-3` |
| Instantaneous power | `insta_power.W`, `insta_power.Wkg` | rel ≤ `1e-3` |
| Per-cycle metrics | `work.J`, `avg_power.W`, `peak_torque.Nm`, `work.Jkg`, `avg_power.Wkg` | rel ≤ `1e-3` |
| Stiffness / damping | `EI1/EIM/EIL.Nm2`, `etaI1/etaIM/etaIL.Nm2s` | rel ≤ `1e-3` |
| Muscle mass | `muscle_csa.m2`, `muscle_mass.kg` | rel ≤ `1e-9` |
| QC — Coughlin | `peak_power_Wkg`, `excluded_trials` set | value rel ≤ `1e-6`; excluded set EXACT |
| QC — sinusoid | `r_harmonic`, `roughness`, `fail_*`, excluded step keys | metrics rel ≤ `1e-6`; flags/keys EXACT |

Relative tolerance `|a−b| ≤ tol·max(|a|,|b|,ε)`, `ε=1e-12`. Any EXACT mismatch or metric beyond
tolerance = FAIL with a per-quantity diff report.

### 4.4 Report location & format
- `R/validate_scup_equivalence.R` (NEW) sources scup (test oracle), runs both, diffs, writes:
  - `.../qualitycontrol/validation/equivalence_report_<date>.csv` (per-quantity max Δ, tol, pass/fail)
  - `.../qualitycontrol/validation/equivalence_report_<date>.md` (human summary + any FAIL rows)
  - overlay PNGs (scup vs Bender3) per trial for eyeball confirmation.
- Base output dir constant:
  `/Users/yjimenez/Library/CloudStorage/OneDrive-ProvidenceCollege/01_JimenezLab/02_ResearchHub/proj_crittergripper/qualitycontrol`.

---

## 5. Milestone-sequenced steps

Each step: goal · files · scup source · in/out · done-when.

- **M1.0 — Dispatcher.** Goal: route all trials by type. File: `R/dispatch.R` (`route_trial`).
  scup: (new; replaces `Type`-string inference). In: file path / loaded `td` (`protocol_sampling_mode`,
  `protocol_type`). Out: `list(mode, protocol, handler)` where `single_finite`→dynamic,
  `segmented_finite`→segmented stub. Done-when: a mixed folder routes every file; dynamic →
  analysis handler, segmented → `analyze_isometric`/`analyze_isovelocity` stubs; unknown →
  logged, no crash.

- **M1.1 — Complete `calc_muscle_torque` parity.** Goal: full faithful port. File: `R/03_analyze.R`.
  scup: `calc_muscle_torque.R`. In: loaded dynamic `td`. Out: adds scup-parity cols
  (`velo.rads`, `sideward`, `contraction`) + existing `muscle_torque.Nm`. Done-when: columns
  match scup names/values on a converted twin within §4.3 tol.

- **M1.2 — Port work-time / phase-power.** File: `R/muscle_work_time.R`. scup:
  `calc_muscle_work_time.R`. In: muscle time series. Out: `prepare_muscle_work_time`,
  `*_levels`, `prepare_muscle_cycle_phase_power`, `*_levels` (+ `side`, `muscle_work.J`,
  signed power). Done-when: functions reproduce scup outputs within tol on a twin.

- **M1.3 — Port torque QC checks.** File: `R/qc_torque.R`. scup: `estimate_body_torque.R`
  QC block. In: loaded `td`. Out: `check_good_torque` (`badness`, pre/post L/R), `check_symmetric_torque`,
  `check_torque_consistency`, `exclude_bad_files`. Done-when: `badness` and flags match scup on twins.

- **M1.4 — Complete Coughlin QC.** File: `R/qc_coughlin.R`. scup: `filter_coughlin_power.R`.
  In: muscle df (+ optional `cycle_metrics`). Out: `summarize_step_peak_power`,
  `audit_coughlin_power_trials` (`by_step`/`flagged_steps`/`excluded_trials`),
  `filter_trials_coughlin_power` (trial drop + CSV). Done-when: excluded-trial set EXACT vs scup;
  `hi=152` one-sided; CSV at hub.

- **M1.5 — Complete sinusoid cycle-quality QC.** File: `R/qc_cycle_quality.R`. scup:
  `filter_cycle_quality.R`. In: muscle df. Out: harmonic fit (`r_sin/r_cos/r_harmonic`),
  `roughness`, 4 fail rules, `audit_sinusoid_quality_steps`, `filter_low_quality_steps`
  (step drop + CSV). Done-when: `fail_*`/excluded step keys EXACT vs scup; metrics within tol.

- **M1.6 — Power-smoothness audit.** File: `R/qc_power_smoothness.R`. scup:
  `power_smoothness_audit.R`. In: muscle df. Out: `audit_trial_power_roughness`,
  `rank_trials_by_power_smoothness` (objective) + `prepare_*`/`plot_trial_power_*` (subjective).
  Done-when: ranking table matches scup; QC panels render.

- **M1.7 — QC orchestrator + single-trial driver.** Files: `R/qc_apply.R`
  (`apply_figure_quality_filters`), `R/run_dynamic_trial.R` (`run_dynamic_trial` = scup
  `process_muscle_trial`). scup: `apply_analysis_filters.R`, `process_muscle_trial.R`.
  In: dynamic file. Out: `list(raw, muscle, cycle_metrics, trial_info, coughlin_audit, quality_audit)`
  + written audit record (pass/fail per trial+step) to
  `.../qualitycontrol/dynamic_qc_audit_<date>.csv`. Done-when: end-to-end dynamic QC produces an
  auditable pass/fail record + subjective figures for a real dynamic file.

- **M1.8 — Converter.** File: `tools/convert_legacy_to_flat.R`. scup: reader mirror of
  `load_bender_data.R` + `nominal_stimulus.R`. In: legacy scup `.h5`. Out: flat twin `.h5`.
  Done-when: `load_bender_flat()` loads the twin and yields the same column set/values as scup's
  loader within tol.

- **M1.9 — Equivalence harness (GATING).** File: `R/validate_scup_equivalence.R`. In: legacy
  original + flat twin. Out: CSV+MD report + overlay PNGs at hub. Done-when: all §4.3 quantities
  PASS on ≥3 dynamic twins (combo, constant-freq, sweep); report committed-path written.

- **M1.10 — Geometry test region (best-effort).** File: `R/geometry.R`. scup:
  `geometry_test_region.R`. Done-when: `compute_test_region` ported; if `bendlocation` absent in
  flat schema, function returns `NA` + logs the OPEN QUESTION (§7) rather than guessing.

---

## 6. QC spec (dynamic / `single_finite`)

### 6.1 Objective screens (ported verbatim)
- **Coughlin steady-state power (trial-level).** Limits `133 ± 19` W/kg (Coughlin et al. 1996).
  Step = `trial_id × freq.Hz × curvature.invm`; `peak_power_Wkg = max(|insta_power.Wkg|)`.
  Fail step if `> 152` (`mean+1SD`); fail TRIAL if ANY step fails → drop trial. `114` (`mean−1SD`)
  computed but one-sided (not used to exclude). Active-only when `cycletype=="act"`.
- **Sinusoid / cycle quality (step-level).** On middle active half-cycles, `min_n=20`:
  - `r_torque_curvature = cor(muscle_torque.Nm, curve.invm)`; `r_harmonic = sqrt(r_sin²+r_cos²)`
    against `sin/cos(2π·phase)`; `roughness = sd(diff(tq)) / max(sd(tq), 1e-9)`; `zc_rate`
    diagnostic-only.
  - Fail if ANY: decoupled (`|r_tc|<0.25 & r_harmonic≥0.55`); non-sinusoidal (`r_harmonic<0.40`);
    rough (`roughness>0.12 & |r_tc|<0.30`); both-weak (`|r_tc|<0.20 & r_harmonic<0.35`).
    Drop the offending `trial|round(freq,2)|round(curve,1)` step (not whole trial).
- **Power-smoothness audit (trial-level, ranking not exclusion).**
  `roughness = median(|diff²(power)|)/peak`; `smoothness_score = peak/(roughness+1e-6)`; ranking
  pool requires `peak_power_Wkg ≥ 20`.
- **Torque sanity (`qc_torque.R`).** `check_good_torque` badness scoring (missing pre-L/R `+10`,
  wrong/ swapped signs `+10/+5`, L/R ratio `>2` log-scaled, act/pass sign agreement `<0.8` `+5`);
  `exclude_bad_files` NA-outs torque when `badness > 2`; `check_torque_consistency` needs
  xt/zt/yf correlations ≥ `0.5`.

### 6.2 Subjective figures (visual inspection)
- Reuse `plot_torque_vs_time_batch.R` (torque-vs-time, stim overlay) + ported
  `plot_trial_power_*_panels` (per-frequency power-vs-phase/time panels, rolling mean `k=5`).
- Figures written alongside the audit under `.../qualitycontrol/figures/`.

### 6.3 Auditable pass/fail record
- Single combined audit CSV per run at
  `.../qualitycontrol/dynamic_qc_audit_<date>.csv` with, per trial (+ per step where relevant):
  `trial_id`, `protocol_type`, Coughlin `peak_power_Wkg`+pass, sinusoid per-step `fail_*`+`quality_fail`,
  smoothness `roughness`/`smoothness_score`, torque `badness`, and an overall `qc_pass` verdict
  with a `reason` string. Excluded-trial and excluded-step CSVs mirror scup names.

---

## 7. Risks & open questions
1. **`protocol_type` string values.** scup keys on Title-Case (`"Constant Frequency"`), Bender3
   `detect_trial_type` uses snake_case. Confirm the ACTUAL `metadata/protocol_type` values written
   for dynamic sub-types (dynamic vs sweep vs combo vs constant-frequency) so routing/layout is correct.
2. **Missing `bendlocation`.** No flat equivalent of legacy `BendLocation_mm` →
   `geometry_test_region` cannot be ported faithfully. Options: add a `measurement_*` field, or
   scope test-region out. (Muscle-mass estimate itself does NOT need it.)
3. **Stim: command vs monitor.** scup used `/RawInput/activation_monitor` (measured). Flat loader
   uses `stim_channel1_command_volt` (commanded) + categorical `stim_side`. For faithful stim
   timing this should be equivalent for `stim` factor derivation, but numeric stim voltage differs.
4. **Species-specific constants.** Muscle-mass CSA regression and Coughlin limits are scup red
   muscle. Applying them to non-scup specimens is a scientific decision — port faithfully now, flag
   for PI before cohort use.
5. **Calibration path difference.** Legacy files ship pre-calibrated `/Calibrated`; Bender3
   re-derives from raw+matrix. Tier-0 parity (§4.2) surfaces any matrix/channel-order discrepancy;
   M1 sidesteps it by carrying calibrated torque into the twin.
6. **LP filter / cutoff alignment.** scup batch uses `cutoff=15`; loader default is
   `cutoffmult*max(freq)`. Validation must pin identical cutoff or torque parity fails spuriously.
7. **`calc_muscle_torque` `t.perc` semantics.** scup `sideward` thresholds (`0.25/0.75`) act on
   within-cycle elapsed seconds `t.perc`, not normalized phase — verify behavior across frequencies
   when completing M1.1 parity.
8. **Converter faithfulness = validation trust anchor.** A converter bug would mask a port bug.
   Mitigate with Tier-0/loader cross-checks and per-quantity overlays.

---

## 8. Future phases (brief)
- **Segmented analysis (M2+).** Flesh out `analyze_isometric` (force-length curve, `F0`/`L0`) and
  `analyze_isovelocity` (force-velocity, `P0`/`Vmax`, peak power), Hill fits, fatigue maps —
  currently stubs.
- **Cohort / notebook layer.** Port scup `load_trial_sheet`/`filter_muscle_to_trial_grid`/cohort
  gates + `batch_muscle_processing` if/when a Bender3 trial manifest beyond the minimal QC audit is
  needed. Out of scope now (depends on the scup notebook CSV, not the H5).
- **Retroactive GUI-field export audit** and apparatus-inertia matrix workflow (schema Deferred
  flag 4) — unrelated to this port.
