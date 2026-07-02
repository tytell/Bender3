"""
Canonical Bender -> HDF5 routing spec (schema v2.1).

THIS IS A SPEC / CONTRACT MODULE. It does not write files. It declares, for every
attribute that can exist on a :class:`bender_functions.Bender` instance at export
time, exactly where that attribute goes in the schema-v2 HDF5 file -- or that it is
deliberately excluded.

Three public objects:

* ``BENDER_ROUTING`` : dict[str, Route]
    Maps a Bender attribute name -> routing record. The future exporter iterates
    ``bender.__dict__`` and, for every public attribute, looks it up here.

* ``EXCLUDED`` : set[str]
    Attributes that are intentionally NOT written to the raw file (internal caches,
    machine-specific paths, objects, legacy aliases, and derived series that belong in
    the research-hub ``derived`` group). All config-sourced hardware/wiring/timing
    provenance is now routed canonically under daq_/calibration_/protocol_ in
    ``BENDER_ROUTING`` (the verbatim config-provenance writer is retired); the only
    config fields intentionally dropped are the legacy ``units`` / ``unit_rules`` dicts.

* ``MISSING_REQUIRED`` : list[dict]
    Schema fields with NO Bender source today. These need a GUI/config wire-up (or
    are populated downstream). Listed so the gap is explicit, not silent.

DESIGN CONTRACT (fail-on-unmapped):
    Every public (non-underscore) attribute on a Bender instance must be either a key
    in ``BENDER_ROUTING`` or a member of ``EXCLUDED``. Any attribute that is neither
    must raise in the exporter. Underscore-prefixed attributes are implicitly excluded.
    Use :func:`check_coverage` to verify an instance against this spec.

Route record fields:
    tier      : 'metadata' | 'timeseries'
    key       : canonical schema key (spelled-out units, see context_jlab_cg_h5schema.md §1)
    required  : bool  -- exporter should warn/raise if absent or empty
    source    : 'GUI' | 'config' | 'DAQ' | 'computed'
    note      : optional clarification (many-to-one merges, pending decisions, etc.)

Schema reference: /Users/yjimenez/code/Bender3/context_jlab_cg_h5schema.md (v2.1)

PROTOCOL VALUE / UNIT CONVENTION (confirmed):
    The file declares its protocol once in ``protocol_type`` (from ``test_type``), e.g.
    "isometric". For protocol parameters whose unit is chosen at runtime by a companion
    ``*_mode`` attribute, the routing uses a value + unit pair:
      * the VALUE attr -> ``protocol_..._<param>``         (the numbers/quantity)
      * the ``*_mode`` attr -> ``protocol_..._<group>_unit`` (resolved spelled-out unit,
        e.g. degree/percent/per_meter, computed at write time) AND, for provenance, the
        raw convention is additionally written as ``protocol_..._<group>_input_mode``.
    One ``_unit`` field covers every value sharing the same mode (e.g. isometric
    initial/final/mirror targets all share ``protocol_isometric_target_unit``).

All routing decisions are confirmed; see bender_routing_decision_log.md for the
recorded answers and rationale.
"""

from __future__ import annotations

from typing import Any, Dict, List, Set, TypedDict


class Route(TypedDict, total=False):
    tier: str
    key: str
    required: bool
    source: str
    note: str


# ---------------------------------------------------------------------------
# BENDER_ROUTING
# ---------------------------------------------------------------------------
BENDER_ROUTING: Dict[str, Route] = {

    # === metadata / specimen_ (non-spatial identity) =======================
    "specimen_id":              {"tier": "metadata", "key": "specimen_id",            "required": True,  "source": "GUI"},
    "specimen_genusspecies":    {"tier": "metadata", "key": "specimen_genusspecies",  "required": False, "source": "GUI"},
    "specimen_prep_condition":  {"tier": "metadata", "key": "specimen_prep_condition", "required": False, "source": "GUI"},
    "specimen_sex":         {"tier": "metadata", "key": "specimen_sex",           "required": False, "source": "GUI"},
    "specimen_muscle_type": {"tier": "metadata", "key": "specimen_muscle_type",   "required": False, "source": "GUI"},
    "segment":              {"tier": "metadata", "key": "specimen_segment",       "required": False, "source": "GUI",
                             "note": "ASSUMED specimen_segment (MC #5 default)"},

    # === metadata / measurement_ (ALL spatial geometry, animal + apparatus) =
    "fishlen_TL":            {"tier": "metadata", "key": "measurement_specimen_bodylength_millimeter",        "required": False, "source": "GUI"},
    "fishlen_SL":            {"tier": "metadata", "key": "measurement_specimen_standardlength_millimeter",    "required": False, "source": "GUI"},
    "fishmass":              {"tier": "metadata", "key": "measurement_specimen_body_mass_gram",               "required": False, "source": "GUI"},
    "xsec_width":            {"tier": "metadata", "key": "measurement_xsec_width_millimeter",                 "required": True,  "source": "GUI"},
    "xsec_height":           {"tier": "metadata", "key": "measurement_xsec_height_millimeter",                "required": False, "source": "GUI"},
    "dclamp":                {"tier": "metadata", "key": "measurement_clamp_separation_millimeter",           "required": True,  "source": "GUI"},
    "dvert":                 {"tier": "metadata", "key": "measurement_clamp_offset_vertical_millimeter",      "required": False, "source": "GUI"},
    "dhoriz":                {"tier": "metadata", "key": "measurement_clamp_offset_horizontal_millimeter",    "required": False, "source": "GUI"},
    "clamp_plate_extension_mm": {"tier": "metadata", "key": "calibration_inertia_apparatus_plate_to_plate_millimeter", "required": False, "source": "GUI",
                                   "note": "plate-to-plate span between body plates within the moving clamp (mm); GUI field 'Inter-clamp span'; "
                                           "feeds apparatus-calibration-matrix lookup (Deferred flag 4); "
                                           "replaces measurement_clamp_plate_extension_millimeter (v2.7)"},
    "target_muscle_depth_mm":   {"tier": "metadata", "key": "measurement_target_muscle_depth_millimeter",     "required": False, "source": "GUI"},
    "specimen_geometry_heights_mm":   {"tier": "metadata", "key": "measurement_specimen_frustum_heights_millimeter",   "required": False, "source": "GUI"},
    "specimen_geometry_depths_mm":    {"tier": "metadata", "key": "measurement_specimen_frustum_depths_millimeter",    "required": False, "source": "GUI"},
    "specimen_geometry_positions_mm": {"tier": "metadata", "key": "measurement_specimen_frustum_positions_millimeter", "required": False, "source": "GUI"},
    "specimen_profile_length_mm":     {"tier": "metadata", "key": "measurement_profile_length_millimeter",     "required": False, "source": "GUI",
                                       "note": "inertial profile model input (spatial -> measurement_ per Rule 4)"},
    "specimen_profile_clamp_offset_mm": {"tier": "metadata", "key": "calibration_inertia_apparatus_aor_to_clamp_millimeter", "required": False, "source": "GUI",
                                       "note": "axial dist from sensor face (AOR in standard inline config) to near edge of rotating clamp (mm); "
                                               "GUI field 'Distance from rotation axis to clamps'. Many-to-one: frustum_clamp_offset_mm maps to same key "
                                               "(first positive value wins). Omit if <= 0; feeds apparatus-cal-matrix lookup (Deferred flag 4). "
                                               "Replaces measurement_profile_clamp_offset_millimeter (v2.7)."},
    # --- frustum model raw inputs (v2.7 amended — previously dropped inside frustum_inputs dict) ---
    # Source: values written into bender.frustum_inputs dict by set_frustum_inertial_model; exporter
    # reads from the dict directly via a dedicated block (ledger also fires if attrs set individually
    # via update_metadata).
    "frustum_height_mm":        {"tier": "metadata", "key": "measurement_specimen_inertia_frustum_height_millimeter",            "required": False, "source": "computed",
                                  "note": "back-end (max) cross-section height for the elliptical frustum MOI model (mm); frustum path only"},
    "frustum_width_mm":         {"tier": "metadata", "key": "measurement_specimen_inertia_frustum_width_millimeter",             "required": False, "source": "computed",
                                  "note": "back-end (max) cross-section width for the elliptical frustum MOI model (mm); frustum path only"},
    "frustum_length_mm":        {"tier": "metadata", "key": "measurement_specimen_inertia_frustum_length_millimeter",            "required": False, "source": "computed",
                                  "note": "axial span of the frustum model (mm); frustum path only"},
    "frustum_clamp_offset_mm":  {"tier": "metadata", "key": "calibration_inertia_apparatus_aor_to_clamp_millimeter",            "required": False, "source": "computed",
                                  "note": "apparatus AOR-to-clamp distance (mm); same physical quantity as specimen_profile_clamp_offset_mm; "
                                          "many-to-one (first positive value wins); read from frustum_inputs dict by exporter. "
                                          "Replaces draft measurement_specimen_frustum_clamp_offset_millimeter (v2.7)."},
    "frustum_density_g_per_mm3": {"tier": "metadata", "key": "measurement_specimen_inertia_frustum_density_gram_per_cubic_millimeter", "required": False, "source": "computed",
                                   "note": "material density for the frustum model (g/mm^3); frustum path only"},
    "dbend":                  {"tier": "metadata", "key": "measurement_bend_position_millimeter", "required": False, "source": "GUI",
                              "note": "MC #8: bend/segment position -> measurement_"},
    "test_segment_position_mm": {"tier": "metadata", "key": "measurement_bend_position_millimeter", "required": False, "source": "GUI",
                              "note": "alias of dbend; same key (first present)"},

    # === metadata / session_ (recording-session logistics) =================
    "session_date":       {"tier": "metadata", "key": "session_date",            "required": True,  "source": "computed",
                           "note": "set at run_experiment start (ISO-8601)"},
    "session_apparatus_id": {"tier": "metadata", "key": "session_apparatus_id",  "required": False, "source": "config"},
    "session_analyst":    {"tier": "metadata", "key": "session_analyst",         "required": False, "source": "GUI"},
    "session_simulated":  {"tier": "metadata", "key": "session_simulated",       "required": True,  "source": "computed",
                           "note": "bool on EVERY file; sim_ filename prefix mirrors this"},
    "simulation_material": {"tier": "metadata", "key": "session_simulation_material", "required": False, "source": "GUI",
                            "note": "meaningful only when session_simulated True"},
    "config_name":        {"tier": "metadata", "key": "session_protocol_template", "required": False, "source": "config"},
    "acquisition_start":  {"tier": "metadata", "key": "session_acquisition_start", "required": False, "source": "computed",
                           "note": "per-run acquisition timestamp; distinct from session_date"},
    "temp_C_room":        {"tier": "metadata", "key": "session_temperature_room_celsius", "required": False, "source": "GUI",
                           "note": "ASSUMED session_ placement (MC #4 default)"},
    "temp_C_tank":        {"tier": "metadata", "key": "session_temperature_tank_celsius", "required": False, "source": "GUI",
                           "note": "ASSUMED session_ placement (MC #4 default)"},

    # === metadata / calibration_ ===========================================
    "calibration":   {"tier": "metadata", "key": "calibration_forcetorque_matrix", "required": True,  "source": "config",
                      "note": "ATI 6x6, raw volts -> newton / newton_meter"},
    "forcetorque_calibration_file": {"tier": "metadata", "key": "calibration_forcetorque_file", "required": False, "source": "config",
                                     "note": "referenced .cal filename; pairs with calibration_forcetorque_matrix"},
    "sono_cal_left":  {"tier": "metadata", "key": "calibration_sono_left_millimeter_per_volt",  "required": False, "source": "config",
                       "note": "V->mm breakpoint table [Low_V,High_V,Low_mm,High_mm] (or longer multi-point)"},
    "sono_cal_right": {"tier": "metadata", "key": "calibration_sono_right_millimeter_per_volt", "required": False, "source": "config",
                       "note": "V->mm breakpoint table [Low_V,High_V,Low_mm,High_mm] (or longer multi-point)"},

    # === metadata / daq_ (acquisition + rig wiring/hardware) ================
    # Per PI decision: all config-sourced hardware/wiring provenance is routed canonically
    # under daq_ (no separate hardware_ group). Nothing from the config is excluded except
    # the legacy units / unit_rules dicts (redundant with spelled-out unit suffixes).
    "daq_ai_sample_rate_hz":    {"tier": "metadata", "key": "daq_ai_sample_rate_hertz",    "required": True, "source": "config"},
    "daq_ao_do_sample_rate_hz": {"tier": "metadata", "key": "daq_ao_do_sample_rate_hertz", "required": True, "source": "config"},
    "input_channels":      {"tier": "metadata", "key": "daq_ai_channel_map",   "required": True,  "source": "computed",
                            "note": "zipped with input_channel_names -> index:identity map"},
    "input_channel_names": {"tier": "metadata", "key": "daq_instrumentation",  "required": False, "source": "computed",
                            "note": "active sensor list; also feeds daq_ai_channel_map"},
    # --- DAQ wiring (NI device + channels + ports) ---
    "device_name":     {"tier": "metadata", "key": "daq_device_name",            "required": False, "source": "config"},
    "motor_port":      {"tier": "metadata", "key": "daq_motor_port",             "required": False, "source": "config"},
    "encoder_chan":    {"tier": "metadata", "key": "daq_encoder_channel",        "required": False, "source": "config"},
    "stim_channels":   {"tier": "metadata", "key": "daq_stim_channels",          "required": False, "source": "config", "note": "list (ao0, ao1)"},
    "S1stim_chan":     {"tier": "metadata", "key": "daq_stim_channel1",          "required": False, "source": "config", "note": "derived from stim_channels[0]"},
    "S2stim_chan":     {"tier": "metadata", "key": "daq_stim_channel2",          "required": False, "source": "config", "note": "derived from stim_channels[1]"},
    "sg_channels":     {"tier": "metadata", "key": "daq_forcetorque_channels",       "required": False, "source": "config"},
    "sg_names":        {"tier": "metadata", "key": "daq_forcetorque_channel_names",  "required": False, "source": "config"},
    "sono_channels":   {"tier": "metadata", "key": "daq_sono_channels",          "required": False, "source": "config"},
    "sono_names":      {"tier": "metadata", "key": "daq_sono_channel_names",     "required": False, "source": "config"},
    "stim_monitor_chan": {"tier": "metadata", "key": "daq_stim_monitor_channel", "required": False, "source": "config",
                          "note": "config-file-only; loaded onto Bender in __init__ so the ledger can route it"},
    "stim_monitor_name": {"tier": "metadata", "key": "daq_stim_monitor_name",    "required": False, "source": "config",
                          "note": "config-file-only; loaded onto Bender in __init__ so the ledger can route it"},
    "use_sono":        {"tier": "metadata", "key": "daq_sono_enabled",           "required": False, "source": "config", "note": "bool"},
    "sono_internal_rate": {"tier": "metadata", "key": "daq_sono_internal_sample_rate_hertz", "required": False, "source": "config",
                           "note": "DS3 internal update rate (~241-242 Hz); distinct from daq_ai_sample_rate_hertz"},
    # --- rig hardware / mechanics ---
    "motor_gear_ratio":        {"tier": "metadata", "key": "daq_motor_gear_ratio",             "required": False, "source": "config", "note": "dimensionless ratio"},
    "motor_full_steps_per_rev": {"tier": "metadata", "key": "daq_motor_full_steps_per_revolution", "required": False, "source": "config"},
    "encoder_pulses_per_rev":  {"tier": "metadata", "key": "daq_encoder_pulses_per_revolution", "required": False, "source": "config"},
    "motor_axis":              {"tier": "metadata", "key": "daq_motor_axis_sensor",            "required": False, "source": "config", "note": "categorical"},
    "bending_axis_sensor":     {"tier": "metadata", "key": "daq_bending_axis_sensor",          "required": False, "source": "config", "note": "categorical"},
    "bending_axis_specimen":   {"tier": "metadata", "key": "daq_bending_axis_specimen",        "required": False, "source": "config",
                                "note": "config-file-only; loaded onto Bender in __init__ so the ledger can route it"},
    "primary_bending_axis":    {"tier": "metadata", "key": "daq_primary_bending_axis",         "required": False, "source": "config", "note": "categorical"},
    "positive_motor_direction": {"tier": "metadata", "key": "daq_positive_motor_direction",    "required": False, "source": "config", "note": "categorical (left/right)"},
    "S1side":                  {"tier": "metadata", "key": "daq_stim_channel1_side",           "required": False, "source": "config"},
    "S2side":                  {"tier": "metadata", "key": "daq_stim_channel2_side",           "required": False, "source": "config"},
    "specimen_lateral_index_on_positive_motor_side": {"tier": "metadata", "key": "daq_specimen_lateral_index_on_positive_motor_side", "required": False, "source": "config",
                                "note": "signed index, dimensionless"},
    "specimen_side_index_left":  {"tier": "metadata", "key": "daq_specimen_side_index_left",   "required": False, "source": "config", "note": "derived"},
    "specimen_side_index_right": {"tier": "metadata", "key": "daq_specimen_side_index_right",  "required": False, "source": "config", "note": "derived"},
    "daq_motor_positive_bend_lateral_index": {"tier": "metadata", "key": "daq_motor_positive_bend_lateral_index", "required": False, "source": "computed",
                                "note": "signed lateral index the positive motor direction bends toward; set in isometric/isovelocity run paths"},

    # === metadata / calibration_inertia_ (computed MOI + calibration scalars) =
    # inertial_calibration_profile is handled by a direct flat-key write block in the
    # exporter (not the ledger pass); it lives in EXCLUDED to suppress fail-on-unmapped.
    # apparatus/bias/axis keys are written by that block, not the ledger (see MISSING_REQUIRED).
    "i_total_system": {"tier": "metadata",
                       "key": "calibration_inertia_total_moi_gram_millimeter_squared",
                       "required": False, "source": "computed",
                       "note": "theoretical apparatus baseline (i_shaft+i_clamps, hardcoded) + "
                               "specimen MOI; NOT apparatus_moi (I_est) + specimen_moi "
                               "(see Deferred flag 5 in migration plan)"},
    "specimen_moi_specimen": {"tier": "metadata",
                               "key": "calibration_inertia_specimen_moi_gram_millimeter_squared",
                               "required": False, "source": "computed",
                               "note": "many-to-one: exporter writes first present of "
                                       "specimen_moi_specimen / profile / frustum"},
    "specimen_moi_profile":  {"tier": "metadata",
                               "key": "calibration_inertia_specimen_moi_gram_millimeter_squared",
                               "required": False, "source": "computed",
                               "note": "see specimen_moi_specimen"},
    "specimen_moi_frustum":  {"tier": "metadata",
                               "key": "calibration_inertia_specimen_moi_gram_millimeter_squared",
                               "required": False, "source": "computed",
                               "note": "see specimen_moi_specimen"},
    "inertial_calibration_file": {"tier": "metadata",
                                  "key": "calibration_inertia_file",
                                  "required": False, "source": "GUI"},
    "specimen_geometry_density_g_per_mm3": {
        "tier": "metadata",
        "key": "measurement_specimen_density_gram_per_cubic_millimeter",
        "required": False, "source": "GUI",
        "note": "measured input for MOI computation; moved from inertial_ to measurement_ (Point 1)"},
    "inertial_specimen_from_geometry": {
        "tier": "metadata",
        "key": "calibration_inertia_specimen_from_geometry",
        "required": False, "source": "computed",
        "note": "bool: analytic specimen MOI available (set in run_experiment FT block)"},

    # === metadata / note_ ==================================================
    "stim_clamp_notices": {"tier": "metadata", "key": "note_bench",    "required": False, "source": "computed",
                           "note": "bench-side stim/clamp warnings (schema §6)"},
    "post_trial_notes":   {"tier": "metadata", "key": "note_posthoc",  "required": False, "source": "GUI",
                           "note": "ASSUMED note_posthoc to avoid collision with stim_clamp_notices->note_bench (overrides MC #1 default)"},

    # === metadata / protocol_ (global defaults; per-step realized -> index_) =
    # File-level protocol selector: one value per file (whole session is one test_type).
    "test_type":            {"tier": "metadata", "key": "protocol_type",                  "required": True,  "source": "GUI",
                             "note": "category: isometric/isovelocity/dynamic/frequency_sweep"},
    "starting_angle_deg":   {"tier": "metadata", "key": "protocol_starting_angle_degree", "required": False, "source": "GUI"},
    "duration":             {"tier": "metadata", "key": "protocol_duration_second",       "required": False, "source": "GUI"},
    "amplitude_frequency_exponent": {"tier": "metadata", "key": "protocol_amplitude_frequency_exponent", "required": False, "source": "GUI"},
    "velocity_exponent":    {"tier": "metadata", "key": "protocol_velocity_exponent",     "required": False, "source": "GUI"},
    "is_stim":              {"tier": "metadata", "key": "protocol_stim_enabled",          "required": False, "source": "GUI", "note": "bool"},
    "recruitment":          {"tier": "metadata", "key": "protocol_recruitment",           "required": False, "source": "GUI"},
    "lateral_mode":         {"tier": "metadata", "key": "protocol_lateral_mode",          "required": False, "source": "GUI"},
    "bilateral_mirror_motor":          {"tier": "metadata", "key": "protocol_bilateral_mirror_motor",          "required": False, "source": "GUI"},
    "bilateral_sequential_left_frac":  {"tier": "metadata", "key": "protocol_bilateral_sequential_left_fraction", "required": False, "source": "GUI"},
    # isometric mirror targets are VALUES; their unit is protocol_isometric_target_unit (from isometric_mode)
    "isometric_mirror_target_left":    {"tier": "metadata", "key": "protocol_isometric_mirror_target_left",    "required": False, "source": "GUI", "note": "value; unit in protocol_isometric_target_unit"},
    "isometric_mirror_target_right":   {"tier": "metadata", "key": "protocol_isometric_mirror_target_right",   "required": False, "source": "GUI", "note": "value; unit in protocol_isometric_target_unit"},
    "strain_shortening_positive_display_sign": {"tier": "metadata", "key": "protocol_strain_shortening_positive_display_sign", "required": False, "source": "GUI"},
    "rest_between_steps_s":  {"tier": "metadata", "key": "protocol_rest_between_steps_second", "required": False, "source": "GUI"},
    "randomize_step_order":  {"tier": "metadata", "key": "protocol_randomize_step_order",      "required": False, "source": "GUI"},
    "reset_between_steps":    {"tier": "metadata", "key": "protocol_reset_between_steps",       "required": False, "source": "computed",
                               "note": "forced True for segmented protocols; no GUI control"},
    "hold_motor_between_steps": {"tier": "metadata", "key": "protocol_hold_motor_between_steps", "required": False, "source": "computed",
                                 "note": "forced True for segmented protocols; no GUI control"},
    "motor_enable_dwell_s":   {"tier": "metadata", "key": "protocol_motor_enable_dwell_second", "required": False, "source": "GUI"},
    "reset_max_speed_deg_per_s": {"tier": "metadata", "key": "protocol_reset_max_speed_degree_per_second", "required": False, "source": "GUI"},
    "pulse_width_ms":         {"tier": "metadata", "key": "protocol_pulse_width_millisecond",  "required": False, "source": "GUI"},
    "left_stim_voltage":      {"tier": "metadata", "key": "protocol_stim_voltage_left_volt",   "required": False, "source": "GUI"},
    "right_stim_voltage":     {"tier": "metadata", "key": "protocol_stim_voltage_right_volt",  "required": False, "source": "GUI"},
    "block_sequence":         {"tier": "metadata", "key": "protocol_block_sequence",           "required": False, "source": "GUI",
                               "note": "list-of-dict block plan; serialized as a JSON string (json.dumps) by the exporter"},
    # run-computed protocol provenance (set during the run path; absent -> skipped)
    "daq_collection_type":        {"tier": "metadata", "key": "daq_collection_type",         "required": False, "source": "computed",
                                   "note": "'continuous' (single_finite) or 'segmented' (segmented_finite), derived from protocol"},
    "protocol_sampling_mode":     {"tier": "metadata", "key": "protocol_sampling_mode",      "required": False, "source": "computed",
                                   "note": "'single_finite' (dynamic/sweep) or 'segmented_finite' (isovelocity/FL/FV; isometric after Step 4)"},
    "protocol_guard_triggered":   {"tier": "metadata", "key": "protocol_guard_triggered",   "required": False, "source": "computed",
                                   "note": "bool: isovelocity angle guard fired"},
    "protocol_guard_angle_degree": {"tier": "metadata", "key": "protocol_guard_angle_degree", "required": False, "source": "computed",
                                   "note": "angle (deg) at which the isovelocity guard fired; NaN if not"},

    # --- protocol_ : config-sourced timing defaults ------------------------
    "waitbefore":        {"tier": "metadata", "key": "protocol_wait_before_second",          "required": False, "source": "config"},
    "waitafter":         {"tier": "metadata", "key": "protocol_wait_after_second",           "required": False, "source": "config"},
    "rampdur":           {"tier": "metadata", "key": "protocol_ramp_duration_second",        "required": False, "source": "config"},
    "prepoststim_dur":   {"tier": "metadata", "key": "protocol_prepoststim_duration_second", "required": False, "source": "config"},
    "prepoststim_sep":   {"tier": "metadata", "key": "protocol_prepoststim_separation_second", "required": False, "source": "config"},
    "prestim_time":      {"tier": "metadata", "key": "protocol_prestim_time_second",         "required": False, "source": "config"},
    "poststim_time":     {"tier": "metadata", "key": "protocol_poststim_time_second",        "required": False, "source": "config"},
    "ramp_mode_default": {"tier": "metadata", "key": "protocol_ramp_mode",                   "required": False, "source": "config", "note": "categorical (linear/exponential)"},
    "amp_step_vel":      {"tier": "metadata", "key": "protocol_amplitude_step_velocity_degree_per_second", "required": False, "source": "config"},

    # --- protocol_ : dynamic / frequency_sweep param families --------------
    # value + _unit pattern: the value field holds the numbers; a sibling _unit holds the
    # resolved spelled-out unit (from the *_mode attr at write time); optional _input_mode
    # preserves the raw convention for provenance (decision: triad option B).
    "all_freqs":          {"tier": "metadata", "key": "protocol_step_frequency_hertz", "required": False, "source": "GUI", "note": "per-step array; fixed unit hertz"},
    "all_amps":           {"tier": "metadata", "key": "protocol_step_amplitude",       "required": False, "source": "GUI", "note": "per-step value; unit in protocol_step_amplitude_unit"},
    "all_amps_mode":      {"tier": "metadata", "key": "protocol_step_amplitude_unit",  "required": False, "source": "computed",
                           "note": "resolved spelled unit (degree/percent/per_meter); raw also kept as protocol_step_amplitude_input_mode"},
    "all_curves":         {"tier": "metadata", "key": "protocol_step_curvature_per_meter", "required": False, "source": "computed", "note": "per-step curvature array; fixed unit"},
    "curve_input_values": {"tier": "metadata", "key": "protocol_curve_input",          "required": False, "source": "GUI", "note": "value; unit in protocol_curve_input_unit"},
    "curve_input_mode":   {"tier": "metadata", "key": "protocol_curve_input_unit",     "required": False, "source": "computed",
                           "note": "resolved spelled unit; raw also kept as protocol_curve_input_input_mode"},
    "all_stimduties":     {"tier": "metadata", "key": "protocol_step_stim_duty_fraction",  "required": False, "source": "GUI"},
    "all_stimphases":     {"tier": "metadata", "key": "protocol_step_stim_phase_fraction", "required": False, "source": "GUI", "note": "phase as fraction of cycle"},
    "stim_cycles_in_step": {"tier": "metadata", "key": "protocol_stim_cycles_in_step", "required": False, "source": "GUI", "note": "count"},
    "stim_pulse_rate":    {"tier": "metadata", "key": "protocol_stim_pulse_rate_hertz", "required": False, "source": "GUI"},
    "cycles_per_step":    {"tier": "metadata", "key": "protocol_cycles_per_step",    "required": False, "source": "GUI", "note": "count"},
    "n_end_cycles":       {"tier": "metadata", "key": "protocol_end_cycles_count",   "required": False, "source": "GUI", "note": "count"},
    "randomize":          {"tier": "metadata", "key": "protocol_randomize",          "required": False, "source": "GUI", "note": "legacy bool; mirrors randomize_step_order"},
    "random_seed":        {"tier": "metadata", "key": "protocol_random_seed",        "required": False, "source": "GUI"},
    "sweep_nominal_frequency": {"tier": "metadata", "key": "protocol_sweep_nominal_frequency_hertz",  "required": False, "source": "GUI"},
    "sweep_nominal_curvature": {"tier": "metadata", "key": "protocol_sweep_nominal_curvature_per_meter", "required": False, "source": "GUI"},

    # --- protocol_ : isometric family --------------------------------------
    # initial/final/mirror targets are VALUES sharing one unit: protocol_isometric_target_unit
    "isometric_initial":   {"tier": "metadata", "key": "protocol_isometric_initial",   "required": False, "source": "GUI", "note": "value; unit in protocol_isometric_target_unit"},
    "isometric_final":     {"tier": "metadata", "key": "protocol_isometric_final",     "required": False, "source": "GUI", "note": "value; unit in protocol_isometric_target_unit"},
    "isometric_num_steps": {"tier": "metadata", "key": "protocol_isometric_num_steps", "required": False, "source": "GUI", "note": "count"},
    "isometric_mode":      {"tier": "metadata", "key": "protocol_isometric_target_unit", "required": False, "source": "computed",
                            "note": "resolved spelled unit (percent/degree) for initial/final/mirror; raw also kept as protocol_isometric_target_input_mode"},
    "isometric_randomize": {"tier": "metadata", "key": "protocol_isometric_randomize", "required": False, "source": "GUI"},
    "isometric_random_seed": {"tier": "metadata", "key": "protocol_isometric_random_seed", "required": False, "source": "GUI"},
    "isometric_inter_step_interval_s": {"tier": "metadata", "key": "protocol_isometric_inter_step_interval_second", "required": False, "source": "GUI", "note": "legacy; mirrors rest_between_steps_s"},
    "isometric_stim_params":    {"tier": "metadata", "key": "protocol_isometric_stim_params",    "required": False, "source": "GUI", "note": "dict"},
    "isometric_stim_overrides": {"tier": "metadata", "key": "protocol_isometric_stim_overrides", "required": False, "source": "GUI", "note": "dormant template back-compat"},

    # --- protocol_ : isovelocity family ------------------------------------
    # min/max velocity are VALUES sharing protocol_isovelocity_velocity_unit; starting_strain
    # is a VALUE with protocol_isovelocity_starting_strain_unit.
    "isovelocity_min_vel":        {"tier": "metadata", "key": "protocol_isovelocity_min_velocity",  "required": False, "source": "GUI", "note": "value; unit in protocol_isovelocity_velocity_unit"},
    "isovelocity_max_vel":        {"tier": "metadata", "key": "protocol_isovelocity_max_velocity",  "required": False, "source": "GUI", "note": "value; unit in protocol_isovelocity_velocity_unit"},
    "isovelocity_velocity_mode":  {"tier": "metadata", "key": "protocol_isovelocity_velocity_unit", "required": False, "source": "computed",
                                   "note": "resolved spelled unit (degree_per_second/percent_per_second); raw also kept as protocol_isovelocity_velocity_input_mode"},
    "isovelocity_starting_strain": {"tier": "metadata", "key": "protocol_isovelocity_starting_strain", "required": False, "source": "GUI", "note": "value; unit in protocol_isovelocity_starting_strain_unit"},
    "isovelocity_starting_strain_mode": {"tier": "metadata", "key": "protocol_isovelocity_starting_strain_unit", "required": False, "source": "computed",
                                         "note": "resolved spelled unit (percent); raw also kept as protocol_isovelocity_starting_strain_input_mode"},
    "isovelocity_num_steps":      {"tier": "metadata", "key": "protocol_isovelocity_num_steps",      "required": False, "source": "GUI", "note": "count"},
    "isovelocity_randomize":      {"tier": "metadata", "key": "protocol_isovelocity_randomize",      "required": False, "source": "GUI"},
    "isovelocity_random_seed":    {"tier": "metadata", "key": "protocol_isovelocity_random_seed",    "required": False, "source": "GUI"},
    "isovelocity_iso_duration_s": {"tier": "metadata", "key": "protocol_isovelocity_iso_duration_second", "required": False, "source": "GUI"},
    "isovelocity_pre_hold_s":     {"tier": "metadata", "key": "protocol_isovelocity_pre_hold_second",     "required": False, "source": "GUI"},
    "isovelocity_stim_params":    {"tier": "metadata", "key": "protocol_isovelocity_stim_params",    "required": False, "source": "GUI", "note": "dict"},
    "isovelocity_stim_overrides": {"tier": "metadata", "key": "protocol_isovelocity_stim_overrides", "required": False, "source": "GUI", "note": "dormant template back-compat"},

    # === metadata / schema version =========================================
    "h5_schema_version": {"tier": "metadata", "key": "schema_version", "required": True, "source": "computed"},

    # === timeseries (continuous, sample-indexed) ===========================
    "t":              {"tier": "timeseries", "key": "time_second",                                   "required": True,  "source": "DAQ"},
    "tnorm":          {"tier": "timeseries", "key": "time_normalized",                               "required": False, "source": "computed"},
    "angle":          {"tier": "timeseries", "key": "angle_commanded_degree",                        "required": True,  "source": "computed"},
    "angle_measured": {"tier": "timeseries", "key": "angle_measured_degree",                         "required": True,  "source": "DAQ"},
    "anglevel":       {"tier": "timeseries", "key": "angular_velocity_commanded_degree_per_second",  "required": False, "source": "computed"},
    "S1stimcmd":      {"tier": "timeseries", "key": "stim_channel1_command_volt",                    "required": False, "source": "computed"},
    "S2stimcmd":      {"tier": "timeseries", "key": "stim_channel2_command_volt",                    "required": False, "source": "computed"},
    "aidata":         {"tier": "timeseries", "key": "aidata",                                        "required": True,  "source": "DAQ",
                       "note": "IMMUTABLE raw 8xN; channel identity in daq_ai_channel_map"},
    "forcetorque_raw": {"tier": "timeseries", "key": "forcetorque_raw",                              "required": True,  "source": "computed",
                        "note": "IMMUTABLE decoded 6xN; NOT split (PI decision). Channel order: "
                                "x_force, y_force, z_force, x_torque, y_torque, z_torque "
                                "(newton / newton_meter). R decodes/names channels downstream."},
    "forcetorque":    {"tier": "timeseries", "key": "forcetorque_raw",
                       "required": False, "source": "computed",
                       "note": "live decoded working array; forcetorque_raw is a copy=True of it "
                               "(bender_functions L1080). Written once as forcetorque_raw; not split."},
    "sono_left_mm":   {"tier": "timeseries", "key": "sono_left_millimeter",   "required": False, "source": "computed", "note": "decoded sono length (raw in aidata)"},
    "sono_right_mm":  {"tier": "timeseries", "key": "sono_right_millimeter",  "required": False, "source": "computed", "note": "decoded sono length (raw in aidata)"},
    "cycle_index_history":   {"tier": "timeseries", "key": "cycle_index",       "required": False, "source": "computed",
                              "note": "per-sample cycle number; -1 = not a numbered cycle (step protocols)"},
    "sweep_instantaneous_freq": {"tier": "timeseries", "key": "instantaneous_frequency_hertz", "required": False, "source": "computed",
                                 "note": "frequency_sweep only"},
    # Per-sample stim streams (PI decision). stim_state/stim_side/stim_type kept per-sample as-is.
    "stim_side":   {"tier": "timeseries", "key": "stim_side",   "required": False, "source": "computed", "note": "per-sample categorical: none/left/right/both"},
    "stim_state":  {"tier": "timeseries", "key": "stim_state",  "required": False, "source": "computed", "note": "per-sample stim phase"},
    "stim_type":   {"tier": "timeseries", "key": "stim_type",   "required": False, "source": "computed", "note": "per-sample activity (active/passive); demote pending PI stim_type flag"},

    # Per-sample index family (NOT Bender attributes -- built per-step inside the protocol helpers
    # and broadcast to length N by the exporter; listed here for documentation only):
    #   step_index, sequence_index, block_index (int);  block_direction, block_stim_sides (str).
    # trial_index is DROPPED (redundant with filename + sequence_index). step_order is DROPPED
    # (reconstructable from per-sample sequence_index + block_index).
}


# ---------------------------------------------------------------------------
# EXCLUDED -- never written to the raw file (fail-on-unmapped must not trip)
# ---------------------------------------------------------------------------
EXCLUDED: Set[str] = {
    # --- internal caches / motor-state trackers (also all _-prefixed attrs) ---
    "_config_daq_ai_sample_rate_hz", "_cfg_input_channels", "_cfg_input_channel_names",
    "_daq_ai_sample_rate_hz", "_last_commanded_angle_deg", "_motor_continuous_step_pos",
    "_last_motor_direction_bit", "_encoder_cumulative_deg", "_reset_encoder_cumulative_deg",
    "_motor_enable_power_up_high", "_motor_driver_energized",

    # --- objects / transport containers (contents routed individually) ---
    "master_logger",            # MasterLogger object
    "h5_protocol_metadata",     # transport dict; its keys are routed as individual fields
    "trial_records",            # raw per-step buffers -> become timeseries + index_step_*

    # --- machine-specific ---
    "outputfile",               # absolute path; portable `filename` written separately
    "timestamp",                # compact run stamp; start_time_iso derived from it

    # --- DAQ scratch buffers ---
    "stimcmdhi", "dig",

    # --- derived series -> research-hub `derived` only (schema §3d / §6) ---
    "primary_torque_raw",       # schema §6: delete from metadata; recomputed in R
    "angledata",                # duplicate of angle_measured (raw encoder buffer)
    "is_stim_cycle",            # transient per-cycle helper

    # --- legacy aliases (mirrored onto canonical names; would duplicate) ---
    "initial", "final", "num_steps", "mode",
    "min_vel", "max_vel", "starting_strain", "starting_strain_mode",
    "test_segment_length_mm",   # alias of dclamp -> measurement_clamp_separation_millimeter (see CONFLICT note in decision log)
    "total_mass",               # MC #2: not a schema field; recomputed in R
    "specimen_mass_specimen",   # run-computed geometry-model mass; not a schema field
    "specimen_mass_frustum",    # run-computed mass from legacy frustum builder; not a schema field
    "specimen_mass_profile",    # run-computed mass from legacy profiled-stations builder; not a schema field
    "specimen_volume_mm3",      # run-computed geometry-model volume; not a schema field
    "specimen_inertial_model",  # internal model tag (elliptical_frustum/profiled_stations/...); provenance downstream

    # --- inertial flags dropped from schema (D12 / Points 1-3) ---
    # use_frustum_inertial_model: notebook auto-build switch only; not a schema field.
    # use_theoretical_inertial_correction: gates nothing in current code; dropped.
    # use_inertial_calibration: controlled calibration_link (dropped per D3/D12).
    # inertial_calibration_profile: handled by direct flat-key block in exporter (not ledger).
    # inertial_system_from_profile: subsumed; dropped per D12.
    "use_frustum_inertial_model", "use_theoretical_inertial_correction",
    "use_inertial_calibration", "inertial_calibration_profile",
    "inertial_system_from_profile",

    # --- legacy profile / frustum model attrs (unused GUI path; not schema fields) ---
    # specimen_profile_* belong to set_profiled_specimen_inertial_model (not the live GUI path).
    # frustum_inputs: the dict object itself is excluded — its individual scalar contents are now
    # routed as first-class measurement_specimen_frustum_* fields (v2.7); the dict is redundant.
    # specimen_profile_num_samples: dropped per Point 1 (no sample count in active code path).
    "frustum_inputs",
    "specimen_profile_stations", "specimen_profile_density_g_per_mm3",
    "specimen_profile_num_samples",

    # --- legacy frustum scalar inputs ---
    # frustum_height/width/length/density are now ROUTED to measurement_specimen_inertia_frustum_* (v2.7 amended).
    # frustum_clamp_offset_mm is ROUTED to calibration_inertia_apparatus_aor_to_clamp_millimeter (many-to-one with
    #   specimen_profile_clamp_offset_mm; apparatus geometry, not specimen frustum geometry).
    # frustum_tip_scale: DROPPED (v2.7) — redundant when N-station geometry captures each cross-section explicitly.
    # frustum_num_samples: DROPPED — no sample count in the active code path (exact per-segment integration).
    "frustum_tip_scale", "frustum_num_samples",

    # --- run-computed scheduling / per-cycle intermediates (DERIVED, recomputable) ---
    # Built during run_experiment from the routed protocol params (all_freqs/all_amps/etc.) to
    # assemble the motor timeline. Not raw data and not schema fields; recomputed downstream.
    "amp_by_cycle", "freq_by_cycle", "phase_by_cycle", "duty_by_cycle",
    "period_by_cycle", "step_by_cycle", "strain_by_cycle", "strainrate_by_cycle",
    "stimburstdur", "Lonoff", "Ronoff",
    "organized_freqs", "organized_curves", "organized_strains", "organized_strainrates",
    "organized_stimduties", "organized_stimphases",
    "all_degs", "all_strains", "all_strainrates",
    "tout", "endTime", "max_commanded_rotation_deg",
    # self.filename: duplicate of the portable name written to 01_Metadata/filename by exporter.
    "filename",
    # stim_monitor: redundant slice of the immutable raw aidata buffer (raw lives in aidata).
    "stim_monitor",
    # derived per-step analysis result objects (recomputed downstream from raw + params).
    "force_length_results", "isovelocity_results",

    # --- legacy transport / convenience containers (unpacked before use) ---
    # protocol_params is a notebook/GUI convenience dict consumed by update_metadata();
    # its contents are unpacked into individual canonical attrs, so the dict itself
    # must never be written to the schema (would double-write every field it contains).
    "protocol_params",

    # --- legacy R-labeling dicts (config-file-only; PI-dropped as redundant) ---
    # Per PI decision: superseded by spelled-out unit suffixes in canonical keys.
    # These are NOT Bender attrs (they live only in the config module), so they never
    # trip fail-on-unmapped; listed here only to document the intentional drop.
    "units", "unit_rules",

    # NOTE: the former hardware-config provenance block (device_name, motor_port,
    # motor_gear_ratio, daq/encoder/stim wiring, sono channels, timing buffers, etc.)
    # is no longer excluded -- every config-sourced field is now routed canonically
    # under daq_ / calibration_ / protocol_ in BENDER_ROUTING above.
}


# ---------------------------------------------------------------------------
# MISSING_REQUIRED -- schema fields with no Bender source today
# ---------------------------------------------------------------------------
MISSING_REQUIRED: List[Dict[str, str]] = [
    {"key": "measurement_specimen_body_width_millimeter", "tier": "metadata",
     "source": "GUI", "note": "Test-section body width; no Bender attr. Needs a GUI widget (distinct from xsec_width)."},
    {"key": "note_posthoc", "tier": "metadata",
     "source": "analysis", "note": "Populated post-hoc in the pipeline; expected absent at acquisition unless post_trial_notes is routed here."},
    # The three flat profile keys are written by the direct block in bender_h5_export.py
    # (lines that replaced the inertial_calibration_profile subgroup writer), NOT by the
    # ledger pass.  No Bender attribute maps 1:1 to these keys; they unpack inertial_calibration_profile.
    {"key": "calibration_inertia_apparatus_moi_gram_millimeter_squared", "tier": "metadata",
     "source": "computed",
     "note": "I_est from inertial_calibration_profile, converted *(180/pi)*1e9 at write time; "
             "written by the direct flat-key block in the exporter, not the ledger."},
    {"key": "calibration_inertia_bias_newton_meter", "tier": "metadata",
     "source": "computed",
     "note": "bias_est from inertial_calibration_profile (N*m, no conversion); "
             "written by the direct flat-key block in the exporter, not the ledger."},
    {"key": "calibration_inertia_axis_sensor", "tier": "metadata",
     "source": "computed",
     "note": "axis_sensor from inertial_calibration_profile (categorical); "
             "written by the direct flat-key block in the exporter, not the ledger."},
    {"key": "index_step_*", "tier": "metadata",
     "source": "computed", "note": "Parallel per-step arrays (sample_start/end/type/target/...); derived from trial_records step boundaries."},
]


# ---------------------------------------------------------------------------
# Coverage check -- for the future fail-on-unmapped exporter
# ---------------------------------------------------------------------------
def check_coverage(bender: Any) -> Dict[str, List[str]]:
    """Classify a Bender instance's public attributes against this spec.

    Returns a dict with:
      'routed'    : public attrs found in BENDER_ROUTING
      'excluded'  : public attrs found in EXCLUDED
      'unmapped'  : public attrs in NEITHER (the exporter must raise on these)

    Underscore-prefixed attrs are ignored (implicitly excluded). Callable
    attributes are ignored (methods).
    """
    routed: List[str] = []
    excluded: List[str] = []
    unmapped: List[str] = []
    for name in sorted(vars(bender)):
        if name.startswith("_"):
            continue
        if callable(getattr(bender, name, None)):
            continue
        if name in BENDER_ROUTING:
            routed.append(name)
        elif name in EXCLUDED:
            excluded.append(name)
        else:
            unmapped.append(name)
    return {"routed": routed, "excluded": excluded, "unmapped": unmapped}


def assert_full_coverage(bender: Any) -> None:
    """Raise ``RuntimeError`` if any public attribute is unmapped (the contract)."""
    result = check_coverage(bender)
    if result["unmapped"]:
        raise RuntimeError(
            "bender_routing_spec: unmapped Bender attributes (add to BENDER_ROUTING "
            f"or EXCLUDED): {result['unmapped']}"
        )


__all__ = ["BENDER_ROUTING", "EXCLUDED", "MISSING_REQUIRED",
           "check_coverage", "assert_full_coverage", "Route"]
