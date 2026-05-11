# Mirrors ``bender_streamlit_gui`` motion / help constants (no Streamlit).
from __future__ import annotations

MOTION_TYPES = frozenset({'dynamic', 'frequency_sweep', 'frequency_step', 'curvature_step', 'step_change'})

MOTION_GUI_FIELDS = [
    ('duration', 'float', 'Duration (s)'),
    ('all_freqs', 'list_float', 'Frequencies Hz (comma-separated, e.g. 1, 5)'),
    (
        'all_amps',
        'list_float',
        'Amplitudes (comma-separated; see all_amps_mode: strain=decimal ε, strain_pct=percent)',
    ),
    ('all_amps_mode', 'select', 'Amplitude interpretation (strain vs strain_pct, etc.)'),
    ('cycles_per_step', 'int', 'Cycles per step'),
    ('n_end_cycles', 'int', 'End cycles'),
    ('randomize', 'bool', 'Randomize order'),
    ('random_seed', 'optional_int', 'Random seed (empty = None)'),
    ('stim_cycles_in_step', 'list_int', 'Stim cycle indices per step (comma-separated, e.g. 2, 3)'),
    ('is_stim', 'bool', 'Stimulation enabled'),
    ('stim_pulse_rate', 'float', 'Stim pulse rate (Hz)'),
    ('S1volts', 'float', 'S1 stim voltage'),
    ('S2volts', 'float', 'S2 stim voltage'),
]

STEP_CHANGE_EXTRA = [
    ('step_change_frequencies', 'list_float', 'step_change frequencies (comma-separated)'),
    ('step_change_curves', 'list_float', 'step_change curves / amplitudes (comma-separated)'),
    ('step_change_cycles_per_step', 'list_int', 'cycles per step (comma-separated)'),
]

FREQUENCY_SWEEP_ONLY_FIELDS = [
    (
        'amplitude_frequency_exponent',
        'float',
        'Amplitude–frequency exponent α (command θ ∝ f^α vs sweep start)',
    ),
]

_MOTION_ROW_BY_NAME = {row[0]: row for row in MOTION_GUI_FIELDS}


def motion_parameter_rows(test_type: str):
    freq_amp_mode = [
        _MOTION_ROW_BY_NAME['all_freqs'],
        _MOTION_ROW_BY_NAME['all_amps'],
        _MOTION_ROW_BY_NAME['all_amps_mode'],
    ]
    organize_block = [
        _MOTION_ROW_BY_NAME['cycles_per_step'],
        _MOTION_ROW_BY_NAME['n_end_cycles'],
        _MOTION_ROW_BY_NAME['randomize'],
        _MOTION_ROW_BY_NAME['random_seed'],
        _MOTION_ROW_BY_NAME['stim_cycles_in_step'],
    ]
    stim_block = [
        _MOTION_ROW_BY_NAME['is_stim'],
        _MOTION_ROW_BY_NAME['stim_pulse_rate'],
        _MOTION_ROW_BY_NAME['S1volts'],
        _MOTION_ROW_BY_NAME['S2volts'],
    ]
    if test_type == 'dynamic':
        return [*freq_amp_mode, *organize_block, *stim_block]
    if test_type == 'frequency_sweep':
        return [
            _MOTION_ROW_BY_NAME['duration'],
            *freq_amp_mode,
            *FREQUENCY_SWEEP_ONLY_FIELDS,
            *stim_block,
        ]
    if test_type in ('frequency_step', 'curvature_step'):
        return [_MOTION_ROW_BY_NAME['duration'], *freq_amp_mode, *stim_block]
    if test_type == 'step_change':
        return [*list(STEP_CHANGE_EXTRA), *stim_block]
    return list(MOTION_GUI_FIELDS)


ALL_AMPS_MODE_OPTIONS = ('strain', 'strain_pct', 'curvature', 'angle')

DATA_FOLDER_HELP = (
    'Choose the folder where experiment files live. Enter the folder path only — do not put the file name here. '
    'Runs, exports, protocol templates, and biometrics files can all use this folder.'
)
DATA_FILE_NAME_HELP = (
    'This is the name of your saved measurements file (HDF5). Enter only the file name, not the full path. '
    'The app joins it with the data folder field above to build where data is saved. You may type .h5 or leave it off. '
    'If that exact file already exists, the app uses a new name like my_run_001.h5 so nothing is overwritten.'
)

BIO_DENSITY_PRESET_LABELS = (
    'Custom — edit the number below',
    'Water-like (~1.00 g/cm³)',
    'Skeletal muscle / soft tissue (~1.06 g/cm³)',
    'Cortical bone (~1.9 g/cm³)',
)
BIO_DENSITY_PRESET_G_PER_MM3 = {
    'Water-like (~1.00 g/cm³)': 1.0e-3,
    'Skeletal muscle / soft tissue (~1.06 g/cm³)': 1.06e-3,
    'Cortical bone (~1.9 g/cm³)': 1.9e-3,
}

BIO_DBEND_FIELD_HELP = (
    'Distance **along the body** (mm) from your length reference (same as TL/SL, often snout or a fixed landmark) to the '
    '**midpoint between the two clamps** — i.e. where the bending test is centered. Use **0** only if that reference is '
    'already at the segment center. Saved as `dbend` / `test_segment_position_mm`.'
)

BIO_PROF_CLAMP_FIELD_HELP = (
    'Used only for **rotating hardware** mass/MOI in the profiled inertial model: offset from the **bend / rotation axis** '
    'to the clamps (mm). The code adds half of its built-in clamp depth to this value when estimating clamp contribution. '
    'Saved as `specimen_profile_clamp_offset_mm`.'
)

ISOMETRIC_STIM_JSON_HELP = (
    'Optional per-step timing and stimulation. **Leave `{}`** unless your protocol needs custom ramps, holds, or stim.'
)
ISOMETRIC_STIM_OVERRIDES_HELP = '**Rare.** Overrides stim **routing** (not timing). **Leave `{}`** unless directed.'
ISOVELOCITY_STIM_JSON_HELP = 'Optional segment timing and stimulation. **Leave `{}`** unless you need overrides.'
ISOVELOCITY_STIM_OVERRIDES_HELP = ISOMETRIC_STIM_OVERRIDES_HELP

RANDOM_SEED_HELP = (
    'Only used when **randomize order** is checked. Enter an integer for a **reproducible** shuffle; leave empty for random.'
)

RECRUITMENT_FIELD_HELP = (
    'How stimulation and/or motor commands are routed across left vs right. '
    'With **Perform test on both sides**, simultaneous or unilateral modes are upgraded to **bilateral sequential** '
    'at run time when needed.'
)

ISOVELOCITY_WIDGET_LABEL = {
    'isovelocity_min_vel': 'Minimum angular velocity (deg/s)',
    'isovelocity_max_vel': 'Maximum angular velocity (deg/s)',
    'isovelocity_starting_strain': 'Starting posture (strain / curvature / angle)',
    'isovelocity_num_steps': 'Number of velocity steps',
    'isovelocity_starting_strain_mode': 'Unit for starting posture value',
    'isovelocity_randomize': 'Randomize order of velocity steps',
    'isovelocity_random_seed': 'Random seed (optional)',
    'isovelocity_iso_duration_s': 'Constant-velocity bend duration (s)',
    'isovelocity_pre_hold_s': 'Pre-hold at starting angle (s)',
}

ISOVELOCITY_FIELD_HELP = {
    'isovelocity_min_vel': 'Lower end of the angular velocity sweep (deg/s).',
    'isovelocity_max_vel': 'Upper end of the angular velocity sweep (deg/s).',
    'isovelocity_starting_strain': 'Initial posture before each velocity step.',
    'isovelocity_num_steps': 'How many commanded angular velocities between min and max (inclusive).',
    'isovelocity_starting_strain_mode': 'Whether **starting posture** is decimal ε, percent strain, curvature κ, or motor angle (deg).',
    'isovelocity_randomize': 'Shuffle the order of velocity steps.',
    'isovelocity_random_seed': RANDOM_SEED_HELP,
    'isovelocity_iso_duration_s': 'After the quiet **pre-hold**, how long (seconds) the motor holds **constant angular velocity**.',
    'isovelocity_pre_hold_s': 'Time (s) at the starting posture before each trial’s constant-velocity segment begins.',
    'recruitment': RECRUITMENT_FIELD_HELP,
    'lateral_mode': 'Expert only. Leave **blank** unless you need a custom stim-router label.',
    'bilateral_mirror_motor': 'Two constant-velocity bends per trial (left then right) with stim on each side; recruitment coerced when needed.',
    'bilateral_sequential_left_frac': 'Share of each step spent on the **left** side before the right (0–1).',
}

ISOMETRIC_FIELD_HELP = {
    'isometric_initial': 'First target in the sweep; units follow **isometric mode**.',
    'isometric_final': 'Last target in the sweep (same units as initial).',
    'isometric_num_steps': 'Number of steps between initial and final (endpoints included).',
    'isometric_mode': 'How **initial** / **final** are interpreted: decimal ε, percent strain, κ, or motor angle.',
    'isometric_randomize': 'Shuffle step order. Use **random seed** for a fixed order across runs.',
    'isometric_random_seed': RANDOM_SEED_HELP,
    'isometric_inter_step_interval_s': 'Idle time after each step’s acquisition before the next ramp (0 = back-to-back).',
    'recruitment': RECRUITMENT_FIELD_HELP,
    'lateral_mode': 'Expert only. Leave **blank** unless you need a custom stim-router label.',
    'bilateral_mirror_motor': 'Ramp/hold toward left then right with matched stimulation; recruitment coerced to sequential when needed.',
    'isometric_mirror_target_left': 'Optional: explicit first-hold magnitude toward LEFT (same units as isometric mode); enable checkbox and set both left and right.',
    'isometric_mirror_target_right': 'Optional: second-hold magnitude toward RIGHT; enable checkbox; use with left target.',
    'bilateral_sequential_left_frac': 'Share of each step on the **left** side before the right (0–1).',
}

MOTION_FIELD_HELP = {
    'random_seed': RANDOM_SEED_HELP,
    'duration': 'Total motion timeline length (s) for protocols that use a single continuous run.',
    'all_freqs': 'Comma-separated list of frequencies (Hz) for each step or segment.',
    'all_amps': 'Comma-separated amplitudes; meaning depends on **amplitude interpretation**.',
    'all_amps_mode': 'How each entry in **amplitudes** is converted to curvature / motor angle.',
    'cycles_per_step': 'Oscillation cycles at each frequency step before advancing.',
    'n_end_cycles': 'Extra cycles at the last frequency at the end of the run.',
    'randomize': 'Shuffle the order of frequency/amplitude steps when applicable.',
    'stim_cycles_in_step': 'Which cycle indices (1-based) within each step receive stimulation.',
    'is_stim': 'Enable patterned stimulation during the motion.',
    'stim_pulse_rate': 'Carrier pulse rate for stimulation (Hz).',
    'S1volts': 'Stimulus channel 1 amplitude (V).',
    'S2volts': 'Stimulus channel 2 amplitude (V).',
    'amplitude_frequency_exponent': 'For frequency sweep: exponent α so amplitude ∝ f^α relative to the sweep start.',
    'step_change_frequencies': 'Comma-separated frequencies for each step-change segment.',
    'step_change_curves': 'Comma-separated curvature / amplitude targets per segment.',
    'step_change_cycles_per_step': 'Comma-separated integers: cycles per step-change segment.',
}


def format_strain_or_amp_mode(opt: str) -> str:
    o = str(opt)
    if o == 'strain':
        return 'strain — decimal ε (0.05 = 5%)'
    if o == 'strain_pct':
        return 'strain_pct — percent scale (5 = 5%)'
    if o == 'curvature':
        return 'curvature — κ (1/m)'
    if o == 'angle':
        return 'angle — motor (deg)'
    return o


RECRUITMENT_OPTIONS = (
    'bilateral_simultaneous',
    'bilateral_sequential',
    'left',
    'right',
)

BILATERAL_MIRROR_LABEL = 'Perform test on both sides (bilateral)'
LATERAL_MODE_LABEL = 'Stim routing override (optional; experts only)'
