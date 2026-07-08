#!/usr/bin/env python3
"""Fit apparatus moment-of-inertia I(aor, width) from empty-apparatus calibration H5 files.

Standalone CLI wrapper around ``fit_apparatus_inertia_calibration`` in ``bender_functions``.
It fits the five candidate forms (F1-F5), selects one, prints a full report, and writes the
versioned artifact JSON that the GUI/exporter later load.

Examples
--------
Fit every ``*.h5`` in a folder, supplying real aor values from a JSON map, excluding a bad file::

    python fit_apparatus_inertia.py /path/to/inertial_calibration_apparatus \\
        --aor-overrides aor.json --exclude _09_ \\
        --out 2026-07-06_apparatus_inertia_calibration.json

Fit an explicit file list::

    python fit_apparatus_inertia.py --files a.h5 b.h5 c.h5 \\
        --out 2026-07-06_apparatus_inertia_calibration.json

Notes
-----
- ``--out`` is REQUIRED and has no generic default: each calibration must be named deliberately
  (there is no auto-generated ``apparatus_inertia_fit.json``). Give it the SAME name the hardware
  config's ``apparatus_inertia_calibration_file`` expects so the run auto-loads it. A date prefix
  (``2026-07-07_...``) records when the empty-apparatus runs were taken.
- ``--aor-overrides`` is a JSON object ``{basename: aor_millimeter}`` used when a file's
  ``calibration_inertia_apparatus_aor_to_clamp_millimeter`` is NaN (the routing bug case).
- ``--exclude`` values are matched as SUBSTRINGS of each file's basename (so ``_09_`` drops
  ``..._bender_09_dynamic.h5``). Excluded files are still reported, never silently dropped.
- All output is plain ASCII (Windows cp1252 safe); prefixes: [info], [warn], [error].
"""
import argparse
import glob
import json
import os
import sys

# stdout/stderr to UTF-8 with replacement as defense-in-depth (matches GUI startup policy).
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        pass

import numpy as np

from bender_functions import fit_apparatus_inertia_calibration, _APPARATUS_FIT_FORMS


def _resolve_files(args):
    """Return a sorted list of H5 paths from either --files or a folder positional."""
    if args.files:
        paths = list(args.files)
    elif args.path:
        if os.path.isdir(args.path):
            paths = sorted(glob.glob(os.path.join(args.path, '*.h5')))
        elif os.path.isfile(args.path):
            paths = [args.path]
        else:
            raise SystemExit('[error] path is neither a folder nor a file: ' + str(args.path))
    else:
        raise SystemExit('[error] provide a folder/file positional path or --files.')
    if not paths:
        raise SystemExit('[error] no .h5 files found to fit.')
    return paths


def _resolve_excludes(paths, exclude_substrings):
    """Map --exclude substrings to the set of matching basenames present in ``paths``."""
    bases = [os.path.basename(p) for p in paths]
    excluded = set()
    for sub in (exclude_substrings or []):
        for b in bases:
            if sub in b:
                excluded.add(b)
    return excluded


def _fmt(x, nd=3):
    """Format a float compactly; ASCII 'nan' for non-finite."""
    try:
        xf = float(x)
    except (TypeError, ValueError):
        return str(x)
    if not np.isfinite(xf):
        return 'nan'
    return ('{:.' + str(nd) + 'g}').format(xf)


def _print_report(art):
    """Print the full human-readable report to stdout (plain ASCII)."""
    print('=' * 78)
    print('APPARATUS INERTIA FIT REPORT')
    print('=' * 78)
    print('[info] build_date    :', art.get('build_date', ''))
    print('[info] source_files  :', ', '.join(art.get('source_files', [])) or '(none)')
    print('[info] excluded      :', ', '.join(art.get('excluded_files', [])) or '(none)')
    print('[info] aor_provenance:', art.get('aor_provenance', 'unknown'),
          '(override = aor taken from lab notes, not recorded in the file)')

    if art.get('blocked'):
        print('[warn] FIT BLOCKED   :', art['blocked'])

    # Per-trial table.
    print('\n-- per-trial ---------------------------------------------------------------')
    print('{:>14} {:>5} {:>7} {:>13} {:>6} {:>5} {:>26}'.format(
        'file', 'aor', 'width', 'I(g*mm^2)', 'R2', 'neg', 'alpha_source'))
    excluded = set(art.get('excluded_files', []))
    outlier_files = {o['file'] for o in art.get('outliers', [])}
    for t in art.get('trials', []):
        if 'error' in t:
            print('{:>14} {:>5} {:>7} {:>13} {:>6} {:>5} {:>26}'.format(
                t.get('file', '?')[-14:], '-', '-', 'EXTRACT_ERR', '-', '-', t['error'][:26]))
            continue
        tag = t['file'][-14:]
        flags = []
        if t['file'] in excluded:
            flags.append('EXCL')
        if t['file'] in outlier_files:
            flags.append('OUTLIER')
        print('{:>14} {:>5} {:>7} {:>13} {:>6} {:>5} {:>26}  {}'.format(
            tag, _fmt(t['aor_millimeter'], 3), _fmt(t['width_millimeter'], 4),
            _fmt(t['i_gram_millimeter_squared'], 6), _fmt(t['r2'], 3),
            'yes' if t.get('sign_negative') else 'no', t.get('alpha_source', ''),
            ' '.join(flags)))

    # All candidate forms.
    print('\n-- candidate forms ---------------------------------------------------------')
    print('{:>4} {:>8} {:>16} {:>14} {:>6}  {}'.format(
        'form', 'R2', 'LOO(g*mm^2)', 'LOO(N*m)', 'valid', 'coefficients'))
    for fid in ('F1', 'F2', 'F3', 'F4', 'F5'):
        c = art.get('candidate_forms', {}).get(fid)
        if not c:
            print('{:>4} {:>8}'.format(fid, '(skipped: too few points)'))
            continue
        coef = {k: round(v, 4) for k, v in c['coefficients'].items()}
        print('{:>4} {:>8} {:>16} {:>14} {:>6}  {}'.format(
            fid, _fmt(c['r2'], 4), _fmt(c['loo_cv_rmse_gram_millimeter_squared'], 6),
            _fmt(c['loo_cv_rmse_newton_meter'], 4),
            'yes' if c['physically_valid'] else 'NO', coef))

    # Selection + geometry check.
    print('\n-- selected form -----------------------------------------------------------')
    sel = art.get('fit_form_id', '')
    if sel:
        print('[info] form          :', sel, '->', _APPARATUS_FIT_FORMS.get(sel, art.get('fit_form', '')))
        print('[info] coefficients  :', {k: round(v, 4) for k, v in art.get('fit_coefficients', {}).items()})
        print('[info] LOO-RMSE      :', _fmt(art.get('loo_cv_rmse_newton_meter'), 4), 'N*m',
              '(' + _fmt(art.get('loo_cv_rmse_gram_millimeter_squared'), 6), 'g*mm^2)')
        print('[info] valid_domain  :', art.get('valid_domain'))
    else:
        print('[warn] no form selected.')

    gc = art.get('geometry_check')
    if gc:
        print('\n-- geometry check (F4 aor^2:width^2 ratio vs expected 4:1) ------------------')
        print('[info] F4 ratio      :', _fmt(gc.get('f4_aor2_over_width2_ratio'), 4),
              '(expected', _fmt(gc.get('expected_ratio'), 3) + ')')
        print('[info] F4 beats F5 by >10%:', gc.get('f4_beats_f5_by_10pct'))
        if gc.get('anomaly'):
            print('[warn] GEOMETRY ANOMALY: F4 clearly beats F5 AND the aor^2:width^2 ratio is far')
            print('[warn]   from 4:1. The single-mass Pythagorean model does not hold -- aor and')
            print('[warn]   width likely move DIFFERENT masses. Use the separable form; investigate')
            print('[warn]   the physical mass distribution before over-interpreting the geometry.')

    # Sign + issues.
    print('\n-- sign convention ---------------------------------------------------------')
    print('[info]', art.get('sign_convention_note', ''))

    outliers = art.get('outliers', [])
    if outliers:
        print('\n[warn] OUTLIERS (reported, not silently dropped):')
        for o in outliers:
            print('[warn]  ', o)
    bad = art.get('bad_fit_files', [])
    if bad:
        print('[warn] BAD-FIT FILES (per-trial R2 below floor):', bad)
    errs = art.get('extraction_errors', [])
    if errs:
        print('[warn] EXTRACTION ERRORS:')
        for e in errs:
            print('[warn]  ', e)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description='Fit apparatus MOI I(aor, width) from empty-apparatus calibration H5 files.')
    ap.add_argument('path', nargs='?', default=None,
                    help='Folder of calibration .h5 files (or a single .h5). Omit if using --files.')
    ap.add_argument('--files', nargs='+', default=None,
                    help='Explicit list of .h5 files (overrides the positional folder).')
    ap.add_argument('--aor-overrides', default=None,
                    help='JSON file: {basename: aor_millimeter} used when a file records NaN aor.')
    ap.add_argument('--exclude', nargs='+', default=['_09_'],
                    help='Basename substrings to drop from the FIT (still reported). Default: _09_.')
    ap.add_argument('--r2-min', type=float, default=0.05,
                    help='Per-trial torque-vs-alpha R2 floor for bad-fit flagging (default 0.05).')
    ap.add_argument('--out', required=True,
                    help='Artifact JSON output path (REQUIRED). Name it deliberately per calibration, '
                         'e.g. 2026-07-06_apparatus_inertia_calibration.json -- there is no generic '
                         'default, and the name must match the config apparatus_inertia_calibration_file '
                         'for auto-load.')
    args = ap.parse_args(argv)

    paths = _resolve_files(args)

    overrides = {}
    if args.aor_overrides:
        if not os.path.isfile(args.aor_overrides):
            raise SystemExit('[error] --aor-overrides file not found: ' + str(args.aor_overrides))
        with open(args.aor_overrides, 'r') as fh:
            overrides = json.load(fh)
        if not isinstance(overrides, dict):
            raise SystemExit('[error] --aor-overrides must be a JSON object {basename: aor_mm}.')
        print('[info] loaded', len(overrides), 'aor override(s) from', args.aor_overrides)

    excluded = _resolve_excludes(paths, args.exclude)
    if excluded:
        print('[info] excluding from fit (still reported):', ', '.join(sorted(excluded)))

    out_path = args.out

    print('[info] fitting', len(paths), 'file(s); writing artifact to', out_path)
    art = fit_apparatus_inertia_calibration(
        paths,
        aor_overrides=overrides,
        exclude_files=excluded,
        r2_min=args.r2_min,
        output_json_path=out_path,
    )

    _print_report(art)
    print('\n[info] artifact written:', out_path)
    if art.get('blocked') or not art.get('fit_form_id'):
        print('[warn] no usable fit was selected -- see the report above.')
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
