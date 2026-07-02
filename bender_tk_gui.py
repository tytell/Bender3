"""
Tkinter operator console for CritterGripper / Bender (replacement for the Streamlit front end).

Design goals (functional, not pretty):
  (a) Hold session state (config + protocol) as plain in-memory objects: the ``Bender`` instance
      is the single source of truth (no st.session_state mirror).
  (b) Load config (.py module) and protocol (.json) templates.
  (c) Append notes to already-written HDF5 files  -- PARKED pending schema Decision A.
  (d) Generate preview plots as PNG files on disk (matplotlib savefig).

This module is GUI-glue only: every action calls an existing Streamlit-free backend callable
(``Bender.run_experiment``, ``bender_h5_export.export_primary_h5`` / ``save_universal_qc_figure``,
``bender_gui_preview.build_protocol_preview``, ``bender_config_builder`` / ``bender_protocol_templates``).
No NI-DAQ or motor logic lives here.

Run with:  python bender_tk_gui.py
"""
from __future__ import annotations

import os
import queue
import sys
import threading
import traceback
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

# Headless-safe matplotlib: we only savefig() preview PNGs, never embed an interactive canvas,
# so force the non-interactive Agg backend before pyplot is imported (avoids Tk/MPL backend clash).
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import numpy as np  # noqa: E402

from bender_functions import Bender  # noqa: E402
from bender_h5_export import export_primary_h5, save_universal_qc_figure  # noqa: E402
from bender_gui_preview import build_protocol_preview  # noqa: E402
from bender_config_builder import default_configs_dir, discover_config_modules  # noqa: E402
from bender_protocol_templates import (  # noqa: E402
    default_templates_dir,
    list_template_files,
    load_protocol_template,
    template_display_label,
)

_PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))

# Test types exposed in this first cut (frequency_sweep deferred per scope decision #8).
TEST_TYPES = ["dynamic", "isometric", "isovelocity"]

_STRAIN_MODES = ["strain", "strain_pct", "curvature", "angle"]
_VELOCITY_MODES = ["angle_vel", "strain_rate", "strain_pct_rate", "curvature_rate"]

# Per-test_type form specs: (label, bender_attr, kind, default, choices?)
#   kind: 'float' | 'int' | 'bool' | 'floatlist' | 'choice'
# Stim voltages use the CANONICAL left/right names; the backend resolves which physical channel
# (S1/S2) is left vs right from the config (S1side/S2side) -- see Bender._deposit_stim_on_side.
# Timing fields tagged TIMING are flagged for the later intuitive-timing cleanup (scope item #8).
PROTOCOL_FIELDS: Dict[str, List[tuple]] = {
    "dynamic": [
        ("Frequencies (Hz, comma list)", "all_freqs", "floatlist", "1.0"),
        ("Amplitudes (comma list)", "all_amps", "floatlist", "0.05"),
        ("Amplitude mode", "all_amps_mode", "choice", "strain", _STRAIN_MODES),
        ("Cycles per step", "cycles_per_step", "int", "5"),
        ("End cycles", "n_end_cycles", "int", "2"),
        ("Stim enabled", "is_stim", "bool", False),
        ("Stim pulse rate (Hz)", "stim_pulse_rate", "float", "75"),
        ("Left stim voltage (V)", "left_stim_voltage", "float", "5"),
        ("Right stim voltage (V)", "right_stim_voltage", "float", "5"),
    ],
    "isometric": [
        ("Initial target", "isometric_initial", "float", "0.0"),
        ("Final target", "isometric_final", "float", "0.1"),
        ("Number of steps", "isometric_num_steps", "int", "5"),
        ("Target mode", "isometric_mode", "choice", "strain", _STRAIN_MODES),
        ("Randomize step order", "randomize_step_order", "bool", False),
        ("Rest between steps (s)  [TIMING]", "rest_between_steps_s", "float", "0.0"),
        ("Left stim voltage (V)", "left_stim_voltage", "float", "5"),
        ("Right stim voltage (V)", "right_stim_voltage", "float", "5"),
        ("Reset max speed (deg/s)", "reset_max_speed_deg_per_s", "float", "15"),
    ],
    "isovelocity": [
        ("Min velocity", "isovelocity_min_vel", "float", "0.5"),
        ("Max velocity", "isovelocity_max_vel", "float", "2.0"),
        ("Starting strain", "isovelocity_starting_strain", "float", "0.0"),
        ("Starting strain mode", "isovelocity_starting_strain_mode", "choice", "strain", _STRAIN_MODES),
        ("Velocity mode", "isovelocity_velocity_mode", "choice", "angle_vel", _VELOCITY_MODES),
        ("Number of steps", "isovelocity_num_steps", "int", "5"),
        ("Iso duration (s)  [TIMING]", "isovelocity_iso_duration_s", "float", "0.2"),
        ("Pre-hold (s)  [TIMING]", "isovelocity_pre_hold_s", "float", "0.3"),
        ("Left stim voltage (V)", "left_stim_voltage", "float", "5"),
        ("Right stim voltage (V)", "right_stim_voltage", "float", "5"),
        ("Reset max speed (deg/s)", "reset_max_speed_deg_per_s", "float", "15"),
    ],
}


def _parse_float_list(text: str) -> List[float]:
    out: List[float] = []
    for tok in str(text or "").replace(";", ",").split(","):
        tok = tok.strip()
        if tok:
            out.append(float(tok))
    if not out:
        raise ValueError("expected at least one number")
    return out


def _coerce(kind: str, raw: Any) -> Any:
    if kind == "float":
        return float(raw)
    if kind == "int":
        return int(round(float(raw)))
    if kind == "bool":
        return bool(raw)
    if kind == "floatlist":
        return _parse_float_list(raw)
    return str(raw)


class BenderTkApp:
    """Tkinter controller. Holds one ``Bender`` instance as the in-memory session state."""

    def __init__(self, root: tk.Tk):
        self.root = root
        self.bender: Optional[Bender] = None
        self.last_h5_path: Optional[str] = None
        self.last_preview: Optional[dict] = None
        # Worker-thread -> UI-thread channel for the (blocking) run.
        self._run_queue: "queue.Queue[tuple]" = queue.Queue()
        self._run_thread: Optional[threading.Thread] = None

        self._field_widgets: Dict[str, Tuple[str, tk.Variable]] = {}

        root.title("CritterGripper / Bender - Tk console")
        root.geometry("760x900")
        self._build_ui()
        self._refresh_protocol_fields()

    # ----------------------------------------------------------------- UI build
    def _build_ui(self) -> None:
        outer = ttk.Frame(self.root, padding=8)
        outer.pack(fill="both", expand=True)

        canvas = tk.Canvas(outer, highlightthickness=0)
        scroll = ttk.Scrollbar(outer, orient="vertical", command=canvas.yview)
        self.body = ttk.Frame(canvas)
        self.body.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=self.body, anchor="nw")
        canvas.configure(yscrollcommand=scroll.set)
        canvas.pack(side="left", fill="both", expand=True)
        scroll.pack(side="right", fill="y")

        self._build_config_section(self.body)
        self._build_protocol_section(self.body)
        self._build_specimen_output_section(self.body)
        self._build_run_section(self.body)
        self._build_preview_section(self.body)
        self._build_parked_sections(self.body)

        self.status = tk.StringVar(value="No config loaded.")
        ttk.Label(self.body, textvariable=self.status, relief="sunken", anchor="w").pack(
            fill="x", pady=(8, 0)
        )

    def _section(self, parent: tk.Widget, title: str) -> ttk.LabelFrame:
        lf = ttk.LabelFrame(parent, text=title, padding=8)
        lf.pack(fill="x", pady=4)
        return lf

    def _build_config_section(self, parent: tk.Widget) -> None:
        lf = self._section(parent, "1 - Hardware configuration")
        try:
            modules = discover_config_modules(_PROJECT_ROOT)
        except Exception:
            modules = []
        self.cfg_var = tk.StringVar(value=modules[0] if modules else "")
        row = ttk.Frame(lf)
        row.pack(fill="x")
        ttk.Label(row, text="Config module:").pack(side="left")
        self.cfg_combo = ttk.Combobox(row, textvariable=self.cfg_var, values=modules, width=42)
        self.cfg_combo.pack(side="left", padx=4)
        ttk.Button(row, text="Load config", command=self._on_load_config).pack(side="left", padx=4)
        ttk.Label(
            lf,
            text=f"Scanning: {default_configs_dir(_PROJECT_ROOT)}",
            foreground="gray",
        ).pack(anchor="w", pady=(4, 0))
        self.cfg_status = tk.StringVar(value="(not loaded)")
        ttk.Label(lf, textvariable=self.cfg_status).pack(anchor="w")

    def _build_protocol_section(self, parent: tk.Widget) -> None:
        lf = self._section(parent, "2 - Protocol")
        row = ttk.Frame(lf)
        row.pack(fill="x")
        ttk.Label(row, text="Test type:").pack(side="left")
        self.tt_var = tk.StringVar(value=TEST_TYPES[0])
        tt_combo = ttk.Combobox(
            row, textvariable=self.tt_var, values=TEST_TYPES, state="readonly", width=18
        )
        tt_combo.pack(side="left", padx=4)
        tt_combo.bind("<<ComboboxSelected>>", lambda e: self._refresh_protocol_fields())
        ttk.Button(row, text="Load template (.json)", command=self._on_load_template).pack(
            side="left", padx=4
        )
        ttk.Button(row, text="Apply protocol", command=self._on_apply_protocol).pack(
            side="left", padx=4
        )

        self.fields_frame = ttk.Frame(lf)
        self.fields_frame.pack(fill="x", pady=(6, 0))

    def _build_specimen_output_section(self, parent: tk.Widget) -> None:
        lf = self._section(parent, "3 - Specimen & output")
        self.specimen_vars: Dict[str, tk.StringVar] = {}
        specs = [
            ("Specimen ID", "specimen_id", ""),
            ("Genus species", "specimen_genusspecies", ""),
            ("Clamp separation dclamp (mm)", "dclamp", "10.0"),
            ("Cross-section width (mm)", "xsec_width", "8.0"),
            ("Body mass (g)", "fishmass", "0.0"),
        ]
        grid = ttk.Frame(lf)
        grid.pack(fill="x")
        for i, (label, attr, default) in enumerate(specs):
            ttk.Label(grid, text=label, width=28).grid(row=i, column=0, sticky="w", pady=1)
            var = tk.StringVar(value=default)
            ttk.Entry(grid, textvariable=var, width=30).grid(row=i, column=1, sticky="w")
            self.specimen_vars[attr] = var
        ttk.Button(lf, text="Apply specimen", command=self._on_apply_specimen).pack(
            anchor="w", pady=(4, 0)
        )

        ttk.Separator(lf, orient="horizontal").pack(fill="x", pady=6)
        orow = ttk.Frame(lf)
        orow.pack(fill="x")
        ttk.Label(orow, text="Data folder:").grid(row=0, column=0, sticky="w")
        self.out_folder_var = tk.StringVar(value=os.path.join(_PROJECT_ROOT, "data"))
        ttk.Entry(orow, textvariable=self.out_folder_var, width=46).grid(row=0, column=1, sticky="w")
        ttk.Button(orow, text="Browse", command=self._on_browse_folder).grid(row=0, column=2, padx=4)
        ttk.Label(orow, text="File name:").grid(row=1, column=0, sticky="w")
        self.out_name_var = tk.StringVar(value="bender_trial_01.h5")
        ttk.Entry(orow, textvariable=self.out_name_var, width=46).grid(row=1, column=1, sticky="w")

        self.sim_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            lf,
            text="Simulation mode (no NI-DAQ hardware)",
            variable=self.sim_var,
        ).pack(anchor="w", pady=(6, 0))

    def _build_run_section(self, parent: tk.Widget) -> None:
        lf = self._section(parent, "4 - Run")
        row = ttk.Frame(lf)
        row.pack(fill="x")
        self.run_btn = ttk.Button(row, text="Run experiment", command=self._on_run)
        self.run_btn.pack(side="left")
        ttk.Button(row, text="KILL DAQ (reset NI device)", command=self._on_kill_daq).pack(
            side="left", padx=8
        )
        self.run_status = tk.StringVar(value="Idle.")
        ttk.Label(lf, textvariable=self.run_status, foreground="navy").pack(anchor="w", pady=(4, 0))
        ttk.Label(
            lf,
            text="Note: acquisition blocks for the full protocol; it runs on a worker thread so "
            "the window stays responsive. KILL DAQ cannot interrupt an in-progress synchronous run.",
            foreground="gray",
            wraplength=700,
            justify="left",
        ).pack(anchor="w")

    def _build_preview_section(self, parent: tk.Widget) -> None:
        lf = self._section(parent, "5 - Preview & QC plots (PNG on disk)")
        row = ttk.Frame(lf)
        row.pack(fill="x")
        ttk.Button(
            row, text="Preview motion -> PNG (no DAQ)", command=self._on_preview_png
        ).pack(side="left")
        ttk.Button(
            row, text="QC PNG from last .h5 (matplotlib)", command=self._on_qc_png_from_h5
        ).pack(side="left", padx=8)
        self.preview_status = tk.StringVar(value="")
        ttk.Label(lf, textvariable=self.preview_status, foreground="gray", wraplength=700).pack(
            anchor="w", pady=(4, 0)
        )

    def _build_parked_sections(self, parent: tk.Widget) -> None:
        lf = self._section(parent, "6 - Parked (pending PI decisions)")
        ttk.Label(
            lf,
            text="Append note to existing .h5  -  PARKED (Decision A resolved: metadata/note_bench).",
            foreground="gray",
            wraplength=700,
            justify="left",
        ).pack(anchor="w")
        ttk.Label(
            lf,
            text="Inertial section (specimen MOI + apparatus MOI / empirical override)  -  PARKED "
            "(Decision B).",
            foreground="gray",
            wraplength=700,
            justify="left",
        ).pack(anchor="w")

    # ----------------------------------------------------------- protocol fields
    def _refresh_protocol_fields(self) -> None:
        for child in self.fields_frame.winfo_children():
            child.destroy()
        self._field_widgets.clear()
        tt = self.tt_var.get()
        for i, spec in enumerate(PROTOCOL_FIELDS.get(tt, [])):
            label, attr, kind, default = spec[0], spec[1], spec[2], spec[3]
            choices = spec[4] if len(spec) > 4 else None
            ttk.Label(self.fields_frame, text=label, width=30).grid(
                row=i, column=0, sticky="w", pady=1
            )
            if kind == "bool":
                var: tk.Variable = tk.BooleanVar(value=bool(default))
                ttk.Checkbutton(self.fields_frame, variable=var).grid(row=i, column=1, sticky="w")
            elif kind == "choice":
                var = tk.StringVar(value=str(default))
                ttk.Combobox(
                    self.fields_frame,
                    textvariable=var,
                    values=choices or [],
                    state="readonly",
                    width=22,
                ).grid(row=i, column=1, sticky="w")
            else:
                var = tk.StringVar(value=str(default))
                ttk.Entry(self.fields_frame, textvariable=var, width=24).grid(
                    row=i, column=1, sticky="w"
                )
            self._field_widgets[attr] = (kind, var)

    # --------------------------------------------------------------- callbacks
    def _require_bender(self) -> bool:
        if self.bender is None:
            messagebox.showerror("No config", "Load a hardware configuration first (section 1).")
            return False
        return True

    def _on_load_config(self) -> None:
        mod = str(self.cfg_var.get() or "").strip()
        if not mod:
            messagebox.showerror("Config", "Select or type a config module name.")
            return
        if mod.endswith(".py"):
            mod = mod[:-3]
        # Ensure the configs dir is importable (discover_config_modules also does this).
        cfg_dir = default_configs_dir(_PROJECT_ROOT)
        for p in (cfg_dir, _PROJECT_ROOT):
            if p and p not in sys.path:
                sys.path.insert(0, p)
        try:
            self.bender = Bender(mod)
        except Exception as e:
            self.bender = None
            self.cfg_status.set(f"FAILED: {type(e).__name__}: {e}")
            messagebox.showerror("Config load failed", f"{type(e).__name__}: {e}")
            return
        self.cfg_status.set(f"Loaded: {mod}")
        self.status.set(f"Config '{mod}' loaded. Bender ready.")
        # Seed output filename/folder from the config's outputfile if present.
        outp = str(getattr(self.bender, "outputfile", "") or "").strip()
        if outp:
            self.out_folder_var.set(os.path.dirname(os.path.normpath(outp)) or self.out_folder_var.get())
            self.out_name_var.set(os.path.basename(os.path.normpath(outp)))

    def _on_load_template(self) -> None:
        if not self._require_bender():
            return
        initial = default_templates_dir(_PROJECT_ROOT)
        if not os.path.isdir(initial):
            initial = _PROJECT_ROOT
        path = filedialog.askopenfilename(
            title="Load protocol template",
            initialdir=initial,
            filetypes=[("Protocol template", "*.json"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            tpl = load_protocol_template(path)
        except Exception as e:
            messagebox.showerror("Template", f"Could not read template: {e}")
            return
        tt = str(tpl.get("test_type") or "").strip()
        if tt not in TEST_TYPES:
            messagebox.showwarning(
                "Template",
                f"Template test_type '{tt}' is not supported in this console "
                f"(supported: {', '.join(TEST_TYPES)}).",
            )
            return
        self.tt_var.set(tt)
        self._refresh_protocol_fields()
        proc = tpl.get("procedure") or {}
        if not isinstance(proc, dict):
            messagebox.showerror("Template", "Template 'procedure' must be a JSON object.")
            return
        # Fill any matching form fields; values still need an explicit Apply.
        filled = 0
        for attr, (kind, var) in self._field_widgets.items():
            if attr in proc and proc[attr] is not None:
                val = proc[attr]
                if kind == "bool":
                    var.set(bool(val))
                elif kind == "floatlist":
                    seq = val if isinstance(val, (list, tuple)) else [val]
                    var.set(", ".join(str(float(x)) for x in seq))
                else:
                    var.set(str(val))
                filled += 1
        self.status.set(
            f"Template '{template_display_label(path)}' loaded into form "
            f"({filled} fields). Click Apply protocol."
        )

    def _collect_protocol_kwargs(self) -> Dict[str, Any]:
        kwargs: Dict[str, Any] = {}
        for attr, (kind, var) in self._field_widgets.items():
            raw = var.get()
            try:
                kwargs[attr] = _coerce(kind, raw)
            except (TypeError, ValueError) as e:
                raise ValueError(f"Field '{attr}': {e}") from e
        return kwargs

    def _on_apply_protocol(self) -> None:
        if not self._require_bender():
            return
        tt = self.tt_var.get()
        try:
            kwargs = self._collect_protocol_kwargs()
        except ValueError as e:
            messagebox.showerror("Protocol", str(e))
            return
        try:
            # update_metadata is the Streamlit-free canonical setter (alias-aware). Stim voltages
            # use canonical left/right names; the backend maps to S1/S2 from config side labels.
            self.bender.update_metadata(test_type=tt, **kwargs)
            self.bender.test_type = tt
        except Exception as e:
            messagebox.showerror("Protocol", f"Apply failed: {type(e).__name__}: {e}")
            return
        ok, missing = self._validate(tt)
        if ok:
            self.status.set(f"Protocol '{tt}' applied and valid.")
        else:
            self.status.set(f"Protocol '{tt}' applied. Missing/invalid: {', '.join(missing)}")

    def _validate(self, tt: str) -> Tuple[bool, List[str]]:
        try:
            rep = self.bender.validate_dispatch_setup(test_type=tt)
            return bool(rep.get("ok", False)), list(rep.get("missing", []) or [])
        except Exception as e:
            return False, [f"{type(e).__name__}: {e}"]

    def _on_apply_specimen(self) -> None:
        if not self._require_bender():
            return
        kwargs: Dict[str, Any] = {}
        for attr, var in self.specimen_vars.items():
            raw = str(var.get()).strip()
            if raw == "":
                continue
            if attr in ("dclamp", "xsec_width", "fishmass"):
                try:
                    kwargs[attr] = float(raw)
                except ValueError:
                    messagebox.showerror("Specimen", f"'{attr}' must be a number.")
                    return
            else:
                kwargs[attr] = raw
        try:
            self.bender.update_metadata(**kwargs)
            # dclamp has a legacy alias used by the strain/curvature math; keep them in sync.
            if "dclamp" in kwargs:
                self.bender.test_segment_length_mm = kwargs["dclamp"]
        except Exception as e:
            messagebox.showerror("Specimen", f"Apply failed: {type(e).__name__}: {e}")
            return
        self.status.set("Specimen/morphometrics applied.")

    def _on_browse_folder(self) -> None:
        d = filedialog.askdirectory(title="Choose data folder", initialdir=self.out_folder_var.get())
        if d:
            self.out_folder_var.set(d)

    def _compose_output_path(self) -> str:
        folder = str(self.out_folder_var.get() or "").strip()
        name = str(self.out_name_var.get() or "").strip()
        if not name:
            raise ValueError("Set a data file name (section 3).")
        if not name.lower().endswith(".h5"):
            name += ".h5"
        return os.path.join(folder, name) if folder else name

    # -------------------------------------------------------------- run (thread)
    def _on_run(self) -> None:
        if not self._require_bender():
            return
        if self._run_thread is not None and self._run_thread.is_alive():
            messagebox.showwarning("Run", "A run is already in progress.")
            return
        tt = self.tt_var.get()
        ok, missing = self._validate(tt)
        if not ok:
            messagebox.showerror("Run", f"Protocol not ready: {', '.join(missing)}")
            return
        try:
            out_path = self._compose_output_path()
        except ValueError as e:
            messagebox.showerror("Run", str(e))
            return
        self.bender.outputfile = out_path
        self.bender.session_simulated = bool(self.sim_var.get())
        self.run_btn.config(state="disabled")
        self.run_status.set("Acquiring... (window stays responsive)")
        self._run_thread = threading.Thread(
            target=self._run_worker, args=(tt, out_path), daemon=True
        )
        self._run_thread.start()
        self.root.after(150, self._poll_run_queue)

    def _run_worker(self, tt: str, out_path: str) -> None:
        try:
            self.bender.run_experiment(test_type=tt)
            rep = export_primary_h5(self.bender, outputfile=out_path, append_post_trial_notes=True)
            qc_path = None
            try:
                base = os.path.splitext(rep.get("outputfile", out_path))[0] + "_qc"
                qc_path, _ = save_universal_qc_figure(self.bender, base_path=base)
            except Exception as qc_err:  # QC is non-fatal; data already saved.
                qc_path = f"(QC plot failed: {qc_err})"
            self._run_queue.put(("done", rep, qc_path))
        except Exception:
            self._run_queue.put(("error", traceback.format_exc()))

    def _poll_run_queue(self) -> None:
        try:
            item = self._run_queue.get_nowait()
        except queue.Empty:
            if self._run_thread is not None and self._run_thread.is_alive():
                self.root.after(150, self._poll_run_queue)
            return
        self.run_btn.config(state="normal")
        if item[0] == "done":
            _, rep, qc_path = item
            self.last_h5_path = rep.get("outputfile")
            self.run_status.set(f"Saved: {self.last_h5_path}")
            self.preview_status.set(f"QC plot: {qc_path}")
            messagebox.showinfo(
                "Run complete",
                f"Data file:\n{self.last_h5_path}\n\nQC plot:\n{qc_path}",
            )
        else:
            self.run_status.set("Run FAILED (see dialog).")
            messagebox.showerror("Run failed", item[1])

    def _on_kill_daq(self) -> None:
        try:
            from bender_daq_kill import daq_emergency_stop

            dev = getattr(self.bender, "device_name", None) if self.bender else None
            if not dev:
                messagebox.showwarning("KILL DAQ", "No device_name (load a config first).")
                return
            motor_port = getattr(self.bender, "motor_port", None) if self.bender else None
            enable_line = f"{motor_port}/line2" if motor_port else None
            daq_emergency_stop(dev, release_motor_enable_line=enable_line)
            self.run_status.set(f"KILL DAQ: reset device '{dev}'.")
        except Exception as e:
            messagebox.showerror("KILL DAQ", f"{type(e).__name__}: {e}")

    # ------------------------------------------------------------------ preview
    def _on_preview_png(self) -> None:
        if not self._require_bender():
            return
        tt = self.tt_var.get()
        try:
            prev = build_protocol_preview(self.bender, requested_test_type=tt)
        except Exception as e:
            messagebox.showerror("Preview", f"{type(e).__name__}: {e}")
            return
        if not prev.get("ok"):
            messagebox.showerror("Preview", prev.get("error") or "Preview failed.")
            return
        self.last_preview = prev
        try:
            folder = str(self.out_folder_var.get() or _PROJECT_ROOT).strip() or _PROJECT_ROOT
            os.makedirs(folder, exist_ok=True)
            stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
            out_png = os.path.join(folder, f"preview_{tt}_{stamp}.png")
            self._render_preview_png(prev, out_png, tt)
        except Exception as e:
            messagebox.showerror("Preview", f"Could not write PNG: {type(e).__name__}: {e}")
            return
        self.preview_status.set(f"Preview PNG: {out_png}")
        messagebox.showinfo("Preview", f"Preview written:\n{out_png}")

    @staticmethod
    def _render_preview_png(prev: dict, out_png: str, tt: str) -> None:
        t = np.asarray(prev.get("t"), dtype=float).reshape(-1)
        angle = np.asarray(prev.get("angle"), dtype=float).reshape(-1)
        anglevel = prev.get("anglevel")
        stim_total = prev.get("stim_total")
        n_rows = 2 if stim_total is not None else 1
        fig, axes = plt.subplots(n_rows, 1, figsize=(11, 3.2 * n_rows), sharex=True, squeeze=False)
        ax0 = axes[0][0]
        ax0.plot(t, angle, color="royalblue", lw=1.0, label="angle_commanded (deg)")
        ax0.set_ylabel("Angle (deg)")
        ax0.grid(True, alpha=0.3)
        if anglevel is not None:
            axv = ax0.twinx()
            axv.plot(t, np.asarray(anglevel, dtype=float).reshape(-1), color="orange", lw=0.7)
            axv.set_ylabel("Velocity (deg/s)", color="orange")
        ax0.legend(loc="upper right", fontsize=8)
        ax0.set_title(f"Motion preview - {tt}")
        if stim_total is not None:
            ax1 = axes[1][0]
            ax1.plot(t, np.asarray(stim_total, dtype=float).reshape(-1), color="firebrick", lw=0.7)
            ax1.set_ylabel("Stim (V)")
            ax1.set_xlabel("Time (s)")
            ax1.grid(True, alpha=0.3)
        else:
            ax0.set_xlabel("Time (s)")
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close(fig)

    def _on_qc_png_from_h5(self) -> None:
        path = self.last_h5_path
        if not path or not os.path.isfile(path):
            path = filedialog.askopenfilename(
                title="Choose .h5 file for QC plot",
                initialdir=self.out_folder_var.get(),
                filetypes=[("HDF5", "*.h5"), ("All files", "*.*")],
            )
        if not path:
            return
        try:
            from validate_plot_h5_batch import plot_qc_from_h5

            out = plot_qc_from_h5(path)
        except Exception as e:
            messagebox.showerror("QC plot", f"{type(e).__name__}: {e}")
            return
        self.preview_status.set(f"QC PNG: {out}")
        messagebox.showinfo("QC plot", f"QC plot written:\n{out}")


def main() -> int:
    root = tk.Tk()
    BenderTkApp(root)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
