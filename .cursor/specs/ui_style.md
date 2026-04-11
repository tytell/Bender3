# CritterGripper / Bender — UI specification

Reference for layout, styling, and UX tone when changing `bender_streamlit_gui.py` or related UI. Implementation lives mainly in `_inject_accessibility_theme`, `_inject_load_save_button_theme`, `_inject_stepwise_compact_layout_css`, and `_landing_style`.

---

## Vibe and audience

- **Industrial / scientific dashboard:** calm, legible, low cognitive load; suitable for lab operators and researchers.
- **Safety and clarity:** obvious primary actions, consistent naming, visible state (disabled steps, progress).
- **Restrained chrome:** white workflow shell, slate typography, bordered cards instead of heavy decoration.
- **Accessibility:** skip link to main content; visible `:focus-visible` rings; optional **large text** mode in Settings.

---

## Information architecture (user-facing language)

**App routes** (see `gui_app_route`): `landing`, `scratch`, `templates`, `stepwise`, `sim_compare` (Simulation & Comparison — offline numpy), `h5_explorer`.

Use this **single triad** consistently in labels and section titles:

| Pillar | Role |
|--------|------|
| **Hardware configuration** | Instrument / DAQ / channels |
| **Biometrics** | Specimen and subject metadata |
| **Protocol** | Experiment procedure (not “procedure” vs “protocol” split elsewhere) |

Stepwise flows use **numbered steps (1–5)** and short subheaders (e.g. `1 · …` / `2 · …`). Prefer **wide, full-width** commit actions over tiny controls for load/save and navigation.

---

## Modes and CSS markers

| Marker | Meaning |
|--------|---------|
| `.bnd-landing-page` | Landing only: hide sidebar, custom shell |
| `.bnd-workflow-active` | Inside a workflow: white app chrome; primary red / secondary white buttons |
| `.bnd-sim-compare-active` | **Simulation & Comparison** route (`sim_compare`): numpy-only teaching view; no NI-DAQ; **workbench** = ~36/64 columns (controls left, plots right) |
| `.bnd-sim-osc-banner` | Oscillatory viscoelastic mode callout on `sim_compare` (navy left accent, slate wash) |
| `.bnd-stepwise-active` | Stepwise layout: tighter vertical rhythm |
| `.bnd-ui-theme` + `bnd-theme-warm` / `bnd-theme-cool` / `bnd-theme-slateivory` | User-chosen background tint (Settings) |
| `.bnd-a11y-large-text` | Slightly larger base font |
| `.bnd-ls-action` | Wrapper marker for **Load/Save** full-width rows (styling hook) |

---

## Color palettes

### Workflow shell (`body:has(.bnd-workflow-active)`)

- **App / header / sidebar background:** `#ffffff`
- **Sidebar border:** `#e2e8f0`
- **Main headings / labels (cool default):** `#334155` (warm theme: `#44403c`)

### Buttons (workflow + landing main)

- **Primary (`type='primary'`):** fill `#dc2626`, border `#b91c1c`, text `#ffffff`, **font-weight 600**, `transition` ~`0.15s` on background/border/shadow.
- **Primary hover (“lights up”):** fill `#f87171`, border `#fca5a5`, soft outer glow `box-shadow: 0 0 0 3px rgba(248, 113, 113, 0.45)`.
- **Secondary (`type='secondary'`):** fill `#ffffff`, border `#cbd5e1`, label `#334155`.
- **Secondary hover:** fill `#fff1f2`, border `#fb7185`, label `#991b1b`, light drop shadow `rgba(220, 38, 38, 0.18)`.
- **Disabled primary:** muted red fill/border, no glow, opacity ~`0.9`.
- **Disabled secondary:** `#f1f5f9` fill, slate border/text, no shadow.
- **Focus-visible:** white inner ring + outer ring (`#dc2626` primary, `#94a3b8` secondary); load/save rows use the same pattern in `:has(.bnd-ls-action)` rules.

**Load/Save** (`_load_save_button`) is **primary** by default (red); secondary load/save stays on the white treatment if used.

### Optional app themes (non-workflow / under global theme)

| Theme | App / sidebar / main text (summary) |
|-------|--------------------------------------|
| **Warm paper** | `#faf8f5`, `#f3eee6`, `#fffefb`; text ~stone |
| **Cool gray** | `#f1f5f9`, `#e2e8f0`, `#f8fafc`; text slate |
| **Slate & ivory** | ivory main, **dark slate sidebar** `#1e293b`; light sidebar text |

### Landing page accents

- **Simulation & Comparison CTA:** deep slate / navy gradient (`#1e3a5f` → `#0f172a`), light text `#f8fafc`, hover cyan accent — scoped via `.bnd-land-sim-btn-marker` (sibling of the Streamlit button row) so live workflow primaries stay red.
- **Headings:** `#334155`
- **Tagline / secondary text:** `#64748b` (tagline ~`1.32rem`, weight ~450)
- **Hero rule (below hero):** solid `4px`, `#c2410c`
- **Bordered cards:** white `#ffffff`, border `#cbd5e1`, radius **14px**, light shadow `rgba(15, 23, 42, 0.06)`; inner body text `#475569`
- **Equation / “learn” callout (EOM box):** teal gradient panel, left accent `#0d9488`

### Bordered panels (`st.container(border=True)`)

- **Default (main):** `2px` border `#94a3b8`, radius **10px**, padding ~`0.65–1.2rem`, white background
- **Workflow override:** always **white** fill `#ffffff`, border `#cbd5e1` (wins over theme-tinted panels)
- **Alerts:** framed with `#cbd5e1` (theme variants adjust border color)

---

## Spacing and layout

- **Landing main column:** `max-width` ~**52rem**, centered; moderate vertical padding on `.block-container`
- **Between major landing blocks:** spacer class ~**1.65rem** height
- **Hero with logo:** reserve right space (`padding-right` ~`min(38%, 15rem)`); logo column ~`min(34%, 13.5rem)`; **hero `min-height`** ~**11rem** so logo and hero rule do not overlap; logo image **`max-height: 100%`** inside the hero box (avoid fixed `11rem` on the img if it exceeds the box)
- **Stepwise (`bnd-stepwise-active`):** reduced block padding, tighter bordered panels (`~0.4–0.55rem` padding, smaller `margin-bottom`), compact `h1–h3` margins

Prefer **consistent vertical rhythm** (smaller gaps in rails, slightly roomier on marketing/landing).

---

## Button behaviors and patterns

1. **Load / Save and other file commits:** `use_container_width=True`, implemented via `_load_save_button` + invisible `.bnd-ls-action` marker for styling.
2. **Mode / branch choice (e.g. Load existing vs Build new):** **two equal columns**, wide buttons; **selected = `primary`**, other = `secondary`** — **not** the load/save marker pattern; avoid horizontal `st.radio` for this binary choice.
3. **Step navigation:** bordered panel with progress and **tab-like step controls**; disabled styling for current step; **Previous / Next** at full width where appropriate.
4. **Landing primary CTAs:** same red / white / min-height **3rem** treatment as workflow (`_landing_style`).

---

## Typography (landing / main)

- **Hero H1:** large clamp (~`2.45rem`–`3.45rem`), tight line height, slate color, negative letter-spacing
- **Section titles:** centered, ~`1.45rem`, semibold, `#475569`
- **Captions:** `#64748b`

---

## Implementation notes for agents

- Global selectors should target **`[data-testid="stMain"]`** (not only `section`) for Streamlit DOM variance.
- **Landing CSS** prefixes with **`:is(#root, body:has(.bnd-landing-page))`** — some Streamlit builds have no `#root` element, so `body:has(.bnd-landing-page)` keeps rules matching while the home marker is present.
- Inject **hero/logo CSS before** the logo `<img>` to avoid layout flash.
- When adding new full-width primary actions, reuse `_load_save_button` **or** document why a separate pattern is needed.
- After visual changes, check **workflow**, **landing**, and **stepwise** together; workflow rules often override theme-tinted panels.

---

*Last aligned with `bender_streamlit_gui.py` patterns (CritterGripper / Bender PC GUI).*
