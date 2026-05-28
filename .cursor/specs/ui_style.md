# CritterGripper / Bender — UI styling specification

## Authority metadata

- Purpose: Streamlit visual style contract for CritterGripper.
- Authority level: Tier 2 (binding for Streamlit UI tasks).
- Scope: Visual style, color, spacing, and button patterns in Streamlit UI.
- Owner: Project maintainer.
- Last reviewed: 2026-05-28.

This file does not override product architecture in `.cursorrules` (Tier 1),
and is the visual companion to `ux_spec.md` (which defines behavior). When the
two specs overlap, `ux_spec.md` governs behavior and this file governs
appearance.

Reference for layout, styling, and UX tone when changing
`bender_streamlit_gui.py`. Implementation lives mainly in
`_inject_accessibility_theme`, `_inject_load_save_button_theme`, and
`_landing_style`.

---

## Vibe and audience

- **Industrial / scientific dashboard:** calm, legible, low cognitive load;
  suitable for lab operators and researchers.
- **Safety and clarity:** obvious primary actions, consistent naming, visible
  state (disabled actions, validation icons).
- **Restrained chrome:** white workflow shell, slate typography, bordered
  cards instead of heavy decoration.
- **Accessibility:** skip link to main content; visible `:focus-visible`
  rings. Text scaling is handled by browser zoom (Cmd/Ctrl +), not an in-app
  toggle.

---

## Information architecture

Two top-level entry points from Home (per `ux_spec.md` Section 2):

| Entry point | Role |
|-------------|------|
| **Run Experiment** | Single scrolling page: the full experiment workflow |
| **Review Data** | H5 explorer: open `.h5`, plot series, append notes |

An offline **Simulate Bending Mechanics** tool may remain accessible from Home.

Within Run Experiment, use this triad consistently in labels and section
titles:

| Pillar | Role |
|--------|------|
| **Hardware configuration** | Instrument / DAQ / channels |
| **Specimen** | Specimen and subject metadata (formerly "Biometrics"; rename to "Morphometrics" pending — see ux_spec Section 10) |
| **Protocol** | Experiment procedure (use "protocol" consistently, never "procedure") |

---

## Modes and CSS markers

| Marker | Meaning |
|--------|---------|
| `.bnd-landing-page` | Home/landing only: hide sidebar, custom shell |
| `.bnd-workflow-active` | Inside Run Experiment: white app chrome; primary red / secondary white buttons |
| `.bnd-sim-compare-active` | Simulate route: numpy-only teaching view; no NI-DAQ; workbench ~36/64 columns (controls left, plots right) |
| `.bnd-sim-osc-banner` | Oscillatory viscoelastic mode callout on the simulate route (navy left accent, slate wash) |
| `.bnd-ls-action` | Wrapper marker for Load/Save full-width rows (styling hook) |

Note: theme-tint markers (`.bnd-ui-theme`, `bnd-theme-*`), the large-text
marker (`.bnd-a11y-large-text`), and the stepwise marker
(`.bnd-stepwise-active`) are **removed**. The app ships with a single fixed
appearance.

---

## Color palette

### Workflow shell (`body:has(.bnd-workflow-active)`)

- **App / header / sidebar background:** `#ffffff`
- **Sidebar border:** `#e2e8f0`
- **Main headings / labels:** `#334155`

### Buttons (workflow + landing main)

- **Primary (`type='primary'`):** fill `#dc2626`, border `#b91c1c`, text
  `#ffffff`, font-weight 600, transition ~`0.15s` on background/border/shadow.
- **Primary hover ("lights up"):** fill `#f87171`, border `#fca5a5`, soft
  outer glow `box-shadow: 0 0 0 3px rgba(248, 113, 113, 0.45)`.
- **Secondary (`type='secondary'`):** fill `#ffffff`, border `#cbd5e1`, label
  `#334155`.
- **Secondary hover:** fill `#fff1f2`, border `#fb7185`, label `#991b1b`,
  light drop shadow `rgba(220, 38, 38, 0.18)`.
- **Disabled primary:** muted red fill/border, no glow, opacity ~`0.9`.
- **Disabled secondary:** `#f1f5f9` fill, slate border/text, no shadow.
- **Focus-visible:** white inner ring + outer ring (`#dc2626` primary,
  `#94a3b8` secondary); load/save rows use the same pattern in
  `:has(.bnd-ls-action)` rules.

**Load/Save** (`_load_save_button`) is **primary** by default (red).

### Landing / Home accents

- **Simulate CTA:** deep slate / navy gradient (`#1e3a5f` → `#0f172a`), light
  text `#f8fafc`, hover cyan accent — scoped via `.bnd-land-sim-btn-marker` so
  live workflow primaries stay red.
- **Headings:** `#334155`
- **Tagline / secondary text:** `#64748b` (tagline ~`1.32rem`, weight ~450)
- **Hero rule (below hero):** solid `4px`, `#c2410c`
- **Bordered cards:** white `#ffffff`, border `#cbd5e1`, radius **14px**, light
  shadow `rgba(15, 23, 42, 0.06)`; inner body text `#475569`
- **Equation / "learn" callout (EOM box):** teal gradient panel, left accent
  `#0d9488`

### Bordered panels (`st.container(border=True)`)

- **Default (main):** `2px` border `#94a3b8`, radius **10px**, padding
  ~`0.65–1.2rem`, white background
- **Workflow override:** always **white** fill `#ffffff`, border `#cbd5e1`
- **Alerts:** framed with `#cbd5e1`

---

## Spacing and layout

- **Landing main column:** `max-width` ~**52rem**, centered; moderate vertical
  padding on `.block-container`
- **Between major landing blocks:** spacer class ~**1.65rem** height
- **Hero with logo:** reserve right space (`padding-right` ~`min(38%, 15rem)`);
  logo column ~`min(34%, 13.5rem)`; hero `min-height` ~**11rem** so logo and
  hero rule do not overlap; logo image `max-height: 100%` inside the hero box.
- **Responsive landing cards/buttons:** below ~`980px`, description cards and
  CTA button rows stack to full width (one per row).

Prefer **consistent vertical rhythm** — smaller gaps in dense areas, slightly
roomier on Home/marketing blocks.

---

## Button behaviors and patterns

1. **Load / Save and other file commits:** `use_container_width=True`, via
   `_load_save_button` + invisible `.bnd-ls-action` marker for styling.
2. **Mode / branch choice (e.g. Load existing vs Build new):** two equal
   columns, wide buttons; selected = `primary`, other = `secondary` — not the
   load/save marker pattern; avoid horizontal `st.radio` for this binary choice.
3. **Apply buttons:** each section's Apply (per `ux_spec.md` Section 3) is a
   full-width primary action at the foot of that section.
4. **Landing primary CTAs:** same red / white / min-height **3rem** treatment
   as workflow (`_landing_style`).
5. **Danger actions:** safety-critical actions (e.g. emergency stop) live in
   their own bordered block for visual separation.

---

## Typography (landing / main)

- **Hero H1:** large clamp (~`2.45rem`–`3.45rem`), tight line height, slate
  color, negative letter-spacing
- **Section titles:** centered, ~`1.45rem`, semibold, `#475569`
- **Captions:** `#64748b`

---

## Implementation notes for agents

- Global selectors should target **`[data-testid="stMain"]`** (not only
  `section`) for Streamlit DOM variance.
- Selector policy is **guarded**: prefer existing hooks/selectors first;
  resilient fallbacks are allowed when Streamlit DOM variance requires them —
  briefly explain the fallback in the change summary.
- **Landing CSS** prefixes with **`:is(#root, body:has(.bnd-landing-page))`** —
  some Streamlit builds have no `#root`, so `body:has(.bnd-landing-page)` keeps
  rules matching while the home marker is present.
- Inject **hero/logo CSS before** the logo `<img>` to avoid layout flash.
- When adding new full-width primary actions, reuse `_load_save_button` or
  document why a separate pattern is needed.
- After visual changes, check **Home** and **Run Experiment** together;
  workflow rules often override panel styling.

---

*bender_ui_style_v2_0 | Jimenez Lab | 2026-05-28 | Aligned with single-page
restructure; theme/stepwise/large-text removed.*
