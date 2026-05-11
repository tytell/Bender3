# CritterGripper New-Agent Bootstrap

Use this when starting a fresh agent session.

## Instruction hierarchy

1. Tier 1 (binding): `.cursorrules`
2. Tier 2 (binding for Streamlit UI tasks): `.cursor/specs/ui_style.md`
3. Tier 3 (hybrid roadmap): `EVOLUTION_PLAN.md` (advisory by default; binding only when explicitly referenced by task/user)

## Operating rules

- Product scope is Streamlit-first and strict unless the task explicitly authorizes an alternative.
- Use consistent user-facing terms: **Hardware configuration**, **Biometrics**, **Protocol**.
- If names are changed by request, apply consistently across affected UI/docs unless instructed otherwise.
- On any instruction conflict or ambiguity, stop and ask before coding.
- If the same conflict pattern repeats, propose a concise rule snippet for future inclusion.
- Prefer existing UI hooks/selectors; resilient fallback selectors are allowed when needed, with a brief rationale.
- Medium-sized refactors are allowed when they solve root causes.
- For DAQ/hardware-affecting changes, include a safety checklist in the summary (stop/reset path, timeout behavior, stale-export prevention, user-visible warnings).
- After substantive edits, run the full test suite by default; if blocked, report why and run the largest feasible subset.
