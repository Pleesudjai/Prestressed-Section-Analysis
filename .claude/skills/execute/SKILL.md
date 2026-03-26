# /execute — Implementation

Clean implementation mode. Read spec, write code, verify.

## Arguments
- `$ARGUMENTS` — spec name (e.g., "ultimate-capacity") or direct task description

## Steps

1. **Load context** — Read only:
   - `CLAUDE.md` (if not already loaded)
   - `specs/$ARGUMENTS.md` (if a spec exists)
   - Relevant `.claude/rules/` file for the domain you're touching

2. **Implement** — Follow the spec steps in order:
   - Check off each step as completed
   - New MATLAB scripts → `08_matlab_program_analysis/`
   - Figures → appropriate `Project*/` subfolder
   - Reports → `06_draft_text/`

3. **Verify** — After implementation:
   - Read back the code and check sign convention (compression = +)
   - Verify coordinate origin (y = 0 at bottom of stems)
   - Check that η is only applied at service, not transfer
   - Confirm all functions have correct signatures matching the call graph in CLAUDE.md

4. **Report blockers** — If stuck for more than 2 attempts on the same issue, stop and report

5. **Run /commit** when done

## Rules
- Do NOT load planning context or research into this session — spec only
- Do NOT refactor unrelated code
- Do NOT change sign conventions in existing files
- Working directory for all MATLAB: `08_matlab_program_analysis/`
