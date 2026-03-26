# /plan-feature — Feature / Task Planning

Plan before coding. Use this for any task that will take more than a few edits.

## Arguments
- `$ARGUMENTS` — brief description of the feature or task (e.g., "ultimate moment capacity check")

## Steps

1. **Clarify scope** — What are the inputs, outputs, and success criteria?
   - For MATLAB scripts: what function signature? what plots? what console output?
   - For analysis: what ACI code sections apply? what equations?
   - For reports: what sections? what figures needed?

2. **Research via sub-agents** — Spin up Explore agents for:
   - Existing code that does something similar (search `08_matlab_program_analysis/`)
   - Relevant theory from `01_sources/` PDFs or NotebookLM
   - ACI 318 provisions (query NotebookLM: `cee-530-prestressed-concrete`)
   - Return **3-5 bullet summaries only** — never dump raw content

3. **Write spec** — Create `specs/$ARGUMENTS.md` with:
   ```
   # [Feature Name]
   ## What & Why
   ## Inputs (data, files, parameters)
   ## Outputs (figures, console, files)
   ## Equations / Theory Reference
   ## Implementation Steps (numbered, checkable)
   ## Open Questions
   ## Out of Scope
   ```

4. **Review with user** — Present the spec, ask for approval before coding

5. **Handoff note** — Tell user: "Run `/execute` with this spec in a fresh session if context is long"

## Rules
- Sign convention: compression = positive, tension = negative
- All new .m files go in `08_matlab_program_analysis/`
- Reference Naaman textbook or ACI 318 for equations — cite section numbers
- Never skip planning for features involving new stress conditions or load cases
