# /prime — Session Initialization

You are resuming work on a CEE 530 Prestressed Concrete project (MATLAB analysis codebase).

## Steps

1. Read `CLAUDE.md` — architecture, rules, current focus
2. Read `docs/decisions.md` — what has been decided and why
3. Read `docs/handoff.md` (if exists) — previous session state
4. Scan `08_matlab_program_analysis/` for .m files: `ls 08_matlab_program_analysis/*.m`
5. Check recent git history: `git log --oneline -10`
6. Check for uncommitted changes: `git status`

## Report Format

After reading, report concisely:

```
SESSION PRIMED
  Done:     [completed items from handoff/decisions]
  Current:  [in-progress work]
  Next:     [priority next action from CLAUDE.md ▶ RESUME HERE]
  Blockers: [any issues found]
  Modified: [uncommitted files, if any]
```

## Rules
- Do NOT start coding yet — just report state
- Do NOT load .claude/rules/ files — only load those when actively working in that domain
- If docs/handoff.md is missing, say so and ask user what they want to work on
