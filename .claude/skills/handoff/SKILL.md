# /handoff — Session Summary

Write a handoff note so the next session (or next Claude instance) can pick up exactly where you left off.

## Steps

1. Write `docs/handoff.md` with this format:

```markdown
# Session Handoff — [date]

## Completed This Session
- [bullet list of what was done]

## Current State
- Working / In-progress / Broken: [which]
- Key files modified: [list]
- Tests/verification status: [pass/fail/not run]

## Next Steps (priority order)
1. [most important next action]
2. [second priority]
3. [third priority]

## Open Questions / Blockers
- [anything unresolved]

## Context for Next Session
- [any non-obvious decisions made, or gotchas discovered]
```

2. Update the `## ▶ RESUME HERE` section in `CLAUDE.md`:
   - Update "Current focus", "Status", "Next action", "Status updated" fields

## Rules
- Keep it concise — next Claude needs facts, not narrative
- Include file paths for everything modified
- If context is >50 messages, strongly recommend starting a fresh session after this
