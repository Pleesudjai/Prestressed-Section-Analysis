# /commit — Document Work & Git Commit

Log what was done, why, and commit to git.

## Steps

1. **Append to `docs/decisions.md`**:
   ```markdown
   ## [date] — [short title]
   **What:** [what was built/changed]
   **Why:** [motivation — assignment requirement, bug fix, design decision]
   **Files:** [key files created or modified]
   **Next:** [what should happen next]
   ```

2. **Update `.claude/rules/`** if a new convention or pattern was discovered this session

3. **Git commit**:
   - Stage only relevant files (never `git add -A` blindly)
   - Commit message format:
     ```
     [type]: [description]

     - [bullet list of changes]
     ```
   - Commit types: `feat`, `fix`, `refactor`, `docs`, `data`, `hw`, `project`

4. **Update Task Checklist** in CLAUDE.md if a checklist item was completed

## Rules
- Never commit files from `01_sources/` (large PDFs)
- Never commit `.env`, credentials, or API keys
- Check `git status` before staging — don't include unintended files
- Ask user before pushing to remote
