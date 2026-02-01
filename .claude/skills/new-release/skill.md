---
name: new-release
description: Create a new release of the piezod package. Bumps version, creates GitHub release, and publishes to PyPI.
allowed-tools: Read, Edit, Bash, AskUserQuestion
---

# New Release

Creates a new release of the piezod Python package.

## Usage

```
/new-release
/new-release 0.2.0
```

If no version is provided, suggests a patch increment.

## Workflow

### Phase 1: Gather Information

1. Read `python/pyproject.toml` to get current version
2. Get latest GitHub release version:
   ```bash
   gh release list --limit 1 --json tagName --jq '.[0].tagName'
   ```
3. Check git status:
   ```bash
   git status --porcelain
   git branch --show-current
   git fetch origin && git status -uno
   ```
4. Get commits since last release:
   ```bash
   git log v{GITHUB_VERSION}..HEAD --oneline --no-merges
   ```
5. Generate release notes: one bullet per commit, excluding:
   - Version bump commits ("Release v...")
   - Trivial commits (typos, formatting)
   - Skill/agent updates (.claude directory changes)
   - CI/workflow updates (.github directory changes)

**Abort conditions:**
- Working tree has uncommitted changes
- Not on master branch
- Local branch is behind remote

### Phase 2: Determine New Version

1. Parse current version as major.minor.patch
2. If user provided version argument, validate format (X.Y.Z) and use it
3. Otherwise, calculate patch increment (e.g., 0.1.0 -> 0.1.1)
4. Use AskUserQuestion to confirm version:
   - Header: "Version"
   - Question: "Release version {SUGGESTED}?"
   - Options: suggested version, minor bump, major bump

### Phase 3: Present Execution Plan

Display this exact format and ask for confirmation:

```
## Release Plan: v{NEW_VERSION}

### Current State
- Local version: {LOCAL_VERSION}
- GitHub release: {GITHUB_VERSION}
- Git: clean, on master, up-to-date

### Release Notes
- {BULLET_POINT_PER_SIGNIFICANT_COMMIT}
- ...

### Steps to Execute
1. Update pyproject.toml version: {OLD} -> {NEW}
2. Run `uv sync` to update lockfile
3. Run `uv run pytest` to verify tests pass
4. Run `uvx ruff check` to verify lint passes
5. Commit: "Release v{NEW_VERSION}"
6. Push to origin/master
7. Trigger GitHub Action "Create Release"
8. Wait for release workflow to complete
9. Wait for publish workflow to complete (PyPI)

Proceed with release?
```

Use AskUserQuestion with options: "Yes, proceed" / "No, abort"

### Phase 4: Execute

Execute each step, reporting progress. Stop immediately on any failure.

#### Step 1: Update version
Use Edit tool to update version in pyproject.toml.

#### Step 2: Sync lockfile
```bash
cd python && uv sync
```

#### Step 3: Run tests
```bash
cd python && uv run pytest
```
Abort if tests fail.

#### Step 4: Lint
```bash
cd python && uvx ruff check src tests --extend-select I,B,SIM,C4,ISC,PIE
```
Abort if lint fails.

#### Step 5: Commit
```bash
git add python/pyproject.toml python/uv.lock
git commit -m "Release v{VERSION}"
```

#### Step 6: Push
```bash
git push origin master
```

#### Step 7: Trigger GitHub Action
```bash
gh workflow run release.yml -f version={VERSION}
```

#### Step 8: Wait for release workflow
```bash
gh run list --workflow=release.yml --limit 1 --json status,conclusion,databaseId
```
Poll every 10 seconds until status is "completed". Abort if conclusion is not "success".

#### Step 9: Wait for publish workflow
```bash
gh run list --workflow=publish.yml --limit 1 --json status,conclusion,databaseId
```
Poll every 10 seconds until status is "completed". Abort if conclusion is not "success".

### Completion

Report success with:
- GitHub release URL: `gh release view v{VERSION} --json url --jq '.url'`
- PyPI URL: https://pypi.org/project/piezod/{VERSION}/

## Error Messages

Provide clear, actionable error messages:

- **Dirty working tree**: "Uncommitted changes detected. Commit or stash changes before releasing."
- **Wrong branch**: "Currently on branch {X}. Switch to master before releasing."
- **Behind remote**: "Local master is behind origin/master. Pull latest changes first."
- **Tests failed**: "Tests failed. Fix test failures before releasing."
- **Lint failed**: "Lint failed. Run `uvx ruff check src tests --fix` to fix issues."
- **Workflow failed**: "GitHub workflow failed. Check Actions tab for details."
