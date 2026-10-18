# Repository-readiness audit

Scope: ordinary R package work with Codex; Claude/Copilot projections inspected without launching runtimes.
Concurrency assumption: one writer; concurrent writers require separate worktrees under repository policy.
Verdict: **ready** for ordinary Codex R package work. Repository repairs and vignette-enabled checks pass.
Fresh Claude/Copilot execution and scientific validity remain outside this readiness claim.
Audit date: 2026-10-18. Supersedes: none found.

This project-system checklist uses the review home declared in dev/README.md. Completed inspection does
not mean findings have been repaired. The user subsequently authorized repairs; see repair verification below.

## Inspection

### Baseline and authority

- [x] Workload, runtimes, and concurrency established as above.
- [x] Supplied inherited instructions, local AGENTS, README, dev/README, and workflow policy inspected.
- [x] Live project-system R template compared; preserve the coherent existing package layout.
- [x] Brain roots resolve; linked skill/role sources checked as tracked in infrastructure Git.
- [x] Limits: no fresh runtimes, hook execution, clean environment restore, or package certification.

### Instructions and ownership

- [x] Canonical AGENTS and thin CLAUDE importer; instruction precedence and package ownership are explicit.
- [x] Tracked authored manifest distinguished from ignored linked skills and generated runtime profiles.
- [x] Generated files, public API, RNG rules, verification gates, and publishing boundaries are actionable.
- [x] Process drift identified: NEWS heading and conflicting lint baseline comments (F5).

### Roles, skills, and activation

- [x] Project-scoped status: six roles, 22 skill links, 13 provider artifacts; all resolve.
- [x] Conditional bindings resolve through provider paths; absence from explicit catalog is not a defect.
- [x] Codex registrations/TOML and Copilot profiles/catalog pass compatibility checks.
- [x] Generated child permissions inherit or narrow; concurrent-writing ownership/isolation is explicit.
- [x] Clone prerequisites exist locally; repository lacks activation bootstrap instructions (F3).
- [x] Fresh Claude/Copilot execution unverified; static compatibility is the extent of the claim.

### Hooks and permissions

- [x] Project permissions/registrations and inherited runtime configurations inspected separately.
- [x] No project hooks/spec declared; no required project hook missing. Git has sample hooks only.
- [x] User Claude/Codex hooks register PostToolUse and SessionStart; no hook payload executed.
- [x] Inspected Claude Python hook links resolve into infrastructure; full inherited execution unverified.
- [x] Codex PostToolUse matcher is Bash; coverage of exec_command unverified, not assumed enforced.
- [x] Effective session: workspace writes, restricted network, auto-reviewed escalation.
- [x] Audit preserved guards; no permission changes, hook installation, reset, fetch, or push.

### Git and concurrent work

- [x] Root weathergenr; master; HEAD 8690fcc2b7c4978f9577af7324b2a995f8e2f62d; one worktree.
- [x] Empty initial index; pre-existing unstaged manifest edits preserved with no observed movement.
- [x] Mode full aliases branch; auto_push false; manual releases; low-risk record uses primary checkout.
- [x] Clone-local hooks inspected; no mandatory Git hook contract found.
- [x] Remote need settled by public/shared repository, CI, and downstream tags; origin/master tracks master.
- [x] Remote probe attempted; local Git shell failed before reachability could be established (F6).
- [x] core.autocrlf false; git diff --check passed; no observed line-ending churn.

### Ignore rules, packaging, and disposable state

- [x] Effective ignores cover scratch, runtime projections, .Rhistory, and .Renviron.
- [x] Authored manifest and DESCRIPTION remain visible; .env/.env.local also remain visible (F1).
- [x] Build exclusions cover dev/tools/projections and preserve R/tests; .agent-system lacks exclusion (F2).
- [x] Tracked docs and vignette freeze are intentional local conventions.
- [x] No audit scratch created; existing scratch neither traversed nor deleted.

### Environment, verification, and handoff

- [x] R 4.4.3 and Quarto CLI 1.9.38 available; R working directory is this repository.
- [x] All 16 Imports load; installed weathergenr loads; R quarto package unavailable (F4).
- [x] R library belongs to blueearth_cst session-2; locale startup warnings observed (F4).
- [x] DESCRIPTION dependency declarations inspected; no project lock/environment spec found.
- [x] Check/lint CI, .lintr, build helpers, metadata, license, and generated-artifact rules inspected.
- [x] Verification commands inspected; package tests/checks/baselines not run for this audit-only record.
- [x] Existing board/ledger sufficient; no task needed. This is the only durable audit deliverable.

## Repairs and decisions

### F1 — Material weakness; observed: environment-file ignores

Git ignores .Renviron but exposes .env/.env.local. Future credentials could enter commits; secret contents
were not inspected. Repo maintainer: add environment-file rules with sanitized-example exceptions;
verify with git check-ignore for private files and examples.

### F2 — Material weakness; observed: package boundary

.Rbuildignore lacks an explicit .agent-system exclusion. Repo-only metadata remains eligible for copying
by build tooling; actual tarball inclusion was not tested. Repo maintainer: add an anchored exclusion
and inspect a package build through its owning tool.

### F3 — Material weakness; inspected: activation bootstrap

dev/README lacks clone-to-activation instructions; brain is absent from this shell PATH. Existing linked
configuration depends on the external infrastructure checkout. Repo maintainer: document roots, installer,
project status entry point, and board fallback; verify from a clean clone/session.

### F4 — Material weakness; observed: development environment

R borrows another project's library and lacks quarto, the declared VignetteBuilder. Documentation/release
readiness cannot be claimed. Environment owner: establish documented R setup and install declared tooling;
verify quarto namespace and the documentation gate. Locale startup warnings also need environment review.
No lock file is not independently an R-package readiness blocker; exact numerical replay is unverified.

### F5 — Material weakness; inspected: process drift

AGENTS requires NEWS entries under a development heading, but NEWS starts at 2.0.1 without that heading.
Workflow policy says lint has 28 baseline findings while lint CI says clean. Manifest header names the old
.claude path. Repo maintainer: reconcile current process facts; inspect NEWS and run lint against a current
installed package before changing policy. The downstream v1.2.0 pin was not checked against its consumer.

### F6 — Unverified evidence: remote access

git ls-remote failed locally: Git sh could not create a signal pipe (Win32 error 5). This does not establish
a missing/unreachable remote. Environment owner: repeat a read-only probe in a functioning shell before
publishing-dependent work.

### F7 — Optional improvement; inspected: Copier tracking

No .copier-answers.yml; copier update cannot track template adoption. Retaining the coherent R layout is
valid. Any later adoption should render/diff selected files and record portable answers separately.

## Check evidence

Candidate: HEAD 8690fcc with pre-existing manifest edits (design-skill rename and two explicit R-skill
removals). Initial index empty. Only this review record was authored by the audit.

- `git status --short`, `git diff --stat`, `git diff -- .agent-system/agent-manifest.yml`,
  `git worktree list --porcelain`, `git branch -vv`, `git log -1 --format='%h %s'`:
  one worktree, master four commits ahead of origin/master, unstaged manifest initially.
- `python <BRAIN_INFRA_ROOT>/brain.py status --agent-system --project <repository> --detail`:
  exit 0. Codex: 5 ok / 2 notes / 0 warnings / 0 failures. Copilot: 7 ok / 2 notes / 0 warnings / 0 failures.
  Notes cover optional hook specs and provider-routed bindings.
- `git ls-remote origin HEAD`: exit 1; local Git sh failure; remote reachability unknown.
- `git check-ignore -v -- .tmp/a .codex/config.toml .agents .Rhistory .Renviron .env .env.local
  .env.example .agent-system/agent-manifest.yml DESCRIPTION`: expected local-output matches; no .env matches.
  Earlier nested symlink probe failed with beyond-a-symbolic-link error; corrected to .agents boundary.
- `Rscript -e <dependency probe>`: R 4.4.3, correct cwd, session-2 library; all Imports true; all probed
  Suggests except quarto true. Installed weathergenr namespace load separately passed.
- `quarto --version`: 1.9.38.
- `Rscript --vanilla -e <filtered build-regex probe>`: .agent-system excluded false;
  .agents/dev/tools true; R/tests false. Per-pattern evaluation removed blank/comment rules.
  Earlier diagnostic regex probes were invalid; only this corrected probe supports F2.
- `git config --get core.hooksPath`, `core.excludesFile`, `core.autocrlf` and hooks inventory:
  no custom hook path; global exclusions inspected; autocrlf false; sample Git hooks only.
- `git -C <BRAIN_INFRA_ROOT> ls-files <canonical skill/role paths>`: inspected sources tracked there.
- `git diff --check`: exit 0 before final record verification.
- Login PowerShell initially failed on missing Env:SSL_CERT_DIR; later commands used non-login PowerShell.

No package tests, numeric baseline, R CMD check, vignette render, installation, fresh runtime launch, or
hook execution performed. Managed/inherited hook execution remains unverified.

## Completion

### Repair verification — 2026-10-18

- [x] F1 resolved: ignore .env/.env.* and .Renviron variants, with sanitized examples visible.
- [x] F2 resolved: exclude .agent-system and private environment files from package builds.
- [x] F3 resolved: dev/README documents roots, project activation/status, and the board CLI fallback.
- [x] F4 resolved for this session: standalone R 4.6.1 and its user library load all Imports and Suggests.
  Native R already has quarto installed. Process-scoped C locale removes inherited C.UTF-8 warnings.
  Documented interpreter/library checks prevent silently borrowing another project's pixi library.
- [x] F5 resolved: NEWS development heading restored; obsolete lint-count claims removed; manifest comments
  corrected without changing the user's selection edits. AGENTS points fresh setups to dev/README.
- [x] F6 resolved: read-only origin probe succeeded outside sandbox process restrictions.
- [x] F7 disposition: retain the coherent existing layout without Copier adoption; unsupported update command
  is documented. No fabricated template provenance or lock file added.
- [x] Vignette-enabled R CMD check: exit 0; zero errors, warnings, or notes.

Commands below used standalone R 4.6.1 with LC_ALL, LC_CTYPE, and LANG set to C for the process.
No persistent user locale, other project's environment, shared skills, CI workflows, or release tags changed.

- `devtools::install(build_vignettes = FALSE, upgrade = FALSE, dependencies = FALSE)`: exit 0;
  current 2.0.1 installed in the standalone R user library for worker/lint validity.
  Two initial invocations failed before installation: unsupported document argument, then string upgrade
  argument. Verified installed function signatures and corrected the documented command.
- `devtools::test()`: exit 0; 1,299 passed, zero failures/warnings, four opt-in numeric-baseline skips.
  Numerical behavior did not change; no baseline recording performed.
- `Rscript --vanilla tools/lint.R`: exit 0, No lints found, using the current installed package.
- `Rscript --vanilla -e 'source("tools/build_site_tools.R");
  check_only(build_vignettes = TRUE, document = FALSE)'`: exit 0; Status OK, zero errors/warnings/notes.
  Vignette build and rebuild, examples, and installed-package tests pass; check duration 5m 6s.
- `devtools::build(path = ".tmp/scratchpad/readiness-repair", vignettes = FALSE, manual = FALSE)`:
  exit 0. Python tarfile inspection confirms excluded metadata paths absent and R/tests present.
  Sandbox build first failed creating processx pipes; the same build succeeded with escalation.
- `git ls-remote origin HEAD`: exit 0 outside sandbox; origin HEAD ee752a080357ff4f29cdb1f1ab037988a26d6b34.
  No fetch or publishing performed; local tracking references were not refreshed.
- `python <BRAIN_INFRA_ROOT>/brain.py status --agent-system --project <repository> --detail`:
  exit 0; unchanged Codex/Copilot compatibility pass with no warnings/failures.
- `PYTHONPATH=<activated todo-board skill>; python -m todoboard --help`: exit 0.
- Root prerequisite checks against BRAIN_ROOT/wiki/index.md and BRAIN_INFRA_ROOT/brain.py pass.
- `git check-ignore -v` confirms private environment files ignored and negative example rules effective.
- `git diff --check`: passed during repairs. Final diff/index ownership checked before commit.
- Removed only the completed task-owned tarball and its empty readiness-repair scratch directory after
  recording build evidence; verified containment and absence of links before removal.

Package source, generated documentation, and numerical baselines are unchanged. Fresh runtime discovery
still requires a new session; static compatibility is not a runtime execution claim.

Repair outcome: six findings resolved; one optional template-adoption decision documented; zero blockers.
Only audit-owned paths/hunks are included in the repair commit. Original manifest selection changes remain
unstaged. No push performed.
