---
title: Release and publish a version
type: workflow
area: release
created: 2026-08-15
updated: 2026-08-15
---

# Release and publish a version

> [!note] Why this note exists
> `v1.3.1` was tagged and pushed to `origin` but never appeared on the
> [releases page](https://github.com/Deltares-research/weathergenr/releases) —
> the tag was there, the GitHub **Release** was not. `git push` creates tags;
> it never creates Releases. That last step (7) is the one that gets forgotten,
> and it is the reason this sequence is written down.

## Scope

The full owner-run gate, from an unbumped development version to a published
release. Every step is deliberate: `.git-workflow.yml` sets `cadence: manual`,
`auto_push: false`, and `release:` null, so nothing here fires automatically.

Command definitions are **not** restated here. The bump and validate
invocations live in `.git-workflow.yml` under the `weathergenr` adapter's
`apply:` and `validate:` keys; read them there so there is one definition, not
two.

## Sequence

1. **Bump the version.** Run the adapter's `apply:` command with the component
   (`patch`, `minor`, `major`) as a trailing argument. It rewrites
   `DESCRIPTION#Version` and the `NEWS.md` "(development version)" heading in
   one call, and aborts on a dirty tree — which is the correct guard here.

2. **Confirm `NEWS.md` is complete.** Entries are written as work lands, not
   reconstructed from `git log` now. Verify the section covers every
   user-visible change since the previous tag:

   ```powershell
   git log --oneline v<previous>..HEAD
   ```

   Internal-only commits (tooling, CI, refactors, roxygen) correctly have no
   entry.

3. **Run the release gate.** The adapter's `validate:` command — `check_only()`
   with vignettes, matching CI. It is a no-write validation path by design.

4. **Run the baseline gate if numeric output could have moved.** Resampling,
   WARM, wavelets, calendar logic, or evaluation statistics — see the *Workflow*
   section of `AGENTS.md`. It is opt-in and local; CI skips it, so a green CI
   does not mean it was checked.

5. **Commit and tag.** Tags in this repo are **annotated** (`git cat-file -t`
   returns `tag`, not `commit`) — keep them that way; `remotes` and the releases
   UI both read the tag message.

   ```powershell
   git commit -m "chore(release): bump version to <X.Y.Z>"
   git tag -a v<X.Y.Z> -m "weathergenr <X.Y.Z>"
   ```

6. **Push the branch, then the tag.** Two separate pushes — `git push` alone
   does not carry tags.

   ```powershell
   git push origin master
   git push origin v<X.Y.Z>
   ```

   At this point the tag is consumable: `remotes::install_github()` resolves
   against tags, not Releases, so downstream can already pin the new version.

7. **Create the GitHub Release.** *The step that was missing for 1.3.1.*
   A Release is a GitHub-side artifact created separately from the tag; without
   it the version is invisible on `/releases` and shows only under `/tags`.

   ```powershell
   gh release create v<X.Y.Z> --title "v<X.Y.Z>" --notes-file <notes.md>
   ```

   House style for the notes body, set by `v1.3.0` and `v1.3.1`: open with a
   link to `NEWS.md` at that tag, then `## Breaking changes` / `## Bug fixes` /
   `## New features` / `## Performance` as applicable, bolding the consequence a
   reader must act on. Verify with `gh release list` — the new version should
   report `Latest`.

## After the release

- **Bump the downstream pin.** `blueearth_cst` installs this package by Git tag
  from `dev/scripts/install_weathergenr.R`, driven by `pixi run install-rdeps`.
  Until that pin moves, the release is invisible downstream. That edit lives in
  the other repo.
- **Return to a development version** when the next change lands — the adapter's
  `apply:` command with the `dev` component.

## Constraints

- Published tags are consumed artifacts: **never move, delete, or re-point one.**
- `v0.1.1`, `v0.2.0`, `v0.9.0`, and `v1.0.0` predate any deliberate scheme. They
  are historical record — do not retrofit or re-tag them, and the gap where a
  `v1.1.0` would sit is not a defect. `.git-workflow.yml` and `NEWS.md` carry
  the same boundary.
- A release is **not** delivered at the commit. With `auto_push: false` the push
  and the Release are owner-run final steps; never report a release complete
  before step 7.
