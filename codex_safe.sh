#!/usr/bin/env bash
set -e

if ! git diff --quiet || ! git diff --cached --quiet; then
  echo "ERROR: Working tree is not clean. Commit or stash before running Codex."
  exit 1
fi

codex edit "$@" --instruction "$(cat codex_rules.txt)"
