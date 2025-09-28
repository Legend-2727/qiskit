#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
git fetch --all --prune
git checkout main && git reset --hard upstream/main
git checkout caes-router-uprev
git rebase main || { echo "Resolve conflicts, then: git add -A && git rebase --continue"; exit 1; }
echo "Rebase complete. Push with: git push -f fork caes-router-uprev"