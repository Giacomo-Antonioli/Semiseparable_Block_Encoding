#!/usr/bin/env bash
# Build the zipped snapshot of the repository for the SIMAX Supplementary
# Materials (code + data, including the git submodules, without .git folders).
#
# Usage:  scripts/make_supplement_zip.sh [TAG_OR_COMMIT]   (default: HEAD)
# Output: $OUT_DIR/semisep-supplement-<ref>.zip  (OUT_DIR defaults to the parent folder)
set -euo pipefail

ref="${1:-HEAD}"
root="$(git rev-parse --show-toplevel)"
cd "$root"
name="semisep-supplement-$(echo "$ref" | tr '/' '_')"
out="$(cd "${OUT_DIR:-$root/..}" && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

mkdir -p "$tmp/$name"
git archive "$ref" | tar -x -C "$tmp/$name"

# Submodules at the commit recorded in $ref
git ls-tree -r "$ref" | awk '$2=="commit"{print $3" "$4}' | while read -r sha path; do
  mkdir -p "$tmp/$name/$path"
  git -C "$path" archive "$sha" | tar -x -C "$tmp/$name/$path"
done

(cd "$tmp" && zip -qr "$out/$name.zip" "$name")
echo "Wrote $out/$name.zip"
