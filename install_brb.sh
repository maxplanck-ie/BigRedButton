#!/usr/bin/env bash
#
# install_brb.sh - create a version-named conda env for BigRedButton and
# install a tagged release into it non-editably. Not part of the Python
# package; a local ops convenience script run after pulling a new release.
#
# Usage: ./install_brb.sh [-t <tag>] [-p <prefix>] [-k <n>] [-f] [-h]
set -euo pipefail
usage() {
  cat <<'EOF'
Usage: install_brb.sh [-t <tag>] [-p <prefix>] [-k <n>] [-f] [-h]
Create a version-named conda env (e.g. BRB_v0.5.0) and install
BigRedButton into it from a git tag, then prune old versioned envs.
Options:
  -t <tag>     Install this tag instead of the latest reachable tag
               (e.g. v0.5.0). Default: latest tag from `git describe`.
  -p <prefix>  Env name prefix. Default: BRB_v
  -k <n>       Number of versioned envs to keep (newest first). Default: 5
  -f           Force: remove the target env first if it already exists.
  -h           Show this help and exit.
EOF
}
prefix="BRB_v"
keep=5
tag=""
force=0
while getopts ":t:p:k:fh" opt; do
  case "$opt" in
    t) tag="$OPTARG" ;;
    p) prefix="$OPTARG" ;;
    k) keep="$OPTARG" ;;
    f) force=1 ;;
    h) usage; exit 0 ;;
    \?) echo "Unknown option: -$OPTARG" >&2; usage; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument" >&2; usage; exit 1 ;;
  esac
done
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"
command -v conda >/dev/null 2>&1 || { echo "conda not found on PATH" >&2; exit 1; }
command -v git   >/dev/null 2>&1 || { echo "git not found on PATH" >&2; exit 1; }
[ -d .git ] || { echo "$script_dir is not a git repository" >&2; exit 1; }
git fetch --tags --quiet
if [ -z "$tag" ]; then
  tag="$(git describe --tags --abbrev=0)"
fi
if ! git rev-parse "refs/tags/$tag" >/dev/null 2>&1; then
  echo "Tag '$tag' not found (did you git fetch --tags?)" >&2
  exit 1
fi
version="${tag#v}"
envname="${prefix}${version}"
# Build from a throwaway worktree at the tag rather than checking out this
# repo in place - avoids touching the caller's working tree (which may hold
# this very script, absent from older tags) or requiring it to be clean.
worktree_dir="$(mktemp -d "${TMPDIR:-/tmp}/brb_install_${version}.XXXXXX")"
cleanup() { git worktree remove --force "$worktree_dir" >/dev/null 2>&1 || rm -rf "$worktree_dir"; }
trap cleanup EXIT
echo "Checking out $tag into $worktree_dir..."
git worktree add --quiet --detach "$worktree_dir" "$tag"
[ -f "$worktree_dir/pyproject.toml" ] || { echo "pyproject.toml not found at $tag" >&2; exit 1; }
# `conda env list` can lag behind the actual env directories on some shared
# filesystems (observed in practice: a directory conda itself later refused
# to recognize as an environment, and so wouldn't remove). Check the actual
# directory on disk too, not just conda's own registry.
envs_base="$(conda info --base)/envs"
env_dir="$envs_base/$envname"
env_exists=0
{ conda env list | awk '{print $1}' | grep -qx "$envname"; } || [ -d "$env_dir" ] && env_exists=1
if [ "$force" -eq 1 ]; then
  echo "Removing existing env $envname (forced, if present)..."
  conda env remove -n "$envname" --yes >/dev/null 2>&1 || true
  # Belt and suspenders: conda env remove is a no-op (not an error) on an
  # unregistered directory conda doesn't recognize as an environment, so
  # fall back to removing it directly if it's still there.
  [ -d "$env_dir" ] && rm -rf "$env_dir"
elif [ "$env_exists" -eq 1 ]; then
  echo "Env '$envname' already exists. Use -f to remove and recreate it." >&2
  exit 1
fi
# BRB's dependencies (see pyproject.toml) are all pure-pip installable --
# unlike dissectBCL, it needs no bioconda-only binaries, so the env just
# needs a Python interpreter, not a full env.yml. git is required too: BRB
# uses setuptools-scm for dynamic versioning, which shells out to `git
# describe` during `pip install`, and pip's isolated build environment only
# sees binaries present in this conda env (see the LookupError note in
# README.md, and the same requirement in dissectBCL's env.yml).
echo "Creating env $envname..."
conda create -n "$envname" "python>=3.10" pip git --yes --quiet
echo "Installing BigRedButton $tag into $envname..."
conda run -n "$envname" pip install "$worktree_dir"
echo "Pruning old '$prefix*' envs, keeping newest $keep..."
mapfile -t old_envs < <(
  conda env list | awk '{print $1}' \
    | grep -E "^${prefix}[0-9]" \
    | sort -Vr \
    | tail -n +"$((keep + 1))"
)
for old in "${old_envs[@]:-}"; do
  [ -z "$old" ] && continue
  echo "Removing old env: $old"
  conda env remove -n "$old" --yes
done
cat <<EOF
Installed BigRedButton $tag into env '$envname'.
To use it: conda activate $envname
Remember to restart your BigRedButton process in the new env.
EOF
