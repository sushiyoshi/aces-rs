#!/usr/bin/env bash
set -euo pipefail

remote="${1:-origin}"
branch="${2:-$(git rev-parse --abbrev-ref HEAD)}"

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "Error: this script must be run inside a Git repository." >&2
  exit 1
fi

if ! git remote get-url "$remote" >/dev/null 2>&1; then
  cat <<MSG
Remote '$remote' is not configured.
Add it with:
  git remote add $remote <REMOTE_URL>
MSG
  exit 2
fi

remote_url=$(git remote get-url "$remote")
echo "Checking Git access for remote '$remote' ($remote_url) and branch '$branch'..."

read_status=0
if git ls-remote "$remote" >/dev/null 2>&1; then
  echo "✔ Read access confirmed (git ls-remote succeeded)."
else
  echo "✖ Unable to read from remote '$remote'."
  read_status=1
fi

push_status=0
if git push --dry-run "$remote" "HEAD:$branch" >/dev/null 2>&1; then
  echo "✔ Push (PR) permissions look good (dry-run push succeeded)."
else
  echo "✖ Push denied. Check your credentials or branch protections."
  push_status=1
fi

if (( read_status != 0 || push_status != 0 )); then
  echo
  echo "Suggested remediation steps:"
  echo "  1. Verify the remote URL and your authentication method (SSH/HTTPS)."
  echo "  2. Ensure your token or SSH key has 'repo' scope and push permissions."
  echo "  3. If using GitHub, confirm you are a collaborator or fork the repository."
  exit 3
fi

echo "All checks passed. You should be able to create a PR from branch '$branch'."
