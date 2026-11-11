#!/usr/bin/env bash
# Check an automated update, then commit it to a new branch and open a pull
# request for review. It never merges.
#
# The scheduled workflow runs this as a Claude Code Stop hook: Claude Code runs
# it automatically when Claude finishes. It runs inside the Claude step, so it
# uses Claude's GitHub login, and the pull request and commit show as
# claude[bot].
#
# Settings come from the workflow's env:
#   PR_BRANCH  branch name prefix; each run uses PR_BRANCH-<date>-<time>
#   PR_PATHS   space-separated paths the update may change
#   PR_CHECK   command that checks the update
#   PR_TITLE   default pull request title (Claude can override it in
#              /tmp/pr-title.txt)
#   PR_BASE    branch the pull request targets
#   START_SHA  commit the run started from
# Claude writes the pull request description to /tmp/pr-body.md.
#
# Exit codes follow Claude Code's hook rules:
#   0  published, or nothing to publish
#   2  the check failed: the output goes back to Claude, which fixes its change
#      and finishes again, which runs this script again (at most 3 tries)
#   1  anything else: Claude stops, and the workflow's last step fails the run
set -uo pipefail
exec 1>&2  # Claude reads hook feedback from stderr
cd "$(dirname "$0")/.."

: "${PR_BRANCH:?}" "${PR_PATHS:?}" "${PR_CHECK:?}" "${PR_TITLE:?}" "${PR_BASE:?}" "${START_SHA:?}"
read -ra PATHS <<< "$PR_PATHS"
DONE=/tmp/publish.done
TRIES_FILE=/tmp/publish.tries
MAX_TRIES=3

git fetch -q origin "$PR_BASE" || exit 1

tries=$(( $(cat "$TRIES_FILE" 2>/dev/null || echo 0) + 1 ))
echo "$tries" > "$TRIES_FILE"
if ! bash -c "$PR_CHECK"; then
  if [ "$tries" -lt "$MAX_TRIES" ]; then
    echo "A check failed (try $tries of $MAX_TRIES). Fix your change, then finish again."
    exit 2
  fi
  echo "Checks still fail after $MAX_TRIES tries. Nothing was published."
  exit 1
fi

set -e

# Keep only the allowed paths. If Claude committed on its own, fold that back
# into the working tree first.
git reset -q "$START_SHA"
if [ -z "$(git status --porcelain -- "${PATHS[@]}")" ]; then
  echo "No changes to publish."
  touch "$DONE"
  exit 0
fi

TITLE=$(head -n1 /tmp/pr-title.txt 2>/dev/null || true)
TITLE=${TITLE:-$PR_TITLE}
BRANCH="$PR_BRANCH-$(date -u +%Y%m%d-%H%M%S)"

git checkout -q -b "$BRANCH"
git add -- "${PATHS[@]}"
git -c user.name="claude[bot]" \
    -c user.email="209825114+claude[bot]@users.noreply.github.com" \
    commit -q -m "$TITLE"
git push -q origin "HEAD:refs/heads/$BRANCH"

[ -s /tmp/pr-body.md ] || echo "Automated update." > /tmp/pr-body.md
gh pr create --base "$PR_BASE" --head "$BRANCH" --title "$TITLE" --body-file /tmp/pr-body.md \
  | tee /tmp/publish.pr  # the PR link, for the workflow's review job
touch "$DONE"
