#!/usr/bin/env bash
# Act on Claude's review of an automated pull request: commit the reviewer's
# edits and removals, post the review as a comment, then merge the pull request
# (when AUTO_MERGE is true), leave it for a person, or close it if nothing is
# left.
#
# The workflow's review job runs this as a Claude Code Stop hook, inside the
# Claude step, so everything it does uses Claude's GitHub login (claude[bot]).
#
# Settings come from the workflow's env:
#   PR_NUMBER   pull request under review
#   PR_HEAD     its branch
#   PR_PATHS    space-separated paths the update may change
#   PR_CHECK    command that checks the update
#   PR_BASE     branch the pull request targets
#   AUTO_MERGE  "true" to merge when entries remain after review
#   START_SHA   commit the review started from (the PR head)
# Claude writes its review to /tmp/review.md.
#
# Exit codes follow Claude Code's hook rules: 0 done, 2 the check failed and
# the output goes back to Claude (at most 3 tries), 1 anything else.
set -uo pipefail
exec 1>&2  # Claude reads hook feedback from stderr
cd "$(dirname "$0")/.."

: "${PR_NUMBER:?}" "${PR_HEAD:?}" "${PR_PATHS:?}" "${PR_CHECK:?}" "${PR_BASE:?}" "${START_SHA:?}"
read -ra PATHS <<< "$PR_PATHS"
DONE=/tmp/review.done
TRIES_FILE=/tmp/review.tries
MAX_TRIES=3

git fetch -q origin "$PR_BASE" || exit 1

tries=$(( $(cat "$TRIES_FILE" 2>/dev/null || echo 0) + 1 ))
echo "$tries" > "$TRIES_FILE"
if ! bash -c "$PR_CHECK"; then
  if [ "$tries" -lt "$MAX_TRIES" ]; then
    echo "A check failed (try $tries of $MAX_TRIES). Fix your change, then finish again."
    exit 2
  fi
  echo "Checks still fail after $MAX_TRIES tries. The review was not applied."
  exit 1
fi

set -e

git reset -q "$START_SHA"
if [ -n "$(git status --porcelain -- "${PATHS[@]}")" ]; then
  git add -- "${PATHS[@]}"
  git -c user.name="claude[bot]" \
      -c user.email="209825114+claude[bot]@users.noreply.github.com" \
      commit -q -m "Apply review changes"
  git push -q origin "HEAD:refs/heads/$PR_HEAD"
fi

[ -s /tmp/review.md ] || echo "Reviewed; no notes." > /tmp/review.md
gh pr comment "$PR_NUMBER" --body-file /tmp/review.md

if git diff --quiet "origin/$PR_BASE" HEAD -- "${PATHS[@]}"; then
  gh pr close "$PR_NUMBER" --delete-branch --comment "Nothing passed review, so this pull request is closed."
elif [ "${AUTO_MERGE:-false}" = true ]; then
  gh pr merge "$PR_NUMBER" --squash --delete-branch
else
  echo "Review posted; leaving the pull request for a person to merge."
fi
touch "$DONE"
