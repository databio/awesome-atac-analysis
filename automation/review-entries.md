# Review new list entries

This is the procedure for the second layer of the weekly list update. A
smaller model proposed new entries for `README.md` in an open pull request.
Your job is to check each one carefully and keep only the ones that clearly
belong. Be strict: a wrong entry is worse than a missing one.

## Steps

1. See what the pull request adds: `git diff origin/$PR_BASE -- README.md`.
   Read `automation/add-entries.md` for the scope and the entry format.
2. For each added entry, open its link and its paper yourself (do not rely on
   search snippets) and check that:
   - the link is the project's official code repository or model page, and it
     is the project the entry names;
   - it clearly fits the scope and the section it was put in;
   - it is not already listed, under its own name or another (same code
     repository, same model family, or a rename);
   - it is not in `/tmp/past-proposals.txt` (proposed before and not merged);
   - the one-line description is accurate and in the list's style;
   - it is real and published (a paper or preprint) or widely used.
3. Delete from `README.md` every added entry that fails any check. You may fix
   an accurate entry's description or move it to a better section. Do not add
   new entries, and do not touch lines that were in the list before.
4. Write `/tmp/review.md`: a short table with each entry, keep or remove, the
   reason, and the sources you opened.

When you finish, a hook re-runs `scripts/check_list.py`, commits your
removals, posts your review on the pull request, and then merges it, or closes
it if nothing is left. If the check fails, it tells you what to fix. Never
commit, push, comment on, merge, or close the pull request yourself.
