# Add new entries to this list

This is the procedure for the weekly automated list update. Each run adds up
to five new, high-quality entries to `README.md` in one pull request. A person
reviews every pull request.

## Scope

Software tools for analyzing chromatin accessibility data: ATAC-seq,
single-cell ATAC-seq, and closely related assays such as DNase-seq and
MNase-seq when the tool is used for accessibility analysis. That covers
quality control, processing pipelines, peak calling, differential
accessibility, motif and footprinting analysis, nucleosome positioning, region
enrichment, and single-cell analysis.

Not in scope: general-purpose genomics tools that have no accessibility-specific
feature, commercial services without a usable tool, tutorials, and datasets.

## Steps

1. Read `README.md`: the intro, the sections, and every existing entry.
2. Read `/tmp/past-proposals.txt`. It lists entries proposed before in pull
   requests that were not merged: either declined or still waiting for review.
   Never propose any of them again, under any name.
3. Look for tools that fit the scope and are not listed yet. Favor recent
   ones with a paper (journal or preprint) and a public code repository or
   model release. Use official sources only: the project's own repository,
   paper, or model page.
4. Keep only candidates that clearly fit the scope and one existing section,
   are not already listed (check names, renamed projects, and links), and have
   a working official link. When in doubt, leave it out. Adding nothing is
   fine.
5. Add at most five entries to `README.md`:
   - Format: `- [Name](link): one-line description`, in the style of the
     neighboring entries.
   - Link to the code repository (GitHub, GitLab, Bioconductor, PyPI or CRAN page).
   - Put each entry in the best existing section, in case-insensitive
     alphabetical order.
   - Do not change, move, or delete existing lines, and do not add sections.
6. Write `/tmp/pr-title.txt` with one line, for example
   `Add Foo, Bar and Baz`.
7. Write `/tmp/pr-body.md` with:
   - a table of the added entries: name, section, why it qualifies, and the
     primary sources (repository, paper) as links;
   - a short list of candidates you considered and left out, each with the
     reason.
8. If nothing qualifies, change no files.

When you finish, a hook runs `scripts/check_list.py`, then commits your change
and opens the pull request. If the check fails, it tells you what to fix: fix
`README.md` and finish again. Never commit, push, or open a pull request
yourself, and never merge one.
