#!/usr/bin/env python3
"""Check the entries the weekly list update adds to README.md.

    python3 scripts/check_list.py                   check README.md
    python3 scripts/check_list.py --past-proposals  print entries from unmerged PRs

The check compares README.md with the base branch (PR_BASE, default master)
and looks only at what changed. It fails if any existing line was changed or
removed, or if a new entry:
  - is not written as "- [Name](https://link): description"
  - repeats an entry already in the list (same name or same link)
  - was proposed before in a pull request that was not merged
    (/tmp/past-proposals.txt, written by --past-proposals)
  - is out of alphabetical order next to its neighbors in its section
  - has a link that does not load

--past-proposals lists every entry added by an open or closed-unmerged pull
request (except PR_NUMBER, the one under review), one per line: name, link
and PR number, separated by tabs.
"""
import json
import os
import re
import subprocess
import sys
import urllib.error
import urllib.request

README = "README.md"
PAST = "/tmp/past-proposals.txt"
ENTRY = re.compile(r"^- \[(?P<name>[^\]]+)\]\((?P<url>https?://[^)\s]+)\)(?::\s+(?P<desc>\S.*))?$")


def norm_name(name):
    return re.sub(r"[^a-z0-9]", "", name.casefold())


def norm_url(url):
    url = re.sub(r"^https?://(www\.)?", "", url.strip().casefold())
    return url.rstrip("/#")


def run(*cmd):
    return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout


def entries(text):
    """(section, name, url, line) for every list entry, in order."""
    out, section = [], None
    for line in text.splitlines():
        if line.startswith("## "):
            section = line[3:].strip()
        m = ENTRY.match(line)
        if section and m:
            out.append((section, m["name"], m["url"], line))
    return out


def past_proposals():
    prs = json.loads(run("gh", "pr", "list", "--state", "all", "--limit", "1000",
                         "--json", "number,state"))
    current = os.environ.get("PR_NUMBER")  # the pull request under review, if any
    for pr in prs:
        if pr["state"] == "MERGED" or str(pr["number"]) == current:
            continue
        diff = run("gh", "pr", "diff", str(pr["number"]))
        for line in diff.splitlines():
            m = ENTRY.match(line[1:]) if line.startswith("+") else None
            if m:
                print(f"{m['name']}\t{m['url']}\t#{pr['number']} {pr['state'].lower()}")


def link_ok(url):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (list link check)"})
    try:
        with urllib.request.urlopen(req, timeout=20):
            return True, ""
    except urllib.error.HTTPError as e:
        if e.code in (401, 403, 429):  # sites that block bots; not proof of a dead link
            return True, f"warning: {url} answered {e.code}"
        return False, f"{url} answered {e.code}"
    except Exception as e:  # DNS failure, timeout, bad certificate
        return False, f"{url} did not load ({e.__class__.__name__})"


def check():
    base = os.environ.get("PR_BASE", "master")
    old = run("git", "show", f"origin/{base}:{README}")
    new = open(README, encoding="utf-8").read()
    errors = []

    diff = run("git", "diff", "--unified=0", f"origin/{base}", "--", README)
    added = []
    for line in diff.splitlines():
        if line.startswith(("+++", "---")):
            continue
        if line.startswith("-"):
            errors.append(f"Existing line changed or removed: {line[1:]!r}. Only add new entries.")
        elif line.startswith("+"):
            text = line[1:]
            if not text.strip():
                continue
            if ENTRY.match(text) and ENTRY.match(text)["desc"]:
                added.append(text)
            else:
                errors.append(f"Not a valid entry: {text!r}. Use '- [Name](https://link): description'.")

    past = []
    if os.path.exists(PAST):
        for row in open(PAST, encoding="utf-8"):
            parts = row.rstrip("\n").split("\t")
            if len(parts) >= 2:
                past.append(parts)

    all_entries = entries(new)
    for i, (section, name, url, line) in enumerate(all_entries):
        if line not in added:
            continue
        for j, (_, other, other_url, _) in enumerate(all_entries):
            if j != i and (norm_name(other) == norm_name(name) or norm_url(other_url) == norm_url(url)):
                errors.append(f"{name} repeats the existing entry {other} ({other_url}).")
        seen = sorted({rest[0].split()[0] if rest else "?" for p_name, p_url, *rest in past
                       if norm_name(p_name) == norm_name(name) or norm_url(p_url) == norm_url(url)})
        if seen:
            errors.append(f"{name} was proposed before and not merged ({', '.join(seen)}). Do not propose it again.")
        neighbors = [e for e in all_entries if e[0] == section]
        k = [e[3] for e in neighbors].index(line)
        if k > 0 and neighbors[k - 1][1].casefold() > name.casefold():
            errors.append(f"{name} should come before {neighbors[k - 1][1]} (alphabetical order in '{section}').")
        if k + 1 < len(neighbors) and name.casefold() > neighbors[k + 1][1].casefold():
            errors.append(f"{name} should come after {neighbors[k + 1][1]} (alphabetical order in '{section}').")
        ok, message = link_ok(url)
        if not ok:
            errors.append(message)
        elif message:
            print(message)

    if old == new:
        print("README.md is unchanged.")
    for e in errors:
        print(f"ERROR: {e}")
    if not errors:
        print(f"OK: {len(added)} new entr{'y' if len(added) == 1 else 'ies'} checked.")
    return 1 if errors else 0


if __name__ == "__main__":
    if sys.argv[1:] == ["--past-proposals"]:
        past_proposals()
    else:
        sys.exit(check())
