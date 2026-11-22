# For tool developers: tell us how to run your tool better

**We want to hear from you.** This benchmark is only fair if every tool is run the way its developers would run it.
We try, but we are not the experts on your tool. If anything below could be better, please say so, even if it
makes your tool look different from what we currently show:

- the command line and parameters we use (see `lib/TrCommands.groovy`, one block per tool);
- the catalog format we convert to for your tool, and how we read its output (`evaluate/calls.py`);
- the container or build we use, and whether the version is the right one to test;
- every other change we made to get your tool running, listed with the reasons in
  [ADAPTATIONS.md](ADAPTATIONS.md);
- how we report your tool (motif-length classes, what counts as a call, what we leave out).

## How
Open an issue in this repository ([new issue](../../issues/new/choose): "How should my tool be run?" or "Comment on an
adaptation"), or a pull request against `lib/TrCommands.groovy` / `conf/versions_current.config`. Short is fine:
"use `--foo 5`, because ..." is enough. An email or a message to the maintainer works too.

## What happens with it
- Your suggestion is discussed in the open, in the issue.
- If it is accepted, the tool is **re-run** with the change, the result is published next to the earlier one (earlier
  results stay visible with the settings they were run with), and the change is recorded in ADAPTATIONS.md.
- If we do not follow a suggestion, we say why in the issue.
- New releases of your tool: tell us and we add them (the rule is in FORK.md; a patch release is re-run when its
  developer says it changes results, as for STRdust 1.4.2).
- Nothing about how a tool is run is changed silently; everything is in the git history and the run log.

## Conflict of interest
The maintainer of this fork develops STRdust and inquiSTR. Those two tools are run under exactly the same rules as all
others, and their developers' suggestions (including the maintainer's own) go through the same public issues, so that anyone
can check them.
