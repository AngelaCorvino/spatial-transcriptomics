# Project instructions

Before creating a module, look for an existing workflow with the same inputs,
processing pattern, and output style. Extend it with minimal code and reuse its
helpers and conventions; create a standalone module only when necessary.

Keep code simple and readable: use clear names, straightforward control flow,
and small, focused functions. Avoid unnecessary abstractions and clever shortcuts.

Develop analysis code locally; run analyses of the experimental data on the
cluster. Keep reports and guidelines concise and avoid repeating information.

Keep experimental data, generated results, and reports outside the code repository.
Resolve local paths from the gitignored `local_config.yaml`; preserve cluster paths.

## Experimental report updates

Update the experimental report only when the user explicitly requests it.
Scientific discussions, code changes, and new images do not trigger an update.

For requested updates, use the `experimental_report` agent defined in
`.codex/agents/experimental_report.toml`. Pass the relevant conversation context,
confirmed decisions, evidence paths, and any unresolved points. When delegation
is unavailable or has no independent work to run alongside it, follow that file's
instructions directly. Keep one report editor active at a time; the report editor
must not delegate another report update.

Read the current conversation; do not search unrelated chats. Mention a report update
briefly in the final response. Do not commit or publish it unless requested.
