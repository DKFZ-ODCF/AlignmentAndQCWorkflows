# CLAUDE.md

## External-Interface Change Sign-Off

This plugin's *external interface* is: cvalues external project configs can
set/override, `<tool>` names, environment-script names, documented
file-naming conventions, and public Groovy/Java API (classes/methods other
plugins could extend or call).

Before closing any task or change in this plugin, check whether it touches
the external interface. If yes:

1. Add/update a row in `docs/developers/ExternalInterfaceChangeLog.adoc`
   describing the change and its risk if unconfirmed.
2. Get explicit user sign-off before merging that change — do not silently
   assume a breaking change is acceptable, even if it looks like dead or
   already-broken code.

Keep that table current for all changes, or whenever you notice it's
drifted from the ground truth (the code) — not just when a task happens to
touch it.

This is a standing, cross-cutting process rule, not a one-shot task — apply
it for the lifetime of this plugin's development, not just the current
branch's task batch.

### Version numbers in change documentation

Never assert the version/tag a not-yet-released change will ship under
(e.g. "new in 1.2.73-205") unless that number is already decided and fixed
— branch-local "next version" labels are placeholders and the actual
release tag (e.g. it could become 1.4.0 instead) is often undecided until
release time. When documenting a still-unreleased addition, either omit the
version or say it's "on this branch, not yet released under a tag". When
documenting a removed/renamed/deprecated item, anchor on the last tag that
is certain to have shipped it (e.g. "present up to and including
1.2.73-204") rather than the uncertain next tag.