#!/usr/bin/env bash
set -ueo pipefail

input="$(cat)"
filePath="$(jq -r '.tool_input.file_path // .tool_response.filePath // empty' <<<"$input")"

[[ "$filePath" =~ \.[Rr]$ ]] || exit 0
[[ -f "$filePath" ]] || exit 0
command -v Rscript >/dev/null 2>&1 || exit 0

output="$(Rscript -e 'lints <- lintr::lint(commandArgs(trailingOnly = TRUE)[1]); if (length(lints) > 0) print(lints)' "$filePath" 2>&1)"

[[ -z "$output" ]] && exit 0

jq -n --arg ctx "lintr found issues in $filePath:
$output" '{hookSpecificOutput: {hookEventName: "PostToolUse", additionalContext: $ctx}}'