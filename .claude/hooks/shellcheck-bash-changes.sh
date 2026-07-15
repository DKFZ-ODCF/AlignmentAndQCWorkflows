#!/usr/bin/env bash
set -ueo pipefail

input="$(cat)"
filePath="$(jq -r '.tool_input.file_path // .tool_response.filePath // empty' <<<"$input")"

[[ "$filePath" == *.sh ]] || exit 0
[[ -f "$filePath" ]] || exit 0

if output="$(shellcheck "$filePath" 2>&1)"; then
    exit 0
fi

jq -n --arg ctx "shellcheck found issues in $filePath:
$output" '{hookSpecificOutput: {hookEventName: "PostToolUse", additionalContext: $ctx}}'