#!/bin/bash
set -e

SPATULA_BIN="/app/bin/spatula"

if [[ ! -f "$SPATULA_BIN" ]]; then
  echo "Error: spatula binary was not found at $SPATULA_BIN"
  exit 1
fi

# If no arguments are passed, show help
if [[ $# -eq 0 ]]; then
  echo " No arguments passed. Showing spatula help:"
  exec "$SPATULA_BIN" --help
else
  # Forward all arguments to spatula
  exec "$SPATULA_BIN" "$@"
fi
