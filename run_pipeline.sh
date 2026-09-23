#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

IMAGE="${PKS_DOCK_IMAGE:-pks-dock:0.2.0}"

mkdir -p "$SCRIPT_DIR/results" "$SCRIPT_DIR/logs"

docker run --rm \
    -v "$SCRIPT_DIR:/opt/PKS_DOCK" \
    "$IMAGE" \
    "$@"
