#!/usr/bin/env bash
# Proves the Phase 9 Vina log parsing fix against a synthetic Vina log,
# without needing Vina itself installed. Run this any time you touch
# the parsing logic in scripts/09_run_docking.sh.
set -euo pipefail

TMP_LOG=$(mktemp)
cat > "$TMP_LOG" << 'EOF'
AutoDock Vina v1.2.5
Estimated Free Energy of Binding    :   -7.23 (kcal/mol) [=(1)+(2)+(3)-(4)]

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1       -7.234      0.000      0.000
   2       -6.987      1.234      2.345
   3       -6.500      2.100      3.400
EOF

RESULT=$(awk '/^-----/{found=1; next} found && NF>=2 {print $2; exit}' "$TMP_LOG")
rm -f "$TMP_LOG"

echo "Parsed value: $RESULT"

if [[ "$RESULT" == "-7.234" ]]; then
    echo "PASS: Vina log parsing correctly extracted -7.234"
    exit 0
else
    echo "FAIL: expected -7.234, got '$RESULT'"
    exit 1
fi
