#!/bin/bash
# =============================================================================
# Transition from coarse mesh (level 4, ~507k cells) to fine mesh (level 6,
# ~2.8M cells) using mapFields to interpolate the solution.
#
# Run this script from the FINE MESH case directory.
# The coarse case (source) path is passed as the first argument.
#
# Usage:
#   bash transition_to_fine.sh <path_to_coarse_case>
#
# Example:
#   cd /scratch/job_FINE/run
#   bash /scratch/job1950709/mapFields/transition_to_fine.sh \
#        /scratch/job1950709/run
# =============================================================================

set -e

COARSE_CASE=${1:?"Usage: $0 <path_to_coarse_case>"}
FINE_CASE=$(pwd)

echo "Source (coarse): $COARSE_CASE"
echo "Target (fine):   $FINE_CASE"
echo ""

# --- Step 1: Reconstruct the latest time on the coarse case -----------------
echo "=== Step 1: Reconstructing coarse case ==="
cd "$COARSE_CASE"
reconstructPar -latestTime > log.reconstructPar_mapFields 2>&1
LATEST=$(foamListTimes -latestTime 2>/dev/null | tail -1)
echo "Latest coarse time: $LATEST"
cd "$FINE_CASE"

# --- Step 2: Copy mapFieldsDict into fine case system/ ----------------------
echo "=== Step 2: Installing mapFieldsDict ==="
cp "$(dirname "$0")/mapFieldsDict" system/mapFieldsDict

# --- Step 3: Run mapFields --------------------------------------------------
echo "=== Step 3: Running mapFields ==="
mapFields "$COARSE_CASE" -sourceTime latestTime -consistent > log.mapFields 2>&1
echo "mapFields done. Check log.mapFields for details."

# --- Step 4: Decompose the mapped time directory ----------------------------
echo "=== Step 4: Decomposing mapped fields ==="
decomposePar -latestTime > log.decomposePar_mapFields 2>&1
echo "decomposePar done."

echo ""
echo "=== Transition complete ==="
echo "Update system/controlDict: set startFrom latestTime"
echo "Then submit the fine-mesh solver job."
