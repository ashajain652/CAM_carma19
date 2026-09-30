#!/bin/bash
# Build the offline CARMA optics driver from the CARMA source in this CAM tree.
#
# Uses CAM's own carma_precision_mod / carma_constants_mod (cam/), all of CARMA's
# base/ code (except the base/ versions of those two modules, which CAM replaces,
# and bhmie.F90, which is replaced by bhmie_patched.F90 from this directory),
# CESM's shr_kind_mod / shr_const_mod, and cam_stubs.F90 for the few CAM modules
# those need. Compiler flags follow the CESM gnu build (-O, free line length).
#
# Usage: ./build.sh            (builds ./carma_optics_driver)
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
# This directory is components/cam/tools/carma_optics_offline in a CESM checkout.
CAM=$(cd "$HERE/../.." && pwd)
CESM=${CESM:-$(cd "$CAM/../.." && pwd)}
CARMA=$CAM/src/physics/carma
SHARE=$CESM/cime/src/share/util
FC=${FC:-gfortran}
FFLAGS="-O -ffree-line-length-none -fallow-argument-mismatch"

OBJ=$HERE/obj
rm -rf "$OBJ"; mkdir -p "$OBJ"
cd "$OBJ"

srcs=("$SHARE/shr_kind_mod.F90" "$SHARE/shr_const_mod.F90" "$HERE/cam_stubs.F90"
      "$CARMA/cam/carma_precision_mod.F90" "$CARMA/cam/carma_constants_mod.F90"
      "$HERE/bhmie_patched.F90")
for f in "$CARMA"/base/*.F90; do
  case $(basename "$f") in
    carma_precision_mod.F90|carma_constants_mod.F90|bhmie.F90) ;;
    *) srcs+=("$f") ;;
  esac
done

# Compile in passes until every file builds (resolves module dependencies).
pending=("${srcs[@]}")
for pass in $(seq 1 30); do
  left=()
  for f in "${pending[@]}"; do
    if ! $FC $FFLAGS -c "$f" -o "$(basename "${f%.F90}").o" > "$(basename "$f").log" 2>&1; then
      left+=("$f")
    fi
  done
  [ ${#left[@]} -eq 0 ] && break
  if [ ${#left[@]} -eq ${#pending[@]} ]; then
    echo "Cannot compile:"; printf '  %s\n' "${left[@]}"
    cat "$(basename "${left[0]}").log"; exit 1
  fi
  pending=("${left[@]}")
done
echo "compiled ${#srcs[@]} CARMA/support files in $pass passes"

$FC $FFLAGS -c "$HERE/carma_optics_driver.F90" -o carma_optics_driver.o
$FC -o "$HERE/carma_optics_driver" *.o
echo "built $HERE/carma_optics_driver"
