#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python3 verify_certificate_match.py
CXX="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:--O3 -DNDEBUG -fno-fast-math -std=c++17 -fopenmp -Wall -Wextra}"
MPFR_LIB="${MPFR_LIB:-}"
if [[ -z "$MPFR_LIB" ]] && command -v ldconfig >/dev/null 2>&1; then MPFR_LIB="$(ldconfig -p 2>/dev/null | awk '/libmpfr\.so(\.6)? /{print $NF; exit}')"; fi
if [[ -n "$MPFR_LIB" ]]; then "$CXX" $CXXFLAGS e_star_lower_bound_verifier.cpp -I. "$MPFR_LIB" -lgmp -o e_star_lower_bound_verifier; else "$CXX" $CXXFLAGS e_star_lower_bound_verifier.cpp -I. -lmpfr -lgmp -o e_star_lower_bound_verifier; fi
OMP_NUM_THREADS="${OMP_NUM_THREADS:-5}" ./e_star_lower_bound_verifier | tee verifier_output_reproduced.log
