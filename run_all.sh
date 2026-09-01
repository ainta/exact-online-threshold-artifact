#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

if command -v sha256sum >/dev/null 2>&1; then
    sha256sum -c SHA256SUMS
elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 -c SHA256SUMS
else
    echo "ERROR: sha256sum or shasum is required" >&2
    exit 1
fi

python3 verify_release_consistency.py
python3 upper/verify_upper_1p944.py upper/estar_upper_1p944_certificate.json
lower/run_verifier.sh
