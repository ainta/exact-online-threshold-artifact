#!/usr/bin/env python3
"""Check the endpoints encoded by the two frozen certificates."""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent
upper = json.loads((ROOT / "upper/estar_upper_1p944_certificate.json").read_text())
lower = json.loads((ROOT / "lower/certificate.json").read_text())

if Fraction(upper["caps"]["B"]) != Fraction(243, 125):
    raise RuntimeError("upper value mismatch")
statement = lower["statement"]["E_star_strict_lower_bound"]
if Fraction(statement["numerator"], statement["denominator"]) != Fraction(173639, 100000):
    raise RuntimeError("lower value mismatch")

print("PASS: frozen upper and lower certificate endpoints are internally consistent")
