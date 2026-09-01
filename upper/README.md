# Certificate for `E_star(0) <= 1.944`

Run `python3 verify_upper_1p944.py estar_upper_1p944_certificate.json`. The verifier uses exact rational arithmetic and integer square root. It checks the spline core, both analytic tails, interface signs, and the rational lower sum for the barrier integral. Successful execution ends with `VERIFIED: E_star <= 243/125`. See Appendix E.1 of the manuscript for continuum soundness.
