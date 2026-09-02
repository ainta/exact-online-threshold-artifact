# Computer-assisted certificates

This repository is the artifact described in Appendix E.3 of *The Exact
Online Threshold for the Asymmetric Binary Perceptron*. It contains frozen
certificate data and separate verifiers for

```text
E_star(0) <= 1.944
E_star(0) > 1.73639
```

The mathematical implications of verifier acceptance are proved in Appendix
E.1 (upper bound) and Appendix E.2 (lower bound) of the manuscript. Search and
optimization programs are not part of the trusted computation.

The manuscript PDF is included as [`submission.pdf`](submission.pdf).

## Reproduce

Requirements: Python 3, a C++17 compiler with OpenMP, GMP, and MPFR 4.x
(`libmpfr` ABI 6). On Debian or Ubuntu, install them with:

```bash
sudo apt-get install python3 g++ libgmp-dev libmpfr-dev
```

Run both checks from the repository root:

```bash
./run_all.sh
```

The upper check uses exact rational arithmetic and normally finishes in a few
seconds. The lower check uses directed MPFR rounding and exact NTT/CRT
convolution and may take several minutes. Successful execution ends with:

```text
VERIFIED: E_star <= 243/125
PASS: E_star > 173639/100000
```

`SHA256SUMS` freezes the released files and is checked by `run_all.sh` before
either verifier runs. Reference output from each verifier is included for
comparison; timing lines are machine-dependent.

## License

The verification code, certificates, and supporting artifact files in this
repository are released under the [MIT License](LICENSE). The manuscript PDF
`submission.pdf` is not covered by that license; its reuse is governed by the
license attached to the corresponding arXiv version.
