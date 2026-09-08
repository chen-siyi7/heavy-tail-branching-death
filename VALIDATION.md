# Validation record

The packaged workflow was tested on 8 September 2026 with R 4.4.1 on macOS. The full run completed successfully in approximately 4.5 seconds on the test machine. Runtime will vary with the machine and R installation.

## Source identity

The supplied `simulation_audit.R` is byte-for-byte identical to the script used for the revised manuscript. Its SHA256 checksum is:

```text
7e77c6d2ac6e29c41718d2fa00600a2f01a05dc05fd1752ad13361614e68affd
```

The new runner does not change the simulation code. It sets the reference random-number generator kinds, runs the original script in a separate output directory, checks output counts, and records the R session and checksums. The master seed remains 20260902.

## Full reproduction check

The runner was invoked by absolute path from outside the repository, with a fresh temporary output directory. It produced six figure PDFs, 15 numerical CSV files, eight LaTeX table fragments, and the audit report. All 15 numerical CSV files and all eight LaTeX table fragments matched the existing manuscript audit outputs byte for byte. PDF byte identity was not used as a test because the export includes creation metadata.

The recorded maximum absolute standardized discrepancies were 2.826 for eventual survival, 2.622 for embedded-walk threshold hitting, and 2.249 for bounded-support means. The maximum threshold approximation error was 3.018e-5. These are Monte Carlo diagnostics on the stated grid, not proofs or guarantees for other parameter settings.

The finite-count tail check had a maximum relative discrepancy of 0.04218 at a count of 10,000. This compares a finite-count result with an asymptotic limit. The leading small-death approximation had relative error 0.9927 at death rate 0.5, documenting its limited range of accuracy.

## Entry-point checks

Both R files parsed successfully. The help option worked. Unsupported arguments, including `--publish`, were rejected by the runner. A second invocation targeting the existing output directory stopped without overwriting it. The run wrote session information, generator settings, elapsed time, and output checksums.

The citation file parsed as YAML and includes the core citation fields. No repository URL, publication DOI, release DOI, or license was fabricated. Generated outputs, reviewer correspondence, manuscript drafts, and legacy simulation files are excluded from the distributable archive.

## Scope

This record covers the supplied code and fixed manuscript settings on the tested R installation. It does not establish cross-platform byte identity, validate arbitrary parameter changes, or certify journal acceptance. The figure-font embedding limitation described in the README remains unchanged. No manuscript files were modified or published during this packaging step.
