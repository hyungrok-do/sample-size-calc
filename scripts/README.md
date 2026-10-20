# Scripts

Tools for working on the source code. None of them is needed to build,
install, test or run the calculator (`setup_venv.py` only makes doing so
easier), and none ships in the wheels or the web page.

| Script | Purpose |
|---|---|
| `setup_venv.py` | Makes the repository's virtual environment (`.venv`) and installs the calculator into it; run it again after a `git pull` to bring it up to date. |
| `gen_temme_coefficients.py` | Generates `core/src/stats/temme_coefficients.inc`, the coefficients of Temme's expansion of the incomplete gamma function, in exact rational arithmetic. |
| `stress_inputs.py` | Types a battery of texts (blank, letters, decimal commas, 1e-300 and 1e300, inf, bad ranges, 3,000 values ...) into every input of every method in every solve mode, as an interface runs them, and reports runs that hang (and whether Stop ends them), failures without a clear message, answers that are no number, and infinite values. |
| `stress_simulations.py` | Stress-tests every simulation procedure: fuzzed and extreme inputs, thread counts, concurrent runs, cancellation, type I error and a million simulated studies (the long version of `tests/test_stress.py`). |
| `update_readme_procedures.py` | Rewrites the list of procedures in `README.md` from the engine's catalog. |

The reference values that the tests compare with (`core/tests/fixtures/`) are
stored as data. The scripts that computed them, in R and with mpmath, are not
part of the source distribution; they remain in the project history (commit
`5e2b407`, folder `scripts/reference/`) for anyone who wants to extend them.
