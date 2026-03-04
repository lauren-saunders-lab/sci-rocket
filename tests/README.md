# Tests

This folder contains:

- `unit/`: fast Python unit tests for workflow scripts.
- `integration/`: end-to-end workflow tests (Snakemake outputs + pytest checks + Playwright UI checks).

These tests all run on every pull request and commit to the main & dev branches as part of the continuous integration setup on github.

Below are instructions for how to run them locally.

## Prerequisites

You will need a Linux shell environment with:

- [Conda](https://conda-forge.org/download/)
- [pnpm](https://pnpm.io/installation#on-posix-systems)

and an active conda environment where you can run sci-rocket, e.g.

```bash
conda create -n sci-rocket-tests -c conda-forge -c bioconda python=3.13 snakemake pytest pandas
conda activate sci-rocket-tests
```

Alternatively you can use an existing conda environment you use to run sci-rocket, just ensure pytest is installed:

```bash
conda install pytest
```

## Run unit tests locally

From repo root:

```bash
pytest tests/unit -v
```

# Run integration tests locally

From the tests/integration folder:

```bash
./run.sh
```

`tests/integration/run.sh` will:

1. Install JS test dependencies with `pnpm install --frozen-lockfile`.
2. Ensure Playwright Chromium is installed.
3. Download integration test data (if `tests/integration/data` is missing).
4. For each integration suite (`bcl`, `aviti`, `velocity`):
   - run the sci-rocket workflow
   - run the Python integration tests in `tests`
   - run Playwright UI tests in `ui_tests`

Steps 1-3 may take some time to run the first time you run the script.