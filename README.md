# CoRaLS

The `CoRaLS.jl` (Cosmic Ray Lunar Sounder) Monte Carlo model computes detection rates of Askaryan emission from cosmic ray interactions in lunar regolith.

## Quick start

CoRaLS currently requires Julia 1.11.x: its `Logging` compatibility bound does
not resolve with Julia 1.10. You also need a Python executable that can import
`matplotlib`; CoRaLS uses it through `PyCall`/`PyPlot`.

```sh
git clone git@github.com:afromero/CoRaLS.jl.git
cd CoRaLS.jl

# Use the Julia 1.11 executable on your system.
PYTHON=/path/to/python \
  /path/to/julia-1.11/bin/julia --project=. scripts/setup.jl

/path/to/julia-1.11/bin/julia --project=. scripts/smoke_test.jl
```

For example, check Python and matplotlib before setup with:

```sh
/path/to/python -c 'import sys, matplotlib; print(sys.executable, matplotlib.__version__)'
```

`scripts/setup.jl` is safe to rerun. It resolves the project packages, builds
`PyCall` against `PYTHON`, and precompiles CoRaLS. The smoke test executes a
small seeded acceptance calculation; a zero acceptance is expected at this
small sample size and is not an error.

To use a Julia REPL after setup:

```sh
/path/to/julia-1.11/bin/julia --project=. -t auto
```

```julia
using CoRaLS
using CoRaLS: km

A = acceptance(10_000, 20;
    region=create_region("psr:south"),
    spacecraft=CircularOrbit(50.0km))
plot_acceptance(A)
```

### Jupyter notebooks

Install IJulia once with the same Julia 1.11 installation, then register a
kernel that activates this project. Select **Julia (CoRaLS) 1.11** in Jupyter.

```sh
/path/to/julia-1.11/bin/julia -e 'using Pkg; Pkg.add("IJulia")'
/path/to/julia-1.11/bin/julia -e \
  'using IJulia; installkernel("Julia (CoRaLS)", "--project=/absolute/path/to/CoRaLS.jl")'
```

[`notebooks/generic_acceptance_diagnostic.ipynb`](notebooks/generic_acceptance_diagnostic.ipynb)
is a small interactive acceptance diagnostic based on `slurm/generic_acceptance.jl`.

### Validation and troubleshooting

Run the package-aware test suite with:

```sh
/path/to/julia-1.11/bin/julia --project=. -e 'using Pkg; Pkg.test()'
```

At the current revision, `test/ut_acceptance.jl` compares a stochastic
acceptance calculation to the checked-in `test/ut_acceptance_test.txt` snapshot
and reports that the reference is stale. This is a model-reference maintenance
issue, not an installation failure; use the smoke test above to validate a new
installation until that snapshot is regenerated.

On Linux, a Conda Python with OpenSSL 3.0 can conflict with Julia plotting
artifacts if a process imports both `PyPlot` and Julia's `Plots`. Upgrade that
Conda environment's OpenSSL (3.3 or newer), or use a separate Python
environment with `matplotlib`. CoRaLS itself uses `PyPlot`.

## Calculating rates with CoRaLS

See the REPL example in the quick-start section. Full documentation is coming
soon.

## Developers

To make and preview docs locally run the following from the root CoRaLS directory:

```bash
julia --project=docs -e 'include("docs/make.jl"); using LiveServer; serve(dir="docs/build")'
```

To run tests:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Citing

Paper coming soon.
