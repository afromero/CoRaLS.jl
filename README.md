# CoRaLS

The `CoRaLS.jl` (Cosmic Ray Lunar Sounder) Monte Carlo model computes detection rates of Askaryan emission from cosmic ray interactions in lunar regolith.

## Installation

1. [Install Julia](https://julialang.org/install/) 
2. Install Python and the `matplotlib` package, the easiest way is with [Anaconda](https://www.anaconda.com/) (e.g. `conda list matplotlib`).
3. Fork or clone this repository:

```sh
git clone git@github.com:cjtu/CoRaLS.jl.git
```

4. (First time only): Setup path to python/matplotlib. In the Julia REPL, supply the path to your python environment from #2 in quotes to `ENV["PYTHON"]=""` (leave blank to use system default python). Then add and build the `PyCall` package:

```bash
$ julia

julia> ENV["PYTHON"]=""
julia> using Pkg
julia> Pkg.add("PyCall")
julia> Pkg.build("PyCall")
```

5. Run the test suite to test the installation. From the CoRaLS.jl directory, run:

```bash
julia --project=. test/runtests.jl
```

If the tests were successful, `CoRaLS.jl` is compiled and ready to use!

## Calculating rates with CoRaLS

Start julia in the CoRaLS.jl project (if you are in the directory use `--project=.`). For multithreaded mode, use `-t` to specify the number of threads (default is 1, "auto" chooses for you).

```bash
julia --project=. -t "auto"
julia> using CoRaLS
```

Compute and plot an acceptance:

```julia
julia> A = acceptance(10000, 20; region=create_region("psr:south"), spacecraft=CircularOrbit(50.0km))
julia> plot_acceptance(A)
```

See full documentation online at... (coming soon)

## Developers

To make and preview docs locally run the following from the root CoRaLS directory:

```bash
julia --project=docs -e 'include("docs/make.jl"); using LiveServer; serve(dir="docs/build")'
```

To run tests:

```bash
julia --project=. test/runtests.jl
```

Running `test/runtests.jl` also runs visual diff cleanup for plot PNGs. Files in `test/plots/test_plots` that do not have visible changes are automatically restored, and changed files are exported to a tmp comparison folder for review.

To re-run only the visual diff step manually:

```bash
julia --project=. test/plots/visual_diff_test_plots.jl
```

## Citing

Paper coming soon.
