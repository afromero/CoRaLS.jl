"""
Set up CoRaLS's Julia environment and build PyCall against a selected Python.

Run from the repository root, for example:

    PYTHON=/path/to/python julia --project=. scripts/setup.jl

The selected Python must be able to import matplotlib.
"""

using Pkg

python = get(ENV, "PYTHON", "")
isempty(python) && error("Set PYTHON to the Python executable that has matplotlib installed.")
isfile(python) || error("PYTHON does not point to a file: $python")

ENV["PYTHON"] = abspath(python)
Pkg.instantiate()
Pkg.build("PyCall")
Pkg.precompile()

println("CoRaLS setup complete.")
println("PyCall is configured to use: $(ENV["PYTHON"])")
