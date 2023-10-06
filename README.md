[![DOI](https://zenodo.org/badge/701301428.svg)](https://zenodo.org/badge/latestdoi/701301428)

# Arbitrary Transport

 This is a [julia](https://julialang.org) code-base for simulating trajectories of (para/dia)-magnetic colloids above magnetic patterns.
 The scripts in the `scripts`-directory produce the trajectories found in the plots.

# Getting Started
## Installing and running julia
Download and install [julia](https://julialang.org/downloads).
Then start a julia REPL from this directory by running `julia` in your terminal.
## Installing DrWatson
DrWatson will handle all the dependencies for this project.
Install the [DrWatson](https://juliadynamics.github.io/DrWatson.jl/dev/) package by typing `] add DrWatson` in the julia REPL.
Hit `backspace` after the installation finished to get back to a julia REPL.
Then type `using DrWatson; @quickactivate`.
In order to install the packages required for this simulation software type `] instantiate`.
Hit `backspace` once more after the installation is done to get back to the julia REPL.
Now everything should be installed and you should be ready to run some simulations.
# Simulating
See the scripts in the `scripts`-directory.
There are scripts for simultaneously drawing some well known symbols on rotated square patterns with `ps-symbols.jl`, drawing the letters `A` through `R` above rotated square patterns with `alphabet.jl`, drawing the letters `A` through `D` above rotated square patterns with `abcd.jl` and for drawing the letter `B` on a specially crafted pattern in `b.jl`. The simple spiraling motion is achieved with `spiral.jl`.

In the julia REPL include the file you want by typing `include("scrips/the-file-you-want.jl")`.
Please start a new julia REPL for each simulation script.
You do this by typing `exit()` after you are done with one simulation.
Then start a new julia REPL by typing `julia` in a terminal.
They are not designed to work in the same REPL.
The simulation scripts include some basic plotting functionality.

# Speeding things up
Many of the simulations are designed to run in parallel.
In order to allow julia to run multiple threads, start julia with the command-line argument `--threads n` with `n` being the number of threads you want to use.
When using multiple threads the progress bars will not be displayed properly.

# Possible Pitfalls
We use the package `PyPlot` for plotting.
Julia does not install the required `python` packages for `PyPlot` itself.
You need to do this manually.
Julia will tell you what to do.
