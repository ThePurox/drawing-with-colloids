using DrWatson
@quickactivate

using Test
testdir(args...) = projectdir("tests", args...)

include(testdir("test-loops.jl"))
include(testdir("test-pattern.jl"))
include(testdir("test-forces.jl"))
include(testdir("test-fromimage.jl"))
