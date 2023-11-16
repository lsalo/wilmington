# General
using Revise
using DelimitedFiles
using StatsBase
using Base.Iterators

# Plots
using GLMakie
GLMakie.activate!()
using Colors

# Mine
lib_dir = "/Users/lluis/Documents/julia/jl-tools/"
includet(joinpath(lib_dir, "modules", "polyconnectivity.jl"))
using .PolyConnectivity