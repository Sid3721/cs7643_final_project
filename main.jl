using Pkg
Pkg.add("PowerModels")
Pkg.add("PGLib")
Pkg.add("Distributions")
Pkg.add("Ipopt")
Pkg.add("TickTock")
Pkg.add("JuMP")
Pkg.add("Random")

using PGLib
using PowerModels
using Distributions
using Ipopt
using Random
using TickTock

include("AC_OPF_Project.jl")
include("utils.jl")
include("generate_data.jl")

(data_str, lowerBound, upperBound, stDev, samples, newfile, x_filename, y_filename) = ARGS

lowerBound = parse(Float16, lowerBound)
upperBound = parse(Float16, upperBound)
stDev   = parse(Float16, stDev)
samples = parse(Int32, samples)

if newfile == "true"
    newfile = true
else
    newfile = false
end

# newfile = 

# network = pglib(data_str)

# generate_data(data_str, network, lowerBound, upperBound, stDev, samples, x_filename, y_filename)
# generate_data(data_str, lowerBound, upperBound, stDev, samples, x_filename, y_filename)
generate_data(data_str, lowerBound, upperBound, stDev, samples, newfile, x_filename, y_filename)

## Possible inputs
# pglib_opf_case793_goc.m 0.8 1.2 0.05 1000 true X_data2.csv y_data2.csv
# pglib_opf_case793_goc.m 0.8 1.2 0.05 10 true X_data2.csv y_data2.csv

# data_str = "pglib_opf_case793_goc.m"
# # network = network_data
# lowerBound = 0.8
# upperBound = 1.2
# stDev = 0.05
# samples = 10
# newfile = true
# x_filename = "X_data2.csv"
# y_filename = "y_data2.csv"
# # generate_networks(network, lowerBound, upperBound, stDev, iterations)

# # d = pglib(data_str)
# # network = data_prep(d)


# generate_data(data_str, lowerBound, upperBound, stDev, samples, newfile, x_filename, y_filename)