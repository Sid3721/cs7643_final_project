using PGLib
using PowerModels
using Distributions
using Ipopt
using Random
using TickTock

include("AC_OPF_Project.jl")
include("utils.jl")


function generate_networks(network, lowerBound, upperBound, stDev, x_filename, y_filename, iterations = 5)

    d = LogNormal(-stDev^2/2, stDev)
    u = Uniform(lowerBound, upperBound)

    network = data_prep(network)

    for iter in 1:iterations

        if iter % 10 == 0
            println("Iteration " * string(iter))
        end

        # determines what factor we multiply the whole system by
        scalar = rand(u)

        new_network = deepcopy(network)
        for (b_id, b_data) in new_network[:load]

            rand1 = rand(d)
            rand2 = rand(d)
            
            # if b_id == 100
            #     println(new_network[:load][b_id]["pd"])
            #     println(b_data["pd"] * scalar * rand1)
            # end

            new_network[:load][b_id]["pd"] = b_data["pd"] * scalar * rand1
            new_network[:load][b_id]["qd"] = b_data["qd"] * scalar * rand2


        end

        model = ac_opf_centralized(new_network)
        optimize!(model)

        if string(termination_status(model)) == "LOCALLY_SOLVED"

            # print("I terminated")
            result = format_results(model, new_network)["solution"]

            # append_output_file(network_data, result) 
            append_output_file(new_network, result, x_filename, y_filename)
        else
            println("Iteration " * string(iter) * " didn't converge. (Universal Mulitplier: ")
            print(string(scalar) * ")")
        end


    end

end



# function generate_data(data_str, network, lowerBound, upperBound, stDev, iterations, newfile = true, x_filename = "X_data.csv", y_filename = "y_data.csv")
function generate_data(data_str, lowerBound, upperBound, stDev, iterations, newfile = true, x_filename = "X_data.csv", y_filename = "y_data.csv")

    tick()

    ## Read data
    data = pglib(data_str)
    network_data = data_prep(data)

    if newfile
        Random.seed!(1)
        ## Create output file
        create_output_file(network_data, x_filename, y_filename)
    end
    
    ## Solve model (ie collect other results)
    model = ac_opf_centralized(network_data)
    optimize!(model)

    # 
    generate_networks(data, lowerBound, upperBound, stDev, x_filename, y_filename, iterations)

    println("Time Elapsed:")
    tock()
end


# ## Solve model (ie collect other results)
# model = ac_opf_centralized(network_data)
# optimize!(model)

## Params
# data_str = "pglib_opf_case793_goc.m"
# lowerBound = 0.8
# upperBound = 1.2
# stDev = 0.05
# samples = 10
# newfile = true
# x_filename = "X_data2.csv"
# y_filename = "y_data2.csv"
# generate_networks(network, lowerBound, upperBound, stDev, iterations)

# d = pglib(data_str)
# network = data_prep(d)


# generate_data(data_str, lowerBound, upperBound, stDev, samples, newfile, x_filename, y_filename)

