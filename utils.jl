using PGLib
using PowerModels
using Distributions
using Ipopt
using Random

include("AC_OPF_Project.jl")

# This function assumes you're passing it in as PGLib
# Reformats data for ease of use
function data_prep(data)
    # Add zeros to turn linear objective functions into quadratic ones
    # so that additional parameter checks are not required
    PowerModels.standardize_cost_terms!(data, order=2)

    # Adds reasonable rate_a values to branches without them
    PowerModels.calc_thermal_limits!(data)

    # use build_ref to filter out inactive components
    data = PowerModels.build_ref(data)[:it][:pm][:nw][0]

    return data
end

# Sorts the results data. Ensures everything is in the correct order
function sort_y_data(results)
    gens = results["gen"]
    out_list = []

    for (g_id, gen) in gens
        qg_res = (g_id * "_" * "qg", gen["qg"])
        pg_res = (g_id * "_" * "pg", gen["pg"])
        push!(out_list, pg_res, qg_res)
    end

    return sort(out_list)
end

function sort_x_data(network_data)
    loads = network_data[:load]
    input_list = []

    for (l_id, load) in loads
        qd_res = (string(l_id) * "_" * "qd", load["qd"])
        pd_res = (string(l_id) * "_" * "pd", load["pd"])
        push!(input_list, pd_res, qd_res)
    end

    return sort(input_list)
end

function rot90(data)
    title_list = []
    value_list = []
    for (title, value) in data
        push!(title_list, title)
        push!(value_list, value)
    end
    return (title_list, value_list)
end

function writecsv(filename, data_list, open_type = "a")
    conn = open(filename, open_type)
    res = ""
    for d in data_list
        res = res * "," * string(d)
    end
    res = res[2:length(res)] * "\n"
    write(conn, res)
    close(conn)
end

function create_output_file(network_data, X_filename = "X_data.csv", y_filename = "y_data.csv")
    model = ac_opf_centralized(network_data)
    optimize!(model)

    result = format_results(model, network_data)["solution"]

    x_data = sort_x_data(network_data)
    y_data = sort_y_data(result)

    X_titles, X_values = rot90(x_data)
    Y_titles, Y_values = rot90(y_data)

    writecsv(X_filename, X_titles, "w")
    writecsv(X_filename, X_values)

    writecsv(y_filename, Y_titles, "w")
    writecsv(y_filename, Y_values)
end

function append_output_file(network_data, result, X_filename = "X_data.csv", y_filename = "y_data.csv")

    x_data = sort_x_data(network_data)
    y_data = sort_y_data(result)

    X_titles, X_values = rot90(x_data)
    Y_titles, Y_values = rot90(y_data)

    writecsv(X_filename, X_values)

    writecsv(y_filename, Y_values)
end

#################
## Example use ##
#################

# ## Read data
# data = pglib("pglib_opf_case793_goc.m")
# network_data = data_prep(data)

# ## Create output file
# create_output_file(network_data)

# ## Solve model (ie collect other results)
# model = ac_opf_centralized(network_data)
# optimize!(model)

# ## 
# result = format_results(model, network_data)["solution"]
# append_output_file(network_data, result)
