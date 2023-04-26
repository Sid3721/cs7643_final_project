# Packages
using JuMP
using PowerModels
using PGLib
using Ipopt
# using NBInclude

function format_results(model, data)

    solutions = Dict("branch"              => Dict()
                    ,"gen"                 => Dict()
    #                 ,"multiinfrastructure" => false
    #                 ,"multinetwork"        => false
                    ,"bus"                 => Dict()
    #                 ,"per_unit"            => true
    )

    for gen in data[:gen]
        i = gen[1]
        pg = value(model[:pg][i])
        qg = value(model[:qg][i])
        solutions["gen"][string(i)] = Dict("pg" => pg, "qg" => qg)
    end

    for bus in data[:bus]
        i = bus[1]
        va = value(model[:va][i])
        vm = value(model[:vm][i])
        solutions["bus"][string(i)] = Dict("va" => va, "vm" => vm)
    end

    for branch in data[:arcs_from]
        l, i, j = branch
        pf = value(model[:pf][(l, i, j)])
        qf = value(model[:qf][(l, i, j)])
        pt = value(model[:pf][(l, j, i)])
        qt = value(model[:qf][(l, j, i)])
        solutions["branch"][string(l)] = Dict("pf" => pf, "qf" => qf, "pt" => pt, "qt" => qt)    
    end

    results = Dict("solve_time" => solve_time(model)
                    ,"optimizer" => "Ipopt"
                    ,"termination_status" => termination_status(model)
                    ,"dual_status" => dual_status(model)
                    ,"primal_status" => primal_status(model)
                    ,"objective" => objective_value(model)
                    ,"solution" => solutions
                )

    return results
end


#############
# Model
#############
function ac_opf_centralized(data, warm_start=false, print_level = 0)
    model = Model(Ipopt.Optimizer)

    set_optimizer_attribute(model, "print_level", print_level)

    G = data[:gen]
    N = data[:bus]
    E_union_Er = data[:arcs]
    
    
    # Power Generated
    @variable(model, G[i]["pmin"] <= pg[i in keys(G)] <= G[i]["pmax"])
    @variable(model, G[i]["qmin"] <= qg[i in keys(G)] <= G[i]["qmax"])

    # Bus Voltage Magnitude and Angles
    @variable(model, N[i]["vmin"] <= vm[i in keys(N)] <= N[i]["vmax"])
    @variable(model, va[i in keys(N)])

    # Line Power Flow
    @variable(model, -data[:branch][l]["rate_a"] <= pf[(l, i, j) in E_union_Er] <= data[:branch][l]["rate_a"])
    @variable(model, -data[:branch][l]["rate_a"] <= qf[(l, i, j) in E_union_Er] <= data[:branch][l]["rate_a"])

    @objective(model, Min, 
        sum(G[a]["cost"][1] * pg[a]^2 
            + G[a]["cost"][2] * pg[a] 
            + G[a]["cost"][3] for a in keys(G))
    )

    if warm_start != false
        gens  = warm_start["gen"]
        buses = warm_start["bus"]
        for pair in gens
            gen_idx, gen_data = pair

            g_id = parse(Int64, gen_idx)

            pg_warm = gen_data["pg"]
            qg_warm = gen_data["qg"]
            
            pg_max = G[g_id]["pmax"]
            pg_min = G[g_id]["pmin"]
            qg_max = G[g_id]["qmax"]
            qg_min = G[g_id]["qmin"]

            # Assign value for pg and qg start values respectively

            # PG
            if pg_warm > pg_max
                pg_start = pg_max
            elseif pg_warm < pg_min
                pg_start = pg_min
            else
                pg_start = pg_warm
            end

            # QG
            if qg_warm > qg_max
                qg_start = qg_max
            elseif pg_warm < pg_min
                qg_start = qg_min
            else
                qg_start = qg_warm
            end

            JuMP.set_start_value(pg[g_id], pg_warm)
            JuMP.set_start_value(qg[g_id], qg_warm)

        end

        # warm start for vm (may remove)
        for pair in buses
            bus_idx, bus_data = pair

            b_id = parse(Int64, bus_idx)

            # va_warm = bus_data["va"]
            vm_warm = bus_data["vm"]

            vm_max = N[b_id]["vmax"]
            vm_min = N[b_id]["vmin"]

            if vm_warm > vm_max
                vm_start = vm_max
            elseif vm_warm < vm_min
                vm_start = vm_min
            else
                vm_start = vm_warm
            end

            JuMP.set_start_value(vm[b_id], vm_start)
        end

        # for triplet in E_union_Er
        #     br_id, i_bus, j_bus = triplet


        #     # TODO: cry...

        # end

    end

    for ref in keys(data[:ref_buses])
        @constraint(model, va[ref] == 0.0)
    end

    for (l, branch) in data[:branch]
        i = branch["f_bus"]
        j = branch["t_bus"]

        g_fr = branch["g_fr"]
        g_to = branch["g_to"]
        b_fr = branch["b_fr"]
        b_to = branch["b_to"]
        tm = branch["tap"]^2

        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)

        # Line Power Flow Constraints
        @constraint(model, pf[(l, i, j)]^2 + qf[(l, i, j)]^2 <= branch["rate_a"]^2) # okay
        @constraint(model, pf[(l, j, i)]^2 + qf[(l, j, i)]^2 <= branch["rate_a"]^2) # okay

        # Angle Differences
        @constraint(model, branch["angmin"] <= va[i] - va[j] <= branch["angmax"]) # okay

        # Ohm's Law, real, forward
        @NLconstraint(model, pf[(l, i, j)] 
                        == (g + g_fr)/tm * vm[i]^2 
                            + (-g * tr + b * ti)/tm * (vm[i]*vm[j]*cos(va[i]-va[j])) 
                            + (-b * tr - g * ti)/tm * (vm[i]*vm[j]*sin(va[i]-va[j]))
        )
        # Ohm's Law real, backward
        @NLconstraint(model, pf[(l, j, i)] 
                        == (g + g_to) * vm[j]^2 
                            + (-g * tr - b * ti)/tm * (vm[i]*vm[j]*cos(va[j]-va[i])) 
                            + (-b * tr + g * ti)/tm * (vm[i]*vm[j]*sin(va[j]-va[i]))
        )
        # Ohm's Law Imaginary, forward
        @NLconstraint(model, qf[(l, i, j)] 
                        == -(b + b_fr)/tm * vm[i]^2 
                            - (-b * tr - g * ti)/tm * (vm[i]*vm[j]*cos(va[i]-va[j])) 
                            + (-g * tr + b * ti)/tm * (vm[i]*vm[j]*sin(va[i]-va[j]))
        )
        # Ohm's Law Imaginary, backward
        @NLconstraint(model, qf[(l, j, i)] 
                        == -(b + b_to) * vm[j]^2 
                            - (-b * tr + g * ti)/tm * (vm[i]*vm[j]*cos(va[j]-va[i])) 
                            + (-g * tr - b * ti)/tm * (vm[i]*vm[j]*sin(va[j]-va[i]))
        )
    end

    for (i, bus) in data[:bus]
        # Real Power, Kirchoff's Law
        p_gen   = sum(pg[a] for a in data[:bus_gens][i]; init=0)
        p_load  = sum(data[:load][a]["pd"] for a in data[:bus_loads][i]; init=0)
        p_shunt = sum(data[:shunt][a]["gs"] for a in data[:bus_shunts][i]; init=0) * vm[i]^2
        p_line  = sum(pf[a] for a in data[:bus_arcs][i]; init=0)

        @constraint(model, p_gen - p_load - p_shunt == p_line)

        # Imaginary Power, Kirchoff's Law
        q_gen   = sum(qg[a] for a in data[:bus_gens][i]; init=0)
        q_load  = sum(data[:load][a]["qd"] for a in data[:bus_loads][i]; init=0)
        q_shunt = sum(data[:shunt][a]["bs"] for a in data[:bus_shunts][i]; init=0) * vm[i]^2
        q_line  = sum(qf[a] for a in data[:bus_arcs][i]; init=0)
        @constraint(model, q_gen - q_load + q_shunt == q_line)
    end 
    return model
end

