using PowerModels
using MAT

dirname = ARGS[1]
filename = ARGS[2]

network_data = PowerModels.parse_file(filename)

num_bus = length(network_data["bus"])
num_branch = length(network_data["branch"])

# ind_bus = [parse(Int,key) for (key,b) in network_data["bus"]]
# ind_branch = [parse(Int,key) for (key,b) in network_data["branch"]]
# ind_gen = [parse(Int,key) for (key,b) in network_data["gen"]]
# ################################################################################
# slack = [bus["index"] for (i,bus) in network_data["bus"] if bus["bus_type"] == 3][1]
# ############## Mapping generators to buses #####################################
# gen_buses = [network_data["gen"][i]["gen_bus"] for i in keys(network_data["gen"])]
# gen_buses = unique(gen_buses)
# load_buses = [parse(Int64,i) for i in keys(network_data["bus"]) if abs(network_data["bus"][i]["pd"]) + abs(network_data["bus"][i]["qd"]) > 1e-2]
# active_buses = union(gen_buses,load_buses)


# inflation_consts = [0.01]
inflation_consts = [0.05,0.1,0.2,0.3,0.4]
algorithms = [0,1,2]

for inflation in inflation_consts, algo in algorithms

    coeff_const = zeros(num_branch)
    coeff_p = zeros(num_branch,num_bus)
    coeff_q = zeros(num_branch,num_bus)
    approx_error = zeros(num_branch)

    for i in keys(network_data["branch"])
        linenum = parse(Int64,i)
        vars = matread(string(dirname)"/linear_approximations_real"string(inflation)"_line_"string(linenum)"_algorithm_"string(algo)".mat")
        coeff_const[linenum] = vars["coeff_const"]
        coeff_p[linenum,:] = vars["coeff_p"]
        coeff_q[linenum,:] = vars["coeff_q"]
        approx_error[linenum] = vars["approx_error"]
    end

    # write aproximations for real power
    matwrite(string(dirname)"/matrix_forms/linear_approximations_real"string(inflation)"_algorithm_"string(algo)".mat",Dict("coeff_const"=>coeff_const,"coeff_p"=>coeff_p,"coeff_q"=>coeff_q,"approx_error"=>approx_error))


    for i in keys(network_data["branch"])
        linenum = parse(Int64,i)
        vars = matread(string(dirname)"/linear_approximations_reactive"string(inflation)"_line_"string(linenum)"_algorithm_"string(algo)".mat")
        coeff_const[linenum] = vars["coeff_const"]
        coeff_p[linenum,:] = vars["coeff_p"]
        coeff_q[linenum,:] = vars["coeff_q"]
        approx_error[linenum] = vars["approx_error"]
    end

    # write aproximations for real power
    matwrite(string(dirname)"/matrix_forms/linear_approximations_reactive"string(inflation)"_algorithm_"string(algo)".mat",Dict("coeff_const"=>coeff_const,"coeff_p"=>coeff_p,"coeff_q"=>coeff_q,"approx_error"=>approx_error))

end
