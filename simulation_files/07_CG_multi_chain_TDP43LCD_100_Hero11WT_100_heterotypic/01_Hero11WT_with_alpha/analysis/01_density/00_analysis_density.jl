#!/usr/bin/env julia

include("/home/ctan/Workspace/genesis_CG_julia/src/lib/conformation.jl")
include("/home/ctan/Workspace/genesis_CG_julia/src/lib/parser_dcd.jl")

using Printf

# EDIT THIS: temperatures
T_list = [i * 10 + 240 for i in 1:10]
n_run = 1

# EDIT THIS: protein names
NAME_PRO1 = "hero11"
NAME_PRO2 = "tdp43"
# EDIT THIS: number of residues in one chain
NUM_ATOM_PRO1 = 99
NUM_ATOM_PRO2 = 154
# EDIT THIS: number of chains
num_pro1 = 100
num_pro2 = 100
num_all = num_pro1 + num_pro2

# EDIT THIS: box size (x, y: nm; z: 100nm)
V = 18 * 18 * 3 * 1.0e-24
NA = 6.02e23

# EDIT THIS: box size z
box_size_z = 3000
box_size_bin = box_size_z / 100


for (iT, T) in enumerate(T_list)
    for i in 1 : n_run
        # EDIT THIS: trajectory file name
        if_name = @sprintf("../../%s_%s_md4_T_%02d.dcd", NAME_PRO1, NAME_PRO2, iT)
        println("Processing ", if_name)

        mytraj = read_dcd(if_name)

        # EDIT THIS: output data file name
        of_pro1_name = @sprintf("%s_md4_T_%02d.density.dat", NAME_PRO1, iT)
        of_pro2_name = @sprintf("%s_md4_T_%02d.density.dat", NAME_PRO2, iT)

        data_pro1_output = open(of_pro1_name, "w")
        data_pro2_output = open(of_pro2_name, "w")

        data_pro1_all = zeros(Float64, (100, length( mytraj.conformations )))
        data_pro2_all = zeros(Float64, (100, length( mytraj.conformations )))
        Threads.@threads for t in 1 : length( mytraj.conformations )
            box_size_z = mytraj.boundary_box_size[6,t]
            mass_pro1_histogram = [0 for j in 1:100]
            mass_pro2_histogram = [0 for j in 1:100]

            for i in 1 : num_pro1 * NUM_ATOM_PRO1
                z = mytraj.conformations[t].coors[3,i]
                z_int = mod( Int(floor(mod(z, box_size_z) / box_size_bin)), 100)
                mass_pro1_histogram[z_int + 1] += 1
            end
            mass_pro1_distribution = mass_pro1_histogram ./ ( NUM_ATOM_PRO1 * V * NA )

            for i in num_pro1 * NUM_ATOM_PRO1 + 1 : mytraj.num_atoms
                z = mytraj.conformations[t].coors[3,i]
                z_int = mod( Int(floor(mod(z, box_size_z) / box_size_bin)), 100)
                mass_pro2_histogram[z_int + 1] += 1
            end
            mass_pro2_distribution = mass_pro2_histogram ./ ( NUM_ATOM_PRO2 * V * NA )

            data_pro1_all[:, t] = mass_pro1_distribution[:]
            data_pro2_all[:, t] = mass_pro2_distribution[:]
        end

        for t in 1 : length( mytraj.conformations )
            for d in data_pro1_all[:, t]
                @printf(data_pro1_output, "%8.6f ", d)
            end
            @printf(data_pro1_output, " \n")

            for d in data_pro2_all[:, t]
                @printf(data_pro2_output, "%8.6f ", d)
            end
            @printf(data_pro2_output, " \n")
        end

        close(data_pro1_output)
        close(data_pro2_output)
    end
end

