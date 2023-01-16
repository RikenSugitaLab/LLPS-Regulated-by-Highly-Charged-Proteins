#!/usr/bin/env julia

# @everywhere begin
include("/home/ctan/Workspace/genesis_CG_julia/src/lib/conformation.jl")
include("/home/ctan/Workspace/genesis_CG_julia/src/lib/parser_dcd.jl")
using Printf
# end

# @everywhere function calculate_density_along_z(T, n_hero, i_md, i_run)
function calculate_density_along_z(T, n_hero, i_md, i_run)
    # EDIT THIS: protein names
    NAME_PRO1 = "tdp43"
    NAME_PRO2 = "hero11"
    # EDIT THIS: number of residues in one chain
    NUM_ATOM_PRO1 = 154
    NUM_ATOM_PRO2 = 99
    # EDIT THIS: number of chains
    num_pro1 = 100
    num_pro2 = n_hero
    num_all = num_pro1 + num_pro2

    # EDIT THIS: box size (x, y: nm; z: 100nm)
    V = 18 * 18 * 3 * 1.0e-24
    NA = 6.02e23

    # EDIT THIS: box size z
    box_size_z = 3000
    box_size_bin = box_size_z / 100

    # EDIT THIS: trajectory file name
    if_name = @sprintf("../../%s_%d_%s_%d_T%03d_md%d_r%02d.dcd", NAME_PRO1, num_pro1, NAME_PRO2, num_pro2, T, i_md, i_run)
    println("Processing ", if_name)

    mytraj = read_dcd(if_name)

    # EDIT THIS: output data file name
    of_pro1_name = @sprintf("%s_%d_%s_%d_T%03d_md%d_r%02d.%s.density.dat", NAME_PRO1, num_pro1, NAME_PRO2, num_pro2, T, i_md, i_run,  NAME_PRO1)
    of_pro2_name = @sprintf("%s_%d_%s_%d_T%03d_md%d_r%02d.%s.density.dat", NAME_PRO1, num_pro1, NAME_PRO2, num_pro2, T, i_md, i_run,  NAME_PRO2)

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

if abspath(PROGRAM_FILE) == @__FILE__
    # @sync @distributed for irun in 1:5
    for irun in 1:5
        for ic in 1:9
            for iT  in 1:2
                T = 285 + 5 * iT
                n_hero = ic * 10
                calculate_density_along_z(T, n_hero, 3, irun)
            end
        end
    end
end




