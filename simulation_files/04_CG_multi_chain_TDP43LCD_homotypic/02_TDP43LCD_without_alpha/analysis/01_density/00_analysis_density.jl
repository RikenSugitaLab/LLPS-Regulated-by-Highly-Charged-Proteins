#!/usr/bin/env julia

include("/home/ctan/Workspace/genesis_CG_julia/src/lib/conformation.jl")
include("/home/ctan/Workspace/genesis_CG_julia/src/lib/parser_dcd.jl")

using Printf

T_list = [i * 10 + 250 for i in 1:10]

NUM_ATOM_HERO = 99
NUM_ATOM_TDP43 = 154

num_tdp43 = 100
num_all = 100
num_hero = 0

V = 18 * 18 * 3 * 1.0e-24
NA = 6.02e23

for (iT, T) in enumerate(T_list)
    # if_name = @sprintf("../../tdp43_100_md1_T_%02d.dcd", iT)
    if_name = @sprintf("../../tdp43_100_md2_T_%02d.dcd", iT)
    println("Processing ", if_name )

    mytraj = read_dcd(if_name)
    box_size_z = mytraj.boundary_box_size[6,end]
    box_size_bin = box_size_z / 100

    # of_tdp43_name = @sprintf("tdp43_md1_T_%02d.density.dat", iT)
    of_tdp43_name = @sprintf("tdp43_md2_T_%02d.density.dat", iT)
    data_tdp43_output = open(of_tdp43_name, "w")
    # of_hero_name = @sprintf("hero_md1_t%3d_%02d.density.dat", T, i)
    # data_hero_output = open(of_hero_name, "w")

    data_tdp43_all = zeros(Float64, (100, length( mytraj.conformations )))
    # data_hero_all = zeros(Float64, (100, length( mytraj.conformations )))
    Threads.@threads for t in 1 : length( mytraj.conformations )
        mass_tdp43_histogram = [0 for j in 1:100]
        # mass_hero_histogram = [0 for j in 1:100]

        for i in 1 : num_tdp43 * NUM_ATOM_TDP43
            z = mytraj.conformations[t].coors[3,i]
            z_int = Int(floor(mod(z, box_size_z) / box_size_bin))
            try
                mass_tdp43_histogram[z_int + 1] += 1
            catch error
                println(z_int)
            end
        end
        mass_tdp43_distribution = mass_tdp43_histogram ./ ( NUM_ATOM_TDP43 * V * NA )

        # for i in num_tdp43 * NUM_ATOM_TDP43 + 1 : mytraj.num_atoms
        #     z = mytraj.conformations[t].coors[3,i]
        #     z_int = Int(floor(mod(z, box_size_z) / box_size_bin))
        #     try
        #         mass_hero_histogram[z_int + 1] += 1
        #     catch error
        #         println(z_int)
        #     end
        # end
        # mass_hero_distribution = mass_hero_histogram ./ ( NUM_ATOM_HERO * V * NA )

        data_tdp43_all[:, t] = mass_tdp43_distribution[:]
        # data_hero_all[:, t] = mass_hero_distribution[:]
    end

    for t in 1 : length( mytraj.conformations )
        for d in data_tdp43_all[:, t]
            @printf(data_tdp43_output, "%8.6f ", d)
        end
        @printf(data_tdp43_output, " \n")

        # for d in data_hero_all[:, t]
        #     @printf(data_hero_output, "%8.6f ", d)
        # end
        # @printf(data_hero_output, " \n")
    end

    # close(data_hero_output)
    close(data_tdp43_output)
end
