#!/usr/bin/env julia

include("/home/ctan/Workspace/genesis_CG_julia/src/lib/conformation.jl")
include("/home/ctan/Workspace/genesis_CG_julia/src/lib/parser_dcd.jl")

using Printf

function main()
    # EDIT THIS: temperatures
    T_list = [i * 10 + 50 for i in 1:30]
    n_run = 5

    # EDIT THIS: number of residues in one chain
    NUM_ATOM_PRO1 = 99

    # EDIT THIS: number of chains
    num_pro1 = 100

    # EDIT THIS: box size (x, y: nm; z: 100nm)
    V = 18 * 18 * 1 * 1.0e-24
    NA = 6.02e23

    for (iT, T) in enumerate(T_list)
        for i_run in 1:n_run
            # EDIT THIS: trajectory name
            if_name = @sprintf("../../hero11_alpha_100_md1_T_%02d_%02d.dcd", iT, i_run)
            println("Processing ", if_name)

            mytraj = read_dcd(if_name)
            box_size_z = mytraj.boundary_box_size[6,end]
            box_size_bin = box_size_z / 100

            # EDIT THIS: output data file name
            of_pro1_name = @sprintf("hero11_alpha_100_md1_T_%02d_%02d.density.dat", iT, i_run)
            data_pro1_output = open(of_pro1_name, "w")

            data_pro1_all = zeros(Float64, (100, length( mytraj.conformations )))
            Threads.@threads for t in 1 : length( mytraj.conformations )
                mass_pro1_histogram = [0 for j in 1:100]

                for i in 1 : num_pro1 * NUM_ATOM_PRO1
                    z = mytraj.conformations[t].coors[3,i]
                    z_int = Int(floor(mod(z, box_size_z) / box_size_bin))
                    try
                        mass_pro1_histogram[z_int + 1] += 1
                    catch error
                        println(z_int)
                    end
                end
                mass_pro1_distribution = mass_pro1_histogram ./ ( NUM_ATOM_PRO1 * V * NA )

                data_pro1_all[:, t] = mass_pro1_distribution[:]
            end

            for t in 1 : length( mytraj.conformations )
                for d in data_pro1_all[:, t]
                    @printf(data_pro1_output, "%8.6f ", d)
                end
                @printf(data_pro1_output, " \n")
            end
            close(data_pro1_output)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
