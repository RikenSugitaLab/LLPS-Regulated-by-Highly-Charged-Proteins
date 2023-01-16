#!/usr/bin/env julia

include("/home/ctan/Workspace/genesis_CG_julia/src/lib/conformation.jl")
include("/home/ctan/Workspace/genesis_CG_julia/src/lib/parser_dcd.jl")

using Printf

function main()
    # EDIT THIS: temperatures and # runs
    T_list = [i * 10 + 50 for i in 1:30]
    n_run = 5

    # EDIT THIS: number of residues in one chain
    NUM_ATOM_PRO1 = 99

    # EDIT THIS: number of chains
    num_pro1 = 100

    # EDIT THIS: box size (x, y, z: nm)
    V = 18 * 18 * 200 * 1.0e-24
    NA = 6.02e23

    # EDIT THIS: num of bins
    n_bin_z = 100

    for (iT, T) in enumerate(T_list)
        for i_run in 1:n_run
            # EDIT THIS: trajectory name
            if_name = @sprintf("../../hero11_alpha_100_md1_T_%02d_%02d.dcd", iT, i_run)
            println("Processing ", if_name)

            mytraj = read_dcd(if_name)
            box_size_x = mytraj.boundary_box_size[1,end]
            box_size_y = mytraj.boundary_box_size[3,end]
            box_size_z = mytraj.boundary_box_size[6,end]
            bin_size_z = box_size_z / n_bin_z
            n_bin_x = Int(ceil(box_size_x / bin_size_z))
            n_bin_y = Int(ceil(box_size_y / bin_size_z))
            bin_size_x = box_size_x / n_bin_x
            bin_size_y = box_size_y / n_bin_y

            v0 = V / (bin_size_z * bin_size_y * bin_size_z)

            # EDIT THIS: output data file name
            of_pro1_name = @sprintf("hero11_alpha_100_md1_T_%02d_%02d.3d.density.dat", iT, i_run)
            data_pro1_output = open(of_pro1_name, "w")

            data_pro1_all = zeros(Float64, (n_bin_z, length( mytraj.conformations )))
            Threads.@threads for t in 1 : length( mytraj.conformations )
                mass_pro1_histogram = zeros(Int, (n_bin_x, n_bin_y, n_bin_z))

                for i in 1 : num_pro1 * NUM_ATOM_PRO1
                    x = mytraj.conformations[t].coors[1,i]
                    y = mytraj.conformations[t].coors[2,i]
                    z = mytraj.conformations[t].coors[3,i]
                    x_int = Int(floor(mod(x, box_size_x) / bin_size_x))
                    y_int = Int(floor(mod(y, box_size_y) / bin_size_y))
                    z_int = Int(floor(mod(z, box_size_z) / bin_size_z))
                    try
                        mass_pro1_histogram[x_int + 1, y_int + 1, z_int + 1] += 1
                    catch error
                        println(x_int, "  ", y_int, "  ", z_int)
                    end
                end

                mass_pro1_distribution_1d = zeros(Float64, n_bin_z)
                for k_z in 1:n_bin_z
                    n_nonempty_cell = 0
                    count_total = 0
                    for k_x in 1:n_bin_x
                        for k_y in 1:n_bin_y
                            h = mass_pro1_histogram[k_x, k_y, k_z]
                            if h > 0
                                count_total += h
                                n_nonempty_cell += 1
                            end
                        end
                    end
                    if count_total > NUM_ATOM_PRO1 * 2
                        d = count_total / (n_nonempty_cell * v0 * NUM_ATOM_PRO1 * NA)
                    elseif n_nonempty_cell > 0
                        d = count_total / (n_bin_x * n_bin_y * v0 * NUM_ATOM_PRO1 * NA)
                    else
                        d = 0
                    end
                    mass_pro1_distribution_1d[k_z] = d
                end

                data_pro1_all[:, t] = mass_pro1_distribution_1d[:]
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
