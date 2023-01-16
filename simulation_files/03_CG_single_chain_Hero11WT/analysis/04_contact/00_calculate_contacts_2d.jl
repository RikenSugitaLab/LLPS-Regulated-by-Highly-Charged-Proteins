#!/usr/bin/env julia

include("/home/ctan/Workspace/genesis_CG_julia/src/lib/gcj.jl")

function main()
    n_res = 99

    n_temperature = 25
    n_run = 5

    for i_temperature in 1:n_temperature
        contact_matrix = zeros(Float64, (n_res, n_res))

        n_struct = 0
        for i_run in 1:n_run
            traj_name = @sprintf("../../hero11_alpha_md1_T_%02d_%02d.dcd", i_temperature, i_run)
            println("Analyzing traj: ", traj_name)
            my_traj = read_dcd(traj_name)

            for step in 1:my_traj.traj_frames
                my_conf = my_traj.conformations[step]
                for i_res in 1:n_res
                    i_coor = my_conf.coors[:, i_res]
                    for j_res in i_res + 3:n_res
                        j_coor = my_conf.coors[:, j_res]
                        d = compute_distance(i_coor, j_coor)
                        c = 1 / (1 + exp(d - 15))
                        contact_matrix[i_res, j_res] += c
                        contact_matrix[j_res, i_res] += c
                    end
                end
                n_struct += 1
            end
        end
        contact_matrix ./= n_struct

        data_out_fname = @sprintf("hero11_alpha_single_chain_contacts_%02d.dat", i_temperature)
        data_out = open(data_out_fname, "w")
        for i_res in 1:n_res
            for j_res in 1:n_res
                @printf(data_out, " %8.3f ", contact_matrix[i_res, j_res])
            end
            println(data_out, " ")
        end

        close(data_out)
    end

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
