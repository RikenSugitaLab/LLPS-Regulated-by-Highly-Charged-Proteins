#!/usr/bin/env julia

include("/home/ctan/Workspace/genesis_CG_julia/src/lib/gcj.jl")

function main()

    n_temperature = 25
    n_run = 5

    for i_temperature in 1:n_temperature
        for i_run in 1:n_run
            traj_name = @sprintf("../../tdp43_md1_T_%02d_%02d.dcd", i_temperature, i_run)
            println("Analyzing traj: ", traj_name)
            my_traj = read_dcd(traj_name)

            data_out_fname = @sprintf("tdp43_md1_T_%02d_r_%02d_rg.dat", i_temperature, i_run)
            data_out = open(data_out_fname, "w")

            for step in 1:my_traj.traj_frames
                my_conf = my_traj.conformations[step]
                rg = radius_of_gyration(my_conf)
                @printf(data_out, "%10d    %8.3f \n", step, rg)
            end

            close(data_out)
        end
    end

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
