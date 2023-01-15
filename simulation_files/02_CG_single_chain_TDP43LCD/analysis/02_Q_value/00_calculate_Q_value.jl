#!/usr/bin/env julia

include("/home/ctan/Workspace/genesis_CG_julia/src/lib/gcj.jl")

function main()
    mytop = read_grotop("../../tdp43.top")

    n_T = 25
    n_run = 5

    for i in 1:n_T
        for j in 1:n_run
            traj_name = @sprintf("../../tdp43_md1_T_%02d_%02d.dcd", i, j)
            println("analyzing qval from traj: ", traj_name)
            my_traj = read_dcd(traj_name)

            fout_name = @sprintf("tdp43_md1_T_%02d_r_%02d_qval.dat", i, j)
            fout = open(fout_name, "w")

            for step in 1:my_traj.traj_frames
                my_conf = my_traj.conformations[step]
                qval = compute_nativeness(mytop, my_conf)
                @printf(fout, "%10d     %8.3f  \n", step, qval)
            end
            close(fout)
        end
    end

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
