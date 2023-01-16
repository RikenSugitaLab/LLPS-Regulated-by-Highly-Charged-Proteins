#!/usr/bin/env julia

include("/home/ctan/Workspace/genesis_CG_julia/src/lib/gcj.jl")

function main()
    rmsd_calc_region_1 = [38:72...]
    rmsd_calc_region_2 = [1:99...]

    my_ref = read_grocrd("../../crd/hero11_alpha.gro")

    my_ref_coor_1 = my_ref.coors[:, rmsd_calc_region_1]
    my_ref_coor_2 = my_ref.coors[:, rmsd_calc_region_2]

    n_temperature = 25
    n_run = 5

    for i_temperature in 1:n_temperature
        for i_run in 1:n_run
            traj_name = @sprintf("../../hero11_alpha_md1_T_%02d_%02d.dcd", i_temperature, i_run)
            println("Analyzing traj: ", traj_name)
            my_traj = read_dcd(traj_name)

            data_out_fname = @sprintf("hero11_alpha_md1_T_%02d_r_%02d.dat", i_temperature, i_run)
            data_out = open(data_out_fname, "w")

            for step in 1:my_traj.traj_frames
                my_conf = my_traj.conformations[step]
                my_coor_1 = my_conf.coors[:, rmsd_calc_region_1]
                my_coor_2 = my_conf.coors[:, rmsd_calc_region_2]
                rmsd_1 = compute_rmsd(my_ref_coor_1, my_coor_1)
                rmsd_2 = compute_rmsd(my_ref_coor_2, my_coor_2)
                @printf(data_out, "%10d    %8.3f     %8.3f  \n", step, rmsd_1, rmsd_2)
            end

            close(data_out)
        end
    end

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
