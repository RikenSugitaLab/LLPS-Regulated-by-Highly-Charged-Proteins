#!/usr/bin/env julia

include("/home/ctan/Workspace/genesis_CG_julia/src/lib/gcj.jl")

using Printf
using DelimitedFiles

function calculate_rg(temperature, n_hero, i_md, i_run)

    # ----------
    # parameters
    # ----------
    name_pro1 = "tdp43"
    name_pro2 = "hero11"

    NUM_ATOM_PRO1 = 154
    NUM_ATOM_PRO2 = 99

    num_pro1 = 100
    num_pro2 = n_hero

    system_name = @sprintf("%s_%d_%s_%d", name_pro1, num_pro1, name_pro2, num_pro2)

    # EDIT THIS: num-of-frames in traj and num-of-frames to calculate
    n_traj_frames = 3000
    # i_begin_frame = div(n_traj_frames, 2)
    i_begin_frame = 1
    i_end_frame = n_traj_frames
    n_structures  = 3000
    # n_interval = div(i_end_frame - i_begin_frame, n_structures)
    n_interval = 1

    # EDIT THIS: trajectory name
    traj_name = @sprintf("../../%s_T%03d_md%1d_r%02d.dcd", system_name, temperature, i_md, i_run)

    # ---------------
    # data structures
    # ---------------
    rg_pro1_list_all = zeros(Float64, (num_pro1, n_structures))
    rg_pro2_list_all = zeros(Float64, (num_pro2, n_structures))

    println("Analyzing trajectory: ", traj_name )
    args = Dict("verbose"=>false,
                "begin" => i_begin_frame,
                "end" => i_end_frame,
                "step" => n_interval)
    # mytraj = read_dcd(traj_name, args)
    mytraj = read_dcd(traj_name)

    Threads.@threads for i_frame in 1:length(mytraj.conformations)
        myconf = mytraj.conformations[i_frame]

        n_atom_count = 0
        for i_chain in 1:num_pro1
            my_coors = myconf.coors[:, n_atom_count + 1:n_atom_count + NUM_ATOM_PRO1]
            rg = radius_of_gyration(my_coors)
            rg_pro1_list_all[i_chain, i_frame] = rg
            n_atom_count += NUM_ATOM_PRO1
        end

        for i_chain in 1:num_pro2
            my_coors = myconf.coors[:, n_atom_count + 1:n_atom_count + NUM_ATOM_PRO2]
            rg = radius_of_gyration(my_coors)
            rg_pro2_list_all[i_chain, i_frame] = rg
            n_atom_count += NUM_ATOM_PRO2
        end
    end

    of1_name = @sprintf("%s_T%03d_md%d_r%02d.%s_rg.dat", system_name, temperature, i_md, i_run, name_pro1)
    of2_name = @sprintf("%s_T%03d_md%d_r%02d.%s_rg.dat", system_name, temperature, i_md, i_run, name_pro2)
    of1 = open(of1_name, "w")
    of2 = open(of2_name, "w")
    for i_frame in 1:n_structures
        for i_chain in 1:num_pro1
            @printf(of1, " %8.3f ", rg_pro1_list_all[i_chain, i_frame])
        end
        @printf(of1, " \n")
    end
    close(of1)
    for i_frame in 1:n_structures
        for i_chain in 1:num_pro2
            @printf(of2, " %8.3f ", rg_pro2_list_all[i_chain, i_frame])
        end
        @printf(of2, " \n")
    end
    close(of2)

end

if abspath(PROGRAM_FILE) == @__FILE__
    for iT in 1:2
        for ic in 1:9
            T = iT * 5 + 285
            for irun in 1:5
                calculate_rg(T, ic * 10, 1, irun)
                calculate_rg(T, ic * 10, 2, irun)
                calculate_rg(T, ic * 10, 3, irun)
            end
        end
    end

end
