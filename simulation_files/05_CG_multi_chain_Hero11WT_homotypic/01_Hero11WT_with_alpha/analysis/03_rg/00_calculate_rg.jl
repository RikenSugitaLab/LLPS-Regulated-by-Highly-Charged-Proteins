#!/usr/bin/env julia

include("/home/ctan/Workspace/genesis_CG_julia/src/lib/gcj.jl")

using Printf

function calculate_contact_map(i_temperature, i_run)

    # ----------
    # parameters
    # ----------
    name_pro1 = "hero11_alpha"

    NUM_ATOM_PRO1 = 99

    num_pro1 = 100

    # EDIT THIS: num-of-frames in traj and num-of-frames to calculate cm
    n_traj_frames = 5000
    i_begin_frame = div(n_traj_frames, 2)
    i_end_frame = n_traj_frames
    # n_structures  = 2000
    # n_interval = div(i_end_frame - i_begin_frame, n_structures)
    n_interval = 1

    # EDIT THIS: trajectory name
    traj_name = @sprintf("../../%s_%d_md1_T_%02d_%02d.dcd", name_pro1, num_pro1, i_temperature, i_run)

    # ---------------
    # data structures
    # ---------------
    rg_list_all = []

    println("Analyzing trajectory: ", traj_name )
    args = Dict("verbose"=>false,
                "begin" => i_begin_frame,
                "end" => i_end_frame,
                "step" => n_interval)
    mytraj = read_dcd(traj_name, args)
    println(mytraj.traj_frames)

    for i_frame in 1:mytraj.traj_frames
        rg_list = []
        myconf = mytraj.conformations[i_frame]

        # ======================================
        # loop over cells and calculate contacts
        # ======================================
        n_atom_count = 0
        for i_chain in 1:num_pro1
            my_coors = myconf.coors[:, n_atom_count + 1:n_atom_count + NUM_ATOM_PRO1]
            rg = radius_of_gyration(my_coors)
            push!(rg_list, rg)
            n_atom_count += NUM_ATOM_PRO1
        end
        push!(rg_list_all, rg_list)
    end

    of_name = @sprintf("%s_rg_t%02d_%02d.dat", name_pro1, i_temperature, i_run)
    offile = open(of_name, "w")
    for rg_list in rg_list_all
        for i_chain in 1:num_pro1
            @printf(offile, " %8.3f ", rg_list[i_chain])
        end
        @printf(offile, " \n")
    end
    close(offile)

end


if abspath(PROGRAM_FILE) == @__FILE__
    for iT in 1:30
        for ir in 1:5
            calculate_contact_map(iT, ir)
        end
    end

end
