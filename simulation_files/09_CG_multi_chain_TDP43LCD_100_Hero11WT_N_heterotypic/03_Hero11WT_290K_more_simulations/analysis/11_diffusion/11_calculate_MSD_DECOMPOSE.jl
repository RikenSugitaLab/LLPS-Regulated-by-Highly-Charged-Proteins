#!/usr/bin/env julia

@everywhere begin
include("/home/ctan/Workspace/genesis_CG_julia/src/lib/gcj.jl")
using DelimitedFiles
using Printf
end

@everywhere function calculate_MSD(temperature, i_conc, i_md, i_run)
    # -------------
    #    parameters
    # -------------
    name_pro1 = "tdp43"
    name_pro2 = "hero11"

    NUM_ATOM_PRO1 = 154
    NUM_ATOM_PRO2 = 99

    num_pro1 = 100
    num_pro2 = i_conc * 10

    box_size_z = 3000
    num_bin = 100
    box_size_bin = box_size_z / num_bin


    # -------------
    # load top file
    # -------------
    top_name = @sprintf( "../../demix_tdp43_ctf_100_hero11_0%2d.top", num_pro2 )
    my_top = read_grotop(top_name)

    # -------------
    # load dcd file
    # -------------
    system_name = @sprintf("%s_%d_%s_%d", name_pro1, num_pro1, name_pro2, num_pro2)
    # EDIT THIS: num-of-frames in traj and num-of-frames to calculate
    n_traj_frames = 3000
    i_begin_frame = 1
    i_end_frame = n_traj_frames
    n_interval = 1
    # EDIT THIS: trajectory name
    traj_name = @sprintf("../../%s_T%3d_md%1d_r%02d.dcd", system_name, temperature, i_md, i_run)

    println("Analyzing trajectory: ", traj_name )
    args = Dict("verbose"=>false,
                "begin" => i_begin_frame,
                "end" => i_end_frame,
                "step" => n_interval)
    mytraj = read_dcd(traj_name)

    # ------------------
    # read position data
    # ------------------
    pro1_phase_data_fname = @sprintf("../03_molecular_position/md%d/%s_T%3d_md%d_r%02d.%s_pos.dat", i_md, system_name, temperature, i_md, i_run, name_pro1)
    pro1_phase_data = readdlm(pro1_phase_data_fname)
    pro2_phase_data_fname = @sprintf("../03_molecular_position/md%d/%s_T%3d_md%d_r%02d.%s_pos.dat", i_md, system_name, temperature, i_md, i_run, name_pro2)
    pro2_phase_data = readdlm(pro2_phase_data_fname)

    # ========================
    # get coordinates for coms
    # ========================
    println(" Calculating COM coordinates...")
    n_real_frames = length(mytraj.conformations)
    pro1_com_coors = zeros(Float64, (3, num_pro1, n_real_frames))
    pro2_com_coors = zeros(Float64, (3, num_pro2, n_real_frames))
    for i_frame in 1:n_real_frames
        my_conf = mytraj.conformations[i_frame]

        n_atom_count = 0
        for i_chain in 1:num_pro1
            com = compute_center_of_mass([n_atom_count+1:n_atom_count+NUM_ATOM_PRO1...], my_top, my_conf)
            pro1_com_coors[:, i_chain, i_frame] = com[:]
            n_atom_count += NUM_ATOM_PRO1
        end
        for i_chain in 1:num_pro2
            com = compute_center_of_mass([n_atom_count+1:n_atom_count+NUM_ATOM_PRO2...], my_top, my_conf)
            pro2_com_coors[:, i_chain, i_frame] = com[:]
            n_atom_count += NUM_ATOM_PRO2
        end
    end

    # ========================================================
    # extract contineous traj pieces and count exchange events
    # ========================================================
    pro1_traj_piece_dense  = [[] for i in 1:num_pro1]
    pro1_traj_piece_dilute = [[] for i in 1:num_pro1]
    pro2_traj_piece_dense  = [[] for i in 1:num_pro2]
    pro2_traj_piece_dilute = [[] for i in 1:num_pro2]
    pro1_count_leave = 0
    pro1_count_enter = 0
    pro2_count_leave = 0
    pro2_count_enter = 0
    for i_chain in 1:num_pro1
        my_state = pro1_phase_data[1, i_chain]
        my_frame_start = 1
        for i_frame in 2:n_real_frames
            i_state = pro1_phase_data[i_frame, i_chain]
            if i_state != my_state
                my_frame_end = i_frame - 1
                if my_state == 1
                    push!( pro1_traj_piece_dense[i_chain] , [my_frame_start, my_frame_end])
                    pro1_count_leave += 1
                else
                    push!( pro1_traj_piece_dilute[i_chain], [my_frame_start, my_frame_end])
                    pro1_count_enter += 1
                end
                my_frame_start = i_frame
                my_state = i_state
            end
        end
    end
    for i_chain in 1:num_pro2
        my_state = pro2_phase_data[1, i_chain]
        my_frame_start = 1
        for i_frame in 2:n_real_frames
            i_state = pro2_phase_data[i_frame, i_chain]
            if i_state != my_state
                my_frame_end = i_frame - 1
                if my_state == 1
                    push!( pro2_traj_piece_dense[i_chain] , [my_frame_start, my_frame_end])
                    pro2_count_leave += 1
                else
                    push!( pro2_traj_piece_dilute[i_chain], [my_frame_start, my_frame_end])
                    pro2_count_enter += 1
                end
                my_frame_start = i_frame
                my_state = i_state
            end
        end
    end
    pro1_event_count_fname = @sprintf("%s_T%3d_md%d_r%02d.%s_diffuse_events.dat", system_name, temperature, i_md, i_run, name_pro1)
    pro1_event_count_file  = open(pro1_event_count_fname, "w")
    @printf(pro1_event_count_file, " # Event of entering condensate: %8d \n", pro1_count_enter)
    @printf(pro1_event_count_file, " # Event of leaving  condensate: %8d \n", pro1_count_leave)
    close(pro1_event_count_file)
    pro2_event_count_fname = @sprintf("%s_T%3d_md%d_r%02d.%s_diffuse_events.dat", system_name, temperature, i_md, i_run, name_pro2)
    pro2_event_count_file  = open(pro2_event_count_fname, "w")
    @printf(pro2_event_count_file, " # Event of entering condensate: %8d \n", pro2_count_enter)
    @printf(pro2_event_count_file, " # Event of leaving  condensate: %8d \n", pro2_count_leave)
    close(pro2_event_count_file)

    # =============
    # calculate MSD
    # =============
    println(" Calculating MSD values...")
    Δt_list = [1, 2, 5, 10, 20, 50]
    n_Δt = length(Δt_list)
    # pro1_MSD_dense_list = zeros(Float64, n_Δt)
    # pro1_MSD_dilute_list = zeros(Float64, n_Δt)
    pro1_MSD_dense_x = zeros(Float64, n_Δt)
    pro1_MSD_dense_y = zeros(Float64, n_Δt)
    pro1_MSD_dense_z = zeros(Float64, n_Δt)
    pro1_MSD_dilute_x = zeros(Float64, n_Δt)
    pro1_MSD_dilute_y = zeros(Float64, n_Δt)
    pro1_MSD_dilute_z = zeros(Float64, n_Δt)
    # pro2_MSD_dense_list = zeros(Float64, n_Δt)
    # pro2_MSD_dilute_list = zeros(Float64, n_Δt)
    pro2_MSD_dense_x = zeros(Float64, n_Δt)
    pro2_MSD_dense_y = zeros(Float64, n_Δt)
    pro2_MSD_dense_z = zeros(Float64, n_Δt)
    pro2_MSD_dilute_x = zeros(Float64, n_Δt)
    pro2_MSD_dilute_y = zeros(Float64, n_Δt)
    pro2_MSD_dilute_z = zeros(Float64, n_Δt)
    for i_Δt in 1:n_Δt
        Δt = Δt_list[i_Δt]
        # ---------------
        # dense diffusion
        # ---------------
        # pro1
        square_displacement_x = []
        square_displacement_y = []
        square_displacement_z = []
        for i_chain in 1:num_pro1
            for j_tp in pro1_traj_piece_dense[i_chain]
                t_start, t_end = j_tp[1], j_tp[2]
                for t1 in t_start:t_end
                    t2 = t1 + Δt
                    if t2 > t_end
                        break
                    end
                    com_t1 = pro1_com_coors[:, i_chain, t1]
                    com_t2 = pro1_com_coors[:, i_chain, t2]
                    com_displacement = com_t1 .- com_t2
                    # com_dist2 = com_displacement' * com_displacement
                    (com_x2, com_y2, com_z2) = com_displacement .* com_displacement
                    # push!(square_displacement_list, com_dist2)
                    push!(square_displacement_x, com_x2)
                    push!(square_displacement_y, com_y2)
                    push!(square_displacement_z, com_z2)
                end
            end
        end
        if length(square_displacement_x) > 0
            # pro1_MSD_dense_list[i_Δt] = sum(square_displacement_list) / length(square_displacement_list)
            pro1_MSD_dense_x[i_Δt] = sum(square_displacement_x) / length(square_displacement_x)
            pro1_MSD_dense_y[i_Δt] = sum(square_displacement_y) / length(square_displacement_y)
            pro1_MSD_dense_z[i_Δt] = sum(square_displacement_z) / length(square_displacement_z)
        else
            # pro1_MSD_dense_list[i_Δt] = -1
            pro1_MSD_dense_x[i_Δt] = -1
            pro1_MSD_dense_y[i_Δt] = -1
            pro1_MSD_dense_z[i_Δt] = -1
        end
        # pro2
        square_displacement_x = []
        square_displacement_y = []
        square_displacement_z = []
        for i_chain in 1:num_pro2
            for j_tp in pro2_traj_piece_dense[i_chain]
                t_start, t_end = j_tp[1], j_tp[2]
                for t1 in t_start:t_end
                    t2 = t1 + Δt
                    if t2 > t_end
                        break
                    end
                    com_t1 = pro2_com_coors[:, i_chain, t1]
                    com_t2 = pro2_com_coors[:, i_chain, t2]
                    com_displacement = com_t1 .- com_t2
                    # com_dist2 = com_displacement' * com_displacement
                    (com_x2, com_y2, com_z2) = com_displacement .* com_displacement
                    # push!(square_displacement_list, com_dist2)
                    push!(square_displacement_x, com_x2)
                    push!(square_displacement_y, com_y2)
                    push!(square_displacement_z, com_z2)
                end
            end
        end
        if length(square_displacement_x) > 0
            pro2_MSD_dense_x[i_Δt] = sum(square_displacement_x) / length(square_displacement_x)
            pro2_MSD_dense_y[i_Δt] = sum(square_displacement_y) / length(square_displacement_y)
            pro2_MSD_dense_z[i_Δt] = sum(square_displacement_z) / length(square_displacement_z)
        else
            pro2_MSD_dense_x[i_Δt] = -1
            pro2_MSD_dense_y[i_Δt] = -1
            pro2_MSD_dense_z[i_Δt] = -1
        end
        # ----------------
        # dilute diffusion
        # ----------------
        # pro1
        square_displacement_x = []
        square_displacement_y = []
        square_displacement_z = []
        for i_chain in 1:num_pro1
            for j_tp in pro1_traj_piece_dilute[i_chain]
                t_start, t_end = j_tp[1], j_tp[2]
                for t1 in t_start:t_end
                    t2 = t1 + Δt
                    if t2 > t_end
                        break
                    end
                    com_t1 = pro1_com_coors[:, i_chain, t1]
                    com_t2 = pro1_com_coors[:, i_chain, t2]
                    com_displacement = com_t1 .- com_t2
                    # com_dist2 = com_displacement' * com_displacement
                    (com_x2, com_y2, com_z2) = com_displacement .* com_displacement
                    # push!(square_displacement_list, com_dist2)
                    push!(square_displacement_x, com_x2)
                    push!(square_displacement_y, com_y2)
                    push!(square_displacement_z, com_z2)
                end
            end
        end
        if length(square_displacement_x) > 0
            pro1_MSD_dilute_x[i_Δt] = sum(square_displacement_x) / length(square_displacement_x)
            pro1_MSD_dilute_y[i_Δt] = sum(square_displacement_y) / length(square_displacement_y)
            pro1_MSD_dilute_z[i_Δt] = sum(square_displacement_z) / length(square_displacement_z)
        else
            pro1_MSD_dilute_x[i_Δt] = -1
            pro1_MSD_dilute_y[i_Δt] = -1
            pro1_MSD_dilute_z[i_Δt] = -1
        end
        # pro2
        square_displacement_x = []
        square_displacement_y = []
        square_displacement_z = []
        for i_chain in 1:num_pro1
            for j_tp in pro1_traj_piece_dilute[i_chain]
                t_start, t_end = j_tp[1], j_tp[2]
                for t1 in t_start:t_end
                    t2 = t1 + Δt
                    if t2 > t_end
                        break
                    end
                    com_t1 = pro1_com_coors[:, i_chain, t1]
                    com_t2 = pro1_com_coors[:, i_chain, t2]
                    com_displacement = com_t1 .- com_t2
                    # com_dist2 = com_displacement' * com_displacement
                    (com_x2, com_y2, com_z2) = com_displacement .* com_displacement
                    # push!(square_displacement_list, com_dist2)
                    push!(square_displacement_x, com_x2)
                    push!(square_displacement_y, com_y2)
                    push!(square_displacement_z, com_z2)
                end
            end
        end
        if length(square_displacement_x) > 0
            # pro2_MSD_dilute_list[i_Δt] = sum(square_displacement_list) / length(square_displacement_list)
            pro2_MSD_dilute_x[i_Δt] = sum(square_displacement_x) / length(square_displacement_x)
            pro2_MSD_dilute_y[i_Δt] = sum(square_displacement_y) / length(square_displacement_y)
            pro2_MSD_dilute_z[i_Δt] = sum(square_displacement_z) / length(square_displacement_z)
        else
            # pro2_MSD_dilute_list[i_Δt] = -1
            pro2_MSD_dilute_x[i_Δt] = -1
            pro2_MSD_dilute_y[i_Δt] = -1
            pro2_MSD_dilute_z[i_Δt] = -1
        end
    end

    pro1_MSD_fname = @sprintf("%s_T%3d_md%d_r%02d.%s_MSD_DECOMPOSED.dat", system_name, temperature, i_md, i_run, name_pro1)
    pro1_MSD_file  = open(pro1_MSD_fname, "w")
    for i_Δt in 1:n_Δt
        # @printf(pro1_MSD_file, "%8d     %18.12f     %18.12f \n", Δt_list[i_Δt], pro1_MSD_dense_list[i_Δt], pro1_MSD_dilute_list[i_Δt])
        @printf(pro1_MSD_file, "%8d  %18.12f  %18.12f  %18.12f    %18.12f  %18.12f  %18.12f \n", Δt_list[i_Δt],
                pro1_MSD_dense_x[i_Δt],
                pro1_MSD_dense_y[i_Δt],
                pro1_MSD_dense_z[i_Δt],
                pro1_MSD_dilute_x[i_Δt],
                pro1_MSD_dilute_y[i_Δt],
                pro1_MSD_dilute_z[i_Δt])
    end
    close(pro1_MSD_file)

    pro2_MSD_fname = @sprintf("%s_T%3d_md%d_r%02d.%s_MSD_DECOMPOSED.dat", system_name, temperature, i_md, i_run, name_pro2)
    pro2_MSD_file  = open(pro2_MSD_fname, "w")
    for i_Δt in 1:n_Δt
        # @printf(pro2_MSD_file, "%8d     %18.12f     %18.12f \n", Δt_list[i_Δt], pro2_MSD_dense_list[i_Δt], pro2_MSD_dilute_list[i_Δt])
        @printf(pro2_MSD_file, "%8d  %18.12f  %18.12f  %18.12f    %18.12f  %18.12f  %18.12f \n", Δt_list[i_Δt],
                pro2_MSD_dense_x[i_Δt],
                pro2_MSD_dense_y[i_Δt],
                pro2_MSD_dense_z[i_Δt],
                pro2_MSD_dilute_x[i_Δt],
                pro2_MSD_dilute_y[i_Δt],
                pro2_MSD_dilute_z[i_Δt])
    end
    close(pro2_MSD_file)

end

if abspath(PROGRAM_FILE) == @__FILE__
    @sync @distributed for i_run in 1:5
        for i_conc in 1:9
            calculate_MSD(290, i_conc, 3, i_run)
        end
    end
end
