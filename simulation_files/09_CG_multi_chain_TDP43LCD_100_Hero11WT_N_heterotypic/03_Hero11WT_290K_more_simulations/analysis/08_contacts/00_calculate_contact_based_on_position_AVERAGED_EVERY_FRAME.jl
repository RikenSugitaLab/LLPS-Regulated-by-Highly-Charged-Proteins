#!/usr/bin/env julia

@everywhere begin
include("/home/ctan/Workspace/genesis_CG_julia/src/lib/gcj.jl")
using Printf
using DelimitedFiles
using Statistics
end


@everywhere function calculate_contact_map(temperature, n_hero, i_md, i_run)
# function calculate_contact_map(temperature, n_hero, i_run)

    # ----------
    # parameters
    # ----------
    contact_r0     = 15.0
    contact_cutoff = 25.0
    com_cutoff     = 250.0

    name_pro1 = "tdp43"
    name_pro2 = "hero11"

    NUM_ATOM_PRO1 = 154
    NUM_ATOM_PRO2 = 99

    num_pro1 = 100
    num_pro2 = n_hero

    system_name = @sprintf("%s_%d_%s_%d", name_pro1, num_pro1, name_pro2, num_pro2)

    # EDIT THIS: trajectory name
    traj_name = @sprintf("../../%s_T%03d_md%1d_r%02d.dcd", system_name, temperature, i_md, i_run)

    top_name = @sprintf("../../demix_tdp43_ctf_100_hero11_%03d.top", n_hero)
    mytop = read_grotop(top_name)


    # --------------------------
    # read position-density data
    # --------------------------
    pro1_phase_data_fname = @sprintf("../03_molecular_position/md%d/%s_T%3d_md%d_r%02d.%s_pos.dat", i_md, system_name, temperature, i_md, i_run, name_pro1)
    pro2_phase_data_fname = @sprintf("../03_molecular_position/md%d/%s_T%3d_md%d_r%02d.%s_pos.dat", i_md, system_name, temperature, i_md, i_run, name_pro2)
    pro1_phase_data = readdlm(pro1_phase_data_fname)
    pro2_phase_data = readdlm(pro2_phase_data_fname)

    # ---------------------------
    # prepare and read trajectory
    # ---------------------------
    println("Analyzing trajectory: ", traj_name )

    # EDIT THIS: num-of-frames in traj and num-of-frames to calculate
    n_traj_frames = length(pro1_phase_data[:, 1])
    println("Length of original traj:", n_traj_frames)
    # i_begin_frame = div(n_traj_frames, 2)
    i_begin_frame = 1
    i_end_frame = n_traj_frames
    n_structures  = 50
    n_interval = div(i_end_frame - i_begin_frame, n_structures)
    # n_interval = 1

    args = Dict("verbose"=>false,
                "begin" => i_begin_frame,
                "end" => i_end_frame,
                "step" => n_interval)
    mytraj = read_dcd(traj_name, args)
    # mytraj = read_dcd(traj_name)
    # println(mytraj.traj_frames)

    # ---------------
    # data structures
    # ---------------
    # parallelization
    num_frames  = length(mytraj.conformations)
    cmtx_pro1_pro1_dense  = zeros(Float64, (NUM_ATOM_PRO1, NUM_ATOM_PRO1, num_frames))
    cmtx_pro2_pro2_dense  = zeros(Float64, (NUM_ATOM_PRO2, NUM_ATOM_PRO2, num_frames))
    cmtx_pro1_pro2_dense  = zeros(Float64, (NUM_ATOM_PRO1, NUM_ATOM_PRO2, num_frames))
    cmtx_pro2_pro1_dense  = zeros(Float64, (NUM_ATOM_PRO2, NUM_ATOM_PRO1, num_frames))
    cmtx_pro1_pro1_dilute = zeros(Float64, (NUM_ATOM_PRO1, NUM_ATOM_PRO1, num_frames))
    cmtx_pro2_pro2_dilute = zeros(Float64, (NUM_ATOM_PRO2, NUM_ATOM_PRO2, num_frames))
    cmtx_pro1_pro2_dilute = zeros(Float64, (NUM_ATOM_PRO1, NUM_ATOM_PRO2, num_frames))
    cmtx_pro2_pro1_dilute = zeros(Float64, (NUM_ATOM_PRO2, NUM_ATOM_PRO1, num_frames))
    cc_pro1_11_dense =  zeros(Float64, (num_pro1, num_frames))
    cc_pro1_12_dense =  zeros(Float64, (num_pro1, num_frames))
    cc_pro2_21_dense =  zeros(Float64, (num_pro2, num_frames))
    cc_pro2_22_dense =  zeros(Float64, (num_pro2, num_frames))
    cc_pro1_11_dilute = zeros(Float64, (num_pro1, num_frames))
    cc_pro1_12_dilute = zeros(Float64, (num_pro1, num_frames))
    cc_pro2_21_dilute = zeros(Float64, (num_pro2, num_frames))
    cc_pro2_22_dilute = zeros(Float64, (num_pro2, num_frames))

    com_pro1_coors = zeros(Float64, (3, num_pro1, num_frames))
    com_pro2_coors = zeros(Float64, (3, num_pro2, num_frames))

    # println(" > > calculating COM of chains...")
    # -------------------
    # analyze every frame
    # -------------------
    # my_pid = myid()
    # my_tid = Threads.threadid()
    for i_frame in 1:num_frames
        # println(" > frame: ", i_frame, " of T ", temperature, " #hero ", n_hero, " on process: ", my_pid, " on thread: ", my_tid)
        println(" > frame: ", i_frame, " of T ", temperature, " #hero ", n_hero)
        myconf = mytraj.conformations[i_frame]
        i_step_real = (i_frame - 1) * n_interval + 1

        # ---------------------------------
        # count num of chains in two phases
        # ---------------------------------
        n1_dense = sum(pro1_phase_data[i_step_real, :] .> 0.5)
        n2_dense = sum(pro2_phase_data[i_step_real, :] .> 0.5)
        n1_dilute = num_pro1 - n1_dense
        n2_dilute = num_pro2 - n2_dense
        n1_dilute = n1_dilute < 1 ? 1 : n1_dilute
        n2_dilute = n2_dilute < 1 ? 1 : n2_dilute

        # Get box size
        box_size_x = mytraj.boundary_box_size[1, i_frame]
        box_size_y = mytraj.boundary_box_size[3, i_frame]
        box_size_z = mytraj.boundary_box_size[6, i_frame]
        box_size = [box_size_x, box_size_y, box_size_z]

        # calculate com
        n_res_count = 0
        for i_chain in 1:num_pro1
            com_pro1_coors[:, i_chain, i_frame] = mean(myconf.coors[:, n_res_count + 1:n_res_count + NUM_ATOM_PRO1], dims=2)
            n_res_count += NUM_ATOM_PRO1
        end
        for j_chain in 1:num_pro2
            com_pro2_coors[:, j_chain, i_frame] = mean(myconf.coors[:, n_res_count + 1:n_res_count + NUM_ATOM_PRO2], dims=2)
            n_res_count += NUM_ATOM_PRO2
        end

        # =====================
        # loop over chain pairs
        # =====================
        # println(" > > calculating pairwise residue-residue contacts...")
        # -----------
        # pro1 - pro1
        # -----------
        for i_chain in 1:num_pro1
            i_chain_phase_state = pro1_phase_data[i_step_real, i_chain]
            i_shift = (i_chain - 1) * NUM_ATOM_PRO1
            com_i = com_pro1_coors[:, i_chain, i_frame]
            for j_chain in i_chain + 1:num_pro1
                j_chain_phase_state = pro1_phase_data[i_step_real, j_chain]

                if abs(i_chain_phase_state - j_chain_phase_state) > 0.5
                    continue
                end

                j_shift = (j_chain - 1) * NUM_ATOM_PRO1
                com_j = com_pro1_coors[:, j_chain, i_frame]
                d_r = com_i - com_j
                com_dr_in_cell = d_r - round.(d_r ./ box_size) .* box_size
                com_dist = sqrt(com_dr_in_cell' * com_dr_in_cell)
                if com_dist > com_cutoff
                    continue
                end

                # begin residue pairwise computation
                for i_res in 1:NUM_ATOM_PRO1
                    coor_i = myconf.coors[:, i_res + i_shift]
                    for j_res in 1:NUM_ATOM_PRO1
                        coor_j = myconf.coors[:, j_res + j_shift]
                        d_r = coor_i - coor_j
                        d_r_in_cell = d_r - round.(d_r ./ box_size) .* box_size
                        dist = sqrt(d_r_in_cell' * d_r_in_cell)
                        if dist > contact_cutoff
                            continue
                        end
                        cnt = 1 / (1 + exp(dist - contact_r0))
                        if i_chain_phase_state * j_chain_phase_state > 0.5
                            cmtx_pro1_pro1_dense[i_res, j_res, i_frame] += cnt
                            cmtx_pro1_pro1_dense[j_res, i_res, i_frame] += cnt
                            cc_pro1_11_dense[i_chain, i_frame] += cnt
                            cc_pro1_11_dense[j_chain, i_frame] += cnt
                        else
                            cmtx_pro1_pro1_dilute[i_res, j_res, i_frame] += cnt
                            cmtx_pro1_pro1_dilute[j_res, i_res, i_frame] += cnt
                            cc_pro1_11_dilute[i_chain, i_frame] += cnt
                            cc_pro1_11_dilute[j_chain, i_frame] += cnt
                        end
                    end
                end
                # end residue pairwise computation
            end
        end

        # -----------
        # pro2 - pro2
        # -----------
        for i_chain in 1:num_pro2
            i_chain_phase_state = pro2_phase_data[i_step_real, i_chain]
            i_shift = (i_chain - 1) * NUM_ATOM_PRO2 + num_pro1 * NUM_ATOM_PRO1
            com_i = com_pro2_coors[:, i_chain, i_frame]
            for j_chain in i_chain + 1:num_pro2
                j_chain_phase_state = pro2_phase_data[i_step_real, j_chain]

                if abs(i_chain_phase_state - j_chain_phase_state) > 0.5
                    continue
                end

                j_shift = (j_chain - 1) * NUM_ATOM_PRO2 + num_pro1 * NUM_ATOM_PRO1
                com_j = com_pro2_coors[:, j_chain, i_frame]
                d_r = com_i - com_j
                com_dr_in_cell = d_r - round.(d_r ./ box_size) .* box_size
                com_dist = sqrt(com_dr_in_cell' * com_dr_in_cell)
                if com_dist > com_cutoff
                    continue
                end

                # begin residue pairwise computation
                for i_res in 1:NUM_ATOM_PRO2
                    coor_i = myconf.coors[:, i_res + i_shift]
                    for j_res in 1:NUM_ATOM_PRO2
                        coor_j = myconf.coors[:, j_res + j_shift]
                        d_r = coor_i - coor_j
                        d_r_in_cell = d_r - round.(d_r ./ box_size) .* box_size
                        dist = sqrt(d_r_in_cell' * d_r_in_cell)
                        if dist > contact_cutoff
                            continue
                        end
                        cnt = 1 / (1 + exp(dist - contact_r0))
                        if i_chain_phase_state * j_chain_phase_state > 0.5
                            cmtx_pro2_pro2_dense[i_res, j_res, i_frame] += cnt
                            cmtx_pro2_pro2_dense[j_res, i_res, i_frame] += cnt
                            cc_pro2_22_dense[i_chain, i_frame] += cnt
                            cc_pro2_22_dense[j_chain, i_frame] += cnt
                        else
                            cmtx_pro2_pro2_dilute[i_res, j_res, i_frame] += cnt
                            cmtx_pro2_pro2_dilute[j_res, i_res, i_frame] += cnt
                            cc_pro2_22_dilute[i_chain, i_frame] += cnt
                            cc_pro2_22_dilute[j_chain, i_frame] += cnt
                        end
                    end
                end
                # end residue pairwise computation
            end
        end

        # -----------
        # pro1 - pro2
        # -----------
        for i_chain in 1:num_pro1
            i_chain_phase_state = pro1_phase_data[i_step_real, i_chain]
            i_shift = (i_chain - 1) * NUM_ATOM_PRO1
            com_i = com_pro1_coors[:, i_chain, i_frame]
            for j_chain in 1:num_pro2
                j_chain_phase_state = pro2_phase_data[i_step_real, j_chain]

                if abs(i_chain_phase_state - j_chain_phase_state) > 0.5
                    continue
                end

                j_shift = (j_chain - 1) * NUM_ATOM_PRO2 + num_pro1 * NUM_ATOM_PRO1
                com_j = com_pro2_coors[:, j_chain, i_frame]
                d_r = com_i - com_j
                com_dr_in_cell = d_r - round.(d_r ./ box_size) .* box_size
                com_dist = sqrt(com_dr_in_cell' * com_dr_in_cell)
                if com_dist > com_cutoff
                    continue
                end

                # begin residue pairwise computation
                for i_res in 1:NUM_ATOM_PRO1
                    coor_i = myconf.coors[:, i_res + i_shift]
                    for j_res in 1:NUM_ATOM_PRO2
                        coor_j = myconf.coors[:, j_res + j_shift]
                        d_r = coor_i - coor_j
                        d_r_in_cell = d_r - round.(d_r ./ box_size) .* box_size
                        dist = sqrt(d_r_in_cell' * d_r_in_cell)
                        if dist > contact_cutoff
                            continue
                        end
                        cnt = 1 / (1 + exp(dist - contact_r0))
                        if i_chain_phase_state * j_chain_phase_state > 0.5
                            cmtx_pro1_pro2_dense[i_res, j_res, i_frame] += cnt
                            cmtx_pro2_pro1_dense[j_res, i_res, i_frame] += cnt
                            cc_pro1_12_dense[i_chain, i_frame] += cnt
                            cc_pro2_21_dense[j_chain, i_frame] += cnt
                        else
                            cmtx_pro1_pro2_dilute[i_res, j_res, i_frame] += cnt
                            cmtx_pro2_pro1_dilute[j_res, i_res, i_frame] += cnt
                            cc_pro1_12_dilute[i_chain, i_frame] += cnt
                            cc_pro2_21_dilute[j_chain, i_frame] += cnt
                        end
                    end
                end
                # end residue pairwise computation
            end
        end

        cmtx_pro1_pro1_dense[:, :,  i_frame] ./= n1_dense
        cmtx_pro1_pro1_dilute[:, :, i_frame] ./= n1_dilute
        cmtx_pro2_pro2_dense[:, :,  i_frame] ./= n2_dense
        cmtx_pro2_pro2_dilute[:, :, i_frame] ./= n2_dilute
        cmtx_pro1_pro2_dense[:, :,  i_frame] ./= n1_dense
        cmtx_pro1_pro2_dilute[:, :, i_frame] ./= n1_dilute
        cmtx_pro2_pro1_dense[:, :,  i_frame] ./= n2_dense
        cmtx_pro2_pro1_dilute[:, :, i_frame] ./= n2_dilute

        cc_pro1_11_dense[ findall(pro1_phase_data[i_step_real, :] .< 0.5), i_frame] .= -1.0
        cc_pro1_12_dense[ findall(pro1_phase_data[i_step_real, :] .< 0.5), i_frame] .= -1.0
        cc_pro2_21_dense[ findall(pro2_phase_data[i_step_real, :] .< 0.5), i_frame] .= -1.0
        cc_pro2_22_dense[ findall(pro2_phase_data[i_step_real, :] .< 0.5), i_frame] .= -1.0
        cc_pro1_11_dilute[findall(pro1_phase_data[i_step_real, :] .> 0.5), i_frame] .= -1.0
        cc_pro1_12_dilute[findall(pro1_phase_data[i_step_real, :] .> 0.5), i_frame] .= -1.0
        cc_pro2_21_dilute[findall(pro2_phase_data[i_step_real, :] .> 0.5), i_frame] .= -1.0
        cc_pro2_22_dilute[findall(pro2_phase_data[i_step_real, :] .> 0.5), i_frame] .= -1.0

    end                         # loop over frame

    # ----------------------
    # calculate the averages
    # ----------------------

    mean_cmtx_pro1_pro1_dense  = mean(cmtx_pro1_pro1_dense,  dims=3)[:, :, 1]
    mean_cmtx_pro1_pro1_dilute = mean(cmtx_pro1_pro1_dilute, dims=3)[:, :, 1]
    mean_cmtx_pro2_pro2_dense  = mean(cmtx_pro2_pro2_dense,  dims=3)[:, :, 1]
    mean_cmtx_pro2_pro2_dilute = mean(cmtx_pro2_pro2_dilute, dims=3)[:, :, 1]
    mean_cmtx_pro1_pro2_dense  = mean(cmtx_pro1_pro2_dense,  dims=3)[:, :, 1]
    mean_cmtx_pro1_pro2_dilute = mean(cmtx_pro1_pro2_dilute, dims=3)[:, :, 1]
    mean_cmtx_pro2_pro1_dense  = mean(cmtx_pro2_pro1_dense,  dims=3)[:, :, 1]
    mean_cmtx_pro2_pro1_dilute = mean(cmtx_pro2_pro1_dilute, dims=3)[:, :, 1]


    # =========
    # output...
    # =========
    pro1_pro1_dense_of_name  = @sprintf("%s_T%03d_md%d_r%02d.%s-%s.%s.contact_matrix.dat", system_name, temperature, i_md, i_run, name_pro1, name_pro1, "dense")
    pro2_pro2_dense_of_name  = @sprintf("%s_T%03d_md%d_r%02d.%s-%s.%s.contact_matrix.dat", system_name, temperature, i_md, i_run, name_pro2, name_pro2, "dense")
    pro1_pro2_dense_of_name  = @sprintf("%s_T%03d_md%d_r%02d.%s-%s.%s.contact_matrix.dat", system_name, temperature, i_md, i_run, name_pro1, name_pro2, "dense")
    pro2_pro1_dense_of_name  = @sprintf("%s_T%03d_md%d_r%02d.%s-%s.%s.contact_matrix.dat", system_name, temperature, i_md, i_run, name_pro2, name_pro1, "dense")
    pro1_pro1_dilute_of_name = @sprintf("%s_T%03d_md%d_r%02d.%s-%s.%s.contact_matrix.dat", system_name, temperature, i_md, i_run, name_pro1, name_pro1, "dilute")
    pro2_pro2_dilute_of_name = @sprintf("%s_T%03d_md%d_r%02d.%s-%s.%s.contact_matrix.dat", system_name, temperature, i_md, i_run, name_pro2, name_pro2, "dilute")
    pro1_pro2_dilute_of_name = @sprintf("%s_T%03d_md%d_r%02d.%s-%s.%s.contact_matrix.dat", system_name, temperature, i_md, i_run, name_pro1, name_pro2, "dilute")
    pro2_pro1_dilute_of_name = @sprintf("%s_T%03d_md%d_r%02d.%s-%s.%s.contact_matrix.dat", system_name, temperature, i_md, i_run, name_pro2, name_pro1, "dilute")

    pro1_pro1_dense_of = open(pro1_pro1_dense_of_name, "w")
    pro2_pro2_dense_of = open(pro2_pro2_dense_of_name, "w")
    pro1_pro2_dense_of = open(pro1_pro2_dense_of_name, "w")
    pro2_pro1_dense_of = open(pro2_pro1_dense_of_name, "w")
    pro1_pro1_dilute_of = open(pro1_pro1_dilute_of_name, "w")
    pro2_pro2_dilute_of = open(pro2_pro2_dilute_of_name, "w")
    pro1_pro2_dilute_of = open(pro1_pro2_dilute_of_name, "w")
    pro2_pro1_dilute_of = open(pro2_pro1_dilute_of_name, "w")

    for i in 1:NUM_ATOM_PRO1
        for j in 1:NUM_ATOM_PRO1
            @printf(pro1_pro1_dense_of,  " %8.3f ", mean_cmtx_pro1_pro1_dense[i, j])
            @printf(pro1_pro1_dilute_of, " %8.3f ", mean_cmtx_pro1_pro1_dilute[i, j])
        end
        @printf(pro1_pro1_dense_of, " \n")
        @printf(pro1_pro1_dilute_of, " \n")
    end
    close(pro1_pro1_dense_of)
    close(pro1_pro1_dilute_of)

    for i in 1:NUM_ATOM_PRO2
        for j in 1:NUM_ATOM_PRO2
            @printf(pro2_pro2_dense_of,  " %8.3f ", mean_cmtx_pro2_pro2_dense[i, j])
            @printf(pro2_pro2_dilute_of, " %8.3f ", mean_cmtx_pro2_pro2_dilute[i, j])
        end
        @printf(pro2_pro2_dense_of, " \n")
        @printf(pro2_pro2_dilute_of, " \n")
    end
    close(pro2_pro2_dense_of)
    close(pro2_pro2_dilute_of)

    for i in 1:NUM_ATOM_PRO1
        for j in 1:NUM_ATOM_PRO2
            @printf(pro1_pro2_dense_of,  " %8.3f ", mean_cmtx_pro1_pro2_dense[i, j])
            @printf(pro1_pro2_dilute_of, " %8.3f ", mean_cmtx_pro1_pro2_dilute[i, j])
        end
        @printf(pro1_pro2_dense_of, " \n")
        @printf(pro1_pro2_dilute_of, " \n")
    end
    close(pro1_pro2_dense_of)
    close(pro1_pro2_dilute_of)

    for i in 1:NUM_ATOM_PRO2
        for j in 1:NUM_ATOM_PRO1
            @printf(pro2_pro1_dense_of,  " %8.3f ", mean_cmtx_pro2_pro1_dense[i, j])
            @printf(pro2_pro1_dilute_of, " %8.3f ", mean_cmtx_pro2_pro1_dilute[i, j])
        end
        @printf(pro2_pro1_dense_of, " \n")
        @printf(pro2_pro1_dilute_of, " \n")
    end
    close(pro2_pro1_dense_of)
    close(pro2_pro1_dilute_of)


    cc_p1_11_dense_of_name   = @sprintf("%s_T%03d_md%d_r%02d.%s-%s.%s.contact_count_ts.dat", system_name, temperature, i_md, i_run, name_pro1, name_pro1, "dense")
    cc_p1_12_dense_of_name   = @sprintf("%s_T%03d_md%d_r%02d.%s-%s.%s.contact_count_ts.dat", system_name, temperature, i_md, i_run, name_pro1, name_pro2, "dense")
    cc_p2_21_dense_of_name   = @sprintf("%s_T%03d_md%d_r%02d.%s-%s.%s.contact_count_ts.dat", system_name, temperature, i_md, i_run, name_pro2, name_pro1, "dense")
    cc_p2_22_dense_of_name   = @sprintf("%s_T%03d_md%d_r%02d.%s-%s.%s.contact_count_ts.dat", system_name, temperature, i_md, i_run, name_pro2, name_pro2, "dense")
    cc_p1_11_dilute_of_name  = @sprintf("%s_T%03d_md%d_r%02d.%s-%s.%s.contact_count_ts.dat", system_name, temperature, i_md, i_run, name_pro1, name_pro1, "dilute")
    cc_p1_12_dilute_of_name  = @sprintf("%s_T%03d_md%d_r%02d.%s-%s.%s.contact_count_ts.dat", system_name, temperature, i_md, i_run, name_pro1, name_pro2, "dilute")
    cc_p2_21_dilute_of_name  = @sprintf("%s_T%03d_md%d_r%02d.%s-%s.%s.contact_count_ts.dat", system_name, temperature, i_md, i_run, name_pro2, name_pro1, "dilute")
    cc_p2_22_dilute_of_name  = @sprintf("%s_T%03d_md%d_r%02d.%s-%s.%s.contact_count_ts.dat", system_name, temperature, i_md, i_run, name_pro2, name_pro2, "dilute")
    cc_11_dense_of  = open(cc_p1_11_dense_of_name, "w")
    cc_22_dense_of  = open(cc_p2_22_dense_of_name, "w")
    cc_12_dense_of  = open(cc_p1_12_dense_of_name, "w")
    cc_21_dense_of  = open(cc_p2_21_dense_of_name, "w")
    cc_11_dilute_of = open(cc_p1_11_dilute_of_name, "w")
    cc_22_dilute_of = open(cc_p2_22_dilute_of_name, "w")
    cc_12_dilute_of = open(cc_p1_12_dilute_of_name, "w")
    cc_21_dilute_of = open(cc_p2_21_dilute_of_name, "w")

    for i_frame in 1:num_frames
        for v in cc_pro1_11_dense[:, i_frame]
            @printf(cc_11_dense_of,  " %8.3f ", v)
        end
        for v in cc_pro1_12_dense[:, i_frame]
            @printf(cc_12_dense_of,  " %8.3f ", v)
        end
        for v in cc_pro2_22_dense[:, i_frame]
            @printf(cc_22_dense_of,  " %8.3f ", v)
        end
        for v in cc_pro2_21_dense[:, i_frame]
            @printf(cc_21_dense_of,  " %8.3f ", v)
        end
        for v in cc_pro1_11_dilute[:, i_frame]
            @printf(cc_11_dilute_of,  " %8.3f ", v)
        end
        for v in cc_pro1_12_dilute[:, i_frame]
            @printf(cc_12_dilute_of,  " %8.3f ", v)
        end
        for v in cc_pro2_22_dilute[:, i_frame]
            @printf(cc_22_dilute_of,  " %8.3f ", v)
        end
        for v in cc_pro2_21_dilute[:, i_frame]
            @printf(cc_21_dilute_of,  " %8.3f ", v)
        end
        @printf(cc_11_dense_of,  "\n")
        @printf(cc_12_dense_of,  "\n")
        @printf(cc_21_dense_of,  "\n")
        @printf(cc_22_dense_of,  "\n")
        @printf(cc_11_dilute_of, "\n")
        @printf(cc_12_dilute_of, "\n")
        @printf(cc_21_dilute_of, "\n")
        @printf(cc_22_dilute_of, "\n")
    end

end

if abspath(PROGRAM_FILE) == @__FILE__
    @sync @distributed for ic in 1:9
        for i_run in 1:5
            calculate_contact_map(290, ic * 10, 2, i_run)
        end
    end
end
