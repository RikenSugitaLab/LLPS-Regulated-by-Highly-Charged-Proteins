#!/usr/bin/env julia

include("/home/ctan/Workspace/genesis_CG_julia/src/lib/gcj.jl")

using Printf
using DelimitedFiles

function calculate_position(temperature, n_hero, i_md, i_run)

    # ----------
    # parameters
    # ----------
    name_pro1 = "tdp43"
    name_pro2 = "hero11"

    NUM_ATOM_PRO1 = 154
    NUM_ATOM_PRO2 = 99

    num_pro1 = 100
    num_pro2 = n_hero

    box_size_z = 3000
    num_bin = 100
    box_size_bin = box_size_z / num_bin

    # -------------
    # load top file
    # -------------
    top_name = @sprintf("../../N_%2d_T_%3d_r_%02d.top", num_pro2, temperature, i_run)
    my_top = read_grotop(top_name)

    # -------------
    # load dcd file
    # -------------
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
    pos_pro1_list_all = zeros(Int, (num_pro1, n_structures))
    pos_pro2_list_all = zeros(Int, (num_pro2, n_structures))

    println("Analyzing trajectory: ", traj_name )
    args = Dict("verbose"=>false,
                "begin" => i_begin_frame,
                "end" => i_end_frame,
                "step" => n_interval)
    # mytraj = read_dcd(traj_name, args)
    mytraj = read_dcd(traj_name)

    # ------------
    # density data
    # ------------
    pro1_density_threshold = 0.0034
    pro1_density_data_fname = @sprintf("../01_density_in_box/md%d/%s_T%d_md%d_r%02d.%s.density.dat", i_md, system_name, temperature, i_md, i_run, name_pro1)
    pro1_density_data = readdlm(pro1_density_data_fname)

    Threads.@threads for i_frame in 1:length(mytraj.conformations)
        my_conf = mytraj.conformations[i_frame]

        n_atom_count = 0
        for i_chain in 1:num_pro1
            com = compute_center_of_mass([n_atom_count+1:n_atom_count+NUM_ATOM_PRO1...], my_top, my_conf)
            z_int = mod( Int(floor(mod(com[3], box_size_z) / box_size_bin)), num_bin)
            if pro1_density_data[i_frame, z_int + 1] > pro1_density_threshold
                pos_pro1_list_all[i_chain, i_frame] = 1
            else
                pos_pro1_list_all[i_chain, i_frame] = 0
            end
            n_atom_count += NUM_ATOM_PRO1
        end

        for i_chain in 1:num_pro2
            com = compute_center_of_mass([n_atom_count+1:n_atom_count+NUM_ATOM_PRO2...], my_top, my_conf)
            z_int = mod( Int(floor(mod(com[3], box_size_z) / box_size_bin)), num_bin)
            if pro1_density_data[i_frame, z_int + 1] > pro1_density_threshold
                pos_pro2_list_all[i_chain, i_frame] = 1
            else
                pos_pro2_list_all[i_chain, i_frame] = 0
            end
            n_atom_count += NUM_ATOM_PRO2
        end
    end

    of1_name = @sprintf("%s_T%03d_md%d_r%02d.%s_pos.dat", system_name, temperature, i_md, i_run, name_pro1)
    of2_name = @sprintf("%s_T%03d_md%d_r%02d.%s_pos.dat", system_name, temperature, i_md, i_run, name_pro2)
    of1 = open(of1_name, "w")
    of2 = open(of2_name, "w")
    for i_frame in 1:n_structures
        for i_chain in 1:num_pro1
            @printf(of1, " %8.3f ", pos_pro1_list_all[i_chain, i_frame])
        end
        @printf(of1, " \n")
    end
    close(of1)
    for i_frame in 1:n_structures
        for i_chain in 1:num_pro2
            @printf(of2, " %8.3f ", pos_pro2_list_all[i_chain, i_frame])
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
                calculate_position(T, ic * 10, 1, irun)
                calculate_position(T, ic * 10, 2, irun)
                calculate_position(T, ic * 10, 3, irun)
            end
        end
    end
end
