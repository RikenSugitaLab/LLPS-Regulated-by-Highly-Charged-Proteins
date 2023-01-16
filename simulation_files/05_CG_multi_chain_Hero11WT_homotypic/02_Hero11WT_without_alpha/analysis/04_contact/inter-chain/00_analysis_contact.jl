#!/usr/bin/env julia

include("/home/ctan/Workspace/genesis_CG_julia/src/lib/gcj.jl")

using Printf

function calculate_contact_map(i_temperature, i_run)

    # ----------
    # parameters
    # ----------
    contact_r0 = 15.0
    contact_cutoff = 25.0

    name_pro1 = "hero11_alpha"

    NUM_ATOM_PRO1 = 99

    num_pro1 = 100

    # EDIT THIS: num-of-frames in traj and num-of-frames to calculate cm
    n_traj_frames = 5000
    i_begin_frame = div(n_traj_frames, 2)
    i_end_frame = n_traj_frames
    n_structures = 20
    n_interval = div(i_end_frame - i_begin_frame, n_structures)

    # EDIT THIS: trajectory name
    traj_name = @sprintf("../../../%s_%d_md1_T_%02d_%02d.dcd", name_pro1, num_pro1, i_temperature, i_run)

    # ---------------
    # data structures
    # ---------------
    data_pro1_all = zeros(Float64, (NUM_ATOM_PRO1, NUM_ATOM_PRO1))

    println("Analyzing trajectory: ", traj_name )
    args = Dict("verbose"=>true,
                "begin" => i_begin_frame,
                "end" => i_end_frame,
                "step" => n_interval)
    mytraj = read_dcd(traj_name, args)


    n_frame_calc = 0
    for i_frame in 1:n_structures
        println(" > Frame: ", i_frame)
        myconf = mytraj.conformations[i_frame]

        n_frame_calc += 1

        # Get box size
        box_size_x = mytraj.boundary_box_size[1, i_frame]
        box_size_y = mytraj.boundary_box_size[3, i_frame]
        box_size_z = mytraj.boundary_box_size[6, i_frame]
        println(" > Simulation box size: ", box_size_x, " x ", box_size_y, " x ", box_size_z)
        box_size = [box_size_x, box_size_y, box_size_z]

        # -------------
        # cell division
        # -------------
        num_cell_x = Int(round(box_size_x / contact_cutoff))
        num_cell_y = Int(round(box_size_y / contact_cutoff))
        num_cell_z = Int(round(box_size_z / contact_cutoff))
        num_cell_all = num_cell_x * num_cell_y * num_cell_z
        println(" > Cells: ", num_cell_x, " x ", num_cell_y, " x ", num_cell_z, " = ", num_cell_all)
        cell_size_x = box_size_x / num_cell_x
        cell_size_y = box_size_y / num_cell_y
        cell_size_z = box_size_z / num_cell_z
        cell_sizes = [cell_size_x, cell_size_y, cell_size_z]

        # -----------------------
        # put particles into cels
        # -----------------------
        println(" > Constructing particle cell list: ")
        cell_particle_list = [[] for j in 1:num_cell_all]
        for j_atom in 1:num_pro1 * NUM_ATOM_PRO1
            coor_j = myconf.coors[:, j_atom]
            cell_vec = Int.(floor.(mod.(coor_j, box_size) ./ cell_sizes))
            cell_idx = 1 + cell_vec[1] + cell_vec[2] * num_cell_x + cell_vec[3] * num_cell_x * num_cell_y
            push!(cell_particle_list[cell_idx], j_atom)
        end

        # ======================================
        # loop over cells and calculate contacts
        # ======================================
        println(" > Calculating contact map: ")
        Threads.@threads for i_cell in 1:num_cell_all
            i_z, i_xy = divrem(i_cell - 1, num_cell_x * num_cell_y)
            i_y, i_x  = divrem(i_xy,   num_cell_x)
            for dx = -1:1
                j_x = mod(i_x + dx, num_cell_x)
                for dy = -1:1
                    j_y = mod(i_y + dy, num_cell_y)
                    for dz = -1:1
                        j_z = mod(i_z + dz, num_cell_z)
                        j_cell = 1 + j_x + j_y * num_cell_x + j_z * num_cell_x * num_cell_y
                        for i_atom in cell_particle_list[i_cell]
                            chain_i = div(i_atom - 1, NUM_ATOM_PRO1)
                            for j_atom in cell_particle_list[j_cell]
                                chain_j = div(j_atom - 1, NUM_ATOM_PRO1)
                                if chain_j <= chain_i
                                    continue
                                end
                                coor_i = myconf.coors[:, i_atom]
                                coor_j = myconf.coors[:, j_atom]
                                d_r = coor_i - coor_j
                                d_r_in_cell = d_r - round.(d_r ./ box_size) .* box_size
                                dist = sqrt(d_r_in_cell' * d_r_in_cell)
                                if dist > contact_cutoff
                                    continue
                                end
                                cnt = 1 / (1 + exp(dist - contact_r0))
                                idx_i = mod(i_atom - 1, NUM_ATOM_PRO1) + 1
                                idx_j = mod(j_atom - 1, NUM_ATOM_PRO1) + 1
                                data_pro1_all[idx_i, idx_j] += cnt
                                data_pro1_all[idx_j, idx_i] += cnt
                            end # j_atom
                        end     # i_atom
                    end         # dz
                end             # dy
            end                 # dx
        end                     # i_cell

    end

    data_pro1_all_ave = data_pro1_all ./ (n_frame_calc * num_pro1)

    of_name = @sprintf("%s_contacts_t%02d_%02d.dat", name_pro1, i_temperature, i_run)
    offile = open(of_name, "w")
    for idx_i in 1:NUM_ATOM_PRO1
        for idx_j in 1:NUM_ATOM_PRO1
            @printf(offile, " %8.3f ", data_pro1_all_ave[idx_i, idx_j])
        end
        @printf(offile, " \n")
    end
    close(offile)

end


if abspath(PROGRAM_FILE) == @__FILE__
    for iT in 1:1
        for irun in 1:5
            calculate_contact_map(iT, irun)
        end
    end

end
