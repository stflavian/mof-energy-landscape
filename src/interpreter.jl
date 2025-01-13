# TODO: The functions are too big and are very messy. They should be split into smaller
# functions. The print statements should also be removed or contained, since now it makes
# reading very difficult.

using Plots

include("constants.jl")
include("force_fields.jl")
include("probe.jl")
include("framework.jl")
include("writer.jl")

"""
    compute_potential_landscape(atom_properties::Dict{SubString{String}, AtomProperties}, 
    framework::Framework, probe::Probe, sizea::Integer, sizeb::Integer, sizec::Integer,
    cutoff::Real, save::SubString{String})

Compute the potential landscape inside the framework for the given probe molecule. 

The unit cell is divided in `sizea` x `sizeb` x `sizec` smaller boxes, at the center of
which the probe is inserted. The Lennard-Jones and Coloumb interactions between the 
probe and the atoms in the framework are then computed and added to the potential value
in that position. If the distance between the probe and one of the atoms in the 
framework is closer than `0.5 sigma`, the potential value is set to 1. If the potential
value in a box is positive, it is set to 0. If the molecule contains more than one atom
this procedure is repeated for different rotations of the molecule, and the results
are averaged. When taking into account the probe-framework interactions, periodic
boundary conditions are used. 

The results are stored in a 4-dimensional array, where the first 3 entries represent 
the x, y, ans z coordinates, and the last entry is the value of the potential

# Arguments
- `atom_properties::Dict{SubString{String}, AtomProperties}`: the dictionary containing
the properties of each atomic species.
- `framework::Framework`: the data structure containing the framework unit cell 
parameters and the constituting atoms.
- `probe::Probe`: the vector containing the constituting atoms of the probe.
- `sizea::Integer`: the number of units in which the a lattice vector is divided. 
- `sizeb::Integer`: the number of units in which the b lattice vector is divided.
- `sizec::Integer`: the number of units in which the c lattice vector is divided.
- `cutoff::Real`: the potential cutoff used for the energy calculations.
- `rotations::Integer: the number of rotations used for the molecule.
- `output_file::IO: the file to which the results are written real-time.
- `save::SubString{String}`: "yes" if the potential in each point should be plotted and
saved.
"""
function compute_potential_landscape(atom_properties::Dict{SubString{String}, AtomProperties}, 
    framework::Framework, probe::Probe, sizea::Int64, sizeb::Int64, sizec::Int64,
    cutoff::Float64, rotations::Int64, output_file::IO, save::SubString{String})
    
    write_section(output_file, "Energy landscape calculations")
    A = compute_conversion_matrix(framework)
   
    # Compute the offset displacements needed for periodic boundary conditions
    write_subsection(output_file, "Periodic boundary conditions")
    
    box_string = rpad("Box [i, j, k]", 30, " ")
    x_string = rpad("X [angstrom]", 16, " ")
    y_string = rpad("Y [angstrom]", 16, " ")
    z_string = rpad("Z [angstrom]", 16, " ")
    
    pbc_framework = generate_pbc(framework)
    
    write(output_file, "$box_string $x_string $y_string $z_string\n")
    pbc_offsets = []
    for i in [-1, 0, 1], j in [-1, 0, 1], k in [-1, 0, 1]
        offset = A * [i, j, k]
        push!(pbc_offsets, offset)
        
        i_string = rpad("$i", 3, " ")
        j_string = rpad("$j", 3, " ")
        k_string = rpad("$k", 24, " ")
        x_string = rpad("$(round(offset[1], digits=8))", 16, " ")
        y_string = rpad("$(round(offset[2], digits=8))", 16, " ")
        z_string = rpad("$(round(offset[3], digits=8))", 16, " ")
        write(output_file, "$i_string$j_string$k_string $x_string $y_string $z_string\n")
    end
    write(output_file, "\n")
    
    # Initialize arrays and assign parameters for probe
    sx = range(start=0, step=1/sizea, length=sizea) .+ 1/(sizea * 2)
    sy = range(start=0, step=1/sizeb, length=sizeb) .+ 1/(sizeb * 2)
    sz = range(start=0, step=1/sizeb, length=sizec) .+ 1/(sizec * 2)
    potential = zeros(sizea, sizeb, sizec)
     
    # Main loop
    total_stats = @timed for (index, _) in enumerate(1:1:rotations)
        
        write_subsection(output_file, "Run $index")
       
        angle = rand() * 2 * pi
        axis_i = 1 - 2 * rand()
        axis_j = 1 - 2 * rand()
        axis_k = 1 - 2 * rand()
        rotated_probe = rotate_probe(probe, angle, axis_i, axis_j, axis_k)
        
        write_subsection(output_file, "Rotated probe")
        write_xyz(output_file, rotated_probe.atoms)

        run_stats = @timed for pos in eachindex(IndexCartesian(), potential)
            
            # Convert fractional coordinates to Cartesian coordinates
            x, y, z = fractional_to_cartesian(A, sx[pos[1]], sy[pos[2]], sz[pos[3]])

            for pbc_atom in pbc_framework.atoms
                
                # Store framework atom properties
                sig2 = atom_properties[pbc_atom.species].sigma
                eps2 = atom_properties[pbc_atom.species].epsilon
                
                # Compute distance between probe center of mass and framework atom
                cmr = compute_distance(x - pbc_atom.x, y - pbc_atom.y, z - pbc_atom.z)
                
                # Check if the probe center of mass falls within the atomic radius
                # of a framework atom 
                if cmr <= 0.5 * sig2 
                    @inbounds potential[pos] = 1.0
                    @goto next_point
                end
                
                for probe_atom in rotated_probe.atoms
                    
                    # Store probe atom properties
                    sig1 = atom_properties[probe_atom.species].sigma
                    eps1 = atom_properties[probe_atom.species].epsilon
                    
                    # Lorentz-Berthelot mixing rules and charge product
                    sig = (sig1 + sig2) * 0.5
                    eps = sqrt(eps1 * eps2)

                    # Compute distance between probe atom and framework atom
                    dx = probe_atom.x + x - pbc_atom.x
                    dy = probe_atom.y + y - pbc_atom.y
                    dz = probe_atom.z + z - pbc_atom.z
                    r = compute_distance(dx, dy, dz)
                    
                    # Check if distance is longer than cutoff
                    if r >= cutoff
                        continue
                    elseif r < cutoff
                        @inbounds potential[pos] += lennard_jones_energy(sig, eps, r) -
                                                 lennard_jones_energy(sig, eps, cutoff)
                    end
                end
            end
            @label next_point
        end
        
        write_subsection(output_file, "Run results")

        total_boxes = sizea * sizeb * sizec
        inaccessible_boxes = length(potential[potential .== 1])
        
        positive_boxes = length(potential[potential .> 0])
        total_pos_pot = sum(potential[potential .> 0]) * NA * 1e-3 / index

        negative_boxes = length(potential[potential .< 0])
        total_neg_pot = sum(potential[potential .< 0]) * NA * 1e-3 / index
        
        write_result(output_file, "Total boxes evaluated [count]", total_boxes)
        write_result(output_file, "Inaccessible box ratio [-]", inaccessible_boxes/total_boxes)
        write_result(output_file, "Negative potential box ratio [-]", negative_boxes/total_boxes)
        write_result(output_file, "Positive potential box ratio [-]", positive_boxes/total_boxes)
        write_result(output_file, "Average potential [kJ/mol]", (total_pos_pot + total_neg_pot) / total_boxes)
        write_result(output_file, "Average positive potential [kJ/mol]", total_pos_pot/positive_boxes)
        write_result(output_file, "Average negative potential [kJ/mol]", total_neg_pot/negative_boxes)
        println(output_file, " ")

        write_subsection(output_file, "Run performance statistics")
        write_result(output_file, "Time elapsed [s]", run_stats.time) 
        write_result(output_file, "GC time elapsed [s]", run_stats.gctime)
        write_result(output_file, "Memory allocated [bytes]", run_stats.bytes)
        println(output_file, " ")
    end
    
    write_section(output_file, "Total results")
    
    total_boxes = sizea * sizeb * sizec
    inaccessible_boxes = length(potential[potential .== 1])
        
    positive_boxes = length(potential[potential .> 0])
    total_pos_pot = sum(potential[potential .> 0]) * NA * 1e-3 / rotations

    negative_boxes = length(potential[potential .< 0])
    total_neg_pot = sum(potential[potential .< 0]) * NA * 1e-3 / rotations
    
    write_result(output_file, "Total boxes evaluated [count]", total_boxes)
    write_result(output_file, "Inaccessible box ratio [-]", inaccessible_boxes/total_boxes)
    write_result(output_file, "Negative potential box ratio [-]", negative_boxes/total_boxes)
    write_result(output_file, "Positive potential box ratio [-]", positive_boxes/total_boxes)
    write_result(output_file, "Average potential [kJ/mol]", (total_pos_pot + total_neg_pot) / total_boxes)
    write_result(output_file, "Average positive potential [kJ/mol]", total_pos_pot/positive_boxes)
    write_result(output_file, "Average negative potential [kJ/mol]", total_neg_pot/negative_boxes)
    println(output_file, " ")

    write_subsection(output_file, "Run performance statistics")
    write_result(output_file, "Time elapsed [s]", total_stats.time) 
    write_result(output_file, "GC time elapsed [s]", total_stats.gctime)
    write_result(output_file, "Memory allocated [bytes]", total_stats.bytes)
    println(output_file, " ")
    
    # if save == "yes"
    #     
    #     mkpath("Output")
    #     
    #     num = sizea * sizeb * sizec
    #     x = zeros(num)
    #     y = zeros(num)
    #     z = zeros(num)
    #     pot = zeros(num)
    #     for i in 1:1:sizea, j in 1:1:sizeb, k in 1:1:sizec
    #         index = i + (j-1) * sizea + (k-1) * sizeb^2
    #         x[index] = potential[i, j, 1, 1]
    #         y[index] = potential[i, j, 1, 2]
    #         z[index] = potential[i, j, 1, 3]
    #         pot[index] = potential[i, j, 1, 4]
    #     end

    #     p = scatter(x, y, marker_z=pot, markersize=2, camera=(0, -90),
    #     showaxis=false, legend=false, colorbar=true, markerstrokewidth=0, 
    #     right_margin=12Plots.mm)
    #     savefig(p, "Output/potential_landscape.png")
    # end
    return potential
end


"""
    compute_characteristic(atom_properties::Dict{SubString{String}, AtomProperties}, 
    framework::Framework, potential::Array{Float64, 4}, sizea::Integer, sizeb::Integer,
    sizec::Integer, npoints::Integer, save::SubString{String})

Compute the characteristic curve from the potential landscape.

Create a range between the minimum potential value recorded and 0 containing `npoints`.
For each energy value in the range, count the number of boxes that have a lower 
potential energy value. Multiply the number of boxes by the volume of a box to obtain
the total volume enclosed by each equipotential line.

# Arguments
- `atom_properties::Dict{SubString{String}, AtomProperties}`: the dictionary containing
the properties of each atomic species.
- `framework::Framework`: the data structure containing the framework unit cell 
parameters and the constituting atoms.
- `potential::Array{Float64, 4}`: the 4-dimensional array containing the value of the
potential at different points in the framework.
- `sizea::Integer`: the number of units in which the a lattice vector is divided. 
- `sizeb::Integer`: the number of units in which the b lattice vector is divided.
- `sizec::Integer`: the number of units in which the c lattice vector is divided.
- `npoints::Integer`: the number of points used for the characteristic curve.
- `save::SubString{String}`: "yes" if the characteristic curve should be plotted and
saved.
"""
function compute_characteristic(atom_properties::Dict{SubString{String}, AtomProperties}, 
    framework::Framework, potential::Array{Float64, 3}, sizea::Int64, sizeb::Int64,
    sizec::Int64, npoints::Int64, save::SubString{String}) 

    framework_mass = compute_framework_mass(atom_properties, framework)
    framework_volume = compute_framework_volume(framework)

    # Compute the volume of a sample point in ml
    sample_volume = framework_volume * 1e-24 / sizea / sizeb / sizec
    
    minimum_potential = minimum(potential)
    potential_range = range(start=minimum_potential, stop=-1e-7, length=npoints)
    
    mkpath("Output")
    output_file = open("Output/characteristic.dat", "w+")
    write(output_file, "# Potential [kJ/mol] \t Volume [ml/g] \n")

    measured_potential = zeros(npoints)
    measured_volumes = zeros(npoints)
    for (index, ads_potential) in enumerate(potential_range)
        
        counter = length(potential[potential .<= ads_potential])

        # Store the volume in ml/g
        volume = counter * sample_volume / framework_mass / 1000
        measured_volumes[index] = volume

        # Store the positive value of potential in kJ/mol
        measured_potential[index] = -ads_potential

        println(output_file, -ads_potential, "\t", volume) 
    end

    close(output_file)
    
    if save == "yes"
        characteristic_plot = plot(measured_potential, measured_volumes)
        xlabel!(characteristic_plot, "Potential [kJ/mol]")
        ylabel!(characteristic_plot, "Volume [ml/g]")
        savefig(characteristic_plot, "Output/characteristic.png") 
    end
end
