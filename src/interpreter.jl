using Plots

include("constants.jl")
include("force_fields.jl")
include("probe.jl")
include("framework.jl")


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
    
    message = rpad("==== Energy landscape calculations ", 81, "=")
    write(output_file, "$message\n")
    write(output_file, "\n")

    A = compute_conversion_matrix(framework)
   
    # Compute the offset displacements needed for periodic boundary conditions
    message = rpad("---- Neighbouring unit cells ", 81, "-")
    write(output_file, "$message\n")
    
    box_string = rpad("Box [i, j, k]", 30, " ")
    x_string = rpad("X [angstrom]", 16, " ")
    y_string = rpad("Y [angstrom]", 16, " ")
    z_string = rpad("Z [angstrom]", 16, " ")

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
    potential = zeros(sizea, sizeb, sizec, 4)
    
    # Main loop
    total_stats = @timed for (index, _) in enumerate(1:1:rotations)

        message = rpad("---- Run $index ", 81, "-")
        write(output_file, "$message\n")
        write(output_file, "\n")
        
        angle = rand() * 2 * pi
        axis_i = 1 - 2 * rand()
        axis_j = 1 - 2 * rand()
        axis_k = 1 - 2 * rand()
        rotated_probe = rotate_probe(probe, angle, axis_i, axis_j, axis_k)
    
        message = rpad("---- Rotated probe ", 81, "-")
        write(output_file, "$message\n")

        index_string = rpad("Index", 10, " ")
        species_string = rpad("Species", 19, " ")
        x_string = rpad("X [angstrom]", 16, " ")
        y_string = rpad("Y [angstrom]", 16, " ")
        z_string = rpad("Z [angstrom]", 16, " ")

        write(output_file, "$index_string $species_string $x_string $y_string $z_string\n")
        for (index, atom) in enumerate(rotated_probe.atoms)

            index_string = rpad(index, 10, " ")
            species_string = rpad(atom.species, 19, " ")
            x_string = rpad(atom.x, 16, " ")
            y_string = rpad(atom.y, 16, " ")
            z_string = rpad(atom.z, 16, " ")

            write(output_file, "$index_string $species_string $x_string $y_string $z_string\n")
        end
        write(output_file, "\n")
        flush(output_file)

        run_stats = @timed for (i, fa) in enumerate(sx), (j, fb) in enumerate(sy), (k, fc) in enumerate(sz)
            
            # Convert fractional coordinates to Cartesian coordinates
            x, y, z = fractional_to_cartesian(A, fa, fb, fc)
            potential[i, j, k, 1] = x
            potential[i, j, k, 2] = y
            potential[i, j, k, 3] = z

            for framework_atom in framework.atoms
                
                # Store framework atom properties
                sig2 = atom_properties[framework_atom.species].sigma
                eps2 = atom_properties[framework_atom.species].epsilon
                #q2 = atom_properties[framework_atom.species].charge
                
                # Compute distance between probe center of mass and framework atom
                dx = x - framework_atom.x
                dy = y - framework_atom.y
                dz = z - framework_atom.z
                cmr = compute_distance(dx, dy, dz)
                
                # Check if the probe center of mass falls within the atomic radius
                # of a framework atom 
                if cmr < 0.5 * sig2 
                    potential[i, j, k, 4] = 1
                    @goto next_point
                end
                
                for probe_atom in rotated_probe.atoms
                    
                    # Store probe atom properties
                    sig1 = atom_properties[probe_atom.species].sigma
                    eps1 = atom_properties[probe_atom.species].epsilon
                    #q1 = atom_properties[probe_atom.species].charge
                    
                    # Lorentz-Berthelot mixing rules and charge product
                    sig = (sig1 + sig2) * 0.5
                    eps = sqrt(eps1 * eps2)
                    #q = q1 * q2 

                    for offset in pbc_offsets
                        
                        # Compute distance between probe atom and framework atom
                        dx = probe_atom.x + x - framework_atom.x - offset[1]
                        dy = probe_atom.y + y - framework_atom.y - offset[2]
                        dz = probe_atom.z + z - framework_atom.z - offset[3]
                        r = compute_distance(dx, dy, dz)
                        
                        # Check if distance is longer than cutoff
                        if r < cutoff
                            potential[i, j, k, 4] += lennard_jones_energy(sig, eps, r)
                            potential[i, j, k, 4] -= lennard_jones_energy(sig, eps, cutoff)
                            #potential[i, j, k, 4] += coloumb_energy(q, r)
                        elseif r > cutoff
                            continue
                        end

                    end
                end
            end
            @label next_point
        end
        
        message = rpad("---- Run results ", 81, "-")
        write(output_file, "$message\n")

        total_boxes = sizea * sizeb * sizec
        inaccessible_boxes = 0

        positive_potential_boxes = 0
        total_positive_potential = 0

        negative_potential_boxes = 0
        total_negative_potential = 0

        for i in 1:1:sizea, j in 1:1:sizeb, k in 1:1:sizec

            if potential[i, j, k, 4] == 1
                inaccessible_boxes += 1
                continue
            elseif potential[i, j, k, 4] > 0
                positive_potential_boxes += 1
                total_positive_potential += potential[i, j, k, 4] * NA * 1e-3 / index
            else
                negative_potential_boxes += 1
                total_negative_potential += potential[i, j, k, 4] * NA * 1e-3 / index
            end
        end

        key_string = rpad("Total boxes evaluated [count]", 40, " ")
        argument_string = lpad(total_boxes, 40, " ")
        write(output_file, "$key_string $argument_string\n")

        key_string = rpad("Inaccessible box ratio [-]", 40, " ")
        argument_string = lpad(inaccessible_boxes/total_boxes, 40, " ")
        write(output_file, "$key_string $argument_string\n")

        key_string = rpad("Negative potential box ratio [-]", 40, " ")
        argument_string = lpad(negative_potential_boxes/total_boxes, 40, " ")
        write(output_file, "$key_string $argument_string\n")

        key_string = rpad("Positive potential box ratio [-]", 40, " ")
        argument_string = lpad(positive_potential_boxes/total_boxes, 40, " ")
        write(output_file, "$key_string $argument_string\n")

        key_string = rpad("Average potential [kJ/mol]", 40, " ")
        val = (total_positive_potential + total_negative_potential) / total_boxes 
        argument_string = lpad(val, 40, " ")
        write(output_file, "$key_string $argument_string\n")

        key_string = rpad("Average positive potential [kJ/mol]", 40, " ")
        argument_string = lpad(total_positive_potential/positive_potential_boxes, 40, " ")
        write(output_file, "$key_string $argument_string\n")

        key_string = rpad("Average negative potential [kJ/mol]", 40, " ")
        argument_string = lpad(total_negative_potential/negative_potential_boxes, 40, " ")
        write(output_file, "$key_string $argument_string\n")

        write(output_file, "\n")
        flush(output_file)
        

        message = rpad("---- Run performance statistics ", 81, "-")
        write(output_file, "$message\n")
        
        key_string = rpad("Time elapsed [s]", 40, " ")
        argument_string = lpad(run_stats.time, 40, " ")
        write(output_file, "$key_string $argument_string\n")

        key_string = rpad("GC time elapsed [s]", 40, " ")
        argument_string = lpad(run_stats.gctime, 40, " ")
        write(output_file, "$key_string $argument_string\n")
        
        key_string = rpad("Memory allocated [bytes]", 40, " ")
        argument_string = lpad(run_stats.bytes, 40, " ")
        write(output_file, "$key_string $argument_string\n")
        
        write(output_file, "\n")
        flush(output_file)
    end
    
    message = rpad("---- Total results ", 81, "-")
    write(output_file, "$message\n")
    
    total_boxes = sizea * sizeb * sizec
    inaccessible_boxes = 0
    
    positive_potential_boxes = 0
    total_positive_potential = 0

    negative_potential_boxes = 0
    total_negative_potential = 0

    for i in 1:1:sizea, j in 1:1:sizeb, k in 1:1:sizec

        if potential[i, j, k, 4] == 1
            inaccessible_boxes += 1
            continue
        elseif potential[i, j, k, 4] > 0
            positive_potential_boxes += 1
            
            # Conversion from J to kJ/mol
            potential[i, j, k, 4] *= NA * 1e-3 / rotations
            total_positive_potential += potential[i, j, k, 4]
            
            # Set potential to 0 for plotting 
            potential[i, j, k, 4] = 0
        else
            negative_potential_boxes += 1
            
            # Conversion from J to kJ/mol
            potential[i, j, k, 4] *= NA * 1e-3 / rotations
            total_negative_potential += potential[i, j, k, 4]
        end
    end
    
    key_string = rpad("Total boxes evaluated [count]", 40, " ")
    argument_string = lpad(total_boxes, 40, " ")
    write(output_file, "$key_string $argument_string\n")
    
    key_string = rpad("Inaccessible box ratio [-]", 40, " ")
    argument_string = lpad(inaccessible_boxes/total_boxes, 40, " ")
    write(output_file, "$key_string $argument_string\n")
    
    key_string = rpad("Negative potential box ratio [-]", 40, " ")
    argument_string = lpad(negative_potential_boxes/total_boxes, 40, " ")
    write(output_file, "$key_string $argument_string\n")
    
    key_string = rpad("Positive potential box ratio [-]", 40, " ")
    argument_string = lpad(positive_potential_boxes/total_boxes, 40, " ")
    write(output_file, "$key_string $argument_string\n")
    
    key_string = rpad("Average potential [kJ/mol]", 40, " ")
    val = (total_positive_potential + total_negative_potential) / total_boxes 
    argument_string = lpad(val, 40, " ")
    write(output_file, "$key_string $argument_string\n")
    
    key_string = rpad("Average positive potential [kJ/mol]", 40, " ")
    argument_string = lpad(total_positive_potential/positive_potential_boxes, 40, " ")
    write(output_file, "$key_string $argument_string\n")
    
    key_string = rpad("Average negative potential [kJ/mol]", 40, " ")
    argument_string = lpad(total_negative_potential/negative_potential_boxes, 40, " ")
    write(output_file, "$key_string $argument_string\n")
    
    write(output_file, "\n")
    flush(output_file)
        

    message = rpad("---- Total performance statistics ", 81, "-")
    write(output_file, "$message\n")
        
    key_string = rpad("Time elapsed [s]", 40, " ")
    argument_string = lpad(total_stats.time, 40, " ")
    write(output_file, "$key_string $argument_string\n")

    key_string = rpad("GC time elapsed [s]", 40, " ")
    argument_string = lpad(total_stats.gctime, 40, " ")
    write(output_file, "$key_string $argument_string\n")
        
    key_string = rpad("Memory allocated [bytes]", 40, " ")
    argument_string = lpad(total_stats.bytes, 40, " ")
    write(output_file, "$key_string $argument_string\n")
        
    write(output_file, "\n")
    flush(output_file)

    if save == "yes"
        
        mkpath("Output")
        
        num = sizea * sizeb * sizec
        x = zeros(num)
        y = zeros(num)
        z = zeros(num)
        pot = zeros(num)
        for i in 1:1:sizea, j in 1:1:sizeb, k in 1:1:sizec
            index = i + (j-1) * sizea + (k-1) * sizeb^2
            x[index] = potential[i, j, 1, 1]
            y[index] = potential[i, j, 1, 2]
            z[index] = potential[i, j, 1, 3]
            pot[index] = potential[i, j, 1, 4]
        end

        p = scatter(x, y, marker_z=pot, markersize=2, camera=(0, -90),
        showaxis=false, legend=false, colorbar=true, markerstrokewidth=0, 
        right_margin=12Plots.mm)
        savefig(p, "Output/potential_landscape.png")
    end
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
    framework::Framework, potential::Array{Float64, 4}, sizea::Int64, sizeb::Int64,
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
        counter = 0
        for i in 1:1:sizea, j in 1:1:sizeb, k in 1:1:sizec
            if potential[i, j, k, 4] <= ads_potential
                counter += 1
            end
        end

        # Store the volume in ml/g
        volume = counter * sample_volume / framework_mass / 1000
        measured_volumes[index] = volume

        # Store the positive value of potential in kJ/mol
        measured_potential[index] = -ads_potential

        write(output_file, "$(-ads_potential) \t $volume \n") 
    end

    close(output_file)
    
    if save == "yes"
        characteristic_plot = plot(measured_potential, measured_volumes)
        xlabel!(characteristic_plot, "Potential [kJ/mol]")
        ylabel!(characteristic_plot, "Volume [ml/g]")
        savefig(characteristic_plot, "Output/characteristic.png") 
    end
end
