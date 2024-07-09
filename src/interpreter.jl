using Plots
using ProgressBars
using LinearAlgebra

include("constants.jl")


"""
Structure for storing the properties of an atomic species.
"""
struct AtomProperties
    epsilon::Real
    sigma::Real
    charge::Real
    mass::Real
end


"""
Structure for storing the species and Cartesian coordinates of atoms. 
"""
struct Atom
    species::String
    x::Real
    y::Real
    z::Real
end


"""
Structure for storing the framework unit cell parameters and position of atoms.
"""
struct Framework
    a::Real
    b::Real
    c::Real
    alpha::Real
    beta::Real
    gamma::Real
    atoms::Vector{Atom}
end


"""
Structure for storing the position of atoms inside the probe.
"""
struct Probe
    atoms::Vector{Atom}
end


"""
    read_probe_file(path::SubString{String})

Parse the probe file and return a vector containing the constituting atoms. 

The probe file should be in .xyz format, with the comment line being ignored.

# Arguments
- `path::SubString{String}`: the path to the probe file.
"""
function read_probe_file(path::SubString{String})
    
    probe_file = readlines(path)

    local probe_atoms
    for (index, line) in enumerate(probe_file)
        
        if index == 1
            line = split(line)
            number_of_atoms = parse(Int64, line[1])
            probe_atoms = Vector{Atom}(undef, number_of_atoms)
        end
        
        if index > 2
            line = split(line)
            species = line[1]
            x = parse(Float64, line[2])
            y = parse(Float64, line[3])
            z = parse(Float64, line[4])

            probe_atoms[index-2] = Atom(species, x, y, z) 
        end
    end

    probe = Probe(probe_atoms)
    
    return probe
end


"""
    read_framework_file(path::SubString{String})

Parse the framework file and return a data structure containing the constituting atoms 
and the unit cell parameters. 

The framework file should be in .xyz format, with the comment line containing the 
lattice vectors a, b, and c and the lattice angles alpha, beta, and gamma, in this 
specific order.

# Arguments
- `path::SubString{String}`: the path to the framework file.
"""
function read_framework_file(path::SubString{String})
    
    framework_file = readlines(path)
    
    local a, b, c, alpha, beta, gamma, framework_atoms
    for (index, line) in enumerate(framework_file)
        
        line = split(line)
        if index == 1
            number_of_atoms = parse(Int64, line[1])
            framework_atoms = Vector{Atom}(undef, number_of_atoms)
        end

        if index == 2
            a = parse(Float64, line[1])
            b = parse(Float64, line[2])
            c = parse(Float64, line[3])
            alpha = parse(Float64, line[4])
            beta = parse(Float64, line[5])
            gamma = parse(Float64, line[6])
        end
        
        if index > 2
            species = line[1]
            x = parse(Float64, line[2])
            y = parse(Float64, line[3])
            z = parse(Float64, line[4])

            framework_atoms[index-2]= Atom(species, x, y, z)
        end
    end
    
    framework = Framework(a, b, c, alpha, beta, gamma, framework_atoms)

    return framework
end


"""
    read_properties_file(path::SubString{String})

Parse the properties file and return a dictionary containing the atomic species as keys
and the AtomProperties data structure as value. 

The file should contain separate lines for each atomic species, each having the name 
of the atom, followed by the epsilon value (K), sigma value (angstrom), charge (e-), 
and mass (m.u.), in this specific order. 

# Arguments
- `path::SubString{String}`: the path to the properties file.
"""
function read_properties_file(path::SubString{String})
    
    properties_file = readlines(path)
    atom_properties = Dict{SubString{String}, AtomProperties}()
    
    for line in properties_file
        
        if length(line) == 0
            continue
        end
        
        line = split(line)
        species = line[1]
        epsilon = parse(Float64, line[2])
        sigma = parse(Float64, line[3])
        charge = parse(Float64, line[4])
        mass = parse(Float64, line[5])
        
        atom_properties[species] = AtomProperties(epsilon, sigma, charge, mass)
    end

    return atom_properties
end


"""
    read_input_file(path::SubString{String})

Parse the input file and return a dictionary containing the simulation settings. 

The file should contain separate lines for each keyword, with the keyword and the 
argument separated by an arbitrary amount of spaces or tabs.

# Arguments
- `path::SubString{String}`: the path to the input file.
"""
function read_input_file(path::String)
    keywords = Dict(
        "FRAMEWORK" => Nothing, 
        "PROPERTIES" => Nothing,
        "CUTOFF" => "15",
        "PROBE" => Nothing,
        "XPOINTS" => "30", 
        "YPOINTS" => "30", 
        "ZPOINTS" => "30",
        "ROTATIONS" => "1",
        "SAVE_POTENTIAL" => "yes", 
        "CHARACTERISTIC_POINTS" => "30",
        "SAVE_CHARACTERISTIC" => "yes")

    input_file = readlines(path)

    for line in input_file
        if length(line) == 0
            continue
        end
        
        line = split(line)
        keyword = line[1]
        if keyword in keys(keywords)
            keywords[keyword] = line[2] 
        end
    end

    return keywords
end


"""
    lennard_jones_energy(sigma::Real, epsilon::Real, distance::Real)

Compute the 6-12 Lennard-Jones energy (J) between two particles, using the sigma value
(angstrom), epsilon value (K), and distance (angstrom).

# Arguments
- `sigma::Real`: the size parameter, usually the sum of the particles' radii.
- `epsilon::Real`: the depth of the potential well, usually the geometric mean of the
particles' epsilon values.
- `distance::Real`: the distance between the center of masses of the two particles.
"""
function lennard_jones_energy(sigma::Real, epsilon::Real, distance::Real)
    frac = sigma / distance    
    return 4 * epsilon * KB * (frac^12 - frac^6)
end


"""
    coloumb_energy(charge::Real, distance::Real)

Compute the electrostatic energy (J) between two particles using the product of the
charges (e-^2) and distance (angstrom).

# Arguments
- `charge::Real`: the product of the charges of the two particles.
- `distance::Real`: the distance between the center of masses of the two particles.
"""
function coloumb_energy(charge::Real, distance::Real)
    return Q^2 * charge * KE * 1e10 / distance
end


"""
    compute_probe_charge(atom_properties::Dict{SubString{String}, AtomProperties},
    probe::Probe)

Compute the total charge (e-) of the probe by taking into account the contribution of 
each constituting atom.

# Arguments
- `atom_properties::Dict{SubString{String}, AtomProperties}`: the dictionary containing
the properties of each atomic species.
- `probe::Probe`: the vector containing the constituting atoms of the probe.
"""
function compute_probe_charge(atom_properties::Dict{SubString{String}, AtomProperties},
    probe::Probe)
    
    charge = 0.0

    for atom in probe.atoms
        charge += atom_properties[atom.species].charge
    end

    return charge
end


"""
    compute_probe_mass(atom_properties::Dict{SubString{String}, AtomProperties},
    probe::Probe)

Compute the total mass (kg) of the probe by taking into account the contribution of 
each constituting atom.

# Arguments
- `atom_properties::Dict{SubString{String}, AtomProperties}`: the dictionary containing
the properties of each atomic species.
- `probe::Probe`: the vector containing the constituting atoms of the probe.
"""
function compute_probe_mass(atom_properties::Dict{SubString{String}, AtomProperties},
    probe::Probe)
    
    mass = 0.0

    for atom in probe.atoms
        mass += atom_properties[atom.species].mass
    end

    return mass * MC
end


"""
    remove_probe_centermass(atom_properties::Dict{SubString{String}, AtomProperties},
    probe::Probe)

Compute the position of the center of mass and subtract the value from the position
of each atom. 

Return a new Probe structure with the center of mass in the origin, together with the 
old position of the center of mass.

# Arguments
- `atom_properties::Dict{SubString{String}, AtomProperties}`: the dictionary containing
the properties of each atomic species.
- `probe::Probe`: the vector containing the constituting atoms of the probe.
"""
function remove_probe_centermass(atom_properties::Dict{SubString{String}, AtomProperties},
    probe::Probe)
    
    centered_probe = deepcopy(probe.atoms)

    probe_mass = compute_probe_mass(atom_properties, probe)
    center_of_mass_x = 0.0
    center_of_mass_y = 0.0
    center_of_mass_z = 0.0

    for atom in probe.atoms
        center_of_mass_x += atom.x * atom_properties[atom.species].mass
        center_of_mass_y += atom.y * atom_properties[atom.species].mass
        center_of_mass_z += atom.z * atom_properties[atom.species].mass
    end

    center_of_mass_x /= probe_mass
    center_of_mass_y /= probe_mass
    center_of_mass_z /= probe_mass
    
    for (index, atom) in enumerate(probe.atoms)
        species = atom.species
        x = atom.x - center_of_mass_x
        y = atom.y - center_of_mass_y
        z = atom.z - center_of_mass_z
        centered_probe[index] = Atom(species, x, y, z)
    end
    
    center_of_mass = [center_of_mass_x, center_of_mass_y, center_of_mass_z]

    return Probe(centered_probe), center_of_mass
end


"""
    rotate_probe(probe::Probe, angle::Real, i::Real, j::Real, k::Real)

Rotate the probe molecule by an arbitrary angle around an arbitrary axis. 

This method is based on quaternion rotation method, where a rotation angle and a 
rotation axis are selected, and then a rotation operator is applied on the each 
constituting atom to compute the new position after the rotation. The rotation axis 
passes through the origin, so the particle needs to have its center of mass in the 
origin.

# Arguments
- `probe::Probe`: the vector containing the constituting atoms of the probe.
- `angle::Real`: the angle of rotation for the molecule.
- `i::Real`: the i (x) component of the rotation axis.
- `j::Real`: the j (y) component of the rotation axis.
- `k::Real`: the k (z) component of the rotation axis.
"""
function rotate_probe(probe::Probe, angle::Real, i::Real, j::Real, k::Real)
    
    rotated_probe = deepcopy(probe.atoms)
    rotation_vector = sin(angle/2) / sqrt(i^2 + j^2 + k^2) * [i; j; k]

    q0 = cos(angle/2)
    qlen2 = rotation_vector[1]^2 + rotation_vector[2]^2 + rotation_vector[3]^2
    
    for (index, atom) in enumerate(probe.atoms)
        species = atom.species
        position_vector = [atom.x; atom.y; atom.z]
       
        rotated_position_vector = (q0^2 - qlen2) * position_vector + 
        2 * dot(rotation_vector, position_vector) * rotation_vector + 
        2 * q0 * cross(rotation_vector, position_vector)
        
        x, y, z = rotated_position_vector
        rotated_probe[index] = Atom(species, x, y, z)
    end
    
    return Probe(rotated_probe)
end


"""
    compute_framework_charge(atom_properties::Dict{SubString{String}, AtomProperties},
    framework::Framework)

Compute the total charge (e-) of the framework unit cell by taking into account the 
contribution of each constituting atom.

# Arguments
- `atom_properties::Dict{SubString{String}, AtomProperties}`: the dictionary containing
the properties of each atomic species.
- `framework::Framework`: the data structure containing the framework unit cell 
parameters and the constituting atoms.
"""
function compute_framework_charge(atom_properties::Dict{SubString{String}, AtomProperties},
    framework::Framework)
    
    charge = 0.0

    for atom in framework.atoms
        charge += atom_properties[atom.species].charge
    end

    return charge
end


"""
    compute_framework_mass(atom_properties::Dict{SubString{String}, AtomProperties},
    framework::Framework)

Compute the total mass (kg) of the framework unit cell by taking into account the 
contribution of each constituting atom.

# Arguments
- `atom_properties::Dict{SubString{String}, AtomProperties}`: the dictionary containing
the properties of each atomic species.
- `framework::Framework`: the data structure containing the framework unit cell 
parameters and the constituting atoms.
"""
function compute_framework_mass(atom_properties::Dict{SubString{String}, AtomProperties},
    framework::Framework)
    
    mass = 0.0

    for atom in framework.atoms
        mass += atom_properties[atom.species].mass
    end

    return mass * MC
end


"""
    compute_framework_volume(atom_properties::Dict{SubString{String}, AtomProperties},
    framework::Framework)

Compute the total volume (angstrom^3) of the framework unit cell from the unit cell
parameters.

# Arguments
- `framework::Framework`: the data structure containing the framework unit cell 
parameters and the constituting atoms.
"""
function compute_framework_volume(framework::Framework)
    
    a = framework.a
    b = framework.b
    c = framework.c

    alpha = framework.alpha * pi / 180
    beta = framework.beta * pi / 180
    gamma = framework.gamma * pi / 180

    volume = a * b * c * sqrt(sin(alpha)^2 + sin(beta)^2 + sin(gamma)^2 + 
            2 * cos(alpha) * cos(beta) * cos(gamma) - 2)

    return volume
end


"""
    compute_framework_density(framework_mass::Real, framework_volume::Real)

Compute the density (kg/m3) of the framework.

# Arguments
- `framework_mass::Real`: the total mass of the atoms inside the unit cell.
- `framework_volume::Real`: the total volume of the unit cell.
"""
function compute_framework_density(framework_mass::Real, framework_volume::Real)
    return framework_mass * 1e30 / framework_volume
end


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
    framework::Framework, probe::Probe, sizea::Integer, sizeb::Integer, sizec::Integer,
    cutoff::Real, rotations::Integer, output_file::IO, save::SubString{String})
    
    message = rpad("==== Energy landscape calculations ", 81, "=")
    write(output_file, "$message\n")
    write(output_file, "\n")
    
    # Compute the transformation matrix for fractional to Cartesian coordinates
    a = framework.a
    b = framework.b
    c = framework.c

    alpha = framework.alpha * pi / 180
    beta = framework.beta * pi / 180
    gamma = framework.gamma * pi / 180

    alphastar = acos((cos(beta) * cos(gamma) - cos(alpha)) / sin(beta) / sin(gamma)) 
    A = [a  b * cos(gamma)  c * cos(beta);
        0  b * sin(gamma)  c * -1 * sin(beta) * cos(alphastar);
        0  0  c * sin(beta) * sin(alphastar)]
   
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
    sx = range(start=0, stop=1, length=sizea)
    sy = range(start=0, stop=1, length=sizeb)
    sz = range(start=0, stop=1, length=sizec)
    potential = zeros(sizea, sizeb, sizec, 4)
    
    for (index, _) in enumerate(1:1:rotations)

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
        
        for (i, fa) in enumerate(sx), (j, fb) in enumerate(sy), (k, fc) in enumerate(sz)

            coordinates = A * [fa, fb, fc]
            x = coordinates[1]
            y = coordinates[2]
            z = coordinates[3]

            potential[i, j, k, 1] = x
            potential[i, j, k, 2] = y
            potential[i, j, k, 3] = z

            for framework_atom in framework.atoms, probe_atom in rotated_probe.atoms
                
                sig1 = atom_properties[probe_atom.species].sigma
                eps1 = atom_properties[probe_atom.species].epsilon
                q1 = atom_properties[probe_atom.species].charge

                sig2 = atom_properties[framework_atom.species].sigma
                eps2 = atom_properties[framework_atom.species].epsilon
                q2 = atom_properties[framework_atom.species].charge
                
                cmr = sqrt((x - framework_atom.x)^2 + (y - framework_atom.y)^2 + 
                (z - framework_atom.z)^2)

                if cmr < 0.5 * sig2
                    potential[i, j, k, 4] = 1
                    @goto next_point
                end

                # Lorentz-Berthelot mixing rules and charge product
                sig = (sig1 + sig2) / 2
                eps = sqrt(eps1 * eps2)
                q = q1 * q2 

                for offset in pbc_offsets

                    f_atom_x = framework_atom.x + offset[1]
                    f_atom_y = framework_atom.y + offset[2]
                    f_atom_z = framework_atom.z + offset[3]

                    p_atom_x = probe_atom.x + x
                    p_atom_y = probe_atom.y + y
                    p_atom_z = probe_atom.z + z
                    
                    r = sqrt((p_atom_x - f_atom_x)^2 + (p_atom_y - f_atom_y)^2 + 
                    (p_atom_z - f_atom_z)^2)

                    if r < cutoff
                        potential[i, j, k, 4] += lennard_jones_energy(sig, eps, r)
                        potential[i, j, k, 4] += coloumb_energy(q, r)
                    elseif r > cutoff
                        continue
                    end

                end
            end
        @label next_point
        end
    end
    
    # Check-up loop
    message = rpad("---- Box counting statistics ", 81, "-")
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
    argument_string = lpad(positive_potential_boxes/total_boxes, 40, " ")
    write(output_file, "$key_string $argument_string\n")
    
    key_string = rpad("Positive potential box ratio [-]", 40, " ")
    argument_string = lpad(negative_potential_boxes/total_boxes, 40, " ")
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

    if save == "yes"
        
        mkpath("Output")
        
        num = Integer(sizea * sizeb * sizec)
        x = zeros(num)
        y = zeros(num)
        z = zeros(num)
        pot = zeros(num)
        for i in 1:1:sizea, j in 1:1:sizeb, k in 1:1:sizec
            index = i + (j-1) * sizea + (k-1) * sizeb^2
            x[index] = potential[i, j, k, 1]
            y[index] = potential[i, j, k, 2]
            z[index] = potential[i, j, k, 3]
            pot[index] = potential[i, j, k, 4]
        end

        p = scatter(x, y, z, marker_z=pot, aspect_ratio=:equal, markersize=2, camera=(0, -90))
        xlabel!(p, "X [\$\\AA\$]")
        ylabel!(p, "Y [\$\\AA\$]")
        zlabel!(p, "Z [\$\\AA\$]")
    
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
    framework::Framework, potential::Array{Float64, 4}, sizea::Integer, sizeb::Integer,
    sizec::Integer, npoints::Integer, save::SubString{String}) 

    framework_mass = compute_framework_mass(atom_properties, framework)
    framework_volume = compute_framework_volume(framework)

    # Compute the volume of a sample point in ml
    sample_volume = framework_volume * 1e-24 / sizea / sizeb / sizec
    
    minimum_potential = minimum(potential)
    potential_range = range(start=minimum_potential, stop=0.0, length=npoints)
    
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
