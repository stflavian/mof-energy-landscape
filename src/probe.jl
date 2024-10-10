using LinearAlgebra

include("constants.jl")
include("reader.jl")

"""
Structure for storing the species and Cartesian coordinates of atoms. 
"""
struct Atom
    species::String
    x::Float64
    y::Float64
    z::Float64
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
function rotate_probe(probe::Probe, angle::Float64, i::Float64, j::Float64, k::Float64)
    
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
