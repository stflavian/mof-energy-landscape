include("reader.jl")
include("constants.jl")

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
Structure for storing the framework unit cell parameters and position of atoms.
"""
struct Framework
    a::Float64
    b::Float64
    c::Float64
    alpha::Float64
    beta::Float64
    gamma::Float64
    atoms::Vector{Atom}
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
function compute_framework_density(framework_mass::Float64, framework_volume::Float64)
    return framework_mass * 1e30 / framework_volume
end


"""
    fractional_to_cartesian(conversion_matrix::Matrix{Float64}, fa::Float64,
    fb::Float64, fc::Float64)

Convert fractional coordinates to cartesian coordinates using the framework's
conversion matrix.

# Arguments
- `conversion_matrix::Matrix{Float64}`: the framework's conversion matrix. 
- `fa::Float64`: fractional position on the a axis.
- `fb::Float64`: fractional position on the b axis.
- `fc::Float64`: fractional position on the c axis.
"""
function fractional_to_cartesian(conversion_matrix::Matrix{Float64}, fa::Float64, 
    fb::Float64, fc::Float64)
    
    return conversion_matrix * [fa, fb, fc]
end


"""
    compute_conversion_matrix(framework::Framework)

Compute the conversion matrix of the framework.

# Arguments
- `framework::Framework`: structure containing the properties of the framework.
"""
function compute_conversion_matrix(framework::Framework)
   
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
    
    return A
end


function generate_pbc(framework::Framework)
    
    number_of_atoms = length(framework.atoms)
    pbc_atoms = Vector{Atom}(undef, 27 * number_of_atoms)
    
    A = compute_conversion_matrix(framework)
    i = 1
    for ox in [0, -1, 1], oy in [0, -1, 1], oz in [0, -1, 1] 
        for atom in framework.atoms
            x, y, z = fractional_to_cartesian(A, atom.x + ox, atom.y + oy, atom.z + oz)
            pbc_atoms[i] = Atom(atom.species, x, y, z)
            i += 1
        end
    end
    
    return Framework(framework.a * 3, framework.b * 3, framework.c * 3, 
                     framework.alpha, framework.beta, framework.gamma, pbc_atoms)
end
