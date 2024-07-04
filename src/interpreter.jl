using Plots
using ProgressBars
using LinearAlgebra

include("constants.jl")


"""
Structure for storing atomic Lennard-Jones parameters and charge.
"""
struct AtomProperties
    epsilon::Real
    sigma::Real
    charge::Real
    mass::Real
end


"""
Structure for storing the Cartesian coordinates of the framework's atoms. 
"""
struct Atom
    species::String
    x::Real
    y::Real
    z::Real
end


"""
Structure for storing the framework's lattice parameters.
"""
struct FrameworkProperties
    a::Real
    b::Real
    c::Real
    alpha::Real
    beta::Real
    gamma::Real
end


"""
TODO
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
TODO
"""
struct Probe
    atoms::Vector{Atom}
end


"""
TODO
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
TODO
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
TODO
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
TODO
"""
function read_input_file_beta(path::String)
    keywords = Dict(
        "FRAMEWORK" => Nothing, 
        "PROPERTIES" => Nothing,
        "CUTOFF" => 15,
        "PROBE" => Nothing,
        "XPOINTS" => "30", 
        "YPOINTS" => "30", 
        "ZPOINTS" => "30", 
        "SAVE_POTENTIAL" => "yes", 
        "CHARACTERISTIC_POINTS" => 30,
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

Compute the 6-12 Lennard-Jones energy between two particles.

# Arguments
- `sigma::Real`: the size parameter, usually the sum of the particle's radii.
- `epsilon::Real`: the depth of the potential well.
- `distance::Real`: the distance between the center of masses of the two particles.
"""
function lennard_jones_energy(sigma::Real, epsilon::Real, distance::Real)
    frac = sigma / distance    
    return 4 * epsilon * KB * (frac^12 - frac^6)
end


"""
    coloumb_energy(charge::Real, distance::Real)

Compute the electrostatic energy between two particles.

# Arguments
- `charge::Real`: the product of the charges of the two particles.
- `distance::Real`: the distance between the center of masses of the two particles.
"""
function coloumb_energy(charge::Real, distance::Real)
    return Q^2 * charge * KE * 1e10 / distance
end


"""
TODO
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
TODO
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
TODO
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
TODO
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
TODO
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
TODO
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
TODO
"""
function compute_framework_density(framework_mass::Real, framework_volume::Real)
    return framework_mass * 1e30 / framework_volume
end


"""
    compute_potential(atom_prop::Dict{SubString{String}, AtomProperties}, 
    frame_prop::FrameworkProperties, framework::Vector{Atom}, probe::SubString{String}, 
    size::Integer)

Compute the potential 

# Arguments
- `atom_prop::Dict{SubString{String}, AtomProperties}`: 
- `frame_prop::FrameworkProperties`:
- `framework::Vector{Atom}`:
- `probe::SubString{String}`:
- `size::Integer`: 
"""
function compute_potential_landscape(atom_prop::Dict{SubString{String}, AtomProperties}, 
    framework::Framework, probe::Probe, sizex::Integer, sizey::Integer, sizez::Integer,
    cutoff::Real, save::SubString{String})
    
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
    pbc_offsets = []
    for i in [-1, 0, 1], j in [-1, 0, 1], k in [-1, 0, 1]
        offset = A * [i, j, k]
        push!(pbc_offsets, offset)
    end
    
    # Initialize arrays and assign parameters for probe
    sx = range(start=0, stop=1, length=sizex)
    sy = range(start=0, stop=1, length=sizey)
    sz = range(start=0, stop=1, length=sizez)
    potential = zeros(sizex, sizey, sizez, 4)

    for (i, fa) in tqdm(enumerate(sx)), (j, fb) in enumerate(sy), (k, fc) in enumerate(sz)
                
        coordinates = A * [fa, fb, fc]
        x = coordinates[1]
        y = coordinates[2]
        z = coordinates[3]
        
        potential[i, j, k, 1] = x
        potential[i, j, k, 2] = y
        potential[i, j, k, 3] = z

        number_of_atoms_in_probe = length(probe.atoms)
        
        if number_of_atoms_in_probe == 1
            rotation_trials = 1
        else
            rotation_trials = 10
        end

        for _ in 1:1:rotation_trials

            angle = rand() * 2 * pi
            axis_i = 1 - 2 * rand()
            axis_j = 1 - 2 * rand()
            axis_k = 1 - 2 * rand()
            rotated_probe = rotate_probe(probe, angle, axis_i, axis_j, axis_k)

            for framework_atom in framework.atoms, probe_atom in rotated_probe.atoms
            
                sig1 = atom_prop[probe_atom.species].sigma
                eps1 = atom_prop[probe_atom.species].epsilon
                q1 = atom_prop[probe_atom.species].charge
        
                sig2 = atom_prop[framework_atom.species].sigma
                eps2 = atom_prop[framework_atom.species].epsilon
                q2 = atom_prop[framework_atom.species].charge
            
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
            
                    if  0.5 * sig < r < cutoff
                        potential[i, j, k, 4] += lennard_jones_energy(sig, eps, r)
                        potential[i, j, k, 4] += coloumb_energy(q, r)
                    elseif r < 0.5 * sig
                        potential[i, j, k, 4] = 0
                        @goto finish_potenial_calculation
                    elseif r > cutoff
                        continue
                    end
                end
            end
        end

        if potential[i, j, k, 4] > 0
            potential[i, j, k, 4] = 0
        end
        
        @label finish_potenial_calculation

        # Conversion from J to kJ/mol
        potential[i, j, k, 4] *= NA * 1e-3 / rotation_trials
    end
   
    if save == "yes"
        
        mkpath("Output")
        
        num = Integer(sizex * sizey * sizez)
        x = zeros(num)
        y = zeros(num)
        z = zeros(num)
        pot = zeros(num)
        for i in 1:1:sizex, j in 1:1:sizey, k in 1:1:sizez
            index = i + (j-1) * sizex + (k-1) * sizey^2
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
    compute_potential(atom_prop::Dict{SubString{String}, AtomProperties}, 
    frame_prop::FrameworkProperties, framework::Vector{Atom}, probe::SubString{String}, 
    size::Integer)

Compute the potential 

# Arguments
- `atom_prop::Dict{SubString{String}, AtomProperties}`: 
- `frame_prop::FrameworkProperties`:
- `framework::Vector{Atom}`:
- `potential::Array{Flat64, 4}`:
- `size::Integer`: 

"""
function compute_characteristic(atom_prop::Dict{SubString{String}, AtomProperties}, 
    framework::Framework, potential::Array{Float64, 4}, sizex::Integer, sizey::Integer,
    sizez::Integer, npoints::Integer, save::SubString{String}) 

    framework_mass = compute_framework_mass(atom_prop, framework)
    framework_volume = compute_framework_volume(framework)

    # Compute the volume of a sample point in ml
    sample_volume = framework_volume * 1e-24 / sizex / sizey / sizez
    
    minimum_potential = minimum(potential)
    potential_range = range(start=minimum_potential, stop=-0.000001, length=npoints)
    
    if save == "yes"
        mkpath("Output")

        output_file = open("Output/characteristic.dat", "w+")
        write(output_file, "# Potential [kJ/mol] \t Volume [ml/g] \n")

        measured_potential = zeros(npoints)
        measured_volumes = zeros(npoints)
        for (index, ads_potential) in enumerate(potential_range)
            counter = 0
            for i in 1:1:sizex, j in 1:1:sizey, k in 1:1:sizez
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

        characteristic_plot = plot(measured_potential, measured_volumes)
        xlabel!(characteristic_plot, "Potential [kJ/mol]")
        ylabel!(characteristic_plot, "Volume [ml/g]")
        savefig(characteristic_plot, "Output/characteristic.png") 
    end
end
