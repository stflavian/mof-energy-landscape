using Plots
using ProgressBars

struct AtomProperties
    epsilon::Real
    sigma::Real
    charge::Real
end

struct FrameworkProperties
    a::Real
    b::Real
    c::Real
    alpha::Real
    beta::Real
    gamma::Real
end

struct Atom
    species::String
    x::Real
    y::Real
    z::Real
end


"""
TODO
"""
function lennard_jones_energy(sigma, epsilon, distance)
    sigma::Real
    epsilon::Real
    distance::Real

    # Boltzmann constant
    kb = 1.3806e-23

    frac = sigma / distance    
    return 4 * epsilon * kb * (frac^12 - frac^6)
end


"""
TODO
"""
function coloumb_energy(charge, distance)
    charge::Real
    distance::Real

    # Vacuum dielectric constant
    q = 1.6022e-19
    ke = 8.9876e9

    return q^2 * charge * ke * 1e10 / distance
end


"""
TODO
"""
function read_input(path)
    path::String

    atom_prop = Dict()
    frame_prop = FrameworkProperties
    framework = []
    probe = String
    
    input_file = readlines(path)

    for line in input_file
        
        if length(line) == 0
            continue
        end

        line = split(line)

        if line[1] == "ATOMPROP"

            species = line[2]
            epsilon = parse(Float64, line[3])
            sigma = parse(Float64, line[4])
            charge = parse(Float64, line[5])
            
            atom_prop[species] = AtomProperties(epsilon, sigma, charge)
        
        elseif line[1] == "FRAMEPROP"
            
            a = parse(Float64, line[2])
            b = parse(Float64, line[3])
            c = parse(Float64, line[4])
            alpha = parse(Float64, line[5])
            beta = parse(Float64, line[6])
            gamma = parse(Float64, line[7])
            
            frame_prop = FrameworkProperties(a, b, c, alpha, beta, gamma)
        
        elseif line[1] == "FRAMEWORK"

            species = line[2]
            x = parse(Float64, line[3])
            y = parse(Float64, line[4])
            z = parse(Float64, line[5])

            frame = Atom(species, x, y, z)
            push!(framework, frame)

        elseif line[1] == "PROBE"
            
            probe = line[2]

        end
    end

    return atom_prop, frame_prop, framework, probe
end


"""
TODO
"""
function compute_potential(atom_prop, frame_prop, framework, probe, size)
    atom_prop::Dict
    frame_prop::FrameworkProperties
    framework::Array
    probe::SubString{String}
    size::Integer
    
    # Compute the transformation matrix for fractional to Cartesian coordinates
    a = frame_prop.a
    b = frame_prop.b
    c = frame_prop.c

    alpha = frame_prop.alpha * pi / 180
    beta = frame_prop.beta * pi / 180
    gamma = frame_prop.gamma * pi / 180

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
    s = range(start=0, stop=1, length=size)
    potential = zeros(size, size, size, 4)

    sig1 = atom_prop[probe].sigma
    eps1 = atom_prop[probe].epsilon
    q1 = atom_prop[probe].charge

    for (i, fa) in tqdm(enumerate(s)), (j, fb) in enumerate(s), (k, fc) in enumerate(s)
                
        coordinates = A * [fa, fb, fc]
        x = coordinates[1]
        y = coordinates[2]
        z = coordinates[3]
        
        potential[i, j, k, 1] = x
        potential[i, j, k, 2] = y
        potential[i, j, k, 3] = z

        for atom in framework
        
            sig2 = atom_prop[atom.species].sigma
            eps2 = atom_prop[atom.species].epsilon
            q2 = atom_prop[atom.species].charge
            
            # Lorentz-Berthelot mixing rules and charge product
            sig = (sig1 + sig2) / 2
            eps = sqrt(eps1 * eps2)
            q = q1 * q2 
            
            for offset in pbc_offsets

                # Compute the position of the atomic image using PBC
                pb_atom_x = atom.x + offset[1]
                pb_atom_y = atom.y + offset[2]
                pb_atom_z = atom.z + offset[3]

                r = sqrt((pb_atom_x - x)^2 + (pb_atom_y - y)^2 + (pb_atom_z - z)^2)
            
                if  sig < r < 14
                    potential[i, j, k, 4] += lennard_jones_energy(sig, eps, r)
                    potential[i, j, k, 4] += coloumb_energy(q, r)
                elseif r < sig
                    potential[i, j, k, 4] = 0
                    @goto finish_potenial_calculation
                elseif r > 14
                    continue
                end

            end
        end
        @label finish_potenial_calculation

        # Conversion from J to kJ/mol
        potential[i, j, k, 4] *= 6.022e23 / 1000
    end
    return potential
end


function compute_characteristic(potential, size, frame_prop)
    potential::Array
    size::Integer
    frame_prop::FrameworkProperties

    maximum_potential = maximum(potential)
    potential_range = range(start=1e-7, stop=maximum_potential, length=30)
    
    a = frame_prop.a
    b = frame_prop.b
    c = frame_prop.c

    alpha = frame_prop.alpha * pi / 180
    beta = frame_prop.beta * pi / 180
    gamma = frame_prop.gamma * pi / 180

    cell_volume = a * b * c * sqrt(sin(alpha)^2 + sin(beta)^2 + sin(gamma)^2 + 
            2 * cos(alpha) * cos(beta) * cos(gamma) - 2)

    sample_volume = cell_volume / size^3
    
    measured_volumes = zeros(30)
    for (index, ads_potential) in enumerate(potential_range)
        counter = 0
        for i in 1:1:size, j in 1:1:size, k in 1:1:size
            if potential[i, j, k, 4] >= ads_potential
                counter += 1
            end
        end
        volume = counter * sample_volume
        measured_volumes[index] = volume
    end
    
    characteristic_plot = plot(potential_range, measured_volumes)
    xlabel!(characteristic_plot, "Adsorption potential [kJ/mol]")
    ylabel!(characteristic_plot, "Filling volume [\$\\AA^{3}\$]")
    savefig(characteristic_plot, "characteristic.png") 
end


"""
TODO
"""
function main()
    
    SIZE = 50

    atom_prop, frame_prop, framework, probe = read_input("geo.in")
    potential = compute_potential(atom_prop, frame_prop, framework, probe, SIZE)
    
    x = zeros(SIZE^3)
    y = zeros(SIZE^3)
    z = zeros(SIZE^3)
    pot = zeros(SIZE^3)
    for i in 1:1:SIZE, j in 1:1:SIZE, k in 1:1:SIZE
        index = i + (j-1) * SIZE + (k-1) * SIZE^2
        x[index] = potential[i, j, k, 1]
        y[index] = potential[i, j, k, 2]
        z[index] = potential[i, j, k, 3]
        pot[index] = potential[i, j, k, 4]
    end

    compute_characteristic(potential, SIZE, frame_prop)

    p = scatter(x, y, z, marker_z=pot, aspect_ratio=:equal, markersize=2, camera=(0, -90))
    xlabel!(p, "X [\$\\AA\$]")
    ylabel!(p, "Y [\$\\AA\$]")
    zlabel!(p, "Z [\$\\AA\$]")
    savefig(p, "random.png")
end

main()
