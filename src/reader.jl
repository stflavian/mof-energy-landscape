"""
Structure for storing the properties of an atomic species.
"""
struct AtomProperties
    epsilon::Float64
    sigma::Float64
    charge::Float64
    mass::Float64
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
