function write_section(file::IO, title::String)
    message = rpad("==== $title ", 81, "=")
    write(file, message)
    println(file, "")
end

function write_subsection(file::IO, title::String)
    message = rpad("---- $title ", 81, "=")
    println(file, message)
    println(file, "")
end

function write_result(file::IO, title::String, result::Real)
    key = rpad(title, 40, " ")
    argument = lpad(result, 40, " ")
    println(file, key, " ", argument)
end

function write_xyz(file::IO, atoms::Vector{Atom})
    
    index = rpad("Index", 10, " ")
    species = rpad("Species", 19, " ")
    x = rpad("X [angstrom]", 16, " ")
    y = rpad("Y [angstrom]", 16, " ")
    z = rpad("Z [angstrom]", 16, " ")
    println(file, index, " ", species, " ", x, " ", y, " ", z)
    
    for (index, atom) in enumerate(atoms)
        index = rpad(index, 10, " ")
        species = rpad(atom.species, 19, " ")
        x = rpad(atom.x, 16, " ")
        y = rpad(atom.y, 16, " ")
        z = rpad(atom.z, 16, " ")
        println(file, index, " ", species, " ", x, " ", y, " ", z)
    end
end
