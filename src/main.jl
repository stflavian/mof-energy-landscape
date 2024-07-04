include("interpreter.jl")


function main()
    output_file = open("escape.out", "w+")
    input_data = read_input_file_beta("config.in")
    probe = read_probe_file(input_data["PROBE"])
    framework = read_framework_file(input_data["FRAMEWORK"])
    properties = read_properties_file(input_data["PROPERTIES"])

    message = rpad("---- Input data ", 81, "-")
    write(output_file, "$message\n")
    for key in keys(input_data)

        key_string = rpad(key, 40, " ")
        argument_string = lpad(input_data[key], 40, " ")

        write(output_file, "$key_string $argument_string\n")
    end
    write(output_file, "\n")

    
    message = rpad("---- Atomic properties ", 81, "-")
    write(output_file, "$message\n")

    species_string = rpad("Species", 10, " ")
    sigma_string = rpad("Sigma [angstrom]", 19, " ")
    epsilon_string = rpad("Epsilon/KB [K]", 16, " ")
    charge_string = rpad("Charge [e-]", 16, " ")
    mass_string = rpad("Mass [a.u.]", 16, " ")

    write(output_file, "$species_string $sigma_string $epsilon_string $charge_string $mass_string\n")
    for species in keys(properties)
    
        species_string = rpad(species, 10, " ")
        sigma_string = rpad(properties[species].sigma, 19, " ")
        epsilon_string = rpad(properties[species].epsilon, 16, " ")
        charge_string = rpad(properties[species].charge, 16, " ")
        mass_string = rpad(properties[species].mass, 16, " ")

        write(output_file, "$species_string $sigma_string $epsilon_string $charge_string $mass_string\n")
    end
    write(output_file, "\n")


    message = rpad("---- Unitcell geometry ", 81, "-")
    write(output_file, "$message\n")
    
    key_string = rpad("A", 40, " ")
    argument_string = lpad(framework.a, 40, " ")
    write(output_file, "$key_string $argument_string\n")
    
    key_string = rpad("B", 40, " ")
    argument_string = lpad(framework.b, 40, " ")
    write(output_file, "$key_string $argument_string\n")
    
    key_string = rpad("C", 40, " ")
    argument_string = lpad(framework.c, 40, " ")
    write(output_file, "$key_string $argument_string\n")
    
    key_string = rpad("ALPHA", 40, " ")
    argument_string = lpad(framework.alpha, 40, " ")
    write(output_file, "$key_string $argument_string\n")
    
    key_string = rpad("BETA", 40, " ")
    argument_string = lpad(framework.beta, 40, " ")
    write(output_file, "$key_string $argument_string\n")
    
    key_string = rpad("GAMMA", 40, " ")
    argument_string = lpad(framework.gamma, 40, " ")
    write(output_file, "$key_string $argument_string\n")
    write(output_file, "\n")


    message = rpad("---- Framework atoms ", 81, "-")
    write(output_file, "$message\n")

    index_string = rpad("Index", 10, " ")
    species_string = rpad("Species", 19, " ")
    x_string = rpad("X [angstrom]", 16, " ")
    y_string = rpad("Y [angstrom]", 16, " ")
    z_string = rpad("Z [angstrom]", 16, " ")

    write(output_file, "$index_string $species_string $x_string $y_string $z_string\n")
    for (index, atom) in enumerate(framework.atoms)
    
        index_string = rpad(index, 10, " ")
        species_string = rpad(atom.species, 19, " ")
        x_string = rpad(atom.x, 16, " ")
        y_string = rpad(atom.y, 16, " ")
        z_string = rpad(atom.z, 16, " ")

        write(output_file, "$index_string $species_string $x_string $y_string $z_string\n")
    end
    write(output_file, "\n")
    

    message = rpad("---- Framework properties ", 81, "-")
    write(output_file, "$message\n")

    framework_charge = compute_framework_charge(properties, framework)
    key_string = rpad("Framework charge [e-]", 40, " ")
    argument_string = lpad(framework_charge, 40, " ")
    write(output_file, "$key_string $argument_string\n")
    
    framework_mass = compute_framework_mass(properties, framework)
    key_string = rpad("Framework mass [kg]", 40, " ")
    argument_string = lpad(framework_mass, 40, " ")
    write(output_file, "$key_string $argument_string\n")
    
    framework_volume = compute_framework_volume(framework)
    key_string = rpad("Framework volume [angstrom^3]", 40, " ")
    argument_string = lpad(framework_volume, 40, " ")
    write(output_file, "$key_string $argument_string\n")
    
    framework_density = compute_framework_density(framework_mass, framework_volume)
    key_string = rpad("Framework density [kg/m^3]", 40, " ")
    argument_string = lpad(framework_density, 40, " ")
    write(output_file, "$key_string $argument_string\n")
    write(output_file, "\n")

    
    message = rpad("---- Probe atoms ", 81, "-")
    write(output_file, "$message\n")

    index_string = rpad("Index", 10, " ")
    species_string = rpad("Species", 19, " ")
    x_string = rpad("X [angstrom]", 16, " ")
    y_string = rpad("Y [angstrom]", 16, " ")
    z_string = rpad("Z [angstrom]", 16, " ")

    write(output_file, "$index_string $species_string $x_string $y_string $z_string\n")
    for (index, atom) in enumerate(probe.atoms)
    
        index_string = rpad(index, 10, " ")
        species_string = rpad(atom.species, 19, " ")
        x_string = rpad(atom.x, 16, " ")
        y_string = rpad(atom.y, 16, " ")
        z_string = rpad(atom.z, 16, " ")

        write(output_file, "$index_string $species_string $x_string $y_string $z_string\n")
    end
    write(output_file, "\n")


    message = rpad("---- Probe properties ", 81, "-")
    write(output_file, "$message\n")

    probe_charge = compute_probe_charge(properties, probe)
    key_string = rpad("Probe charge [e-]", 40, " ")
    argument_string = lpad(probe_charge, 40, " ")
    write(output_file, "$key_string $argument_string\n")
    
    probe_mass = compute_probe_mass(properties, probe)
    key_string = rpad("Probe mass [kg]", 40, " ")
    argument_string = lpad(probe_mass, 40, " ")
    write(output_file, "$key_string $argument_string\n")
    
    probe, center_of_mass = remove_probe_centermass(properties, probe)
    key_string = rpad("Probe center of mass [angstrom]", 36, " ")
    argument_string1 = lpad(center_of_mass[1], 14, " ")
    argument_string2 = lpad(center_of_mass[2], 14, " ")
    argument_string3 = lpad(center_of_mass[3], 14, " ")
    write(output_file, "$key_string $argument_string1 $argument_string2 $argument_string3\n")
    write(output_file, "\n")

    sizex = parse(Int64, input_data["XPOINTS"])
    sizey = parse(Int64, input_data["YPOINTS"])
    sizez = parse(Int64, input_data["ZPOINTS"])
    cutoff = parse(Float64, input_data["CUTOFF"])

    potential = compute_potential_landscape(properties, framework, probe, sizex, sizey,
    sizez, cutoff, input_data["SAVE_POTENTIAL"])

    npoints = parse(Int64, input_data["CHARACTERISTIC_POINTS"])
    compute_characteristic(properties, framework, potential, sizex, sizey, sizez, 
    npoints, input_data["SAVE_CHARACTERISTIC"])
    
    close(output_file)
end


main()
