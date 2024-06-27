include("interpreter.jl")


function main()
    SIZE = 30

    atom_prop, frame_prop, framework, probe = read_input_file("geo.in")
    potential = compute_potential(atom_prop, frame_prop, framework, probe, SIZE)
    compute_characteristic(atom_prop, frame_prop, framework, potential, SIZE)
end


main()
