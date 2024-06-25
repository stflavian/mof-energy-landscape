include("interpreter.jl")


function main()
    SIZE = 10

    atom_prop, frame_prop, framework, probe = read_input_file("geo.in")
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
