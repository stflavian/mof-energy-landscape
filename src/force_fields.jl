include("constants.jl")

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
function lennard_jones_energy(sigma::Float64, epsilon::Float64, distance::Float64)
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
function coloumb_energy(charge::Float64, distance::Float64)
    return Q^2 * charge * KE * 1e10 / distance
end


"""
    compute_distance(dx::Float64, dy::Float64, dz::Float64)

Compute the distance between two points given the difference in positions. 

# Arguments
- `dx::Float64`: position difference on the x axis.
- `dy::Float64`: position difference on the y axis.
- `dz::Float64`: position difference on the z axis.
"""
function compute_distance(dx::Float64, dy::Float64, dz::Float64)
    return sqrt(dx^2 + dy^2 + dz^2)
end
