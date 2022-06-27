# Control function for 2D Fisher-Stefan level-set solutions
# Alex Tam, 02/02/2022

# Load packages
using Parameters
using Printf
using Dierckx
using LinearAlgebra
using DifferentialEquations
using Plots
using Measures
using LaTeXStrings
using DelimitedFiles

# Include external files
include("twic.jl")
include("domain.jl")
include("fkpp.jl")
include("velocity_extension.jl")
include("interface_density.jl")
include("interface_speed.jl")
include("level-set.jl")
include("reinitialisation.jl")

"Data structure for parameters"
@with_kw struct Params
    D::Float64 = 1.0 # [-] Diffusion coefficient
    λ::Float64 = 1.0 # [-] Reaction rate
    κ::Float64 = -0.5 # [-] Inverse Stefan number
    α::Float64 = 1.0 # [-] Maximum initial density
    β::Float64 = 6.0 # [-] Initial interface position
    uf::Float64 = 0.1 # [-] Background density at interface
    θb::Float64 = 0.01 # [-] Threshold for whether a grid point is close to interface (relative to Δx)
    θ::Float64 = 1.99 # [-] Parameter for minmod flux-limiter
    Lx::Float64 = 7.5 # [-] Spatial domain limit (x)
    Ly::Float64 = 10.0 # [-] Spatial domain limit (y)
    T::Float64 = 1.0 # [-] End time
    Nx::Int = 301 # [-] Number of grid points (x)
    Ny::Int = 401 # [-] Number of grid points (y)
    Nt::Int = 1001 # [-] Number of time steps
    Nξ::Int = 2001 # [-] Number of grid points for travelling wave (ξ)
    V_Iterations::Int = 20 # [-] Number of iterations for velocity extrapolation PDE
    ϕ_Iterations::Int = 20 # [-] Number of iterations for reinitialisation PDE
    γ::Float64 = 0.1 # [-] Surface tension coefficient
    ε::Float64 = 0.1 # [-] Small amplitude of perturbations
    q::Float64 = 0.0 # [-] Wave number of perturbations
end

"Interpolate to obtain initial condition"
function ic(par, x, y)
    U = Array{Float64}(undef, par.Nx, par.Ny) # Pre-allocate 2D array of U
    ϕ = Array{Float64}(undef, par.Nx, par.Ny) # Pre-allocate 2D array of ϕ
    # Compute travelling wave solution and correction
    u0, u1, ξ = twic(par)
    # Construct linear splines for interpolation
    spl_0 = Spline1D(ξ, u0) # Generate 1D spline
    spl_1 = Spline1D(ξ, u1) # Generate 1D spline
    # Compute ϕ = ξ at each grid point
    for i = 1:length(x)
        for j = 1:length(y)
            ϕ[i,j] = x[i] - par.β - par.ε*cos(par.q*y[j])
        end
    end
    # Interpolate travelling wave to obtain initial densities
    for i = 1:length(x)
        for j = 1:length(y)
            if ϕ[i,j] < 0 # If grid point is in Ω(0)
                U[i,j] = spl_0(ϕ[i,j]) + par.ε*spl_1(ϕ[i,j])*cos(par.q*y[j]) # Perturbed travelling wave
            else # If grid point is not in Ω(0)
                U[i,j] = par.uf
            end
        end
    end
    return U, ϕ
end

"Build vector from matrix, ordered by entries in D"
function build_vector(U::Array{Float64}, D)
    u = Vector{Float64}() # Pre-allocate empty vector
    for gp in D
        push!(u, U[gp.xInd, gp.yInd])
    end
    return u
end

"Build matrix from vector ordered by entries in D"
function build_u_matrix(u::Vector, y, par, D)
    U = zeros(par.Nx, par.Ny) # Pre-allocate (incorporate a Dirichlet condition on right boundary)
    U[par.Nx,:] .= par.uf # Dirichlet condition on right boundary
    U[1,:] .= par.α # Dirichlet condition on left boundary
    for i = 1:length(D)
        U[D[i].xInd, D[i].yInd] = u[i]
    end
    return U
end

"Build matrix from vector ordered by entries in D"
function build_v_matrix(v::Vector, par, D)
    V = zeros(par.Nx, par.Ny) # Pre-allocate (incorporate a Dirichlet condition on computational boundary)
    for i = 1:length(D)
        V[D[i].xInd, D[i].yInd] = v[i]
    end
    return V
end

"Locate interface position"
function front_position(x, ϕ, par, ny, dx)
    L = 0.0 # Pre-allocate front-position
    Lmax = 0.0 # Pre-allocate max front position
    Lmin = par.Lx # Pre-allocate min front position
    # Find front position using x-direction slice
    for j = 1:par.Ny # Loop over y
        ϕv = ϕ[:,j] # Obtain 1D vector of ϕ
        for i = 1:par.Nx
            if (ϕv[i] < 0) && (ϕv[i+1] >= 0)
                θ = ϕv[i]/(ϕv[i] - ϕv[i+1])
                Lj = x[i] + θ*dx
                if Lj <= Lmin
                    Lmin = Lj
                end
                if Lj >= Lmax
                    Lmax = Lj
                end
                if j == ny
                    L = Lj
                end
            end
        end
    end
    return L, (Lmax-Lmin)/2 # Return amplitude of perturbation
end

"Compute a solution"
function fisher_stefan_2d()
    # Parameters and domain
    par = Params() # Initialise data structure of model parameters
    nx::Int = (par.Nx-1)/2; ny::Int = (par.Ny-1)/2 # Indices for slice plots
    x = range(0, par.Lx, length = par.Nx); dx = x[2] - x[1] # Computational domain (x)
    y = range(0, par.Ly, length = par.Ny); dy = y[2] - y[1] # Computational domain (y)
    t = range(0, par.T, length = par.Nt); dt = t[2] - t[1] # Time domain
    writedlm("x.csv", x); writedlm("y.csv", y); writedlm("t.csv", t) # Write data to files
    # Initial condition
    U, ϕ = ic(par, x, y) # Obtain initial density and ϕ
    writedlm("U-0.csv", U); writedlm("Phi-0.csv", ϕ) # Write data to files
    plot_times = Vector{Int}() # Vector of time-steps at which data is obtained
    writedlm("plot_times.csv", plot_times)
    L = Vector{Float64}() # Preallocate empty vector of interface position
    Amp = Vector{Float64}() # Preallocate empty vector of perturbation amplitude
    Li, amp = front_position(x, ϕ, par, ny, dx)
    push!(L, Li); push!(Amp, amp)
    # Time stepping
    for i = 1:par.Nt-1
        # 1. Find Ω, dΩ, and irregular grid points
        D = find_domain(par, ϕ)
        dΩ = find_interface(par, D, ϕ)
        # 2. Solve FKPP equation on Ω
        uf = interface_density(dΩ, ϕ, par, dx, dy) # Density on interface for BC
        @time U = fkpp(D, dΩ, U, ϕ, uf, y, par, dx, dy, dt, i)
        # 3. Compute extension velocity field
        V = extend_velocity(D, dΩ, U, ϕ, par, dx, dy)
        # 4. Solve level-set equation
        ϕ = level_set(V, ϕ, par, dx, dy, dt)
        # 5. Re-initialise level-set function as a signed-distance
        if mod(i, 100) == 0
            ϕ = reinitialisation(ϕ, par, dx, dy, par.ϕ_Iterations)
        end
        # Optional: Post-processing
        if mod(i, 100) == 0
            writedlm("ux-$i.csv", U[:,ny])
            writedlm("uy-$i.csv", U[nx,:])
            writedlm("U-$i.csv", U)
            writedlm("V-$i.csv", V)
            writedlm("Phi-$i.csv", ϕ)
            push!(plot_times, i)
            writedlm("plot_times.csv", plot_times)
        end
        Li, amp = front_position(x, ϕ, par, ny, dx)
        push!(L, Li)
        push!(Amp, amp)
        writedlm("L.csv", L)
        writedlm("Amp.csv", Amp)
    end
end

@time fisher_stefan_2d()