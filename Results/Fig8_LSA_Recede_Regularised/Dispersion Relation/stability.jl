# Linear stability analysis for Fisher-Stefan model
# Alex Tam, 02/03/2022

# Load packages
using Parameters
using LinearAlgebra
using Printf
using DelimitedFiles
using Plots
using Measures
using LaTeXStrings

"Data structure for model parameters"
@with_kw struct Params
    κ::Float64 = 0.1 # [-] Inverse Stefan number
    γ::Float64 = 0.0 # [-] Surface-tension coefficient
    uf::Float64 = 0.0 # [-] Background density at interface
    Lξ::Float64 = 20.0 # [-] Spatial domain limit
    Nξ::Int = 1001 # [-] Number of grid points
end

"Leading-order solution"
function leading_order(ξ, dξ, par)
    # Shooting method
    c = -1.0 # Initial guess
    for j = 1:20
        # Solve BVP using Newton's method
        c_old = c
        ϵ = 1e-6 # Small perturbation for numerical derivative
        u = bvp_lo(dξ, c, par)
        up = bvp_lo(dξ, c+ϵ, par)
        f = (3*u[par.Nξ] - 4*u[par.Nξ-1] + u[par.Nξ-2])/(2*dξ) + c/par.κ
        fp = (3*up[par.Nξ] - 4*up[par.Nξ-1] + up[par.Nξ-2])/(2*dξ) + (c+ϵ)/par.κ
        d = (fp-f)/ϵ
        c = c - f/d
        if abs(c - c_old) < 1e-6
            return c, bvp_lo(dξ, c, par)
        end
    end
    @printf("Leading-order shooting did not converge.\n")
end

"Solve leading-order two-point BVP"
function bvp_lo(dξ, c, par)
    u = ones(par.Nξ) # Initial guess
    for i = 1:20
        u_old = u
        F = F_lo(dξ, c, u, par) # Construct F
        J = J_lo(dξ, c, u, par) # Construct Jacobian
        u = u - J\F
        if norm(u - u_old) < 1e-6
            return u
        end
    end
    @printf("Leading-order system did not converge.\n")
end

"Vector function for leading-order two-point BVP"
function F_lo(dξ, c, u, par)
    F = Vector{Float64}(undef, par.Nξ) # Pre-allocate F
    F[1] = u[1] - 1
    for i = 2:par.Nξ-1
        F[i] = (u[i+1]-2*u[i]+u[i-1])/(dξ^2) + c*(u[i+1]-u[i-1])/(2*dξ) + u[i]*(1-u[i])
    end
    F[par.Nξ] = u[par.Nξ] - par.uf
    return F
end

"Jacobian for leading-order two-point BVP"
function J_lo(dξ, c, u, par)
    J = zeros(par.Nξ, par.Nξ) # Pre-allocate Jacobian
    J[1,1] = 1.0; J[par.Nξ, par.Nξ] = 1.0
    for i = 2:par.Nξ-1
        J[i,i-1] = 1/(dξ^2) - c/(2*dξ)
        J[i,i] = -2/(dξ^2) + 1 - 2*u[i]
        J[i,i+1] = 1/(dξ^2) + c/(2*dξ)
    end
    return J
end

"Shooting method for stability problem"
function stability_shoot(dξ, par, c, u0, q)
    # Compute derivative of u0
    du0 = Vector{Float64}(undef, par.Nξ)
    du0[1] = (-3*u0[1] + 4*u0[2] - u0[3])/(2*dξ) # One-sided
    du0[par.Nξ] = -c/par.κ
    for i = 2:par.Nξ-1
        du0[i] = (u0[i+1]-u0[i-1])/(2*dξ) # Central
    end
    # Shooting method
    ω = 0.1 # Initial guess
    for j = 1:20
        # Solve BVP using Newton's method
        ω_old = ω
        ϵ = 1e-6 # Small perturbation for numerical derivative
        u = bvp_stability(dξ, ω, c, q, u0, du0, par)
        up = bvp_stability(dξ, ω+ϵ, c, q, u0, du0, par)
        f = (3*u[par.Nξ] - 4*u[par.Nξ-1] + u[par.Nξ-2])/(2*dξ) + ω/par.κ
        fp = (3*up[par.Nξ] - 4*up[par.Nξ-1] + up[par.Nξ-2])/(2*dξ) + (ω+ϵ)/par.κ
        d = (fp-f)/ϵ
        ω = ω - f/d
        if abs(ω - ω_old) < 1e-6
            return ω, bvp_stability(dξ, ω, c, q, u0, du0, par)
        end
    end
    @printf("Stability shooting did not converge.\n")
end

"Solve stability two-point BVP"
function bvp_stability(dξ, ω, c, q, u0, du0, par)
    f = 1 .- 2*u0 .- ω .- q^2 
    g = -(ω+q^2).*du0
    u = 0.1*ones(par.Nξ) # Initial guess
    for i = 1:20
        u_old = u
        F = F_stability(dξ, c, u, f, g, q, par) # Construct F
        J = J_stability(dξ, c, u, f, g, q, par) # Construct Jacobian
        u = u - J\F
        if norm(u - u_old) < 1e-6
            return u
        end
    end
    @printf("Stability system did not converge.\n")
end

"Vector function for stability two-point BVP"
function F_stability(dξ, c, u, f, g, q, par)
    F = Vector{Float64}(undef, par.Nξ) # Pre-allocate F
    F[1] = u[1]
    for i = 2:par.Nξ-1
        F[i] = (u[i+1]-2*u[i]+u[i-1])/(dξ^2) + c*(u[i+1]-u[i-1])/(2*dξ) + f[i]*u[i] - g[i]
    end
    F[par.Nξ] = u[par.Nξ] + par.γ*q^2
    return F
end

"Jacobian for stability two-point BVP"
function J_stability(dξ, c, u, f, g, q, par)
    J = zeros(par.Nξ, par.Nξ) # Pre-allocate Jacobian
    J[1,1] = 1.0; J[par.Nξ, par.Nξ] = 1.0
    for i = 2:par.Nξ-1
        J[i,i-1] = 1/(dξ^2) - c/(2*dξ)
        J[i,i] = -2/(dξ^2) + f[i]
        J[i,i+1] = 1/(dξ^2) + c/(2*dξ)
    end
    return J
end

"Main function for stability analysis"
function main()
    par = Params(κ = -0.8, uf = 0.1, γ = 0.1) # Initialise data structure of model parameters
    ξ = range(-par.Lξ, 0.0, length = par.Nξ); dξ = ξ[2] - ξ[1] # Computational domain (x)
    # Pre-allocate
    Q = range(0.0, 10.0, length = 65)
    Ω = Vector{Float64}()
    # Obtain leading-order travelling-wave solution
    c, u0 = leading_order(ξ, dξ, par)
    # Obtain dispersion relation
    for q in Q
        ω, u = stability_shoot(dξ, par, c, u0, q)
        push!(Ω, ω)
    end
    # Save data to files
    k = par.κ; uf = par.uf; g = par.γ
    writedlm("q-($k,$uf,$g).csv", Q)
    writedlm("omega-($k,$uf,$g).csv", Ω)
    return Q, Ω
end

@time Q, Ω = main()