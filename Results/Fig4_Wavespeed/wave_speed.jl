# Linear stability analysis for Fisher-Stefan model
# Alex Tam, 02/03/2022

# Load packages
using Parameters
using LinearAlgebra
using Printf
using Plots
using Measures
using LaTeXStrings

"Data structure for model parameters"
@with_kw struct Params
    uf::Float64 = 0.1 # [-] Background density at interface
    Lξ::Float64 = 20.0 # [-] Spatial domain limit
    Nξ::Int = 1001 # [-] Number of grid points
end

"Leading-order solution"
function leading_order(ξ, dξ, par, κ)
    # Shooting method
    c = 1.0 # Initial guess
    for j = 1:20
        # Solve BVP using Newton's method
        c_old = c
        ϵ = 1e-6 # Small perturbation for numerical derivative
        u = bvp_lo(dξ, c, par)
        up = bvp_lo(dξ, c+ϵ, par)
        f = (3*u[par.Nξ] - 4*u[par.Nξ-1] + u[par.Nξ-2])/(2*dξ) + c/κ
        fp = (3*up[par.Nξ] - 4*up[par.Nξ-1] + up[par.Nξ-2])/(2*dξ) + (c+ϵ)/κ
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

"Main function for stability analysis"
function main()
    par = Params() # Initialise data structure of model parameters
    ξ = range(-par.Lξ, 0.0, length = par.Nξ); dξ = ξ[2] - ξ[1] # Computational domain (x)
    # Pre-allocate
    Κ = range(-1.05, 5.0, length = 101)
    C = Vector{Float64}()
    # Obtain leading-order solution
    for κ in Κ
        c, u0 = leading_order(ξ, dξ, par, κ)
        push!(C, c)
    end
    # Numerical results
    κ_n = [ -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    c_n = [ -1.265792, -0.935469, -0.712235, -0.542538, -0.410977, -0.301454, -0.209769,  -0.130543, -0.061015, 0.053005, 0.102110, 0.145130, 0.187406, 0.223620, 0.259390, 0.293878, 0.322382, 0.351350, 0.380074, 0.495444, 0.587141, 0.663968, 0.730893, 0.785856, 0.836070, 0.882904, 0.919883, 0.995538, 1.040286, 1.099003, 1.128895, 1.174941 ]
    # Plot
    gr(); plot() # Load GR plotting backend and clear previous plots
    default(fontfamily = "Computer Modern", titlefontsize = 14, guidefontsize = 20, tickfontsize = 14, legendfontsize = 12)
    plot(Κ, C, xlabel = L"$\kappa$", ylabel = L"$c$", linecolor = :black, linewidth = 2, grid = true, margin=2mm, xlims=(-1.5,5.0), ylims=(-3,1.0), xticks=-1:1:10, yticks=-3:0.5:1.0, legend =:bottomright, label = "Leading-Order Shooting")
    scatter!(κ_n, c_n, markersize = 5, markershape = :xcross, markercolor = :red, label = "Full 2D Numerical")
    plot!([-1/(1-par.uf), -1/(1-par.uf)], [-3,1.0], label=false, linewidth = 2, linestyle=:dash, linecolor=:black)
    savefig("Fig4.pdf")
end

@time main()