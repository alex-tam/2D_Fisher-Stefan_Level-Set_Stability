# Find Ω(t), dΩ(t), and computational stencils
# Alex Tam, 02/02/2022

"Data structure for an interior grid point"
struct GridPoint
    xInd::Int # Index of x-coordinate
    yInd::Int # Index of y-coordinate
    Ω::Bool # Boolean for whether grid point is inside Ω
    dΩ::Bool # Boolean for whether grid point is close to dΩ
    θxp::Float64 # Relative distance to plus stencil node (x)
    θyp::Float64 # Relative distance to plus stencil node (y)
    θxm::Float64 # Relative distance to minus stencil node (x)
    θym::Float64 # Relative distance to minus stencil node (y)
end

"Data structure for interface point"
struct InterfacePoint
    m::GridPoint # Grid point to negative side of interface
    p::GridPoint # Grid point to positive side of interface
end

"Classify and find stencils for interior grid points"
function find_domain(par, ϕ)
    D = Vector{GridPoint}() # Pre-allocate empty vector of interior grid points
    for i = 2:par.Nx-1 # Loop over grid points in D
        for j = 1:par.Ny # Loop over grid points in D
            if (j == 1) || (j == par.Ny) # Periodic boundary points
                if ϕ[i,j] < 0 # Grid points inside Ω
                    if ϕ[i+1,j] >= 0
                        θxp = ϕ[i,j]/(ϕ[i,j] - ϕ[i+1,j]) # Linear interpolation to find interface
                    else
                        θxp = 1.0
                    end
                    if ϕ[i-1,j] >= 0
                        θxm = ϕ[i,j]/(ϕ[i,j] - ϕ[i-1,j]) # Linear interpolation to find interface
                    else
                        θxm = 1.0
                    end
                    if ϕ[i,2] >= 0
                        θyp = ϕ[i,j]/(ϕ[i,j] - ϕ[i,2]) # Linear interpolation to find interface
                    else
                        θyp = 1.0
                    end
                    if ϕ[i,par.Ny-1] >= 0
                        θym = ϕ[i,j]/(ϕ[i,j] - ϕ[i,par.Ny-1]) # Linear interpolation to find interface
                    else
                        θym = 1.0
                    end
                    if (θxp < par.θb) || (θxm < par.θb) || (θyp < par.θb) || (θym < par.θb) # If point is close to interface
                        push!(D, GridPoint(i, j, true, true, θxp, θyp, θxm, θym))
                    else # If point is far enough from interface
                        push!(D, GridPoint(i, j, true, false, θxp, θyp, θxm, θym))
                    end
                else # Grid points outside Ω
                    if ϕ[i+1,j] < 0
                        θxp = ϕ[i,j]/(ϕ[i,j] - ϕ[i+1,j]) # Linear interpolation to find interface
                    else
                        θxp = 1.0
                    end
                    if ϕ[i-1,j] < 0
                        θxm = ϕ[i,j]/(ϕ[i,j] - ϕ[i-1,j]) # Linear interpolation to find interface
                    else
                        θxm = 1.0
                    end
                    if ϕ[i,2] < 0
                        θyp = ϕ[i,j]/(ϕ[i,j] - ϕ[i,2]) # Linear interpolation to find interface
                    else
                        θyp = 1.0
                    end
                    if ϕ[i,par.Ny-1] < 0
                        θym = ϕ[i,j]/(ϕ[i,j] - ϕ[i,par.Ny-1]) # Linear interpolation to find interface
                    else
                        θym = 1.0
                    end
                    if (θxp < par.θb) || (θxm < par.θb) || (θyp < par.θb) || (θym < par.θb) # If point is close to interface
                        push!(D, GridPoint(i, j, false, true, θxp, θyp, θxm, θym))
                    else # If point is far enough from interface
                        push!(D, GridPoint(i, j, false, false, θxp, θyp, θxm, θym))
                    end
                end
            else # Interior grid points
                if ϕ[i,j] < 0 # Grid points inside Ω
                    if ϕ[i+1,j] >= 0
                        θxp = ϕ[i,j]/(ϕ[i,j] - ϕ[i+1,j]) # Linear interpolation to find interface
                    else
                        θxp = 1.0
                    end
                    if ϕ[i-1,j] >= 0
                        θxm = ϕ[i,j]/(ϕ[i,j] - ϕ[i-1,j]) # Linear interpolation to find interface
                    else
                        θxm = 1.0
                    end
                    if ϕ[i,j+1] >= 0
                        θyp = ϕ[i,j]/(ϕ[i,j] - ϕ[i,j+1]) # Linear interpolation to find interface
                    else
                        θyp = 1.0
                    end
                    if ϕ[i,j-1] >= 0
                        θym = ϕ[i,j]/(ϕ[i,j] - ϕ[i,j-1]) # Linear interpolation to find interface
                    else
                        θym = 1.0
                    end
                    if (θxp < par.θb) || (θxm < par.θb) || (θyp < par.θb) || (θym < par.θb) # If point is close to interface
                        push!(D, GridPoint(i, j, true, true, θxp, θyp, θxm, θym))
                    else # If point is far enough from interface
                        push!(D, GridPoint(i, j, true, false, θxp, θyp, θxm, θym))
                    end
                else # Grid points outside Ω
                    if ϕ[i+1,j] < 0
                        θxp = ϕ[i,j]/(ϕ[i,j] - ϕ[i+1,j]) # Linear interpolation to find interface
                    else
                        θxp = 1.0
                    end
                    if ϕ[i-1,j] < 0
                        θxm = ϕ[i,j]/(ϕ[i,j] - ϕ[i-1,j]) # Linear interpolation to find interface
                    else
                        θxm = 1.0
                    end
                    if ϕ[i,j+1] < 0
                        θyp = ϕ[i,j]/(ϕ[i,j] - ϕ[i,j+1]) # Linear interpolation to find interface
                    else
                        θyp = 1.0
                    end
                    if ϕ[i,j-1] < 0
                        θym = ϕ[i,j]/(ϕ[i,j] - ϕ[i,j-1]) # Linear interpolation to find interface
                    else
                        θym = 1.0
                    end
                    if (θxp < par.θb) || (θxm < par.θb) || (θyp < par.θb) || (θym < par.θb) # If point is close to interface
                        push!(D, GridPoint(i, j, false, true, θxp, θyp, θxm, θym))
                    else # If point is far enough from interface
                        push!(D, GridPoint(i, j, false, false, θxp, θyp, θxm, θym))
                    end
                end
            end
        end
    end
    return D
end

"Locate interface and associated grid points"
function find_interface(par, D, ϕ)
    dΩ = Vector{InterfacePoint}() # Pre-allocate interface
    for gp in D # Loop over all grid points
        if gp.Ω == true # Check points inside domain
            if ϕ[gp.xInd+1, gp.yInd] >= 0 # Check point to right
                for p in D # Find grid point opposite interface
                    if (p.xInd == gp.xInd+1) && (p.yInd == gp.yInd)
                        push!(dΩ, InterfacePoint(gp, p))
                    end
                end 
            end
            if ϕ[gp.xInd-1, gp.yInd] >= 0 # Check point to left
                for m in D # Find grid point opposite interface
                    if (m.xInd == gp.xInd-1) && (m.yInd == gp.yInd)
                        push!(dΩ, InterfacePoint(m, gp))
                    end
                end
            end
            if gp.yInd != par.Ny
                if ϕ[gp.xInd, gp.yInd+1] >= 0 # Check point above
                    for p in D # Find grid point opposite interface
                        if (p.xInd == gp.xInd) && (p.yInd == gp.yInd+1)
                            push!(dΩ, InterfacePoint(gp, p))
                        end
                    end
                end
            end
            if gp.yInd != 1
                if ϕ[gp.xInd, gp.yInd-1] >= 0 # Check point below
                    for m in D # Find grid point opposite interface
                        if (m.xInd == gp.xInd) && (m.yInd == gp.yInd-1)
                            push!(dΩ, InterfacePoint(m, gp))
                        end
                    end
                end
            end
        end
    end
    return dΩ
end