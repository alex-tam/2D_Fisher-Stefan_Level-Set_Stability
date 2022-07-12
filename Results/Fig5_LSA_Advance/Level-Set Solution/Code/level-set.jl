# Solve level-set equations in 2D using the Jiang-Lin-Tadmor method
# Alex Tam, 02/02/2022

"Function to solve level-set equation"
function level_set(V, ϕ, par, dx, dy, dt)
    ϕn = Array{Float64}(undef, par.Nx, par.Ny) # Pre-allocate array of ϕ
    # Compute solutions on staggered grid using Lin-Tadmor
    ϕs = lt_staggered(V, ϕ, par, dx, dy, dt)
    # Compute derivatives on staggered grid
    dxϕs = Array{Float64}(undef, par.Nx-1, par.Ny-1) # Pre-allocate x-derivatives at staggered grid points
    dyϕs = Array{Float64}(undef, par.Nx-1, par.Ny-1) # Pre-allocate y-derivatives at staggered grid points
    for i = 1:par.Nx-1
        for j = 1:par.Ny-1 # Loop over staggered grid points
            # dϕ/dx
            if i == 1
                dxϕs[i,j] = (ϕs[i+1,j]-ϕs[i,j])/dx # ϕ'_{i+1/2,j+1/2}/Δx
            elseif i == par.Nx-1
                dxϕs[i,j] = (ϕs[i,j]-ϕs[i-1,j])/dx # ϕ'_{i+1/2,j+1/2}/Δx
            else
                dxϕs[i,j] = minmod(par.θ*(ϕs[i+1,j]-ϕs[i,j])/dx, (ϕs[i+1,j]-ϕs[i-1,j])/(2*dx), par.θ*(ϕs[i,j]-ϕs[i-1,j])/dx) # ϕ'_{i+1/2,j+1/2}/Δx
            end
            # dϕ/dy
            # Check these BCs, because they are different for the staggered grid
            if j == 1
                # dyϕs[i,j] = (ϕs[i,j+1]-ϕs[i,j])/dy # ϕ`_{i+1/2,j+1/2}/Δy
                dyϕs[i,j] = minmod(par.θ*(ϕs[i,2]-ϕs[i,j])/dy, (ϕs[i,2]-ϕs[i,par.Ny-1])/(2*dy), par.θ*(ϕs[i,j]-ϕs[i,par.Ny-1])/dy) # ϕ`_{i+1/2,j+1/2}/Δy
            elseif j == par.Ny-1
                # dyϕs[i,j] = (ϕs[i,j]-ϕs[i,j-1])/dy # ϕ`_{i+1/2,j+1/2}/Δy
                dyϕs[i,j] = minmod(par.θ*(ϕs[i,1]-ϕs[i,j])/dy, (ϕs[i,1]-ϕs[i,j-1])/(2*dy), par.θ*(ϕs[i,j]-ϕs[i,j-1])/dy) # ϕ`_{i+1/2,j+1/2}/Δy
            else
                dyϕs[i,j] = minmod(par.θ*(ϕs[i,j+1]-ϕs[i,j])/dy, (ϕs[i,j+1]-ϕs[i,j-1])/(2*dy), par.θ*(ϕs[i,j]-ϕs[i,j-1])/dy) # ϕ`_{i+1/2,j+1/2}/Δy
            end
        end
    end
    # Perform cell-averaging at non-staggered points (Jiang 1998)
    for i = 1:par.Nx
        for j = 1:par.Ny # Loop over non-staggered grid points
            if (i == 1) || (i == par.Nx) # Ignore x-boundaries
                ϕn[i,j] = ϕ[i,j]
            elseif (i != 1) && (i != par.Nx) && ((j == 1) || (j == par.Ny)) # Apply periodic conditions on y-boundaries
                # This rule needs to be carefully checked, I think it's correct
                ϕn[i,j] = (ϕs[i,1] + ϕs[i-1,1] + ϕs[i,par.Ny-1] + ϕs[i-1,par.Ny-1])/4 + 
                (dx*(dxϕs[i-1,par.Ny-1]-dxϕs[i,par.Ny-1]) + dx*(dxϕs[i-1,1]-dxϕs[i,1]) + dy*(dyϕs[i-1,par.Ny-1]-dyϕs[i-1,1]) + dy*(dyϕs[i,par.Ny-1]-dyϕs[i,1]))/16 # Compute non-staggered cell average
            else # Interior grid points
                ϕn[i,j] = (ϕs[i,j] + ϕs[i-1,j] + ϕs[i,j-1] + ϕs[i-1,j-1])/4 + 
                (dx*(dxϕs[i-1,j-1]-dxϕs[i,j-1]) + dx*(dxϕs[i-1,j]-dxϕs[i,j]) + dy*(dyϕs[i-1,j-1]-dyϕs[i-1,j]) + dy*(dyϕs[i,j-1]-dyϕs[i,j]))/16 # Compute non-staggered cell average
            end
        end
    end
    return ϕn
end

"Staggered Lin-Tadmor scheme"
function lt_staggered(V, ϕ, par, dx, dy, dt)
    ϕ_stag = Array{Float64}(undef, par.Nx-1, par.Ny-1) # Pre-allocate solution on staggered grid ϕ_stag[i,j] = ϕ[i+1/2, j+1/2], one-based indexing
    V_stag = Array{Float64}(undef, par.Nx-1, par.Ny-1) # Pre-allocate velocity on staggered grid V_stag[i,j] = ϕ[i+1/2, j+1/2], one-based indexing
    dxϕ = Array{Float64}(undef, par.Nx, par.Ny) # Pre-allocate x-derivatives at non-staggered grid points
    dyϕ = Array{Float64}(undef, par.Nx, par.Ny) # Pre-allocate y-derivatives at non-staggered grid points
    ϕ_mid = Array{Float64}(undef, par.Nx, par.Ny) # Pre-allocate mid-values at non-staggered grid points
    # Evaluate derivatives and missing mid-values on full non-staggered grid
    for i = 1:par.Nx
        for j = 1:par.Ny
            # Compute dϕ/dx
            if i == 1 # Boundary condition
                dxϕ[i,j] = (ϕ[i+1,j] - ϕ[i,j])/dx
            elseif i == par.Nx # Boundary condition
                dxϕ[i,j] = (ϕ[i,j] - ϕ[i-1,j])/dx
            else # Apply min-mod at interior grid points
                dxϕ[i,j] = minmod(par.θ*(ϕ[i+1,j]-ϕ[i,j])/dx, (ϕ[i+1,j]-ϕ[i-1,j])/(2*dx), par.θ*(ϕ[i,j]-ϕ[i-1,j])/dx)
            end
            # Compute dϕ/dy
            if j == 1 # Boundary condition
                # dyϕ[i,j] = (ϕ[i,j+1] - ϕ[i,j])/dy
                dyϕ[i,j] = minmod(par.θ*(ϕ[i,j+1]-ϕ[i,j])/dy, (ϕ[i,j+1]-ϕ[i,par.Ny-1])/(2*dy), par.θ*(ϕ[i,j]-ϕ[i,par.Ny-1])/dy)
            elseif j == par.Ny # Boundary condition
                # dyϕ[i,j] = (ϕ[i,j] - ϕ[i,j-1])/dy
                dyϕ[i,j] = minmod(par.θ*(ϕ[i,2]-ϕ[i,j])/dy, (ϕ[i,2]-ϕ[i,j-1])/(2*dy), par.θ*(ϕ[i,j]-ϕ[i,j-1])/dy)
            else # Apply min-mod at interior grid points
                dyϕ[i,j] = minmod(par.θ*(ϕ[i,j+1]-ϕ[i,j])/dy, (ϕ[i,j+1]-ϕ[i,j-1])/(2*dy), par.θ*(ϕ[i,j]-ϕ[i,j-1])/dy)
            end
            ϕ_mid[i,j] = ϕ[i,j] - dt/2*hamiltonian_ls(V[i,j], dxϕ[i,j], dyϕ[i,j]) # Missing mid-values
        end
    end
    # Obtain velocity field on staggered grid
    for i = 1:par.Nx-1
        for j = 1:par.Ny-1 # Loop over staggered grid points
            V_stag[i,j] = (V[i,j] + V[i+1,j] + V[i,j+1] + V[i+1,j+1])/4
        end
    end
    # Evaluate solution on staggered grid
    for i = 1:par.Nx-1
        for j = 1:par.Ny-1 # Loop over staggered grid points
            ϕ_stag[i,j] = (ϕ[i,j] + ϕ[i+1,j] + ϕ[i,j+1] + ϕ[i+1,j+1])/4 + 
            (dxϕ[i,j]-dxϕ[i+1,j]+dxϕ[i,j+1]-dxϕ[i+1,j+1])*dx/16 + 
            (dyϕ[i,j]-dyϕ[i,j+1]+dyϕ[i+1,j]-dyϕ[i+1,j+1])*dy/16 - 
            dt/2*(hamiltonian_ls(V_stag[i,j], (ϕ_mid[i+1,j]-ϕ_mid[i,j])/dx, (ϕ_mid[i+1,j+1]-ϕ_mid[i+1,j])/dy) + 
            hamiltonian_ls(V_stag[i,j], (ϕ_mid[i+1,j+1]-ϕ_mid[i,j+1])/dx, (ϕ_mid[i,j+1]-ϕ_mid[i,j])/dy))
        end
    end
    return ϕ_stag
end

"Minmod function for flux-limiter"
function minmod(a, b, c)
    if (a < 0) && (b < 0) && (c < 0)
        return max(a, b, c)
    elseif (a > 0) && (b > 0) && (c > 0)
        return min(a, b, c)
    else
        return 0.0
    end
end 

"Evaluate Hamiltonian at one grid point"
function hamiltonian_ls(v, ϕx, ϕy)
    return v*sqrt(ϕx^2 + ϕy^2)
end