# Solve reinitialisation equation to reset ϕ as signed-distance function
# Alex Tam, 02/02/2022

function reinitialisation(ϕ, par, dx, dy, nI)
    dτ = 0.5*(dx/5 + dy/5)
    for i = 1:nI
        ϕ = reinit_lt(ϕ, par, dx, dy, dτ)
    end
    return ϕ
end

"Compute one τ-step of the reinitialisation equation using Lin-Tadmor"
function reinit_lt(ϕ, par, dx, dy, dt)
    ϕn = Array{Float64}(undef, par.Nx, par.Ny) # Pre-allocate array of ϕ
    dxϕ = Array{Float64}(undef, par.Nx, par.Ny) # Pre-allocate x-derivatives at non-staggered grid points
    dyϕ = Array{Float64}(undef, par.Nx, par.Ny) # Pre-allocate y-derivatives at non-staggered grid points
    S = Array{Float64}(undef, par.Nx, par.Ny) # Pre-allocate array of S
    # Compute derivatives at non-staggered points at τ_n
    for i = 1:par.Nx
        for j = 1:par.Ny # Loop over grid points
            # dϕ/dx
            if i == 1 # Boundary condition
                dxϕ[i,j] = (ϕ[i+1,j] - ϕ[i,j])/dx 
            elseif i == par.Nx # Boundary condition
                dxϕ[i,j] = (ϕ[i,j] - ϕ[i-1,j])/dx
            else # Apply min-mod at interior grid points
                dxϕ[i,j] = minmod(par.θ*(ϕ[i+1,j]-ϕ[i,j])/dx, (ϕ[i+1,j]-ϕ[i-1,j])/(2*dx), par.θ*(ϕ[i,j]-ϕ[i-1,j])/dx)
            end
            # dϕ/dy
            if j == 1
                dyϕ[i,j] = minmod(par.θ*(ϕ[i,j+1]-ϕ[i,j])/dy, (ϕ[i,j+1]-ϕ[i,par.Ny-1])/(2*dy), par.θ*(ϕ[i,j]-ϕ[i,par.Ny-1])/dy) # ϕ`_{i+1/2,j+1/2}/Δy
            elseif j == par.Ny
                dyϕ[i,j] = minmod(par.θ*(ϕ[i,2]-ϕ[i,j])/dy, (ϕ[i,2]-ϕ[i,j-1])/(2*dy), par.θ*(ϕ[i,j]-ϕ[i,j-1])/dy) # ϕ`_{i+1/2,j+1/2}/Δy
            else
                dyϕ[i,j] = minmod(par.θ*(ϕ[i,j+1]-ϕ[i,j])/dy, (ϕ[i,j+1]-ϕ[i,j-1])/(2*dy), par.θ*(ϕ[i,j]-ϕ[i,j-1])/dy) # ϕ`_{i+1/2,j+1/2}/Δy
            end
        end
    end
    # Compute sign function on grid points at τ_n
    for i = 1:par.Nx
        for j = 1:par.Ny
            S[i,j] = mod_sign(ϕ[i,j], dxϕ[i,j], dyϕ[i,j], dx)
        end
    end
    # Compute solutions on staggered grid
    ϕs = lt_staggered_reinit(S, ϕ, dxϕ, dyϕ, par, dx, dy, dt)
    # Compute derivatives on staggered grid
    dxϕs = Array{Float64}(undef, par.Nx-1, par.Ny-1) # Pre-allocate x-derivatives at staggered grid points
    dyϕs = Array{Float64}(undef, par.Nx-1, par.Ny-1) # Pre-allocate y-derivatives at staggered grid points
    for i = 1:par.Nx-1
        for j = 1:par.Ny-1
            if i == 1
                dxϕs[i,j] = (ϕs[i+1,j]-ϕs[i,j])/dx # ϕ'_{i+1/2,j+1/2}/Δx
            elseif i == par.Nx-1
                dxϕs[i,j] = (ϕs[i,j]-ϕs[i-1,j])/dx # ϕ'_{i+1/2,j+1/2}/Δx
            else
                dxϕs[i,j] = minmod(par.θ*(ϕs[i+1,j]-ϕs[i,j])/dx, (ϕs[i+1,j]-ϕs[i-1,j])/(2*dx), par.θ*(ϕs[i,j]-ϕs[i-1,j])/dx) # ϕ'_{i+1/2,j+1/2}/Δx
            end
            if j == 1
                dyϕs[i,j] = minmod(par.θ*(ϕs[i,j+1]-ϕs[i,j])/dy, (ϕs[i,j+1]-ϕs[i,par.Ny-1])/(2*dy), par.θ*(ϕs[i,j]-ϕs[i,par.Ny-1])/dy) # ϕ`_{i+1/2,j+1/2}/Δy
            elseif j == par.Ny-1
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

"Implement staggered Lin-Tadmor scheme"
function lt_staggered_reinit(S, ϕ, dxϕ, dyϕ, par, dx, dy, dt)
    ϕ_stag = Array{Float64}(undef, par.Nx-1, par.Ny-1) # Pre-allocate solution on staggered grid ϕ_stag[i,j] = ϕ[i+1/2, j+1/2], one-based indexing
    S_stag = Array{Float64}(undef, par.Nx-1, par.Ny-1) # Pre-allocate velocity on staggered grid V_stag[i,j] = ϕ[i+1/2, j+1/2], one-based indexing
    ϕ_mid = Array{Float64}(undef, par.Nx, par.Ny) # Pre-allocate mid-values at non-staggered grid points
    # Evaluate missing mid-values on full non-staggered grid
    for i = 1:par.Nx
        for j = 1:par.Ny  
            ϕ_mid[i,j] = ϕ[i,j] - dt/2*hamiltonian_reinit(S[i,j], dxϕ[i,j], dyϕ[i,j]) # Missing mid-values
        end
    end
    # Obtain sign function on staggered grid
    for i = 1:par.Nx-1
        for j = 1:par.Ny-1
            S_stag[i,j] = (S[i,j] + S[i+1,j] + S[i,j+1] + S[i+1,j+1])/4
        end
    end
    # Evaluate solution on staggered grid
    for i = 1:par.Nx-1
        for j = 1:par.Ny-1
            ϕ_stag[i,j] = (ϕ[i,j] + ϕ[i+1,j] + ϕ[i,j+1] + ϕ[i+1,j+1])/4 + 
            (dxϕ[i,j]-dxϕ[i+1,j]+dxϕ[i,j+1]-dxϕ[i+1,j+1])*dx/16 + 
            (dyϕ[i,j]-dyϕ[i,j+1]+dyϕ[i+1,j]-dyϕ[i+1,j+1])*dy/16 - 
            dt/2*(hamiltonian_reinit(S_stag[i,j], (ϕ_mid[i+1,j]-ϕ_mid[i,j])/dx, (ϕ_mid[i+1,j+1]-ϕ_mid[i+1,j])/dy) + 
            hamiltonian_reinit(S_stag[i,j], (ϕ_mid[i+1,j+1]-ϕ_mid[i,j+1])/dx, (ϕ_mid[i,j+1]-ϕ_mid[i,j])/dy))
        end
    end
    return ϕ_stag
end

"Implement the modified sign function"
function mod_sign(ϕ, ϕx, ϕy, dx)
    return ϕ/sqrt(ϕ^2 + dx^2*(ϕx^2 + ϕy^2))
end

"Hamiltonian for reinitialisation equation"
function hamiltonian_reinit(s, ϕx, ϕy)
    return s*(sqrt(ϕx^2 + ϕy^2) - 1)
end