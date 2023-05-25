# Apply Stefan condition to calculate interface speed
# Alex Tam, 02/02/2022

"Apply Stefan condition to find speed at interface points"
function interface_speed(dΩ, U, ϕ, par, dx, dy)
    uf = interface_density(dΩ, ϕ, par, dx, dy) # Obtain interface density
    vi = Vector{Float64}() # Pre-allocate empty vector of interface speed
    for i in eachindex(dΩ)
        if dΩ[i].m.xInd != dΩ[i].p.xInd # Interface lies in x-direction
            ux = find_ux_x(dΩ, dΩ[i], uf, uf[i], U, ϕ, par, dx)
            uy = find_uy_x(dΩ[i], uf[i], U, ϕ, par, dx, dy)
            ϕx = find_ϕx_x(dΩ[i], ϕ, dx)
            ϕy = find_ϕy_x(dΩ[i], ϕ, par, dy)
        else # Interface lies in y-direction
            ux = find_ux_y(dΩ[i], uf[i], U, ϕ, par, dx, dy)
            uy = find_uy_y(dΩ, dΩ[i], uf, uf[i], U, ϕ, par, dy)
            ϕx = find_ϕx_y(dΩ[i], ϕ, dx)
            ϕy = find_ϕy_y(dΩ[i], ϕ, par, dy)
        end
        push!(vi, -par.κ*(ux*ϕx + uy*ϕy))
    end
    return vi
end

"Compute du/dx for interface in x-direction"
function find_ux_x(dΩ, ip, ufv, uf, U, ϕ, par, dx)
    if ϕ[ip.m.xInd, ip.m.yInd] < 0 # If left point lies inside Ω
        i = ip.m.xInd; j = ip.m.yInd; θE = ip.m.θxp; θW = ip.m.θxm # Left point is reference index
        if θE >= par.θb # If left point is not close to dΩ, apply non-standard finite-difference
            if θW == 1.0 # If (x_{i-1}, y_j) is inside Ω
                return (1/dx)*( (2*θE+θW)*uf/(θE*(θW+θE)) - (θW + θE)*U[i,j]/(θW*θE) + θE*U[i-1,j]/(θW*(θW+θE)) )
            else # If (x_{i-1}, y_j) is not inside Ω
                ufg = get_interface_density(i-1, j, i, j, dΩ, ufv)
                return (1/dx)*( (2*θE+θW)*uf/(θE*(θW+θE)) - (θW + θE)*U[i,j]/(θW*θE) + θE*ufg/(θW*(θW+θE)) )
            end
        else # If left point is close to dΩ, apply one-sided difference using (x_i,y_j) as right node, U_ij = 0
            if (ϕ[i-1,j] < 0) && (ϕ[i-2,j] < 0) # If two points to left of interface are in Ω
                return 1/(2*dx)*( 3*U[i,j] - 4*U[i-1,j] + U[i-2,j])
            elseif (ϕ[i-1,j] < 0) && (ϕ[i-2,j] >= 0) # If only one point to left of interface are in Ω
                θg = ϕ[i-1,j]/(ϕ[i-1,j]-ϕ[i-2,j])
                ufg = get_interface_density(i-2, j, i-1, j, dΩ, ufv)
                if θg < par.θb # If opposite ghost node is also close to interface
                    return (uf - ufg)/((1+θg+θE)*dx) # First-order approximation
                else 
                    return (1/dx)*( (uf*(2+θg))/(1+θg) - (U[i-1,j]*(1+θg))/θg + ufg/(θg*(1+θg)) )
                end
            else # No points to left of interface are in Ω
                θg = ϕ[i,j]/(ϕ[i,j]-ϕ[i-1,j])
                ufg = get_interface_density(i-1, j, i, j, dΩ, ufv)
                return (uf - ufg)/((θE+θg)*dx) # First-order aproximation
            end
        end
    else # If right point lies inside Ω (should never occur for perturbed planar fronts)
        i = ip.p.xInd; j = ip.p.yInd; θW = ip.p.θxm; θE = ip.p.θxp # Right point is reference index
        if θW >= par.θb # If right point is not close to dΩ, apply non-standard finite-difference
            if θE == 1.0 # If (x_{i+1}, y_j) is inside Ω
                return (1/dx)*( -uf*(2*θW+θE)/(θW*(θW+θE)) + (θW+θE)*U[i,j]/(θW*θE) - θW*U[i+1,j]/(θE*(θW+θE)) )
            else # If (x_{i+1}, y_j) is not inside Ω
                ufg = get_interface_density(i, j, i+1, j, dΩ, ufv)
                return (1/dx)*( -uf*(2*θW+θE)/(θW*(θW+θE)) + (θW+θE)*U[i,j]/(θW*θE) - θW*ufg/(θE*(θW+θE)) )
            end
        else # If right point is close to dΩ, apply one-sided difference using U_{i,j} as left node
            if (ϕ[i+1,j] < 0) && (ϕ[i+2,j] < 0) # If two points to right of interface are in Ω
                return 1/(2*dx)*( -3*U[i,j] + 4*U[i+1,j] - U[i+2,j])
            elseif (ϕ[i+1,j] < 0) && (ϕ[i+2,j] >= 0)  # If only one point to right of interface is in Ω
                θg = ϕ[i+1,j]/(ϕ[i+1,j]-ϕ[i+2,j])
                ufg = get_interface_density(i+1, j, i+2, j, dΩ, ufv)
                if θg < par.θb # If opposite ghost node is also close to interface
                    return (ufg - uf)/((1+θg+θW)*dx) # First-order approximation
                else
                    return (1/dx)*( -(uf*(2+θg))/(1+θg) + (U[i+1,j]*(1+θg))/θg - ufg/(θg*(1+θg)) )
                end
            else # If point to right of reference is not in Ω
                θg = ϕ[i,j]/(ϕ[i,j]-ϕ[i+1,j])
                ufg = get_interface_density(i, j, i+1, j, dΩ, ufv)
                return (ufg - uf)/((θW+θg)*dx) # First-order aproximation
            end 
        end
    end
end

"Compute du/dy for interface in x-direction"
function find_uy_x(ip, uf, U, ϕ, par, dx, dy)
    i = ip.m.xInd; j = ip.m.yInd; θE = ip.m.θxp # Left point is reference index
    # Upper ghost point
    if j != par.Ny # Interior nodes
        ϕrp = θE*ϕ[i+1,j+1] + (1-θE)*ϕ[i,j+1] # ϕ above interface point
        urp = θE*U[i+1,j+1] + (1-θE)*U[i,j+1] # U above interface point
    else # Apply periodic BCs
        ϕrp = θE*ϕ[i+1,2] + (1-θE)*ϕ[i,2] # ϕ above interface point
        urp = θE*U[i+1,2] + (1-θE)*U[i,2] # U above interface point
    end
    # Lower ghost point
    if j != 1 # Interior nodes
        ϕrm = θE*ϕ[i+1,j-1] + (1-θE)*ϕ[i,j-1] # ϕ below interface point
        urm = θE*U[i+1,j-1] + (1-θE)*U[i,j-1] # U below interface point
    else # Apply periodic BCs
        ϕrm = θE*ϕ[i+1,par.Ny-1] + (1-θE)*ϕ[i,par.Ny-1] # ϕ below interface point
        urm = θE*U[i+1,par.Ny-1] + (1-θE)*U[i,par.Ny-1] # U below interface point
    end
    # Calculate derivative du/dy
    if (ϕrp < 0) && (ϕrm < 0) # If ghost nodes directly above and below interface are both in Ω
        return (urp - urm)/(2*dy)
    elseif (ϕrp < 0) && (ϕrm >= 0) # If ghost node above only is in Ω
        # Apply periodic BC if necessary
        if (j != par.Ny) && (j != par.Ny-1) # Interior nodes
            ϕrpp = θE*ϕ[i+1,j+2] + (1-θE)*ϕ[i,j+2] # ϕ two points above interface point
            urpp = θE*U[i+1,j+2] + (1-θE)*U[i,j+2] # U two points above interface point
        elseif j == par.Ny # Apply periodic BC
            ϕrpp = θE*ϕ[i+1,3] + (1-θE)*ϕ[i,3] # ϕ two points above interface point
            urpp = θE*U[i+1,3] + (1-θE)*U[i,3] # U two points above interface point
        else # Apply periodic BC
            ϕrpp = θE*ϕ[i+1,2] + (1-θE)*ϕ[i,2] # ϕ two points above interface point
            urpp = θE*U[i+1,2] + (1-θE)*U[i,2] # U two points above interface point
        end
        # Compute derivative
        if ϕrpp < 0
            return (-3*uf + 4*urp - urpp)/(2*dy) # Three-point one-sided difference
        else
            θg = ϕrp/(ϕrp-ϕrpp) # θ for second interface ghost point
            if (j != par.Ny) && (j!= par.Ny-1)
                K = get_curvature_ghost(i, j+1, θE, θg, ϕ, par, dx, dy)
            elseif j == par.Ny-1 # Apply periodic BC if j = Ny-1
                K = get_curvature_ghost(i, 1, θE, θg, ϕ, par, dx, dy)
            else # Apply periodic BC if j = Ny
                K = get_curvature_ghost(i, 2, θE, θg, ϕ, par, dx, dy)
            end
            ufg = par.uf - par.γ*K # Replace with curvature at ghost point later
            return (1/dy)*( -(2+θg)*uf/(1+θg) + (1+θg)*urp/θg - ufg/(θg*(1+θg)) ) # Non-standard one-sided difference
        end
    elseif (ϕrp >= 0) && (ϕrm < 0) # If ghost node below only is in Ω
        # Apply periodic BC if necessary
        if (j != 1) && (j != 2) # Interior nodes
            ϕrmm = θE*ϕ[i+1,j-2] + (1-θE)*ϕ[i,j-2] # ϕ two points below interface point
            urmm = θE*U[i+1,j-2] + (1-θE)*U[i,j-2] # U two points below interface point
        elseif j == 1 # Apply periodic BC
            ϕrmm = θE*ϕ[i+1,par.Ny-2] + (1-θE)*ϕ[i,par.Ny-2] # ϕ two points below interface point
            urmm = θE*U[i+1,par.Ny-2] + (1-θE)*U[i,par.Ny-2] # U two points below interface point
        else # Apply periodic BC
            ϕrmm = θE*ϕ[i+1,par.Ny-1] + (1-θE)*ϕ[i,par.Ny-1] # ϕ two points below interface point
            urmm = θE*U[i+1,par.Ny-1] + (1-θE)*U[i,par.Ny-1] # U two points below interface point
        end
        # Compute derivative
        if ϕrmm < 0
            return (3*uf - 4*urm + urmm)/(2*dy) # Three-point one-sided difference
        else
            θg = ϕrm/(ϕrm-ϕrmm) # θ for second interface ghost point
            if (j != 1) && (j != 2) # Interior points
                K = get_curvature_ghost(i, j-2, θE, 1-θg, ϕ, par, dx, dy)
            elseif j == 2 # Periodic BC for j = 2
                K = get_curvature_ghost(i, par.Ny-1, θE, 1-θg, ϕ, par, dx, dy)
            else # Periodic BC for j = 1
                K = get_curvature_ghost(i, par.Ny-2, θE, 1-θg, ϕ, par, dx, dy)
            end
            ufg = par.uf - par.γ*K # Replace with curvature at ghost point later
            return (1/dy)*( (2+θg)*uf/(1+θg) - (1+θg)*urm/θg + ufg/(θg*(1+θg)) ) # Non-standard one-sided difference
        end
    else # If ghost nodes above and below interface are both not in Ω
        return 0.0
    end
end

"Compute du/dx for interface in y-direction"
function find_ux_y(ip, uf, U, ϕ, par, dx, dy)
    i = ip.m.xInd; j = ip.m.yInd; θN = ip.m.θyp # Bottom point is reference index
    # Left and right ghost points
    ϕrp = θN*ϕ[i+1,j+1] + (1-θN)*ϕ[i+1,j] # ϕ to right of interface point
    urp = θN*U[i+1,j+1] + (1-θN)*U[i+1,j] # U to right of interface point
    ϕrm = θN*ϕ[i-1,j+1] + (1-θN)*ϕ[i-1,j] # ϕ to left of interface point
    urm = θN*U[i-1,j+1] + (1-θN)*U[i-1,j] # ϕ to left of interface point
    # Calculate derivative du/dx
    if (ϕrp < 0) && (ϕrm < 0) # If ghost nodes to left and right of interface are both in Ω (should never occur for perturbed planar fronts)
        return (urp-urm)/(2*dx) # Central difference
    elseif (ϕrp < 0) && (ϕrm >= 0) # If ghost node to right only is in Ω (should never occur for perturbed planar fronts)
        ϕrpp = θN*ϕ[i+2,j+1] + (1-θN)*ϕ[i+2,j] # ϕ two points to right of interface point
        urpp = θN*U[i+2,j+1] + (1-θN)*U[i+2,j] # U two points to right of interface point
        if ϕrpp < 0 # If point two to the right of interface is in Ω
            return (-3*uf + 4*urp - urpp)/(2*dx) # One-sided difference
        else
            θg = ϕrp/(ϕrp-ϕrpp) # θ for second interface ghost point
            K = get_curvature_ghost(i+1, j, θg, θN, ϕ, par, dx, dy)
            ufg = par.uf - par.γ*K
            return (1/dx)*( -(2+θg)*uf/(1+θg) + (1+θg)*urp/θg - ufg/(θg*(1+θg)) ) # Non-standard one-sided difference
        end
    elseif (ϕrm < 0) && (ϕrp >= 0) # If ghost node to left only is in Ω (should always occur for perturbed planar fronts)
        ϕrmm = θN*ϕ[i-2,j+1] + (1-θN)*ϕ[i-2,j] # ϕ two points to left of interface point
        urmm = θN*U[i-2,j+1] + (1-θN)*U[i-2,j] # U two points to left of interface point
        if ϕrmm < 0 # If point two to the left of interface is in Ω (should always occur for perturbed planar fronts)
            return (3uf - 4*urm + urmm)/(2*dx) # One-sided difference
        else
            θg = ϕrm/(ϕrm-ϕrmm) # θ for second interface ghost point
            K = get_curvature_ghost(i-2, j, 1-θg, θN, ϕ, par, dx, dy)
            ufg = par.uf - par.γ*K # Replace with curvature at ghost point later
            return (1/dx)*( (2+θg)*uf/(1+θg) - (1+θg)*urm/θg + ufg/(θg*(1+θg)) ) # Non-standard one-sided difference
        end
    else # If ghost nodes to left and right of interface are both not in Ω (should never occur for perturbed planar fronts)
        return 0.0
    end
end

"Compute du/dy for interface in y-direction"
function find_uy_y(dΩ, ip, ufv, uf, U, ϕ, par, dy)
    if ϕ[ip.m.xInd, ip.m.yInd] < 0 # If bottom point is inside Ω
        i = ip.m.xInd; j = ip.m.yInd; θN = ip.m.θyp; θS = ip.m.θym # Bottom point is reference index
        if θN >= par.θb # If bottom point is not close to dΩ, apply non-standard finite-difference
            if θS == 1.0 # If (x_i, y_{j-1}) is inside Ω
                if j == 1 # Apply periodic BC if necessary
                    us = U[i-1,par.Ny-1]
                else # Interior grid points
                    us = U[i,j-1]
                end
                return (1/dy)*( (2*θN+θS)*uf/(θN*(θS+θN)) - (θS+θN)*U[i,j]/(θS*θN) + θN*us/(θS*(θS+θN)) )
            else # If (x_i, y_{j-1}) is outside Ω
                if j != 1
                    ufg = get_interface_density(i, j-1, i, j, dΩ, ufv)
                else # Apply periodic BC on j = 1
                    ufg = get_interface_density(i, par.Ny-1, i, par.Ny, dΩ, ufv)
                end
                return (1/dy)*( (2*θN+θS)*uf/(θN*(θS+θN)) - (θS+θN)*U[i,j]/(θS*θN) + θN*ufg/(θS*(θS+θN)) )
            end
        else # If bottom point is close to dΩ, apply one-sided difference using (x_i, y_j) as top node, U_ij = 0
            # Apply periodic BCs if necessary
            if (j != 1) && (j != 2) # Interior grid points
                ϕm = ϕ[i,j-1]; um = U[i,j-1]
                ϕmm = ϕ[i,j-2]; umm = U[i,j-2]
            elseif j == 2 # Apply periodic BC for j = 2
                ϕm = ϕ[i,1]; um = U[i,1]
                ϕmm = ϕ[i,par.Ny-1]; umm = U[i,par.Ny-1]
            else # Apply periodic BC for j = 1
                ϕm = ϕ[i,par.Ny-1]; um = U[i,par.Ny-1]
                ϕmm = ϕ[i,par.Ny-2]; umm = U[i,par.Ny-2]
            end
            # Compute derivative
            if (ϕm < 0) && (ϕmm < 0) # Standard finite difference if both lower points are in Ω
                return 1/(2*dy)*(3*uf - 4*um + umm)
            elseif (ϕm < 0) && (ϕmm >= 0) # Only one lower point is in Ω
                θg = ϕm/(ϕm - ϕmm) # θ for second ghost node
                # Obtain interface density
                if (j != 1) && (j != 2)
                    ufg = get_interface_density(i, j-2, i, j-1, dΩ, ufv)
                elseif j == 2 # Apply periodic BC for j = 2
                    ufg = get_interface_density(i, par.Ny-1, i, par.Ny, dΩ, ufv)
                else # Apply periodic BC for j = 1
                    ufg = get_interface_density(i, par.Ny-2, i, par.Ny-1, dΩ, ufv)
                end
                # Compute derivative
                if θg < par.θb
                    return (uf-ufg)/((1+θg+θN)*dy) # First-order approximation
                else
                    return (1/dy)*( (uf*(2+θg))/(1+θg) - (um*(1+θg))/θg + ufg/(θg*(1+θg)) )
                end
            else # No lower points are in Ω
                θg = ϕ[i,j]/(ϕ[i,j]-ϕm)
                if j != 1
                    ufg = get_interface_density(i, j-1, i, j, dΩ, ufv) # Periodic BC
                else # Apply periodic BC on j = 1
                    ufg = get_interface_density(i, par.Ny-1, i, par.Ny, dΩ, ufv) # Periodic BC
                end
                return (uf - ufg)/((θN+θg)*dy) # First-order aproximation  
            end
        end
    else # If top point is inside Ω
        i = ip.p.xInd; j = ip.p.yInd; θN = ip.p.θyp; θS = ip.p.θym # Top point is reference index
        if θS >= par.θb # If interface is not close to top point, apply non-standard one-sided difference
            if θN == 1.0 # If (x_i, y_{j+1}) is in Ω
                if j != par.Ny # Interior points
                    un = U[i,j+1]
                else # Apply periodic BC if necessary
                    un = U[i,2]
                end
                return (1/dy)*( -uf*(2*θS+θN)/(θS*(θS+θN)) + (θS+θN)*U[i,j]/(θS*θN) - θS*un/(θN*(θN+θS)) )
            else # If (x_i, y_{j+1}) is not in Ω
                if j != par.Ny
                    ufg = get_interface_density(i, j, i, j+1, dΩ, ufv)
                else # Apply periodic BC on j = Ny
                    ufg = get_interface_density(i, 1, i, 2, dΩ, ufv)
                end
                return (1/dy)*( -uf*(2*θS+θN)/(θS*(θS+θN)) + (θS+θN)*U[i,j]/(θS*θN) - θS*ufg/(θN*(θN+θS)) )
            end
        else # If interface is close to top point, apply one-sided difference with (x_i,y_j) as bottom point 
            # Apply periodic BCs if necessary
            if (j != par.Ny) && (j != par.Ny-1) # Interior grid points
                ϕp = ϕ[i,j+1]; up = U[i,j+1]
                ϕpp = ϕ[i,j+2]; upp = U[i,j+2]
            elseif j == par.Ny-1 # Apply periodic BC near boundary
                ϕp = ϕ[i,par.Ny]; up = U[i,par.Ny]
                ϕpp = ϕ[i,2]; upp = U[i,2]
            else # Apply periodic BC at boundary
                ϕp = ϕ[i,2]; up = U[i,2]
                ϕpp = ϕ[i,3]; upp = U[i,3]
            end
            # Compute derivative
            if (ϕp < 0) && (ϕpp < 0) # If both nodes above interface are in Ω
                return (-3*uf + 4*up - upp)/(2*dy)
            elseif (ϕp < 0) && (ϕpp >= 0) # Only one upper node in Ω
                θg = ϕp/(ϕp - ϕpp) # θ for second ghost node
                # Obtain interface density
                if (j != par.Ny) && (j != par.Ny-1)
                    ufg = get_interface_density(i, j+1, i, j+2, dΩ, ufv)
                elseif j == par.Ny-1 # Apply periodic BC on j = Ny-1
                    ufg = get_interface_density(i, 1, i, 2, dΩ, ufv)
                else # Apply periodic BC on j = Ny
                    ufg = get_interface_density(i, 2, i, 3, dΩ, ufv)
                end
                # Compute derivative
                if θg < par.θb
                    return (ufg-uf)/((1+θg+θS)*dy) # First-order approximation
                else
                    return (1/dy)*( -(uf*(2+θg))/(1+θg) + (up*(1+θg))/θg - ufg/(θg*(1+θg)) )
                end
            else # Neither upper node in Ω
                θg = ϕ[i,j]/(ϕ[i,j]-ϕp)
                if j != par.Ny # Interior points
                    ufg = get_interface_density(i, j, i, j+1, dΩ, ufv)
                else # Apply periodic BC on j = Ny
                    ufg = get_interface_density(i, 1, i, 2, dΩ, ufv)
                end
                return (ufg - uf)/((θS+θg)*dy) # First-order aproximation
            end
        end
    end
end

"Compute dϕ/dx for interface in x-direction"
function find_ϕx_x(ip, ϕ, dx)
    i = ip.m.xInd; j = ip.m.yInd; θE = ip.m.θxp # Left point is reference index
    ϕxl = (ϕ[i+1,j] - ϕ[i-1,j])/(2*dx)
    ϕxr = (ϕ[i+2,j] - ϕ[i,j])/(2*dx) # Assume far enough from x-boundary
    return θE*ϕxr + (1-θE)*ϕxl
end

"Compute dϕ/dy for interface in x-direction"
function find_ϕy_x(ip, ϕ, par, dy)
    i = ip.m.xInd; j = ip.m.yInd; θE = ip.m.θxp # Left point is reference index
    if (j == 1) || (j == par.Ny) # Boundary points
        ϕrp = θE*ϕ[i+1,2] + (1-θE)*ϕ[i,2]
        ϕrm = θE*ϕ[i+1,par.Ny-1] + (1-θE)*ϕ[i,par.Ny-1] 
    else # Interior points
        ϕrp = θE*ϕ[i+1,j+1] + (1-θE)*ϕ[i,j+1]
        ϕrm = θE*ϕ[i+1,j-1] + (1-θE)*ϕ[i,j-1]     
    end
    return (ϕrp - ϕrm)/(2*dy)
end

"Compute dϕ/dx for interface in y-direction"
function find_ϕx_y(ip, ϕ, dx)
    i = ip.m.xInd; j = ip.m.yInd; θN = ip.m.θyp # Bottom point is reference index
    ϕrp = θN*ϕ[i+1,j+1] + (1-θN)*ϕ[i+1,j]
    ϕrm = θN*ϕ[i-1,j+1] + (1-θN)*ϕ[i-1,j]
    return (ϕrp - ϕrm)/(2*dx)
end

"Compute dϕ/dy for interface in y-direction"
function find_ϕy_y(ip, ϕ, par, dy)
    i = ip.m.xInd; j = ip.m.yInd; θN = ip.m.θyp # Bottom point is reference index
    # Apply periodic conditions for ϕ at boundaries
    if j == 1 # Interface near bottom boundary
        ϕyt = (ϕ[i,3] - ϕ[i,1])/(2*dy)
        ϕyb = (ϕ[i,2] - ϕ[i,par.Ny-1])/(2*dy)
    elseif j == par.Ny-1 # Interface near top boundary
        ϕyt = (ϕ[i,2] - ϕ[i,par.Ny-1])/(2*dy)
        ϕyb = (ϕ[i,par.Ny] - ϕ[i,par.Ny-2])/(2*dy)
    else # Interior nodes
        ϕyt = (ϕ[i,j+2] - ϕ[i,j])/(2*dy)
        ϕyb = (ϕ[i,j+1] - ϕ[i,j-1])/(2*dy)
    end
    return θN*ϕyt + (1-θN)*ϕyb
end

"Compute curvature at a ghost point using two-dimensional interpolation"
function get_curvature_ghost(i, j, θE, θN, ϕ, par, dx, dy) # Let (i,j) be indices of bottom left point in rectangle
    # ϕx
    ϕxbl = (ϕ[i+1,j] - ϕ[i-1,j])/(2*dx)
    ϕxbr = (ϕ[i+2,j] - ϕ[i,j])/(2*dx)
    ϕxtl = (ϕ[i+1,j+1] - ϕ[i-1,j+1])/(2*dx)
    ϕxtr = (ϕ[i+2,j+1] - ϕ[i,j+1])/(2*dx)
    # ϕy
    if (j != 1) && (j != par.Ny-1) # Interior grid points
        ϕybl = (ϕ[i,j+1] - ϕ[i,j-1])/(2*dy)
        ϕybr = (ϕ[i+1,j+1] - ϕ[i+1,j-1])/(2*dy)
        ϕytl = (ϕ[i,j+2] - ϕ[i,j])/(2*dy)
        ϕytr = (ϕ[i+1,j+2] - ϕ[i+1,j])/(2*dy)
    elseif j == 1 # Bottom boundary point
        ϕybl = (ϕ[i,j+1] - ϕ[i,par.Ny-1])/(2*dy)
        ϕybr = (ϕ[i+1,j+1] - ϕ[i+1,par.Ny-1])/(2*dy)
        ϕytl = (ϕ[i,j+2] - ϕ[i,j])/(2*dy)
        ϕytr = (ϕ[i+1,j+2] - ϕ[i+1,j])/(2*dy)
    else # Top boundary point
        ϕybl = (ϕ[i,j+1] - ϕ[i,j-1])/(2*dy)
        ϕybr = (ϕ[i+1,j+1] - ϕ[i+1,j-1])/(2*dy)
        ϕytl = (ϕ[i,2] - ϕ[i,j])/(2*dy)
        ϕytr = (ϕ[i+1,2] - ϕ[i+1,j])/(2*dy)
    end
    # ϕxx
    ϕxxbl = (ϕ[i+1,j] - 2*ϕ[i,j] + ϕ[i-1,j])/(dx^2)
    ϕxxbr = (ϕ[i+2,j] - 2*ϕ[i+1,j] + ϕ[i,j])/(dx^2)
    ϕxxtl = (ϕ[i+1,j+1] - 2*ϕ[i,j+1] + ϕ[i-1,j+1])/(dx^2)
    ϕxxtr = (ϕ[i+2,j+1] - 2*ϕ[i+1,j+1] + ϕ[i,j+1])/(dx^2)
    # ϕyy
    if (j != 1) && (j != par.Ny-1) # Interior grid points
        ϕyybl = (ϕ[i,j+1] - 2*ϕ[i,j] + ϕ[i,j-1])/(dy^2)
        ϕyybr = (ϕ[i+1,j+1] - 2*ϕ[i+1,j] + ϕ[i+1,j-1])/(dy^2)
        ϕyytl = (ϕ[i,j+2] - 2*ϕ[i,j+1] + ϕ[i,j])/(dy^2)
        ϕyytr = (ϕ[i+1,j+2] - 2*ϕ[i+1,j+1] + ϕ[i+1,j])/(dy^2)
    elseif j == 1 # Bottom boundary point
        ϕyybl = (ϕ[i,j+1] - 2*ϕ[i,j] + ϕ[i,par.Ny-1])/(dy^2)
        ϕyybr = (ϕ[i+1,j+1] - 2*ϕ[i+1,j] + ϕ[i+1,par.Ny-1])/(dy^2)
        ϕyytl = (ϕ[i,j+2] - 2*ϕ[i,j+1] + ϕ[i,j])/(dy^2)
        ϕyytr = (ϕ[i+1,j+2] - 2*ϕ[i+1,j+1] + ϕ[i+1,j])/(dy^2)
    else # Top boundary point
        ϕyybl = (ϕ[i,j+1] - 2*ϕ[i,j] + ϕ[i,j-1])/(dy^2)
        ϕyybr = (ϕ[i+1,j+1] - 2*ϕ[i+1,j] + ϕ[i+1,j-1])/(dy^2)
        ϕyytl = (ϕ[i,2] - 2*ϕ[i,j+1] + ϕ[i,j])/(dy^2)
        ϕyytr = (ϕ[i+1,2] - 2*ϕ[i+1,j+1] + ϕ[i+1,j])/(dy^2)
    end
    # ϕxy
    if (j != 1) && (j != par.Ny-1) # Interior grid points
        ϕxybl = (ϕ[i+1,j+1] - ϕ[i+1,j-1] - ϕ[i-1,j+1] + ϕ[i-1,j-1])/(dx*dy)
        ϕxybr = (ϕ[i+2,j+1] - ϕ[i+2,j-1] - ϕ[i,j+1] + ϕ[i,j-1])/(dx*dy)
        ϕxytl = (ϕ[i+1,j+2] - ϕ[i+1,j] - ϕ[i-1,j+2] + ϕ[i-1,j])/(dx*dy)
        ϕxytr = (ϕ[i+2,j+2] - ϕ[i+2,j] - ϕ[i,j+2] + ϕ[i,j])/(dx*dy)
    elseif j == 1 # Bottom boundary point
        ϕxybl = (ϕ[i+1,j+1] - ϕ[i+1,par.Ny-1] - ϕ[i-1,j+1] + ϕ[i-1,par.Ny-1])/(dx*dy)
        ϕxybr = (ϕ[i+2,j+1] - ϕ[i+2,par.Ny-1] - ϕ[i,j+1] + ϕ[i,par.Ny-1])/(dx*dy)
        ϕxytl = (ϕ[i+1,j+2] - ϕ[i+1,j] - ϕ[i-1,j+2] + ϕ[i-1,j])/(dx*dy)
        ϕxytr = (ϕ[i+2,j+2] - ϕ[i+2,j] - ϕ[i,j+2] + ϕ[i,j])/(dx*dy)
    else # Top boundary point
        ϕxybl = (ϕ[i+1,j+1] - ϕ[i+1,j-1] - ϕ[i-1,j+1] + ϕ[i-1,j-1])/(dx*dy)
        ϕxybr = (ϕ[i+2,j+1] - ϕ[i+2,j-1] - ϕ[i,j+1] + ϕ[i,j-1])/(dx*dy)
        ϕxytl = (ϕ[i+1,2] - ϕ[i+1,j] - ϕ[i-1,2] + ϕ[i-1,j])/(dx*dy)
        ϕxytr = (ϕ[i+2,2] - ϕ[i+2,j] - ϕ[i,2] + ϕ[i,j])/(dx*dy)
    end
    # 2D weighted-mean interpolation
    ϕx = (1-θE)*(1-θN)*ϕxbl + θE*(1-θN)*ϕxbr + (1-θE)*θN*ϕxtl + θE*θN*ϕxtr
    ϕy = (1-θE)*(1-θN)*ϕybl + θE*(1-θN)*ϕybr + (1-θE)*θN*ϕytl + θE*θN*ϕytr
    ϕxx = (1-θE)*(1-θN)*ϕxxbl + θE*(1-θN)*ϕxxbr + (1-θE)*θN*ϕxxtl + θE*θN*ϕxxtr
    ϕyy = (1-θE)*(1-θN)*ϕyybl + θE*(1-θN)*ϕyybr + (1-θE)*θN*ϕyytl + θE*θN*ϕyytr
    ϕxy = (1-θE)*(1-θN)*ϕxybl + θE*(1-θN)*ϕxybr + (1-θE)*θN*ϕxytl + θE*θN*ϕxytr
    return (ϕxx*ϕy^2 - 2*ϕy*ϕx*ϕxy + ϕyy*ϕx^2)/((ϕx^2 + ϕy^2)^1.5)
end
