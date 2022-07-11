# Apply penalisation to interface density
# Alex Tam, 22/2/2022

"Compute density at interface using surface tension regularisation"
function interface_density(dΩ, ϕ, par, dx, dy)
    uf = Vector{Float64}() # Pre-allocate empty vector of interface density
    for ip in dΩ
        # Compute derivatives of level-set function
        if ip.m.xInd != ip.p.xInd # Interface lies in x-direction
            ϕx = find_ϕx_x(ip, ϕ, dx)
            ϕy = find_ϕy_x(ip, ϕ, par, dy)
            ϕxx = find_ϕxx_x(ip, ϕ, dx)
            ϕxy = find_ϕxy_x(ip, ϕ, par, dx, dy)
            ϕyy = find_ϕyy_x(ip, ϕ, par, dy)
        else # Interface lies in y-direction
            ϕx = find_ϕx_y(ip, ϕ, dx)
            ϕy = find_ϕy_y(ip, ϕ, par, dy)
            ϕxx = find_ϕxx_y(ip, ϕ, dx)
            ϕxy = find_ϕxy_y(ip, ϕ, par, dx, dy)
            ϕyy = find_ϕyy_y(ip, ϕ, par, dy)
        end
        # Compute local curvature and apply penalisation
        K = (ϕxx*ϕy^2 - 2*ϕy*ϕx*ϕxy + ϕyy*ϕx^2)/((ϕx^2 + ϕy^2)^1.5)
        push!(uf, par.uf - par.γ*K) # Modify later to add curvature formula
    end
    return uf
end

"Obtain velocity at a specified interface point"
function get_interface_density(xm, ym, xp, yp, dΩ, uf)
    for i = 1:length(dΩ) # Loop over candidate points
        if (xm == dΩ[i].m.xInd) && (ym == dΩ[i].m.yInd) && (xp == dΩ[i].p.xInd) && (yp == dΩ[i].p.yInd) # If all indices are correct
            return uf[i]
        end
    end
end

# Functions for second derivatives of ϕ on interface
"ϕ_xx for interface in x-direction"
function find_ϕxx_x(ip, ϕ, dx)
    i = ip.m.xInd; j = ip.m.yInd; θE = ip.m.θxp # Left point is reference index
    ϕxxl = (ϕ[i+1,j] - 2*ϕ[i,j] + ϕ[i-1,j])/(dx^2)
    ϕxxr = (ϕ[i+2,j] - 2*ϕ[i+1,j] + ϕ[i,j])/(dx^2) # Assume far enough from x-boundary
    return θE*ϕxxr + (1-θE)*ϕxxl
end

"ϕ_xy for interface in x-direction"
function find_ϕxy_x(ip, ϕ, par, dx, dy)
    i = ip.m.xInd; j = ip.m.yInd; θE = ip.m.θxp # Left point is reference index
    if (j != 1) && (j != par.Ny) # Interior grid points
        ϕxlt = (ϕ[i+1,j+1]-ϕ[i-1,j+1])/(2*dx)
        ϕxrt = (ϕ[i+2,j+1]-ϕ[i,j+1])/(2*dx)
        ϕxlb = (ϕ[i+1,j-1]-ϕ[i-1,j-1])/(2*dx)
        ϕxrb = (ϕ[i+2,j-1]-ϕ[i,j-1])/(2*dx)
    elseif j == 1 # Bottom boundary point
        ϕxlt = (ϕ[i+1,j+1]-ϕ[i-1,j+1])/(2*dx)
        ϕxrt = (ϕ[i+2,j+1]-ϕ[i,j+1])/(2*dx)
        ϕxlb = (ϕ[i+1,par.Ny-1]-ϕ[i-1,par.Ny-1])/(2*dx)
        ϕxrb = (ϕ[i+2,par.Ny-1]-ϕ[i,par.Ny-1])/(2*dx)
    else # Top boundary point
        ϕxlt = (ϕ[i+1,2]-ϕ[i-1,2])/(2*dx)
        ϕxrt = (ϕ[i+2,2]-ϕ[i,2])/(2*dx)
        ϕxlb = (ϕ[i+1,j-1]-ϕ[i-1,j-1])/(2*dx)
        ϕxrb = (ϕ[i+2,j-1]-ϕ[i,j-1])/(2*dx)
    end
    ϕxt = θE*ϕxrt + (1-θE)*ϕxlt
    ϕxb = θE*ϕxrb + (1-θE)*ϕxlb
    return (ϕxt - ϕxb)/(2*dy)
end

"ϕ_yy for interface in x-direction"
function find_ϕyy_x(ip, ϕ, par, dy)
    i = ip.m.xInd; j = ip.m.yInd; θE = ip.m.θxp # Left point is reference index
    if (j != 1) && (j != par.Ny) # Interior grid points
        ϕt = θE*ϕ[i+1,j+1] + (1-θE)*ϕ[i,j+1]
        ϕb = θE*ϕ[i+1,j-1] + (1-θE)*ϕ[i,j-1]
    elseif j == 1 # Bottom boundary point
        ϕt = θE*ϕ[i+1,j+1] + (1-θE)*ϕ[i,j+1]
        ϕb = θE*ϕ[i+1,par.Ny-1] + (1-θE)*ϕ[i,par.Ny-1]
    else # Top boundary point
        ϕt = θE*ϕ[i+1,2] + (1-θE)*ϕ[i,2]
        ϕb = θE*ϕ[i+1,j-1] + (1-θE)*ϕ[i,j-1]
    end
    return (ϕt + ϕb)/(dy^2) # Assume ϕ = 0 on interface
end

"ϕ_xx for interface in y-direction"
function find_ϕxx_y(ip, ϕ, dx)
    i = ip.m.xInd; j = ip.m.yInd; θN = ip.m.θyp # Bottom point is reference index
    ϕl = θN*ϕ[i-1,j+1] + (1-θN)*ϕ[i-1,j]
    ϕr = θN*ϕ[i+1,j+1] + (1-θN)*ϕ[i+1,j]
    return (ϕr + ϕl)/(dx^2)
end

"ϕ_xy for interface in y-direction"
function find_ϕxy_y(ip, ϕ, par, dx, dy)
    i = ip.m.xInd; j = ip.m.yInd; θN = ip.m.θyp # Bottom point is reference index
    if (j != 1) && (j != par.Ny-1) # Interior grid points
        ϕylt = ϕ[i-1,j+2] - ϕ[i-1,j]
        ϕylb = ϕ[i-1,j+1] - ϕ[i-1,j-1]
        ϕyrt = ϕ[i+1,j+2] - ϕ[i+1,j]
        ϕyrb = ϕ[i+1,j+1] - ϕ[i+1,j-1]
    elseif (j == 1) # Bottom boundary point
        ϕylt = ϕ[i-1,j+2] - ϕ[i-1,j]
        ϕylb = ϕ[i-1,j+1] - ϕ[i-1,par.Ny-1]
        ϕyrt = ϕ[i+1,j+2] - ϕ[i+1,j]
        ϕyrb = ϕ[i+1,j+1] - ϕ[i+1,par.Ny-1]
    else # Top boundary point
        ϕylt = ϕ[i-1,2] - ϕ[i-1,j]
        ϕylb = ϕ[i-1,j+1] - ϕ[i-1,j-1]
        ϕyrt = ϕ[i+1,2] - ϕ[i+1,j]
        ϕyrb = ϕ[i+1,j+1] - ϕ[i+1,j-1]
    end
    ϕyl = θN*ϕylt + (1-θN)*ϕylb
    ϕyr = θN*ϕyrt + (1-θN)*ϕyrb
    return (ϕyr - ϕyl)/(2*dx)
end

"ϕ_yy for interface in y-direction"
function find_ϕyy_y(ip, ϕ, par, dy)
    i = ip.m.xInd; j = ip.m.yInd; θN = ip.m.θyp # Bottom point is reference index
    if (j != 1) && (j != par.Ny-1) # Interior grid points
        ϕyyt = (ϕ[i,j+2] - 2*ϕ[i,j+1] + ϕ[i,j])/(dy^2)
        ϕyyb = (ϕ[i,j+1] - 2*ϕ[i,j] + ϕ[i,j-1])/(dy^2)
    elseif (j == 1) # Bottom boundary point
        ϕyyt = (ϕ[i,j+2] - 2*ϕ[i,j+1] + ϕ[i,j])/(dy^2)
        ϕyyb = (ϕ[i,j+1] - 2*ϕ[i,j] + ϕ[i,par.Ny-1])/(dy^2)
    else # Top boundary point
        ϕyyt = (ϕ[i,2] - 2*ϕ[i,j+1] + ϕ[i,j])/(dy^2)
        ϕyyb = (ϕ[i,j+1] - 2*ϕ[i,j] + ϕ[i,j-1])/(dy^2)
    end
    return θN*ϕyyt + (1-θN)*ϕyyb
end