# Compute extension velocity field for Fisher-Stefan level-set method
# Alex Tam, 02/02/2022

"Compute velocity extension by orthogonal extrapolation"
function extend_velocity(D, dΩ, U, ϕ, par, dx, dy)
    dτ = 0.5*(dx/5 + dy/5) # "Time" step for extrapolation PDE
    vi = interface_speed(dΩ, U, ϕ, par, dx, dy) # Find speed at interface using Stefan condition
    V = velocity_init(D, dΩ, par, vi) # Obtain "initial" condition for V(x,y,t) using speed at interface
    v = build_vector(V, D) # Create vector using data at interior grid points
    # Solve PDE to extrapolate V outwards from interface
    prob = ODEProblem((du, u, p, t) -> vel_outward_rhs!(du, u, p, t, D, dΩ, vi, par, ϕ, dx, dy), v, (0, par.V_Iterations*dτ))
    sol = solve(prob, Tsit5(), reltol = 1e-3, abstol = 1e-6, saveat = par.V_Iterations*dτ)
    v = sol[:,end] # Update vector of V
    # # Solve PDE to extrapolate V inwards from interface
    prob = ODEProblem((du, u, p, t) -> vel_inward_rhs!(du, u, p, t, D, dΩ, vi, par, ϕ, dx, dy), v, (0, par.V_Iterations*dτ))
    sol = solve(prob, Tsit5(), reltol = 1e-3, abstol = 1e-6, saveat = par.V_Iterations*dτ)
    # Reshape solution to matrix
    V = build_v_matrix(sol[:,end], par, D)
    return V
end

# This could be changed, perhaps
"Obtain initial condition for V(x,y,t) using speed at interface"
function velocity_init(D, dΩ, par, vi)
    V = zeros(par.Nx, par.Ny) # Pre-allocate matrix of V(x,y,t)
    # Manually overwrite speed close to interface
    for gp in D # Loop over interior grid points
        if (gp.dΩ == true) # If grid point is close to interface
            # Locate closest interface point to current grid point
            if (gp.θxp < gp.θxm) && (gp.θxp < gp.θyp) && (gp.θxp < gp.θym) # Closest point in positive x-direction
                xm = gp.xInd; xp = gp.xInd+1; ym = gp.yInd; yp = gp.yInd
            elseif (gp.θxm < gp.θxp) && (gp.θxm < gp.θyp) && (gp.θxm < gp.θym) # Closest point in negative x-direction
                xm = gp.xInd-1; xp = gp.xInd; ym = gp.yInd; yp = gp.yInd
            elseif (gp.θyp < gp.θxm) && (gp.θyp < gp.θxp) && (gp.θyp < gp.θym) # Closest point in positive y-direction
                xm = gp.xInd; xp = gp.xInd; ym = gp.yInd; yp = gp.yInd+1
            else # Closest point in negative y-direction
                xm = gp.xInd; xp = gp.xInd; ym = gp.yInd-1; yp = gp.yInd
            end
            # Assign V to interface speed at closest point
            for i = 1:length(dΩ)
                if (xm == dΩ[i].m.xInd) && (ym == dΩ[i].m.yInd) && (xp == dΩ[i].p.xInd) && (yp == dΩ[i].p.yInd)
                    V[gp.xInd, gp.yInd] = vi[i]
                end
            end
        end
    end
    return V
end

"Obtain velocity at a specified interface point"
function get_interface_speed(xm, ym, xp, yp, dΩ, vi)
    for i = 1:length(dΩ) # Loop over candidate points
        if (xm == dΩ[i].m.xInd) && (ym == dΩ[i].m.yInd) && (xp == dΩ[i].p.xInd) && (yp == dΩ[i].p.yInd) # If all indices are correct
            return vi[i]
        end
    end
end

"Construct right-hand vector of outward velocity extension PDE, for use in DifferentialEquations.jl"
function vel_outward_rhs!(du, u, p, t, D, dΩ, vi, par, ϕ, dx, dy)
    V = build_v_matrix(u, par, D) # Generate matrix
    for i = 1:length(D) # Loop over interior grid points
        if (D[i].Ω == true) || (D[i].dΩ == true) # Grid points in Ω or close to interface
            du[i] = 0.0 # Set RHS to zero
        else
            # Compute advection coefficients
            ϕx = 1/(2*dx)*(ϕ[D[i].xInd+1, D[i].yInd]-ϕ[D[i].xInd-1, D[i].yInd]) # dϕ/dx
            # dϕ/dy
            if D[i].yInd == 1
                ϕy = 1/(2*dy)*(ϕ[D[i].xInd, 2]-ϕ[D[i].xInd, par.Ny-1]) # Periodic BC
            elseif D[i].yInd == par.Ny
                ϕy = 1/(2*dy)*(ϕ[D[i].xInd, 2]-ϕ[D[i].xInd, par.Ny-1]) # Periodic BC
            else # Interior grid point
                ϕy = 1/(2*dy)*(ϕ[D[i].xInd, D[i].yInd+1]-ϕ[D[i].xInd, D[i].yInd-1])
            end
            n = sqrt(ϕx^2 + ϕy^2) # Normalisation
            if n == 0 # Prevent division by zero
                a = 0; b = 0
            else
                a = ϕx/n # x-component
                b = ϕy/n # y-component
            end
            # Discretise dV/dx (first-order upwind)
            if a < 0 # Forward difference
                if D[i].θxp == 1.0 # Regular discretisation
                    Vx = (V[D[i].xInd+1, D[i].yInd] - V[D[i].xInd, D[i].yInd])/dx
                else # Irregular discretisation
                    v = get_interface_speed(D[i].xInd, D[i].yInd, D[i].xInd+1, D[i].yInd, dΩ, vi)
                    Vx = (v - V[D[i].xInd, D[i].yInd])/(D[i].θxp*dx)
                end
            else # Backward difference
                if D[i].θxm == 1.0 # Regular discretisation
                    Vx = (V[D[i].xInd, D[i].yInd] - V[D[i].xInd-1, D[i].yInd])/dx
                else # Irregular discretisation
                    v = get_interface_speed(D[i].xInd-1, D[i].yInd, D[i].xInd, D[i].yInd, dΩ, vi)
                    Vx = (V[D[i].xInd, D[i].yInd]-v)/(D[i].θxm*dx)
                end
            end
            # Discretise dV/dy (first-order upwind)
            if D[i].yInd == 1 # Lower boundary
                if b < 0 # Forward difference
                    if D[i].θyp == 1.0 # Regular discretisation
                        Vy = (V[D[i].xInd, D[i].yInd+1] - V[D[i].xInd, D[i].yInd])/dy
                    else # Irregular discretisation
                        v = get_interface_speed(D[i].xInd, D[i].yInd, D[i].xInd, D[i].yInd+1, dΩ, vi)
                        Vy = (v - V[D[i].xInd, D[i].yInd])/(D[i].θyp*dy)
                    end
                else # Backward difference
                    if D[i].θym == 1.0
                        Vy = (V[D[i].xInd, 1] - V[D[i].xInd, par.Ny-1])/dy
                    else
                        # Obtain interface velocity
                        v = get_interface_speed(D[i].xInd, par.Ny-1, D[i].xInd, D[i].yInd, dΩ, vi)
                        Vy = (V[D[i].xInd, D[i].yInd]-v)/(D[i].θym*dy)
                    end
                end
            elseif D[i].yInd == par.Ny # Upper boundary
                if b < 0 # Forward difference
                    if D[i].θyp == 1.0 # Regular discretisation
                        Vy = (V[D[i].xInd, 2] - V[D[i].xInd, D[i].yInd])/dy
                    else # Irregular discretisation
                        v = get_interface_speed(D[i].xInd, 1, D[i].xInd, 2, dΩ, vi)
                        Vy = (v - V[D[i].xInd, D[i].yInd])/(D[i].θyp*dy)
                    end
                else # Backward difference
                    if D[i].θym == 1.0
                        Vy = (V[D[i].xInd, D[i].yInd] - V[D[i].xInd, D[i].yInd-1])/dy
                    else
                        # Obtain interface velocity
                        v = get_interface_speed(D[i].xInd, D[i].yInd-1, D[i].xInd, D[i].yInd, dΩ, vi)
                        Vy = (V[D[i].xInd, D[i].yInd]-v)/(D[i].θym*dy)
                    end
                end
            else # Interior points
                if b < 0 # Forward difference
                    if D[i].θyp == 1.0 # Regular discretisation
                        Vy = (V[D[i].xInd, D[i].yInd+1] - V[D[i].xInd, D[i].yInd])/dy
                    else # Irregular discretisation
                        v = get_interface_speed(D[i].xInd, D[i].yInd, D[i].xInd, D[i].yInd+1, dΩ, vi)
                        Vy = (v - V[D[i].xInd, D[i].yInd])/(D[i].θyp*dy)
                    end
                else # Backward difference
                    if D[i].θym == 1.0
                        Vy = (V[D[i].xInd, D[i].yInd] - V[D[i].xInd, D[i].yInd-1])/dy
                    else
                        # Obtain interface velocity
                        v = get_interface_speed(D[i].xInd, D[i].yInd-1, D[i].xInd, D[i].yInd, dΩ, vi)
                        Vy = (V[D[i].xInd, D[i].yInd]-v)/(D[i].θym*dy)
                    end
                end
            end
            du[i] = -a*Vx - b*Vy # Update RHS term
        end
    end
end

"Construct right-hand vector of inward velocity extension PDE, for use in DifferentialEquations.jl"
function vel_inward_rhs!(du, u, p, t, D, dΩ, vi, par, ϕ, dx, dy)
    V = build_v_matrix(u, par, D) # Generate matrix
    for i = 1:length(D) # Loop over interior grid points
        if (D[i].Ω == false) || (D[i].dΩ == true) # Grid points outside Ω or close to interface
            du[i] = 0.0 # Set RHS to zero
        else
            # Compute advection coefficients
            ϕx = 1/(2*dx)*(ϕ[D[i].xInd+1, D[i].yInd]-ϕ[D[i].xInd-1, D[i].yInd]) # dϕ/dx
            # dϕ/dy
            if D[i].yInd == 1
                ϕy = 1/(2*dy)*(ϕ[D[i].xInd, 2]-ϕ[D[i].xInd, par.Ny-1]) # Periodic BC
            elseif D[i].yInd == par.Ny
                ϕy = 1/(2*dy)*(ϕ[D[i].xInd, 2]-ϕ[D[i].xInd, par.Ny-1]) # Periodic BC
            else # Interior grid point
                ϕy = 1/(2*dy)*(ϕ[D[i].xInd, D[i].yInd+1]-ϕ[D[i].xInd, D[i].yInd-1])
            end
            n = sqrt(ϕx^2 + ϕy^2) # Normalisation
            if n == 0 # Prevent division by zero
                a = 0; b = 0
            else
                a = -ϕx/n # x-component
                b = -ϕy/n # y-component
            end
            # Discretise dV/dx (first-order upwind)
            if a < 0 # Forward difference
                if D[i].θxp == 1.0 # Regular discretisation
                    Vx = (V[D[i].xInd+1, D[i].yInd] - V[D[i].xInd, D[i].yInd])/dx
                else # Irregular discretisation
                    v = get_interface_speed(D[i].xInd, D[i].yInd, D[i].xInd+1, D[i].yInd, dΩ, vi)
                    Vx = (v - V[D[i].xInd, D[i].yInd])/(D[i].θxp*dx)
                end
            else # Backward difference
                if D[i].θxm == 1.0 # Regular discretisation
                    Vx = (V[D[i].xInd, D[i].yInd] - V[D[i].xInd-1, D[i].yInd])/dx
                else # Irregular discretisation
                    v = get_interface_speed(D[i].xInd-1, D[i].yInd, D[i].xInd, D[i].yInd, dΩ, vi)
                    Vx = (V[D[i].xInd, D[i].yInd]-v)/(D[i].θxm*dx)
                end
            end
            # Discretise dV/dy (first-order upwind)
            if D[i].yInd == 1 # Lower boundary
                if b < 0 # Forward difference
                    if D[i].θyp == 1.0 # Regular discretisation
                        Vy = (V[D[i].xInd, D[i].yInd+1] - V[D[i].xInd, D[i].yInd])/dy
                    else # Irregular discretisation
                        v = get_interface_speed(D[i].xInd, D[i].yInd, D[i].xInd, D[i].yInd+1, dΩ, vi)
                        Vy = (v - V[D[i].xInd, D[i].yInd])/(D[i].θyp*dy)
                    end
                else # Backward difference
                    if D[i].θym == 1.0
                        Vy = (V[D[i].xInd, 1] - V[D[i].xInd, par.Ny-1])/dy
                    else
                        # Obtain interface velocity
                        v = get_interface_speed(D[i].xInd, par.Ny-1, D[i].xInd, par.Ny, dΩ, vi)
                        Vy = (V[D[i].xInd, D[i].yInd]-v)/(D[i].θym*dy)
                    end
                end
            elseif D[i].yInd == par.Ny # Upper boundary
                if b < 0 # Forward difference
                    if D[i].θyp == 1.0 # Regular discretisation
                        Vy = (V[D[i].xInd, 2] - V[D[i].xInd, par.Ny])/dy
                    else # Irregular discretisation
                        v = get_interface_speed(D[i].xInd, 1, D[i].xInd, 2, dΩ, vi)
                        Vy = (v - V[D[i].xInd, D[i].yInd])/(D[i].θyp*dy)
                    end
                else # Backward difference
                    if D[i].θym == 1.0
                        Vy = (V[D[i].xInd, D[i].yInd] - V[D[i].xInd, D[i].yInd-1])/dy
                    else
                        # Obtain interface velocity
                        v = get_interface_speed(D[i].xInd, D[i].yInd-1, D[i].xInd, D[i].yInd, dΩ, vi)
                        Vy = (V[D[i].xInd, D[i].yInd]-v)/(D[i].θym*dy)
                    end
                end
            else # Interior points
                if b < 0 # Forward difference
                    if D[i].θyp == 1.0 # Regular discretisation
                        Vy = (V[D[i].xInd, D[i].yInd+1] - V[D[i].xInd, D[i].yInd])/dy
                    else # Irregular discretisation
                        v = get_interface_speed(D[i].xInd, D[i].yInd, D[i].xInd, D[i].yInd+1, dΩ, vi)
                        Vy = (v - V[D[i].xInd, D[i].yInd])/(D[i].θyp*dy)
                    end
                else # Backward difference
                    if D[i].θym == 1.0
                        Vy = (V[D[i].xInd, D[i].yInd] - V[D[i].xInd, D[i].yInd-1])/dy
                    else
                        # Obtain interface velocity
                        v = get_interface_speed(D[i].xInd, D[i].yInd-1, D[i].xInd, D[i].yInd, dΩ, vi)
                        Vy = (V[D[i].xInd, D[i].yInd]-v)/(D[i].θym*dy)
                    end
                end
            end
            du[i] = -a*Vx - b*Vy # Update RHS term
        end
    end
end