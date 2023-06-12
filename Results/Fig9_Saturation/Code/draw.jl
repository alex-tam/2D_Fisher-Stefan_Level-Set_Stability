# Plot solutions to 2D Fisher-Stefan
# Alex Tam, 02/02/2022

using Plots
using Measures
using LaTeXStrings
using DelimitedFiles
using Printf
using Polynomials

"Plot solutions as 2D heat maps"
function draw_heat(x, y, U, V, ϕ, i, Lx, Ly)
    gr(); plot() # Load GR plotting backend and clear previous plots
    default(fontfamily = "Computer Modern", titlefontsize = 18, guidefontsize = 20, tickfontsize = 12, legendfontsize = 14)
    p1 = heatmap(x,y,transpose(U), xlabel = L"$x$", ylabel = L"$y$", margin=3mm, aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma, clims=(0.0, 1.0))
    contour!(p1, x, y, transpose(ϕ), c=:red, linewidth = 2, levels=[0.000001], clim=(0.0, 1.0))
    savefig("u-$i.pdf")
    heatmap(x,y,transpose(V), xlabel = L"$x$", ylabel = L"$y$", margin=3mm, aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma)
    savefig("V-$i.pdf")
    p2 = heatmap(x,y,transpose(ϕ), xlabel = L"$x$", ylabel = L"$y$", margin=3mm, aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma)
    contour!(p2, x, y, transpose(ϕ), c=:red, linewidth = 2, levels=[0.0])
    savefig("phi-$i.pdf")
end

function draw_heat(x, y, U, ϕ, i, Lx, Ly)
    gr(); plot() # Load GR plotting backend and clear previous plots
    default(fontfamily = "Computer Modern", titlefontsize = 18, guidefontsize = 20, tickfontsize = 12, legendfontsize = 14)
    p1 = heatmap(x,y,transpose(U), xlabel = L"$x$", ylabel = L"$y$", margin=3mm, aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma, clims=(0.0, 1.0))
    contour!(p1, x, y, transpose(ϕ), c=:red, linewidth = 2, levels=[0.000001], clim=(0.0, 1.0))
    savefig("u-$i.pdf")
    p2 = heatmap(x,y,transpose(ϕ), xlabel = L"$x$", ylabel = L"$y$", margin=3mm, aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma)
    contour!(p2, x, y, transpose(ϕ), c=:red, linewidth = 2, levels=[0.0])
    savefig("phi-$i.pdf")
end

"Density slice plots"
function draw_slices(x, y, nx, ny, Lx, Ly, plot_times)
    gr(); plot() # Load GR plotting backend and clear previous plots
    default(fontfamily = "Computer Modern", titlefontsize = 18, guidefontsize = 26, tickfontsize = 18, legendfontsize = 18)
    U = readdlm("U-0.csv")
    plot(x, U[:,ny], xlabel = L"$x$", ylabel = L"$u(x,10,t)$", linecolor = :black, linewidth = 2, aspect_ratio = 20.0, grid = false, margin=3mm, legend = false, xlims=(0,Lx), ylims=(0,1))
    for i in plot_times
        u = readdlm("ux-$i.csv")
        plot!(x, u, linecolor = :black, linestyle = :dash, linewidth = 2)
    end
    savefig("u_slice_x.pdf")
    plot(y, U[nx,:], xlabel = L"$y$", ylabel = L"$u(10,y,t)$", linecolor = :black, linewidth = 2, aspect_ratio = 20.0, grid = false, margin=3mm, legend = false, xlims=(0,Ly), ylims=(0,1))
    for i in plot_times
        u = readdlm("uy-$i.csv")
        plot!(y, u, linecolor = :black, linestyle = :dash, linewidth = 2)
    end
    savefig("u_slice_y.pdf")
end

"Compute growth rate"
function draw_growth(t, Amp, t_min::Float64)
    gr(); plot() # Load GR plotting backend and clear previous plots
    default(fontfamily = "Computer Modern", titlefontsize = 14, guidefontsize = 20, tickfontsize = 14, legendfontsize = 14)
    ind::Int = 0 # Index corresponding to t = 0.1
    ind_end::Int = 1001 # Index for t = 1
    for i = 1:length(t)
        if t[i] <= t_min
            ind += 1
        end
    end
    poly = fit(t[ind:ind_end], log.(Amp)[ind:ind_end], 1) # Fit straight line to data
    plot(t, log.(Amp), xlabel = L"$t$", ylabel = L"$\log(A)$", label = "Numerical Data", linewidth = 2, margin = 5mm, legend=:topleft) # Plot data
    scatter!([t[ind], t[ind_end]], [log(Amp[ind]), log(Amp[ind_end])], markersize = 5, markershape = :xcross, markercolor = :red, label = false) # Scatter plot of t_min
    plot!(t, poly.coeffs[2].*t .+ poly.coeffs[1], label = "Linear Fit", linestyle=:dash, linewidth = 2) # Plot linear trendline
    ω = poly.coeffs[2] # Obtain slope
    @printf("The numerical growth rate is: %f.\n", ω)
    savefig("perturbation_amplitude.pdf")
end

"Plot interface position versus time"
function draw_interface(t, L, ε, x, y, dx, Lx, Ly, Nx, Ny, plot_times)
    gr(); plot() # Load GR plotting backend and clear previous plots
    default(titlefont = (18, "Computer Modern"), guidefont = (26, "Computer Modern"), tickfont = (18, "Computer Modern"))
    plot(t, L, xlabel = L"$t$", ylabel = L"$L(t)$", margin=3mm, xlims=(0,maximum(t)), ylims=(0,maximum(L)), legend = false)
    savefig("L.pdf")
    if ε == 0.0
        @printf("Numerical travelling wave speed is %f.\n", (L[end]-L[1])/par.T)
    end
    # Optional: Draw interface
    plot() # Load GR plotting backend and clear previous plots
    plot_interface(x, y, dx, Lx, Ly, Nx, Ny, plot_times)
    savefig("interface.pdf")
end

"Draw interface as series of points"
function plot_interface(x, y, dx, Lx, Ly, Nx, Ny, plot_times)
    pt = vcat(0, plot_times)
    pt = [0,1000,2000,3000,4000,5000,6000,7000]
    for k in pt
        # Import level-set function at relevant time
        ϕ = readdlm("Phi-$k.csv")
        # Pre-allocate empty vectors
        xi = Vector{Float64}()
        yi = Vector{Float64}()
        # Locate interface grid points
        for j = 1:Ny
            ϕv = ϕ[:,j] # Obtain 1D vector of ϕ
            for i = 1:Nx
                if (ϕv[i] < 0) && (ϕv[i+1] >= 0)
                    θ = ϕv[i]/(ϕv[i] - ϕv[i+1])
                    push!(xi, x[i] + θ*dx)
                    push!(yi, y[j])
                end
            end
        end
        # Plot interface
        scatter!([xi],[yi], xlabel = L"$x$", ylabel = L"$y$", margin=3mm, xlims=(0,Lx), ylims=(0,Ly), color="green", legend = false, aspect_ratio=:equal, markersize=2)
    end
    savefig("interface.pdf")
end

"Control function for plotting"
function draw()
    # Import data
    plot_times = convert(Vector{Int}, vec(readdlm("plot_times.csv")))
    x = vec(readdlm("x.csv")); dx = x[2] - x[1]
    y = vec(readdlm("y.csv"))
    nx::Int = (length(x)-1)/2; ny::Int = (length(y)-1)/2
    Lx = maximum(x); Ly = maximum(y); 
    Nx = length(x); Ny = length(y)
    U = readdlm("U-0.csv")
    ϕ = readdlm("Phi-0.csv")
    t = vec(readdlm("t.csv"))
    L = vec(readdlm("L.csv"))
    Amp = vec(readdlm("Amp.csv"))
    ε = 0.1
    # Plot
    draw_heat(x, y, U, ϕ, 0, Lx, Ly) # Heat maps
    draw_slices(x, y, nx, ny, Lx, Ly, plot_times) # Slice plots
    draw_interface(t, L, ε, x, y, dx, Lx, Ly, Nx, Ny, plot_times)
    for i in plot_times
        U = readdlm("U-$i.csv")
        V = readdlm("V-$i.csv")
        ϕ = readdlm("Phi-$i.csv")
        draw_heat(x, y, U, V, ϕ, i, Lx, Ly)
    end
    # Numerical growth rate
    draw_growth(t, Amp, ε)
end

@time draw()