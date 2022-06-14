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
    default(fontfamily = "Computer Modern", titlefontsize = 18, guidefontsize = 26, tickfontsize = 18, legendfontsize = 18)
    heatmap(x,y,transpose(U), xlabel = L"$x$", ylabel = L"$y$", margin=3mm, aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma, clims=(0.0, 1.0))
    savefig("u-$i.pdf")
    heatmap(x,y,transpose(V), xlabel = L"$x$", ylabel = L"$y$", margin=3mm, aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma)
    savefig("V-$i.pdf")
    heatmap(x,y,transpose(ϕ), xlabel = L"$x$", ylabel = L"$y$", margin=3mm, aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma)
    savefig("phi-$i.pdf")
end

function draw_heat(x, y, U, ϕ, i, Lx, Ly)
    gr(); plot() # Load GR plotting backend and clear previous plots
    default(fontfamily = "Computer Modern", titlefontsize = 18, guidefontsize = 26, tickfontsize = 18, legendfontsize = 18)
    heatmap(x,y,transpose(U), xlabel = L"$x$", ylabel = L"$y$", margin=3mm, aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma, clims=(0.0, 1.0))
    savefig("u-$i.pdf")
    heatmap(x,y,transpose(ϕ), xlabel = L"$x$", ylabel = L"$y$", margin=3mm, aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma)
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
    for i = 1:length(t)
        if t[i] <= t_min
            ind += 1
        end
    end
    poly = fit(t[ind:end], log.(Amp)[ind:end], 1) # Fit straight line to data
    plot(t, log.(Amp), xlabel = L"$t$", ylabel = L"$\log(\textrm{Amplitude})$", label = "Numerical Data", linewidth = 2, margin = 5mm) # Plot data
    scatter!([t[ind], t[end]], [log(Amp[ind]), log(Amp[end])], markersize = 5, markershape = :xcross, markercolor = :red, label = false) # Scatter plot of t_min
    plot!(t, poly.coeffs[2].*t .+ poly.coeffs[1], label = "Linear Fit", linestyle=:dash, linewidth = 2) # Plot linear trendline
    ω = poly.coeffs[2] # Obtain slope
    @printf("The numerical growth rate is: %f.\n", ω)
    savefig("perturbation_amplitude.pdf")
end

"Control function for plotting"
function draw()
    # Import data
    plot_times = convert(Vector{Int}, vec(readdlm("plot_times.csv")))
    x = vec(readdlm("x.csv"))
    y = vec(readdlm("y.csv"))
    nx::Int = (length(x)-1)/2; ny::Int = (length(y)-1)/2
    Lx = maximum(x); Ly = maximum(y)
    U = readdlm("U-0.csv")
    ϕ = readdlm("Phi-0.csv")
    t = vec(readdlm("t.csv"))
    Amp = vec(readdlm("Amp.csv"))
    # Plot
    draw_heat(x, y, U, ϕ, 0, Lx, Ly)
    draw_slices(x, y, nx, ny, Lx, Ly, plot_times)
    for i in plot_times
        U = readdlm("U-$i.csv")
        V = readdlm("V-$i.csv")
        ϕ = readdlm("Phi-$i.csv")
        draw_heat(x, y, U, V, ϕ, i, Lx, Ly)
    end
    # Numerical growth rate
    draw_growth(t, Amp, 0.1)
end

@time draw()