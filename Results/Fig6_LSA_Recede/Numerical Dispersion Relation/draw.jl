# Plot numerical dispersion relations from linear stability analysis
# Alex Tam, 01/06/2022

using DelimitedFiles
using Plots
using LaTeXStrings
using Measures

function draw_numerical()
    # Import data (κ,uf,γ)
    q = vec(readdlm("q-(-0.3,0.0,0.0).csv")); ω = vec(readdlm("omega-(-0.3,0.0,0.0).csv"))
    qn = [0, 2, 4, 6, 8, 10, 12, 14].*(π/10) # Numerical wave numbers
    ωn = [0.0, 0.073005, 0.191092, 0.297989, 0.387592, 0.463523, 0.534146, 0.601290] # Numerical growth rates
    # Plot data
    plot(q, ω, label = L"$\textrm{Theoretical}$", linecolor=:black, linewidth = 2, xlabel = L"$q$", ylabel = L"$\omega$", margin = 3mm, xlims = (0, 5.0), ylims = (-0.5,1.0), grid = true, legend=:bottomleft)
    scatter!(qn, ωn, label = L"$\textrm{Numerical}$", linecolor=:red, markercolor=:red, markersize = 6)
    savefig("Recede_Numerical_Dispersion.pdf")
end

# Plot configuration
gr() # Load GR plotting backend
default(titlefont = (18, "Computer Modern"), guidefont = (18, "Computer Modern"), tickfont = (18, "Computer Modern"), legendfontsize = 14)
# Generate plots
@time draw_numerical()
