# Plot numerical dispersion relations from linear stability analysis
# Alex Tam, 01/06/2022

using DelimitedFiles
using Plots
using LaTeXStrings
using Measures

function draw_kappa()
    # Import data (κ,uf,γ)
    q1 = vec(readdlm("q-(-0.1,0.0,0.0).csv")); ω1 = vec(readdlm("omega-(-0.1,0.0,0.0).csv"))  
    q2 = vec(readdlm("q-(-0.2,0.0,0.0).csv")); ω2 = vec(readdlm("omega-(-0.2,0.0,0.0).csv")) 
    q3 = vec(readdlm("q-(-0.3,0.0,0.0).csv")); ω3 = vec(readdlm("omega-(-0.3,0.0,0.0).csv")) 
    q4 = vec(readdlm("q-(-0.4,0.0,0.0).csv")); ω4 = vec(readdlm("omega-(-0.4,0.0,0.0).csv")) 
    q5 = vec(readdlm("q-(-0.5,0.0,0.0).csv")); ω5 = vec(readdlm("omega-(-0.5,0.0,0.0).csv"))  
    # Plot data
    plot(q1, ω1, label = L"$\kappa = -0.1$", linewidth = 2, xlabel = L"$q$", ylabel = L"$\omega$", margin = 3mm, xlims = (0, 5.0), ylims = (-0.5,3.0), grid = true, legend=:topleft)
    plot!(q2, ω2, label = L"$\kappa = -0.2$", linewidth = 2)
    plot!(q3, ω3, label = L"$\kappa = -0.3$", linewidth = 2)
    plot!(q4, ω4, label = L"$\kappa = -0.4$", linewidth = 2)
    plot!(q5, ω5, label = L"$\kappa = -0.5$", linewidth = 2)
    savefig("Fig6a.pdf")
end

function draw_uf()
    # Import data (κ,uf,γ)
    q1 = vec(readdlm("q-(-0.5,0.0,0.0).csv")); ω1 = vec(readdlm("omega-(-0.5,0.0,0.0).csv"))  
    q2 = vec(readdlm("q-(-0.5,0.1,0.0).csv")); ω2 = vec(readdlm("omega-(-0.5,0.1,0.0).csv")) 
    q3 = vec(readdlm("q-(-0.5,0.2,0.0).csv")); ω3 = vec(readdlm("omega-(-0.5,0.2,0.0).csv")) 
    q4 = vec(readdlm("q-(-0.5,0.3,0.0).csv")); ω4 = vec(readdlm("omega-(-0.5,0.3,0.0).csv")) 
    q5 = vec(readdlm("q-(-0.5,0.4,0.0).csv")); ω5 = vec(readdlm("omega-(-0.5,0.4,0.0).csv"))
    q6 = vec(readdlm("q-(-0.5,0.5,0.0).csv")); ω6 = vec(readdlm("omega-(-0.5,0.5,0.0).csv"))  
    # Plot data
    plot(q1, ω1, label = L"$u_f = 0.0$", linewidth = 2, xlabel = L"$q$", ylabel = L"$\omega$", margin = 3mm, xlims = (0, 5.0), ylims = (-0.5,1.5), grid = true, legend=:bottomright)
    plot!(q2, ω2, label = L"$u_f = 0.1$", linewidth = 2)
    plot!(q3, ω3, label = L"$u_f = 0.2$", linewidth = 2)
    plot!(q4, ω4, label = L"$u_f = 0.3$", linewidth = 2)
    plot!(q5, ω5, label = L"$u_f = 0.4$", linewidth = 2)
    plot!(q6, ω6, label = L"$u_f = 0.5$", linewidth = 2)
    savefig("Fig6b.pdf")
end

# Plot configuration
gr() # Load GR plotting backend
default(titlefont = (18, "Computer Modern"), guidefont = (18, "Computer Modern"), tickfont = (18, "Computer Modern"), legendfontsize = 14)
# Generate plots
@time draw_kappa()
@time draw_uf()