# Plot numerical dispersion relations from linear stability analysis
# Alex Tam, 01/06/2022

using DelimitedFiles
using Plots
using LaTeXStrings
using Measures

function draw_kappa()
    # Import data (κ,uf,γ)
    q01 = vec(readdlm("q-(-0.01,0.1,0.1).csv")); ω01 = vec(readdlm("omega-(-0.01,0.1,0.1).csv"))
    # q05 = vec(readdlm("q-(-0.05,0.1,0.1).csv")); ω05 = vec(readdlm("omega-(-0.05,0.1,0.1).csv"))
    q1 = vec(readdlm("q-(-0.1,0.1,0.1).csv")); ω1 = vec(readdlm("omega-(-0.1,0.1,0.1).csv"))  
    q2 = vec(readdlm("q-(-0.2,0.1,0.1).csv")); ω2 = vec(readdlm("omega-(-0.2,0.1,0.1).csv")) 
    # q3 = vec(readdlm("q-(-0.3,0.1,0.1).csv")); ω3 = vec(readdlm("omega-(-0.3,0.1,0.1).csv")) 
    q4 = vec(readdlm("q-(-0.4,0.1,0.1).csv")); ω4 = vec(readdlm("omega-(-0.4,0.1,0.1).csv")) 
    # q5 = vec(readdlm("q-(-0.5,0.1,0.1).csv")); ω5 = vec(readdlm("omega-(-0.5,0.1,0.1).csv")) 
    q6 = vec(readdlm("q-(-0.6,0.1,0.1).csv")); ω6 = vec(readdlm("omega-(-0.6,0.1,0.1).csv")) 
    # q7 = vec(readdlm("q-(-0.7,0.1,0.1).csv")); ω7 = vec(readdlm("omega-(-0.7,0.1,0.1).csv")) 
    # q8 = vec(readdlm("q-(-0.8,0.1,0.1).csv")); ω8 = vec(readdlm("omega-(-0.8,0.1,0.1).csv"))  
    # Plot data
    plot(q01, ω01, label = L"$\kappa = -0.01$", linewidth = 2)
    # plot!(q05, ω05, label = L"$\kappa = -0.05$", linewidth = 2)
    plot!(q1, ω1, label = L"$\kappa = -0.1$", linewidth = 2, xlims = (0, 5.0), ylims = (-1.0,0.5), legend=:bottomleft, xlabel = L"$q$", ylabel = L"$\omega$", margin = 3mm, grid = true)
    plot!(q2, ω2, label = L"$\kappa = -0.2$", linewidth = 2)
    # plot!(q3, ω3, label = L"$\kappa = -0.3$", linewidth = 2)
    plot!(q4, ω4, label = L"$\kappa = -0.4$", linewidth = 2)
    # plot!(q5, ω5, label = L"$\kappa = -0.5$", linewidth = 2)
    plot!(q6, ω6, label = L"$\kappa = -0.6$", linewidth = 2)
    # plot!(q7, ω7, label = L"$\kappa = -0.7$", linewidth = 2)
    # plot!(q8, ω8, label = L"$\kappa = -0.8$", linewidth = 2)
    savefig("Fig8a.pdf")
end

function draw_uf()
    # Import data (κ,uf,γ)
    q0 = vec(readdlm("q-(-0.5,0.0,0.1).csv")); ω0 = vec(readdlm("omega-(-0.5,0.0,0.1).csv"))
    q1 = vec(readdlm("q-(-0.5,0.1,0.1).csv")); ω1 = vec(readdlm("omega-(-0.5,0.1,0.1).csv"))  
    q2 = vec(readdlm("q-(-0.5,0.2,0.1).csv")); ω2 = vec(readdlm("omega-(-0.5,0.2,0.1).csv")) 
    # q3 = vec(readdlm("q-(-0.5,0.3,0.1).csv")); ω3 = vec(readdlm("omega-(-0.5,0.3,0.1).csv")) 
    q4 = vec(readdlm("q-(-0.5,0.4,0.1).csv")); ω4 = vec(readdlm("omega-(-0.5,0.4,0.1).csv")) 
    # q5 = vec(readdlm("q-(-0.5,0.5,0.1).csv")); ω5 = vec(readdlm("omega-(-0.5,0.5,0.1).csv"))
    q6 = vec(readdlm("q-(-0.5,0.6,0.1).csv")); ω6 = vec(readdlm("omega-(-0.5,0.6,0.1).csv"))
    # q7 = vec(readdlm("q-(-0.5,0.7,0.1).csv")); ω7 = vec(readdlm("omega-(-0.5,0.7,0.1).csv"))
    # q8 = vec(readdlm("q-(-0.5,0.8,0.1).csv")); ω8 = vec(readdlm("omega-(-0.5,0.8,0.1).csv"))
    # Plot data
    plot(q0, ω0, label = L"$u_f = 0$", linewidth = 2)
    plot!(q1, ω1, label = L"$u_f = 0.1$", linewidth = 2, xlims = (0, 5.0), ylims = (-1.0,0.5), legend=:topright, xlabel = L"$q$", ylabel = L"$\omega$", margin = 3mm, grid = true)
    plot!(q2, ω2, label = L"$u_f = 0.2$", linewidth = 2)
    # plot!(q3, ω3, label = L"$u_f = 0.3$", linewidth = 2)
    plot!(q4, ω4, label = L"$u_f = 0.4$", linewidth = 2)
    # plot!(q5, ω5, label = L"$u_f = 0.5$", linewidth = 2)
    plot!(q6, ω6, label = L"$u_f = 0.6$", linewidth = 2)
    # plot!(q7, ω7, label = L"$u_f = 0.7$", linewidth = 2)
    # plot!(q8, ω8, label = L"$u_f = 0.8$", linewidth = 2)
    savefig("Fig8b.pdf")
end

function draw_gamma()
    # Import data (κ,uf,γ)
    q0 = vec(readdlm("q-(-0.5,0.1,0.0).csv")); ω0 = vec(readdlm("omega-(-0.5,0.1,0.0).csv"))
    q01 = vec(readdlm("q-(-0.5,0.1,0.01).csv")); ω01 = vec(readdlm("omega-(-0.5,0.1,0.01).csv")) 
    # q05 = vec(readdlm("q-(-0.5,0.1,0.05).csv")); ω05 = vec(readdlm("omega-(-0.5,0.1,0.05).csv"))
    q1 = vec(readdlm("q-(-0.5,0.1,0.1).csv")); ω1 = vec(readdlm("omega-(-0.5,0.1,0.1).csv"))
    q2 = vec(readdlm("q-(-0.5,0.1,0.2).csv")); ω2 = vec(readdlm("omega-(-0.5,0.1,0.2).csv"))
    # q3 = vec(readdlm("q-(-0.5,0.1,0.3).csv")); ω3 = vec(readdlm("omega-(-0.5,0.1,0.3).csv"))
    q4 = vec(readdlm("q-(-0.5,0.1,0.4).csv")); ω4 = vec(readdlm("omega-(-0.5,0.1,0.4).csv"))
    # q5 = vec(readdlm("q-(-0.5,0.1,0.5).csv")); ω5 = vec(readdlm("omega-(-0.5,0.1,0.5).csv"))
    q6 = vec(readdlm("q-(-0.5,0.1,0.6).csv")); ω6 = vec(readdlm("omega-(-0.5,0.1,0.6).csv"))
    # q7 = vec(readdlm("q-(-0.5,0.1,0.7).csv")); ω7 = vec(readdlm("omega-(-0.5,0.1,0.7).csv"))
    # q8 = vec(readdlm("q-(-0.5,0.1,0.8).csv")); ω8 = vec(readdlm("omega-(-0.5,0.1,0.8).csv"))
    # Plot data
    plot(q0, ω0, label = L"$\gamma = 0$", linewidth = 2, xlims = (0, 5.0), ylims = (-1.5,1.5), legend=:right, xlabel = L"$q$", ylabel = L"$\omega$", margin = 3mm,  grid = true)
    plot!(q01, ω01, label = L"$\gamma = 0.01$", linewidth = 2)
    # plot!(q05, ω05, label = L"$\gamma = 0.05$", linewidth = 2)
    plot!(q1, ω1, label = L"$\gamma = 0.1$", linewidth = 2)
    plot!(q2, ω2, label = L"$\gamma = 0.2$", linewidth = 2)
    # plot!(q3, ω3, label = L"$\gamma = 0.3$", linewidth = 2)
    plot!(q4, ω4, label = L"$\gamma = 0.4$", linewidth = 2)
    # plot!(q5, ω5, label = L"$\gamma = 0.5$", linewidth = 2)
    plot!(q6, ω6, label = L"$\gamma = 0.6$", linewidth = 2)
    # plot!(q7, ω7, label = L"$\gamma = 0.7$", linewidth = 2)
    # plot!(q8, ω8, label = L"$\gamma = 0.8$", linewidth = 2)
    savefig("Fig8c.pdf")
end

# Plot configuration
gr() # Load GR plotting backend
default(titlefont = (18, "Computer Modern"), guidefont = (18, "Computer Modern"), tickfont = (18, "Computer Modern"), legendfontsize = 14) # Fonts

# Generate plots
@time draw_kappa()
@time draw_uf()
@time draw_gamma()