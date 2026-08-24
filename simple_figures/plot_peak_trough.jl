using DelimitedFiles
using DataFrames
using LaTeXStrings
using StatsBase
using MAT
using Measurements
using Peaks
using Plots
using Printf
using Images
using LinearAlgebra


dt_phys = 0.001 # seconds
out_int = 10 # sampling frequency 
dt = dt_phys*out_int

x_oo = readdlm("sample_wave_data.txt", ' '; skipstart=5)
x_oo = x_oo[:,1:3]
x_oo = Float64.(x_oo)
x_oo = x_oo[sortperm(x_oo[:, 1]), :] # reordering array
x_n = x_oo[:,1]
z = x_oo[:,3] .- zero_sea_level # shift back to origin
eta_plot = plot(x_n,z,color=:black,linewidth=2,legend=:bottomleft,label="",legendfontsize=10)

plot!([11.8,11.8],[-0.085,0.23],linewidth=4,color=:red,label="\$H_1\$")
plot!([10.4,11.8],[0.23,0.23],linestyle=:dot,color=:gray,label="")
plot!([10.4,10.4],[-0.148,0.23],linewidth=4,color=:blue,label="\$H_2\$")
plot!([8.35,10.4],[-0.148,-0.148],linestyle=:dot,color=:gray,label="")
plot!([8.35,8.35],[-0.148,0.211],linewidth=4,color=:green,label="\$H_3\$")
plot!([6.4,8.35],[0.211,0.211],linestyle=:dot,color=:gray,label="")
xlabel!("\$x\$ [m]")
ylabel!("\$\\eta(x)\$ [m]")
xlims!(6,13)
savefig(eta_plot,"plots/peaks_and_troughs.pdf")
