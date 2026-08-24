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

# choose case (1, 2, 3)
icase = 3

# physical parameters for simulation
wave_length  = 4. # meters
wave_min     = 0.0 # expected peak
wave_max     = 0.0 # expected trough
wave_height  = 0.0 # (specified later)
wave_period  = 1.6 # seconds
wave_number  = 2*pi/wave_length # meters^-1
wave_frequency = 1/wave_period # Hertz. Note: this is the fundamental frequency: the peaks of the signal will be observed here
angular_frequency = 2*pi*wave_frequency # rad/sec
zero_sea_level = 0.0 # meters
water_depth = 18.0 # meters
dt_phys = 0.001 # seconds
box_length= 24.0 # meters
rel_gen_length = 2*wave_length # 2 wavelengths for longer domain

# getting the data from the gauges
initial_time = 5.2
final_time = 5.5 # seconds
out_int = 10 # sampling frequency 
dt = dt_phys*out_int

path_prefix = "/Users/mkuhn/Library/CloudStorage/OneDrive-NLR/testruns_data/MKuhn_regularwavesdata/"

# plot wave elevation at different times
x_oo = readdlm(path_prefix * "/FifthOrder/post_processing/sampling_fs_domain" * lpad(string(Int(5.5 / dt_phys)), 5, '0') * ".txt", ' '; skipstart=5)
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
