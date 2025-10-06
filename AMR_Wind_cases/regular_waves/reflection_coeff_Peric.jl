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

# linear H = 0.0256
# II order H = 0.16
# V order H = 0.3585

# physical parameters for simulation
wave_length  = 4. # meters
wave_height  = 0.0256 # meters (this is through to crest)
wave_period  = 1.6 # seconds
wave_number  = 2*pi/wave_length # meters^-1
wave_frequency = 1/wave_period # Hertz. Note: this is the fundamental frequency: the peaks of the signal will be observed here
angular_frequency = 2*pi*wave_frequency # rad/sec
zero_sea_level = 0.0 # meters
water_depth = 18.0 # meters
final_time = 18. # seconds
dt_phys = 0.001 # seconds
box_length= 24.0 # meters
rel_gen_length = 2*wave_length # 2 wavelengths for longer domain

# getting the data from the gauges
out_int = 10 # sampling frequency 
dt = dt_phys*out_int
nt = Int(final_time/dt) + 1 # number of timesteps / frequency + 1

x = []
eta = []

time_points = 0
space_points = 0
# parse data
valid_ns = [n for n in 1400:(nt-1) if mod(n * out_int, 1) == 0]
for n in valid_ns # Last 4 seconds of the simulations to avoid transient times (flat waveprofiles)
    global time_points = time_points + 1
    x_oo = readdlm("/PATH/TO/LINEARWAVES/post_processing/sampling_fs_domain" * lpad(string(n * out_int), 5, '0') * ".txt", ' '; skipstart=5)
    #x_oo = readdlm("/PATH/TO/SECONDORDERWAVES/post_processing/sampling_fs_domain" * lpad(string(n * out_int), 5, '0') * ".txt", ' '; skipstart=5)
    #x_oo = readdlm("/PATH/TO/FIFTHORDER/post_processing/sampling_fs_domain" * lpad(string(n * out_int), 5, '0') * ".txt", ' '; skipstart=5)
    x_oo = x_oo[:,1:3]
    x_oo = Float64.(x_oo)
    x_oo = x_oo[sortperm(x_oo[:, 1]), :] # reordering array
    x_n = x_oo[:,1]
    global space_points = length(x_n)
    global x = [x;x_n]
    z = x_oo[:,3] .- zero_sea_level # shift back to origin
    global eta = [eta;z]
end

# ensure data is floats
x = Float64.(x)
eta = Float64.(eta)

x = reshape(x,(space_points,time_points))
eta = reshape(eta,(space_points,time_points))

# get only the data that is 1.5lambdas adjacent to the relaxation zone to compute the reflection coefficient
domain_start = box_length - rel_gen_length - 1.5*wave_length
domain_end = box_length - rel_gen_length
x_Cr = []
eta_Cr = []
space_points_Cr = 0
for ii in axes(x,2) # looking at each saved time
    x_col = x[:,ii]
    eta_col = eta[:,ii]
    idx_cr = findall(x_col->((x_col>=domain_start) && (x_col<=domain_end)),x_col)
    global x_Cr = [x_Cr;x_col[idx_cr]]
    global eta_Cr = [eta_Cr;eta_col[idx_cr]]
    global space_points_Cr = length(idx_cr)
end
x_Cr = reshape(x_Cr,(space_points_Cr,time_points))
eta_Cr = reshape(eta_Cr,(space_points_Cr,time_points))

# plot wave elevation at different times
eta_plot = plot(x[:,1],eta[:,1],label="t = 14 sec",legend=:outertop)
forcing_zone = []
eta_rx_zone = []
rx_points = 0
for ii = 0:time_points
    if mod(ii,100)==0
        if ii>=100 && mod(ii,100)==0
            plot!(x[:,ii],eta[:,ii],label=@sprintf("t = %.0f", (14 + ii*dt)) * " sec")
        end
        x_relax = x[:,ii+1]
        eta_relax = eta[:,ii+1]
        relax_zone_location = findall(x_relax->x_relax>=16.0,x_relax)
        global forcing_zone = [forcing_zone;x_relax[relax_zone_location]]
        global eta_rx_zone = [eta_rx_zone;eta_relax[relax_zone_location]]
        global rx_points = length(relax_zone_location)
    end
end
hline!(ones(length(x))*wave_height*0.5,linestyle=:dash,color=:black,label="Theoretical Amplitude")
hline!(ones(length(x))*wave_height*0.5*-1,linestyle=:dash,color=:black,primary=false)
vline!(ones(length(eta))*16,linestyle=:dash,color=:blue,label="Relaxation Zone")
xlabel!("\$x\$ [m]")
ylabel!("\$\\eta(x)\$ [m]")
savefig(eta_plot,"eta.pdf")

# get all peaks for all of the relevant times
# wave height is the distance between a through (peak) and a crest
# within our domain of interest (1.5lambdas-long) we will either have 1 peak and 2 crests, or 1 crests and 2 peaks 
H_list = []
for ii in axes(eta_Cr,2)
    peaks_idx, _ , = findmaxima(eta_Cr[:,ii])
    crests_idx, _ , = findminima(eta_Cr[:,ii])
    max_length = (length(peaks_idx)>=length(crests_idx)) ?  length(peaks_idx) : length(crests_idx)
    for jj=1:max_length
        if (length(peaks_idx)==length(crests_idx))
            H = eta_Cr[peaks_idx[jj],ii] - eta_Cr[crests_idx[jj],ii]
            global H_list = [H_list;H]
        elseif (length(peaks_idx)>length(crests_idx)) && jj!=max_length
            H_1 = eta_Cr[peaks_idx[jj],ii] - eta_Cr[crests_idx[jj],ii]
            H_2 = eta_Cr[peaks_idx[jj+1],ii] - eta_Cr[crests_idx[jj],ii]
            global H_list = [H_list;H_1;H_2]
        elseif (length(peaks_idx)<length(crests_idx)) && jj!=max_length
            H_1 = eta_Cr[peaks_idx[jj],ii] - eta_Cr[crests_idx[jj],ii]
            H_2 = eta_Cr[peaks_idx[jj],ii] - eta_Cr[crests_idx[jj+1],ii]
            global H_list = [H_list;H_1;H_2]
        end
    end
end

# print maximum and minimum heights found and the computed reflection coefficient
indices = findall(H_list->(H_list>(0.2*wave_height) && H_list<(1.2*wave_height)),H_list) # to get rid of outliers that would skew the computation
H_list = H_list[indices]
H_max = maximum(H_list)
println("H_max = ",H_max)
H_min = minimum(H_list)
println("H_min = ",H_min)
C_r = (H_max-H_min)/(H_max+H_min)
println("Reflection Coefficient = ",C_r)
