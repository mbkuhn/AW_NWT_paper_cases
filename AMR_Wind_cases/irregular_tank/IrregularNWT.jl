using DelimitedFiles
using Interpolations
using DataFrames
using LaTeXStrings
using StatsBase
using MAT
using Measurements
include("helper_scripts/FFTAnalysis.jl")
include("helper_scripts/get_evenly_spaced_time_series.jl")

# Define smoothing type and relevant parameters
# Options are: 
# Ensemble Averaging -> smoothing_type = 1; smoothing_parameter = 10000
# Constant Filter (MARIN's Method) -> smoothing_type = 2; smoothing_parameter = 30
# Gaussian -> smoothing_type = 3; smoothing_parameter = 30
smoothing_type = 2
smoothing_parameter = -1
if smoothing_type==1
    smoothing_parameter = 10000
else
    smoothing_parameter = 30
end

file_amrwind = "/Users/mkuhn/testruns_data/OC5_semisubmersible/HOS-NWT/amr-wind-standalone/same_res/output_same_res.txt"
file_hosnwt = "/Users/mkuhn/testruns_data/OC5_semisubmersible/HOS-NWT/amr-wind-standalone/HOS_NWT_case4_probes.dat"

rows_to_skip = 49
# Data from probes for JONSWAP spectrum from HOS-NWT simulation
HOS_NWT_Case4_data = readdlm(file_hosnwt)
HOS_NWT_Case4_data = HOS_NWT_Case4_data[rows_to_skip:end,1:3]
HOS_NWT_Case4_data = Float64.(HOS_NWT_Case4_data)
HOS_NWT_Case4_time = HOS_NWT_Case4_data[:,1]
HOS_NWT_Case4_wave = HOS_NWT_Case4_data[:,3]

# Rewrite HOS_NWT time due to lack of precision
HOS_NWT_dt = HOS_NWT_Case4_time[2] - HOS_NWT_Case4_time[1]
for ii in 1:length(HOS_NWT_Case4_time)
    HOS_NWT_Case4_time[ii] = HOS_NWT_Case4_time[1] + (ii-1)*HOS_NWT_dt
end

# AMR-Wind
amrwind_data = readdlm(file_amrwind)
amrwind_data = Float64.(amrwind_data)
amrwind_time = amrwind_data[:,1]
amrwind_wave = amrwind_data[:,2]

# Plot Waves Time series (Model Scale)
Time_series_plot = plot(HOS_NWT_Case4_time,HOS_NWT_Case4_wave,label="HOS-NWT", color=:black)
plot!(amrwind_time,amrwind_wave,label="AMR-Wind", color=:red, linestyle=:dot)
xlabel!("Time [sec]")
ylabel!("Wave Elevation [m]")
savefig(Time_series_plot,"plotting_outputs/wave_elevation_time_series.pdf")

### **** PSD COMPUTATIONS AND PLOTS ******

f, P1, PSD, fdouble, P2 = FFTAnalysis(HOS_NWT_Case4_time,HOS_NWT_Case4_wave,true,false,smoothing_type,smoothing_parameter)
wave_plot = plot(f[1:end],abs.(PSD[1:end]),yscale=:log10,ylims=(1e-7,1e-2),xlims=(0,1.5),label="HOS-NWT",color=:black)

dt_array = zeros(length(amrwind_data[:,1]))
for ii=1:length(dt_array)-1
    local dt = amrwind_data[ii+1,1].-amrwind_data[ii,1]
    dt_array[ii] = dt
end
t_spacing = median(dt_array)

AMRWind_WAVE_time_interp = Interpolations.deduplicate_knots!(amrwind_time)
itp_waves_amrwind =  interpolate((AMRWind_WAVE_time_interp,), amrwind_wave, Gridded(Linear()))
amr_time_evenly_spaced = collect(amrwind_data[1,1]:t_spacing:amrwind_data[end,1])
amrwind_wave_evenly_spaced = itp_waves_amrwind.(amr_time_evenly_spaced)
AMRWind_time = amr_time_evenly_spaced
AMRWind_Wave = amrwind_wave_evenly_spaced

f, P1, PSD, fdouble, P2 = FFTAnalysis(AMRWind_time,AMRWind_Wave,true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="AMR-Wind Standalone",color=:red,linestyle=:dash)
xlabel!("Frequency [Hz]")
ylabel!("PSD of Incident Wave [m\$^2\$/Hz]")
savefig(wave_plot,"plotting_outputs/wave_PSD.pdf")