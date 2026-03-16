using DelimitedFiles
using Interpolations
using DataFrames
using LaTeXStrings
using StatsBase
using MAT
using Measurements
include("supplemental_scripts/FFTAnalysis.jl")
include("supplemental_scripts/rangedTrapz.jl")
include("supplemental_scripts/get_evenly_spaced_time_series.jl")

# Set relevant frequency ranges based on MARIN's set-up
lowfrequencies_range = [0.005 0.05]
wavefrequencies_range = [0.05519 0.2]

# Define smoothing type and relevant parameters
# Options are: 
# Ensemble Averaging -> smoothing_type = 1; smoothing_parameter = 10000
# Constant Filter (MARIN's Method) -> smoothing_type = 2; smoothing_parameter = 30
# Gaussian -> smoothing_type = 3; smoothing_parameter = 30
smoothing_type = 3
smoothing_parameter = -1
if smoothing_type==1
    smoothing_parameter = 10000
else
    smoothing_parameter = 30
end

file_exawind_pfx = "/Users/mkuhn/testruns_data/offshore_array/results/run"
file_exawind_sfx1 = "/jacket"
file_exawind_sfx2 = "_forces.dat"
file_exawind_labels = ["0","1","2"]
jacket_labels = ["1","2","3","4"]

# NaluWind DATA

NaluWave_Fx = [Float64[],Float64[],Float64[],Float64[]]
NaluWave_Fy = [Float64[],Float64[],Float64[],Float64[]]
NaluWave_Fz = [Float64[],Float64[],Float64[],Float64[]]
NaluWave_Mx = [Float64[],Float64[],Float64[],Float64[]]
NaluWave_My = [Float64[],Float64[],Float64[],Float64[]]
NaluWave_Mz = [Float64[],Float64[],Float64[],Float64[]]
ji = 0
for jn in jacket_labels
    global ji += 1
    file_exawind = file_exawind_pfx * file_exawind_labels[1] * file_exawind_sfx1 * jn * file_exawind_sfx2
    Nalu_data = readdlm(file_exawind,skipstart=2)
    for label in file_exawind_labels[2:end]
        file_exawind_tmp = file_exawind_pfx * label * file_exawind_sfx1 * jn * file_exawind_sfx2
        Nalu_data_tmp = readdlm(file_exawind_tmp,skipstart=1)
        first_time = Nalu_data_tmp[1,1]
        last_index = length(Nalu_data[:,1])
        for ii=1:length(Nalu_data[:,1])
            last_time = Nalu_data[last_index,1]
            if (last_time >= first_time)
                last_index -= 1
            else
                break
            end
        end
        Nalu_data = vcat(Nalu_data[1:last_index,:],Nalu_data_tmp)
    end

    NaluWAVE_data = Float64.(Nalu_data[2:end,:])
    # Make unevenly space time series into an evenly spaced one
    dt_array = zeros(length(NaluWAVE_data[:,1]))
    for ii=1:length(dt_array)-1
        dt = NaluWAVE_data[ii+1,1].-NaluWAVE_data[ii,1]
        dt_array[ii] = dt
    end
    t_spacing = median(dt_array)
    # Create evenly-spaced time array
    time_evenly_spaced = collect(NaluWAVE_data[1,1]:t_spacing:NaluWAVE_data[end,1])

    # Interpolate the unevely-spaced time series
    # Fx
    Fx_Nalu = NaluWAVE_data[:,2].+NaluWAVE_data[:,5]
    _ , Nalu_EvenlySpaced_Fx = get_evenly_spaced_time_series(NaluWAVE_data[:,1],Fx_Nalu)
    # Fy
    Fy_Nalu = NaluWAVE_data[:,3].+NaluWAVE_data[:,6]
    _ , Nalu_EvenlySpaced_Fy = get_evenly_spaced_time_series(NaluWAVE_data[:,1],Fy_Nalu)
    # Fz
    Fz_Nalu = NaluWAVE_data[:,4].+NaluWAVE_data[:,7]
    _ , Nalu_EvenlySpaced_Fz = get_evenly_spaced_time_series(NaluWAVE_data[:,1],Fz_Nalu)
    # Mx
    _ , Nalu_EvenlySpaced_Mx = get_evenly_spaced_time_series(NaluWAVE_data[:,1],NaluWAVE_data[:,8])
    # My
    _ , Nalu_EvenlySpaced_My = get_evenly_spaced_time_series(NaluWAVE_data[:,1],NaluWAVE_data[:,9])
    # Mz
    _ , Nalu_EvenlySpaced_Mz = get_evenly_spaced_time_series(NaluWAVE_data[:,1],NaluWAVE_data[:,10])

    NaluWave_time = time_evenly_spaced

    # In case one needs to discard any transient time data
    t_skip = 240.
    t_final = t_skip + 1e10
    indices = findall(NaluWave_time->(NaluWave_time>=(t_skip) && NaluWave_time<=(t_final)),NaluWave_time)

    global NaluWave_time = NaluWave_time[indices]
    NaluWave_Fx[ji] = (Nalu_EvenlySpaced_Fx[indices])
    NaluWave_Fy[ji] = (Nalu_EvenlySpaced_Fy[indices])
    NaluWave_Fz[ji] = (Nalu_EvenlySpaced_Fz[indices])
    NaluWave_Mx[ji] = (Nalu_EvenlySpaced_Mx[indices])
    NaluWave_My[ji] = (Nalu_EvenlySpaced_My[indices])
    NaluWave_Mz[ji] = (Nalu_EvenlySpaced_Mz[indices])
end

stepRangeNalu = collect(1:length(NaluWave_time)) # include all of the Nalu time

### **** PSD COMPUTATIONS AND PLOTS ******
# Fx
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Fx[1][stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
F_plot = plot(f[1:end],abs.(PSD[1:end]),yscale=:log10,xlims=(0,0.3),ylims=(1e8,1e15),label="Jacket 1 - x",linestyle=:solid,color=:blue,legend=:topright)
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Fx[2][stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="Jacket 2 - x",linestyle=:solid,color=:red)
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Fx[3][stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="Jacket 3 - x",linestyle=:solid,color=:green)
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Fx[4][stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="Jacket 4 - x",linestyle=:solid,color=:black)

# Fy
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Fy[1][stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,xlims=(0,0.4),label="Jacket 1 - y",linestyle=:dash,color=:blue)
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Fy[2][stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label=:none,linestyle=:dash,color=:red)
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Fy[3][stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label=:none,linestyle=:dash,color=:green)
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Fy[4][stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label=:none,linestyle=:dash,color=:black)

# Fz
# remove buoyancy force, assumed to be average heave
buoy = mean(mean(NaluWave_Fz))
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Fz[1][stepRangeNalu] .- buoy,true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,xlims=(0,0.4),label="Jacket 1 - z",linestyle=:dot,color=:blue)
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Fz[2][stepRangeNalu] .- buoy,true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label=:none,linestyle=:dot,color=:red)
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Fz[3][stepRangeNalu] .- buoy,true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label=:none,linestyle=:dot,color=:green)
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Fz[4][stepRangeNalu] .- buoy,true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label=:none,linestyle=:dot,color=:black)

xlabel!("Frequency [Hz]")
ylabel!("Global \$F [N^2/Hz]\$")
savefig(F_plot,"plotting_outputs/F_PSD.pdf")

# Mx
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Mx[1][stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
M_plot = plot(f[1:end],abs.(PSD[1:end]),yscale=:log10,xlims=(0,0.4),ylims=(1e9,1e17),label="Jacket 1 - x",legend=:topright,linestyle=:solid,color=:blue)
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Mx[2][stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="Jacket 2 - x",linestyle=:solid,color=:red)
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Mx[3][stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="Jacket 3 - x",linestyle=:solid,color=:green)
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Mx[4][stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="Jacket 4 - x",linestyle=:solid,color=:black)

# My
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_My[1][stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,xlims=(0,0.4),label="Jacket 1 - y",linestyle=:dash,color=:blue)
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_My[2][stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="",linestyle=:dash,color=:red)
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_My[3][stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="",linestyle=:dash,color=:green)
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_My[4][stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="",linestyle=:dash,color=:black)

# Mz
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Mz[1][stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,xlims=(0,0.4),label="Jacket 1 - z",linestyle=:dot,color=:blue)
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Mz[2][stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="",linestyle=:dot,color=:red)
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Mz[3][stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="",linestyle=:dot,color=:green)
f, P1, PSD, fdouble, P2 = FFTAnalysis(NaluWave_time[stepRangeNalu],NaluWave_Mz[4][stepRangeNalu],true,false,smoothing_type,smoothing_parameter)
plot!(f[1:end],abs.(PSD[1:end]),yscale=:log10,label="",linestyle=:dot,color=:black)

xlabel!("Frequency [Hz]")
ylabel!("Global \$M [(Nm)^2/Hz]\$")
savefig(M_plot,"plotting_outputs/M_PSD.pdf")