using DelimitedFiles
using Statistics
using Polynomials
using Plots
using DataFrames
using LaTeXStrings


# physical parameters for simulation
wave_length  = 4. # meters
wave_height  = 0.3585 # meters (this is through to crest)
wave_period  = 1.6 # seconds
wave_number  = 2*pi/wave_length # meters^-1
wave_frequency = 1/wave_period # Hertz. Note: this is the fundamental frequency: the peaks of the signal will be observed here
angular_frequency = 2*pi*wave_frequency # rad/sec
zero_sea_level = 0.0 # meters
water_depth = 18.0 # meters
final_time = 18. # seconds
dt_phys = 0.001 # seconds

# getting the data from the gauges
out_int = 10 # sampling frequency # keeping it at 10
nt = Int(final_time*out_int^2) + 1 # number of timesteps / frequency + 1
z1_list = Vector{Real}(undef,nt)
z2_list = Vector{Real}(undef,nt)
z3_list = Vector{Real}(undef,nt)
z4_list = Vector{Real}(undef,nt)
z5_list = Vector{Real}(undef,nt)
z6_list = Vector{Real}(undef,nt)

x1 = 0
x2 = 0
x3 = 0
x4 = 0
x5 = 0
x6 = 0

path_prefix = "/Users/mkuhn/testruns_data/VofConvergenceData/vof/post_processing"
ow_path_prefix = "/Users/mkuhn/testruns_data/VofConvergenceData/ow_vof/post_processing"
local_prefix = "/Users/mkuhn/testruns_data/stokes_waves/5th_order_paper_runs_convergence/lev-"


for n in 0:(nt-1)

    # ow_vof
    x_oo_1 = readdlm(ow_path_prefix*"/sampling_fs_gauge" * lpad(string(n * out_int), 5, '0') * ".txt", ' '; skipstart=5)
    x_oo_1 = x_oo_1[1:3]
    x_oo_1 = Float64.(x_oo_1)
    global x1 = x_oo_1[1]
    z1 = x_oo_1[3] - zero_sea_level # shift back to origin
    global z1_list[n+1] = z1

    # vof, two levels of refinement
    x_oo_2 = readdlm(path_prefix*"_2levs/sampling_fs_gauge" * lpad(string(n * out_int), 5, '0') * ".txt", ' '; skipstart=5)
    x_oo_2 = x_oo_2[1:3]
    x_oo_2 = Float64.(x_oo_2)
    global x2 = x_oo_2[1]
    z2 = x_oo_2[3] - zero_sea_level # shift back to origin
    global z2_list[n+1] = z2

    # vof, one level of ref
    x_oo_3 = readdlm(path_prefix*"_1levs/sampling_fs_gauge" * lpad(string(n * out_int), 5, '0') * ".txt", ' '; skipstart=5)
    x_oo_3 = x_oo_3[1:3]
    x_oo_3 = Float64.(x_oo_3)
    global x3 = x_oo_3[1]
    z3 = x_oo_3[3] - zero_sea_level # shift back to origin
    global z3_list[n+1] = z3

    # vof, zero level of ref
    x_oo_4 = readdlm(path_prefix*"_0levs/sampling_fs_gauge" * lpad(string(n * out_int), 5, '0') * ".txt", ' '; skipstart=5)
    x_oo_4 = x_oo_4[1:3]
    x_oo_4 = Float64.(x_oo_4)
    global x4 = x_oo_4[1]
    z4 = x_oo_4[3] - zero_sea_level # shift back to origin
    global z4_list[n+1] = z4

    # vof, -1 level of ref
    x_oo_5 = readdlm(local_prefix*"1/post_processing/sampling_fs_gauge" * lpad(string(n * out_int), 5, '0') * ".txt", ' '; skipstart=5)
    x_oo_5 = x_oo_5[1:3]
    x_oo_5 = Float64.(x_oo_5)
    global x5 = x_oo_5[1]
    z5 = x_oo_5[3] - zero_sea_level # shift back to origin
    global z5_list[n+1] = z5

    # vof, -2 level of ref
    x_oo_6 = readdlm(local_prefix*"2/post_processing/sampling_fs_gauge" * lpad(string(n * out_int), 5, '0') * ".txt", ' '; skipstart=5)
    x_oo_6 = x_oo_6[1:3]
    x_oo_6 = Float64.(x_oo_6)
    global x6 = x_oo_6[1]
    z6 = x_oo_6[3] - zero_sea_level # shift back to origin
    global z6_list[n+1] = z6

end

# Convert to arrays of floats
z1_list = Float64.(z1_list)
z2_list = Float64.(z2_list)
z3_list = Float64.(z3_list)
z4_list = Float64.(z4_list)
z5_list = Float64.(z5_list)
z6_list = Float64.(z6_list)

dt = dt_phys*out_int # this dt is actually the sampling rate: dt*10 because we are sampling every 10 time steps
t_physical = range(0, step=dt, stop=final_time) # Note that we might not have enough frequency samples.
t_skip = 10.0
times_idx = findall(t_physical->(t_physical>=t_skip),t_physical)

diffs_two_levs = (z1_list[times_idx] .- z2_list[times_idx])
diffs_one_levs = (z1_list[times_idx] .- z3_list[times_idx])
diffs_zero_levs = (z1_list[times_idx] .- z4_list[times_idx])
diffs_minus_1_levs = (z1_list[times_idx] .- z5_list[times_idx])
diffs_minus_2_levs = (z1_list[times_idx] .- z6_list[times_idx])

sum_err_two_levs = sum(abs.(diffs_two_levs))
sum_err_one_levs = sum(abs.(diffs_one_levs))
sum_err_zero_levs = sum(abs.(diffs_zero_levs))
sum_err_minus_1_levs = sum(abs.(diffs_minus_1_levs))
sum_err_minus_2_levs = sum(abs.(diffs_minus_2_levs))
error_vec = [sum_err_minus_2_levs, sum_err_minus_1_levs, sum_err_zero_levs, sum_err_one_levs, sum_err_two_levs]

dx_minus_2 = 4.0*0.0625
dx_minus_1 = 2.0*0.0625
dx_zero = 0.0625
dx_one = 0.5*0.0625
dx_two = 0.25*0.0625
dx_vec = [dx_minus_2, dx_minus_1, dx_zero, dx_one, dx_two]

# Linear fit in log-log space
log_dx = log10.(dx_vec)
log_err = log10.(error_vec)
p = Polynomials.fit(log_dx, log_err, 1)
slope = p.coeffs[2]  # convergence rate
intercept = p.coeffs[1]
println("Convergence Rate = ",slope)

# Evaluate fitted line
fit_line = 1.0 .* log_dx .+ (log_err[1] - log_dx[1]) .- 0.05
fit_line2 = 2.0 .* log_dx .+ (log_err[1] - 2. * log_dx[1])
fit_errors = (10 .^ fit_line)  # back-transform to linear space
fit_errors2 = (10 .^ fit_line2)


wave_profiles = plot(t_physical,z1_list,lw=1.5,color=:black,label="Target",legend=:outerbottom,legend_column=3)
plot!(t_physical,z2_list,lw=1.5,color=:red,linestyle=:solid,label="\$H = 23 \\Delta x\$")
plot!(t_physical,z3_list,lw=1.5,color=:blue,linestyle=:dash,label="\$H = 11.5 \\Delta x\$")
plot!(t_physical,z4_list,lw=1.5,color=:green,linestyle=:dashdot,label="\$H = 5.5 \\Delta x\$")
plot!(t_physical,z5_list,lw=1.5,color=:orange,linestyle=:dashdotdot,label="\$H = 3 \\Delta x\$")
plot!(t_physical,z6_list,lw=1.5,color=:purple,linestyle=:dot,label="\$H = 1.5 \\Delta x\$")
xlabel!("\$t\$ [sec]")
ylabel!("\$\\eta(x_\\mathrm{p}=10,t)\$ [m]")
savefig(wave_profiles,"../plotting_outputs/vof_vs_owvof.pdf")

error_plot = plot(t_physical[times_idx],diffs_two_levs,color=:red,lw=2,linestyle=:solid,label="\$H = 23 \\Delta x\$",legend=:outerbottom,legend_column=3)
plot!(t_physical[times_idx],diffs_one_levs,color=:blue,lw=2,linestyle=:dash,label="\$H = 11.5 \\Delta x\$")
plot!(t_physical[times_idx],diffs_zero_levs,color=:green,lw=2,linestyle=:dashdot,label="\$H = 5.5 \\Delta x\$")
plot!(t_physical[times_idx],diffs_minus_1_levs,color=:orange,lw=2,linestyle=:dashdotdot,label="\$H = 3 \\Delta x\$")
plot!(t_physical[times_idx],diffs_minus_2_levs,color=:purple,lw=2,linestyle=:dot,label="\$H = 1.5 \\Delta x\$")
xlabel!("\$t\$ [sec]")
ylabel!("\$\\eta_{\\mathrm{target}}(x_\\mathrm{p}=10,t) - \\eta(x_\\mathrm{p}=10,t) \$ [m]")
savefig(error_plot,"../plotting_outputs/error.pdf")

norm_plot = plot(dx_vec, error_vec, marker = :circle, lw=2, color=:red, linestyle = :solid, label = "Error", xscale = :log10, yscale = :log10)
plot!(dx_vec,fit_errors,lw=2, linestyle = :dash,color=:black,label="1st order")
#plot!(dx_vec,fit_errors2,lw=2, linestyle = :dot,color=:black,label="2nd order")
xlabel!("\$\\Delta{x}\$ [m]")
ylabel!("\$\\Sigma|\\eta_{\\mathrm{target}}(t) - \\eta(t)|\$ [m]")
# title!("Shape Error \$(L_1)\$ at \$x=10\$ for \$t\\in[10,18]sec\$")
savefig(norm_plot,"../plotting_outputs/norm_plot.pdf")
