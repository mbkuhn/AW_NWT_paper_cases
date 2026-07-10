function get_evenly_spaced_time_series(aw_time,aw_eta)

    dt_array = zeros(length(aw_time))
    for ii=1:length(dt_array)-1
        dt = aw_time[ii+1].-aw_time[ii]
        dt_array[ii] = dt
    end
    t_spacing = median(dt_array)

    time_evenly_spaced = collect(aw_time[1]:t_spacing:aw_time[end])

    aw_time_interp = Interpolations.deduplicate_knots!(aw_time)
    aw_itp = interpolate((aw_time_interp,), aw_eta, Gridded(Linear()))

    eta_aw = aw_itp.(time_evenly_spaced)

    return time_evenly_spaced, eta_aw
end