import marimo

__generated_with = "0.15.2"
app = marimo.App()


@app.cell
def _():
    import numpy as np
    import matplotlib.pyplot as plt 
    plt.rc('text',usetex=True)
    plt.rc('font',family='serif')

    # List of some relevant physical variables
    #Lx = 9.3717632253682877E+03
    #Ly = Lx
    #depth = 1000
    #rho_l = 1000
    #g = 9.81
    return np, plt


@app.cell
def _():
    path = '/Users/mpolimen/Desktop/data_to_process'
    return (path,)


@app.cell
def _(np):
    import re

    def parse_file(filepath):
        data_by_time = {}
        kx_ky_ref = None
        current_time = None
        a_eta_list = []

        zone_pattern = re.compile(r'ZONE T = "\s*([0-9.Ee+-]+)"')

        parsing_first_zone = True
        kx_list, ky_list = [], []

        with open(filepath, 'r') as f:
            for line in f:
                # Check for start of new zone
                zone_match = zone_pattern.match(line)
                if zone_match:
                    # Save data from previous zone
                    if current_time is not None and a_eta_list:
                        a_eta_array = np.array(a_eta_list)

                        # If first zone, build kx_ky_ref
                        if parsing_first_zone:
                            kx_arr = np.array(kx_list)
                            ky_arr = np.array(ky_list)
                            kx_ky_ref = np.column_stack((kx_arr, ky_arr))

                            parsing_first_zone = False  # Done parsing the first zone

                        # Validate and stack data
                        if len(a_eta_array) != kx_ky_ref.shape[0]:
                            raise ValueError(
                                f"Mismatch between kx/ky length ({kx_ky_ref.shape[0]}) "
                                f"and a_eta length ({len(a_eta_array)}) at time {current_time}"
                            )

                        full_data = np.column_stack((kx_ky_ref, a_eta_array))
                        data_by_time[current_time] = full_data

                        a_eta_list = []  # Reset for next zone

                    current_time = float(zone_match.group(1))
                    continue

                parts = line.strip().split()

                # Parse line based on zone
                if parsing_first_zone and len(parts) >= 3:
                    try:
                        kx, ky, a_eta = map(float, parts[:3])
                        kx_list.append(kx)
                        ky_list.append(ky)
                        a_eta_list.append(a_eta)
                    except ValueError:
                        continue

                elif not parsing_first_zone and len(parts) >= 1:
                    try:
                        a_eta = float(parts[0])
                        a_eta_list.append(a_eta)
                    except ValueError:
                        continue

            # Handle the last zone after loop ends
            if current_time is not None and a_eta_list:
                a_eta_array = np.array(a_eta_list)

                if kx_ky_ref is None:
                    raise RuntimeError("No kx/ky reference data found.")

                if len(a_eta_array) != kx_ky_ref.shape[0]:
                    raise ValueError(
                        f"Mismatch between kx/ky length ({kx_ky_ref.shape[0]}) "
                        f"and a_eta length ({len(a_eta_array)}) at time {current_time}"
                    )

                full_data = np.column_stack((kx_ky_ref, a_eta_array))
                data_by_time[current_time] = full_data

        return data_by_time

    # Usage:
    data = parse_file("/Users/mpolimen/Desktop/data_to_process/HOS/a_3d.dat") # full domain modes

    def get_data(time_value):
        return data.get(time_value, None)
    return (get_data,)


@app.cell
def _(np):
    # list of available times
    dt = 1.225
    times = np.arange(0,1599.85+dt,dt) # this array contains the relevant time for the HOS data. Last AMR-Wind time is 1600s
    time_wanted = 1600.0 # change this to update time to compare for spectrum
    time_to_check = np.argmin(np.abs(times-time_wanted))
    my_time = times[time_to_check]
    print(my_time)
    return


@app.cell
def _(get_data):
    _zone_data = get_data(1599.85)
    if _zone_data is not None:
        print(_zone_data.shape)
        print(_zone_data[:5])
    else:
        print('Time not found.')
    return


@app.cell
def _(get_data, plt):
    zone_time = 1599.85
    _zone_data = get_data(zone_time)
    zone_wave_numbers = get_data(zone_time)
    if _zone_data is None:
        print(f'No data found for time {zone_time}')
    else:
        I, J = (129, 256)
        kx = zone_wave_numbers[:, 0].reshape(J, I)
        ky = zone_wave_numbers[:, 1].reshape(J, I)
        a_eta = _zone_data[:, 2].reshape(J, I)
        plt.figure(figsize=(8, 6))
        pcm = plt.pcolormesh(kx, ky, a_eta, shading='auto', cmap='viridis')
        plt.colorbar(pcm, label='E0')
        plt.xlabel('kx')
        plt.ylabel('ky')
        plt.title(f'E0 at t = {zone_time:.5f}')
        plt.tight_layout()
        plt.show()
    return a_eta, kx


@app.cell
def _(np, path):
    Dfinal_three_levels = np.genfromtxt(path + '/post_processing/three_levels_of_ref/sampling16000_fs.txt', skip_header=2, delimiter='')
    Dfinal_wholedom_fivelevs = np.genfromtxt(path + '/post_processing/whole_dom_five_levs/sampling32000_fs.txt', skip_header=2, delimiter='')

    _nx = 256
    _ny = 256

    X_256 = np.zeros((_nx, _ny))
    Y_256 = np.zeros((_nx, _ny))

    X_256_threelevs = np.zeros((_nx, _ny))
    Y_256_threelevs = np.zeros((_nx, _ny))

    Eta_CFD_three_levels = np.zeros((_nx, _ny))
    Eta_CFD_wholedom_fivelevs = np.zeros((_nx, _ny))

    for i in range(_nx):
        for _j in range(_ny):
            X_256[i, _j] = Dfinal_wholedom_fivelevs[i + _j * _nx, 0]
            Y_256[i, _j] = Dfinal_wholedom_fivelevs[i + _j * _nx, 1]

            X_256_threelevs[i, _j] = Dfinal_three_levels[i + _j * _nx, 0]
            Y_256_threelevs[i, _j] = Dfinal_three_levels[i + _j * _nx, 1]

            Eta_CFD_three_levels[i, _j] = Dfinal_three_levels[i + _j * _nx, 2]
            Eta_CFD_wholedom_fivelevs[i, _j] = Dfinal_wholedom_fivelevs[i + _j * _nx, 2]
    return Eta_CFD_three_levels, Eta_CFD_wholedom_fivelevs, X_256


@app.cell
def _(Eta_CFD_three_levels, Eta_CFD_wholedom_fivelevs, X_256, np):
    from scipy.signal import welch, csd, butter, lfilter, freqz
    from scipy.stats import gamma

    def PSD(frequency, TimeSeries, NumbModes, nfft):
        f, psd = welch(TimeSeries, fs=frequency, nperseg=NumbModes, scaling='spectrum', nfft=256)
        return (f, psd)
    _nx = 256
    _ny = _nx

    # works for all the cases with the domain of typical length
    Nperseg_three_levels = 256
    SamplingFreq = 2 * np.pi / (X_256[1, 0] - X_256[0, 0])

    kxCFD_three_levels, PhixCFDTotal_three_levels = PSD(SamplingFreq, Eta_CFD_three_levels[:, 0], Nperseg_three_levels, 256)
    kxCFD_wholedom_fivelevs, PhixCFDTotal_wholedom_fivelevs = PSD(SamplingFreq, Eta_CFD_wholedom_fivelevs[:, 0], Nperseg_three_levels, 256)

    for _j in range(1, _ny):
        # three levels
        kxCFD_three_levels, PhixCFD_three_levels = PSD(SamplingFreq, Eta_CFD_three_levels[:, _j], Nperseg_three_levels, 256)
        PhixCFDTotal_three_levels = PhixCFDTotal_three_levels + PhixCFD_three_levels
        # five levels
        kxCFD_wholedom_fivelevs, PhixCFD_wholedom_fivelevs= PSD(SamplingFreq, Eta_CFD_wholedom_fivelevs[:, _j], Nperseg_three_levels, 256)
        PhixCFDTotal_wholedom_fivelevs = PhixCFDTotal_wholedom_fivelevs + PhixCFD_wholedom_fivelevs


    PhixCFDTotal_256_three_levels = PhixCFDTotal_three_levels / (2 * _ny)
    PhixCFDTotal_wholedom_fivelevs = PhixCFDTotal_wholedom_fivelevs / (2 * _ny)
    return (
        PhixCFDTotal_256_three_levels,
        PhixCFDTotal_wholedom_fivelevs,
        kxCFD_three_levels,
        kxCFD_wholedom_fivelevs,
    )


@app.cell
def _(
    PhixCFDTotal_256_three_levels,
    PhixCFDTotal_wholedom_fivelevs,
    a_eta,
    kx,
    kxCFD_three_levels,
    kxCFD_wholedom_fivelevs,
    plt,
):
    _fig = plt.figure(2)
    plt.title('One dimensional spectra, t = 1600 s')
    plt.plot(kxCFD_three_levels, PhixCFDTotal_256_three_levels, label='CFD solution, $\\Delta x$ = 4.5 m (three levels)')
    plt.plot(kxCFD_wholedom_fivelevs, PhixCFDTotal_wholedom_fivelevs, label='CFD solution, $\\Delta x$ = 2.25 m (five levels)')
    plt.plot(kx[127, :], a_eta[127, :], 'w', label='HOS state at t=1600')
    plt.xlabel('$k_x$', fontsize=14)
    plt.ylabel('E', fontsize=14)
    plt.legend()
    plt.savefig('one_dimensional_spectra_x.png', bbox_inches=None)
    plt.legend(fontsize=14, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.0)
    return


@app.cell
def _(np, path):
    ## Load the HOS data
    HOSEnergy = np.genfromtxt(path+"/HOS/vol_energy.dat",skip_header=67)
    return (HOSEnergy,)


@app.cell
def _(np, path):
    ## Load the CFD data
    CFDE_array256_threelev = np.genfromtxt(path+"/post_processing/three_levels_of_ref/we00000.txt",skip_header=1,delimiter='')

    CFDE_mech256_threelev = CFDE_array256_threelev[:,2]+CFDE_array256_threelev[:,3]

    CFD_time256_threelev = CFDE_array256_threelev[:,1]
    return CFDE_mech256_threelev, CFD_time256_threelev


@app.cell
def _(np, path):
    ## Load the CFD data
    CFDE_array256_wholedom_fivelev = np.genfromtxt(path+"/post_processing/whole_dom_five_levs/we00000.txt",skip_header=1,delimiter='')

    CFDE_mech256_wholedom_fivelev = CFDE_array256_wholedom_fivelev[:,2]+CFDE_array256_wholedom_fivelev[:,3]

    CFD_time256_wholedom_fivelev = CFDE_array256_wholedom_fivelev[:,1]
    return CFDE_mech256_wholedom_fivelev, CFD_time256_wholedom_fivelev


@app.cell
def _(
    CFDE_mech256_threelev,
    CFDE_mech256_wholedom_fivelev,
    CFD_time256_threelev,
    CFD_time256_wholedom_fivelev,
    HOSEnergy,
    plt,
):
    _fig = plt.figure(1)
    ME_ref = HOSEnergy[0, 4]
    plt.plot(HOSEnergy[:, 0], HOSEnergy[:, 4] / ME_ref, 'w', label='HOS reference solution')
    plt.plot(CFD_time256_wholedom_fivelev, CFDE_mech256_wholedom_fivelev / CFDE_mech256_wholedom_fivelev[0], label='AMR-Wind, $\\Delta x$ = 2.28 m (five levels)')
    plt.plot(CFD_time256_threelev, CFDE_mech256_threelev / CFDE_mech256_threelev[0], label='AMR-Wind, $\\Delta x$ = 4.5 m (three levels)')
    plt.xlim(0, 1600)
    plt.xlabel('time [s]', fontsize=16)
    plt.ylabel('$E/E_0$', fontsize=16)
    plt.legend(fontsize=14, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.0)
    plt.title('Mechanical Energy Evolution')
    plt.ylim([0.5, 1.2])
    plt.show()
    plt.savefig('Energy.png', format='png', dpi=100)
    return


if __name__ == "__main__":
    app.run()
