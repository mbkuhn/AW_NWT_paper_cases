import marimo

__generated_with = "0.16.4"
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
    path = '/Users/mkuhn/testruns_data/IrregularWavesDataOcean/'
    return (path,)


@app.cell
def _(np, path):
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
    data = parse_file(path + "../IrregularWavesDataOcean_discard/HOS/a_3d.dat") # full domain modes

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
    return


#@app.cell
#def _(get_data):
#    _zone_data = get_data(1599.85)
#    if _zone_data is None:
#        print('Time not found.')
#    return

@app.cell
def _(get_data, plt, np, path):
    from scipy.signal import welch, csd, butter, lfilter, freqz
    from scipy.stats import gamma

    _nx = 256
    _ny = 256

    def PSD(frequency, TimeSeries, NumbModes, nfft):
        f, psd = welch(TimeSeries, fs=frequency, nperseg=NumbModes, scaling='spectrum', nfft=256)
        return (f, psd)

    def make_spectrum_plot(itime):
        # HOS section - getting data
        zone_time = 999.6 #1599.85
        _zone_data = get_data(zone_time)
        zone_wave_numbers = get_data(zone_time)
        if _zone_data is None:
            print(f'No data found for time {zone_time}')
        else:
            I, J = (129, 256)
            kx = zone_wave_numbers[:, 0].reshape(J, I)
            ky = zone_wave_numbers[:, 1].reshape(J, I)
            a_eta = _zone_data[:, 2].reshape(J, I)
        #return a_eta, kx

        # SGF section - getting data
        Dfinal_three_levels = np.genfromtxt(path + '/SGF/3_ref_levels/sampling05000_fs.txt', skip_header=2, delimiter='')
        Dfinal_four_levels = np.genfromtxt(path + '/SGF/4_ref_levels/sampling10000_fs.txt', skip_header=2, delimiter='')
        Dfinal_four_levels_AR2x = np.genfromtxt(path + '/SGF/4_ref_levels_AR2x/sampling20000_fs.txt', skip_header=2, delimiter='')
        Dfinal_five_levels = np.genfromtxt(path + '/SGF/5_ref_levels/sampling20000_fs.txt', skip_header=2, delimiter='')

        X_256 = np.zeros((_nx, _ny))
        Y_256 = np.zeros((_nx, _ny))

        Eta_CFD_three_levels = np.zeros((_nx, _ny))
        Eta_CFD_four_levels = np.zeros((_nx, _ny))
        Eta_CFD_four_levels_AR2x = np.zeros((_nx, _ny))
        Eta_CFD_five_levels = np.zeros((_nx, _ny))

        for i in range(_nx):
            for _j in range(_ny):
                X_256[i, _j] = Dfinal_five_levels[i + _j * _nx, 0]
                Y_256[i, _j] = Dfinal_five_levels[i + _j * _nx, 1]

                Eta_CFD_three_levels[i, _j] = Dfinal_three_levels[i + _j * _nx, 2]
                Eta_CFD_four_levels[i, _j] = Dfinal_four_levels[i + _j * _nx, 2]
                Eta_CFD_four_levels_AR2x[i, _j] = Dfinal_four_levels_AR2x[i + _j * _nx, 2]
                Eta_CFD_five_levels[i, _j] = Dfinal_five_levels[i + _j * _nx, 2]
        
        # Calculate PSD

        # works for all the cases with the domain of typical length
        Nperseg_three_levels = 256 # plan to remove this
        SamplingFreq = 2 * np.pi / (X_256[1, 0] - X_256[0, 0])

        kxCFD_three_levels, PhixCFDTotal_three_levels = PSD(SamplingFreq, Eta_CFD_three_levels[:, 0], Nperseg_three_levels, 256)
        kxCFD_four_levels, PhixCFDTotal_four_levels = PSD(SamplingFreq, Eta_CFD_four_levels[:, 0], Nperseg_three_levels, 256)
        kxCFD_four_levels_AR2x, PhixCFDTotal_four_levels_AR2x = PSD(SamplingFreq, Eta_CFD_four_levels_AR2x[:, 0], Nperseg_three_levels, 256)
        kxCFD_five_levels, PhixCFDTotal_five_levels = PSD(SamplingFreq, Eta_CFD_five_levels[:, 0], Nperseg_three_levels, 256)

        for _j in range(1, _ny):
            # three levels
            kxCFD_three_levels, PhixCFD_three_levels = PSD(SamplingFreq, Eta_CFD_three_levels[:, _j], Nperseg_three_levels, 256)
            PhixCFDTotal_three_levels = PhixCFDTotal_three_levels + PhixCFD_three_levels
            # four levels
            kxCFD_four_levels, PhixCFD_four_levels = PSD(SamplingFreq, Eta_CFD_four_levels[:, _j], Nperseg_three_levels, 256)
            PhixCFDTotal_four_levels = PhixCFDTotal_four_levels + PhixCFD_four_levels
            # four levels, AR 2x
            kxCFD_four_levels_AR2x, PhixCFD_four_levels_AR2x = PSD(SamplingFreq, Eta_CFD_four_levels_AR2x[:, _j], Nperseg_three_levels, 256)
            PhixCFDTotal_four_levels_AR2x = PhixCFDTotal_four_levels_AR2x + PhixCFD_four_levels_AR2x
            # five levels
            kxCFD_five_levels, PhixCFD_five_levels = PSD(SamplingFreq, Eta_CFD_five_levels[:, _j], Nperseg_three_levels, 256)
            PhixCFDTotal_five_levels = PhixCFDTotal_five_levels + PhixCFD_five_levels

        PhixCFDTotal_three_levels = PhixCFDTotal_three_levels / (_ny)
        PhixCFDTotal_four_levels = PhixCFDTotal_four_levels / (_ny)
        PhixCFDTotal_four_levels_AR2x = PhixCFDTotal_four_levels_AR2x / (_ny)
        PhixCFDTotal_five_levels = PhixCFDTotal_five_levels / (_ny)

        # Make plot
        _fig = plt.figure()
        plt.title('One dimensional spectra, t = 1000 s')
        plt.plot(kxCFD_three_levels, PhixCFDTotal_three_levels, label='$\\Delta z$ = 4.7 m, AR = 2')
        plt.plot(kxCFD_four_levels, PhixCFDTotal_four_levels, label='$\\Delta z$ = 2.3 m, AR = 2')
        plt.plot(kxCFD_four_levels_AR2x, PhixCFDTotal_four_levels_AR2x, label='$\\Delta z$ = 1.2 m, AR = 4')
        plt.plot(kxCFD_five_levels, PhixCFDTotal_five_levels, label='$\\Delta z$ = 1.2 m, AR = 2')
        plt.plot(kx[127, :], a_eta[127, :], 'k', label='HOS-Ocean')
        plt.xlabel('$k_x$', fontsize=14)
        plt.ylabel('E', fontsize=14)
        plt.legend()
        plt.savefig('plotting_outputs/one_dimensional_spectra_x_1000s.png', bbox_inches='tight', dpi=300)
        
        return
    
    return make_spectrum_plot

@app.cell
def _(make_spectrum_plot):
    make_spectrum_plot(0)

@app.cell
def _(np, path):
    ## Load the HOS data
    HOSEnergy = np.genfromtxt(path+"/HOS/vol_energy.dat",skip_header=67)
    return (HOSEnergy,)


@app.cell
def _(np, path):
    ## Load the CFD data
    CFDE_array256_threelev = np.genfromtxt(path+"/SGF/3_ref_levels/we00000.txt",skip_header=1,delimiter='')

    CFDE_mech256_threelev = CFDE_array256_threelev[:,2]+CFDE_array256_threelev[:,3]

    CFD_time256_threelev = CFDE_array256_threelev[:,1]
    return CFDE_mech256_threelev, CFD_time256_threelev

@app.cell
def _(np, path):
    ## Load the CFD data
    CFDE_array256_fourlev = np.genfromtxt(path+"/SGF/4_ref_levels/we00000.txt",skip_header=1,delimiter='')

    CFDE_mech256_fourlev = CFDE_array256_fourlev[:,2]+CFDE_array256_fourlev[:,3]

    CFD_time256_fourlev = CFDE_array256_fourlev[:,1]
    return CFDE_mech256_fourlev, CFD_time256_fourlev

@app.cell
def _(np, path):
    ## Load the CFD data
    CFDE_array256_fourlev_AR2x = np.genfromtxt(path+"/SGF/4_ref_levels_AR2x/we00000.txt",skip_header=1,delimiter='')

    CFDE_mech256_fourlev_AR2x = CFDE_array256_fourlev_AR2x[:,2]+CFDE_array256_fourlev_AR2x[:,3]

    CFD_time256_fourlev_AR2x = CFDE_array256_fourlev_AR2x[:,1]
    return CFDE_mech256_fourlev_AR2x, CFD_time256_fourlev_AR2x

@app.cell
def _(np, path):
    ## Load the CFD data
    CFDE_array256_fivelev = np.genfromtxt(path+"/SGF/5_ref_levels/we00000.txt",skip_header=1,delimiter='')
    CFDE_array256_fivelev = np.append(CFDE_array256_fivelev, np.genfromtxt(path+"/SGF/5_ref_levels/we07950.txt",skip_header=1,delimiter=''), axis=0)
    CFDE_array256_fivelev = np.append(CFDE_array256_fivelev, np.genfromtxt(path+"/SGF/5_ref_levels/we10600.txt",skip_header=1,delimiter=''), axis=0)
    

    CFDE_mech256_fivelev = CFDE_array256_fivelev[:,2]+CFDE_array256_fivelev[:,3]

    CFD_time256_fivelev = CFDE_array256_fivelev[:,1]
    return CFDE_mech256_fivelev, CFD_time256_fivelev

@app.cell
def _(np, path):
    ## Load the CFD data
    CFDE_array256_fivelev_AR2x = np.genfromtxt(path+"/SGF/5_ref_levels_AR2x/we00000.txt",skip_header=1,delimiter='')
    CFDE_array256_fivelev_AR2x = np.append(CFDE_array256_fivelev_AR2x, np.genfromtxt(path+"/SGF/5_ref_levels_AR2x/we20950.txt",skip_header=1,delimiter=''), axis=0)
    

    CFDE_mech256_fivelev_AR2x = CFDE_array256_fivelev_AR2x[:,2]+CFDE_array256_fivelev_AR2x[:,3]

    CFD_time256_fivelev_AR2x = CFDE_array256_fivelev_AR2x[:,1]
    return CFDE_mech256_fivelev_AR2x, CFD_time256_fivelev_AR2x

# @app.cell
# def _(np, path):
#     ## Load the CFD data
#     CFDE_array256_sixlev = np.genfromtxt(path+"/SGF/6_ref_levels/we00000.txt",skip_header=1,delimiter='')
#     CFDE_array256_sixlev = np.append(CFDE_array256_sixlev, np.genfromtxt(path+"/SGF/6_ref_levels/we22000.txt",skip_header=1,delimiter=''), axis=0)
#     CFDE_array256_sixlev = np.append(CFDE_array256_sixlev, np.genfromtxt(path+"/SGF/6_ref_levels/we24000.txt",skip_header=1,delimiter=''), axis=0)
#     
# 
#     CFDE_mech256_sixlev = CFDE_array256_sixlev[:,2]+CFDE_array256_sixlev[:,3]
# 
#     CFD_time256_sixlev = CFDE_array256_sixlev[:,1]
#     return CFDE_mech256_sixlev, CFD_time256_sixlev


@app.cell
def _(
    CFDE_mech256_threelev,
    CFDE_mech256_fourlev,
    CFDE_mech256_fourlev_AR2x,
    CFDE_mech256_fivelev,
    CFDE_mech256_fivelev_AR2x,
    # CFDE_mech256_sixlev,
    CFD_time256_threelev,
    CFD_time256_fourlev,
    CFD_time256_fourlev_AR2x,
    CFD_time256_fivelev,
    CFD_time256_fivelev_AR2x,
    # CFD_time256_sixlev,
    HOSEnergy,
    plt,
):
    _fig = plt.figure(1)
    from scipy.signal import savgol_filter
    ME_ref = HOSEnergy[0, 4]
    ME_ref_CFD = CFDE_mech256_fivelev[0]
    plt.plot(CFD_time256_threelev, savgol_filter(CFDE_mech256_threelev, 11, 3) / ME_ref_CFD, label='$\\Delta z$ = 4.7 m, AR = 2')
    plt.plot(CFD_time256_fourlev, savgol_filter(CFDE_mech256_fourlev, 11, 3) / ME_ref_CFD, label='$\\Delta z$ = 2.3 m, AR = 2')
    plt.plot(CFD_time256_fourlev_AR2x, savgol_filter(CFDE_mech256_fourlev_AR2x, 11, 3) / ME_ref_CFD, label='$\\Delta z$ = 1.2 m, AR = 4')
    plt.plot(CFD_time256_fivelev, savgol_filter(CFDE_mech256_fivelev, 11, 3) / ME_ref_CFD, label='$\\Delta z$ = 1.2 m, AR = 2')
    plt.plot(CFD_time256_fivelev_AR2x, savgol_filter(CFDE_mech256_fivelev_AR2x, 11, 3) / ME_ref_CFD, label='$\\Delta z$ = 0.6 m, AR = 4')
    # plt.plot(CFD_time256_sixlev, savgol_filter(CFDE_mech256_sixlev, 11, 3) / ME_ref_CFD, label='$\\Delta z$ = 0.6 m, AR = 2')
    plt.plot(HOSEnergy[:, 0], HOSEnergy[:, 4] / ME_ref, 'k', label='HOS reference solution')
    plt.xlim(0, 1600)
    plt.xlabel('time [s]', fontsize=16)
    plt.ylabel('$E/E_0$', fontsize=16)
    plt.legend()
    plt.title('Mechanical Energy Evolution')
    plt.ylim([0.1, 1.1])
    plt.savefig('plotting_outputs/Energy.png', format='png', bbox_inches='tight', dpi=300)
    return


if __name__ == "__main__":
    app.run()
