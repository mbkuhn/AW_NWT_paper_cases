using FFTW
using DSP
using ImageFiltering
using Plots
using Statistics

function FFTAnalysis(time,x,bSmoothing,bPlotting,smoothtype,parameter)

    dt = (time[end]-time[1])/(length(time)-1)
    if sum(abs.(time[2:end]-time[1:end-1].-dt)>(1e-7.*time[2:end]))>0
        println("FFTAnalysis: Time Step is not constant!")
    end

    N = length(time)                     # Length of signal
    Fs = 1/dt                            # Sampling frequency
    df = Fs/N                            # Frequency resolution
    f = collect(0:ceil(N/2)-1).*df        # (Non-negative) Frequency array
    f = reshape(f,length(f),1)

    # Make sure time and signal are both column vectors
    time = reshape(time,N,1)
    time = vec(time)
    x = reshape(x,N,1)
    x = vec(x)

    #Compute the Fourier tranform of signal
    x_fft = fft(x)

    #Compute two-sided spectrum 
    P2 = x_fft./N

    #Compute one-sided specturm
    ceilN2 = Int.(ceil(N/2))
    P1 = P2[1:ceilN2]
    P1[2:end] = P1[2:end].*2

    #Rearrange 2-sided spectrum with DC (zero-frequency) component at the center
    fdouble = [-f[end:-1:2]; f]  #Double-sided frequency array with negative frequencies
    P2 = fftshift(P2);            #Rearrange double-sided spectrum
    if mod(length(P2),2)==0
        P2 = P2[2:end] #Drop Nyquist component for full symmetry
    end
    
    #Correct phase of amplitude spectrum if time does not start from 0.
    P1 = P1 .* exp.(-im*2*pi*f*time[1])
    P2 = P2 .* exp.(-im*2*pi*fdouble*time[1])
    
    #Calculate one-sided spectrum (PSD)
    spectrum = (P1.*conj(P1))./df
    spectrum[2:end] = spectrum[2:end].*0.5
    
    #Smooth the spectrum
    if bSmoothing
        if smoothtype == 1
            # Ensemble Averaging
            ni = parameter    # Number of frequencies/length of sub-timeseries
            nens = floor(N/ni) # Number of ensembles
            spectrumSum = zeros(Int.(0.5*ni),1) # this needs to be only half of ni, because the FFT will cut the second half of the frequencies
            for sectionNo=1:nens
                section = collect(ni*(sectionNo-1)+1:ni*sectionNo)
                section = Int.(section)
                f2, P1, spectrum_smoothed_iter, fdouble, P2  = FFT_smoothing(time[section],x[section],false,false,smoothtype,parameter)
                spectrumSum =  spectrumSum + spectrum_smoothed_iter
                f2 = f2 # To repress scope warning
            end
            spectrum_smoothed = spectrumSum./nens
        end
        if smoothtype == 2
            # MARIN's method - Constant filter
            nskip = parameter;
            window = ones(nskip,1)./nskip
            f2 = f
            spectrum_smoothed = imfilter(reshape(spectrum,length(spectrum),1), reflect(centered(window)), Fill(0))
        end
        if smoothtype == 3
            # Gaussian filter
            f2 = f
            spectrum_smoothed = imfilter(reshape(spectrum,length(spectrum),1), Kernel.gaussian(parameter))
        end
    end
    # Plot spectrums
    if bPlotting==true
        #1-sided spectrum on semilog plot

        # full frequency range
        spectrum_comparison_fullrange = plot(f,abs.(spectrum),yscale=:log10,label="Unsmoothed Spectrum")
        #spectrum_comparison_fullrange = plot(f,abs.(spectrum),label="Unsmoothed Spectrum") # linear scale
        xlabel!("Frequency [Hz]")
        ylabel!("Amplitude [m]")
        if bSmoothing==true  #plot smoothed spectrum on top
            plot!(f2,abs.(spectrum_smoothed),yscale=:log10,label="Smoothed Spectrum")
            #plot!(f2,abs.(spectrum_smoothed),label="Smoothed Spectrum")
        end
        savefig(spectrum_comparison_fullrange,"Spectrum_comparison_full_range.pdf")

        # cut-off range
        spectrum_comparison = plot(f,abs.(spectrum),yscale=:log10,label="Unsmoothed Spectrum",xlims=(0,10))
        #spectrum_comparison = plot(f,abs.(spectrum),label="Unsmoothed Spectrum",xlims=(0,10))
        xlabel!("Frequency [Hz]")
        ylabel!("Amplitude [m]")
        if bSmoothing==true  #plot smoothed spectrum on top
            plot!(f2,abs.(spectrum_smoothed),yscale=:log10,label="Smoothed Spectrum")
            #plot!(f2,abs.(spectrum_smoothed),label="Smoothed Spectrum")
        end
        savefig(spectrum_comparison,"Spectrum_comparison.pdf")
    end
    
    # Output smoothed spectrum if requested
    if bSmoothing==true
        f = vec(f2)
        spectrum = vec(spectrum_smoothed)
        P1 = P1
    end

    return f, P1, spectrum, fdouble, P2

end