import numpy as np
from numpy import argmax, asarray, copy, diff, log, mean
from numpy.fft import rfft
from scipy.signal import correlate, decimate
from scipy.signal.windows import kaiser
import math;

class freq_estimate:

    @classmethod
    def freq_from_crossings(signal, fs, interp='linear'):
        """
        Estimate frequency by counting zero crossings
    
        Works well for long low-noise sines, square, triangle, etc.
    
        Pros: Fast, accurate (increasing with signal length).
    
        Cons: Doesn't work if there are multiple zero crossings per cycle,
        low-frequency baseline shift, noise, inharmonicity, etc.
        """
        signal = asarray(signal) + 0.0
    
        # Find all indices right before a rising-edge zero crossing
        indices = find((signal[1:] >= 0) & (signal[:-1] < 0))
    
        if interp == 'linear':
            # More accurate, using linear interpolation to find intersample
            # zero-crossings (Measures 1000.000129 Hz for 1000 Hz, for instance)
            crossings = [i - signal[i] / (signal[i+1] - signal[i])
                         for i in indices]
        elif interp == 'none' or interp is None:
            # Naive (Measures 1000.185 Hz for 1000 Hz, for instance)
            crossings = indices
        else:
            raise ValueError('Interpolation method not understood')
    
            # TODO: Some other interpolation based on neighboring points might be
            # better.  Spline, cubic, whatever  Can pass it a function?
    
        return fs / mean(diff(crossings))
        
    @classmethod
    def freq_from_fft(signal, fs):
        """
        Estimate frequency from peak of FFT
    
        Pros: Accurate, usually even more so than zero crossing counter
        (1000.000004 Hz for 1000 Hz, for instance).  Due to parabolic
        interpolation being a very good fit for windowed log FFT peaks?
        https://ccrma.stanford.edu/~jos/sasp/Quadratic_Interpolation_Spectral_Peaks.html
        Accuracy also increases with signal length
    
        Cons: Doesn't find the right value if harmonics are stronger than
        fundamental, which is common.
        """
        signal = asarray(signal)
    
        N = len(signal)
    
        # Compute Fourier transform of windowed signal
        windowed = signal * kaiser(N, 100)
        f = rfft(windowed)
    
        # Find the peak and interpolate to get a more accurate peak
        i_peak = argmax(abs(f))  # Just use this value for less-accurate result
        i_interp = parabolic(log(abs(f)), i_peak)[0]
    
        # Convert to equivalent frequency
        return fs * i_interp / N  # Hz
        
    @classmethod
    def freq_from_autocorr(signal, fs):
        """
        Estimate frequency using autocorrelation
    
        Pros: Best method for finding the true fundamental of any repeating wave,
        even with strong harmonics or completely missing fundamental
    
        Cons: Not as accurate, doesn't find fundamental for inharmonic things like
        musical instruments, this implementation has trouble with finding the true
        peak
        """
        signal = asarray(signal) + 0.0
    
        # Calculate autocorrelation, and throw away the negative lags
        signal -= mean(signal)  # Remove DC offset
        corr = correlate(signal, signal, mode='full')
        corr = corr[len(corr)//2:]
    
        # Find the first valley in the autocorrelation
        d = diff(corr)
        start = find(d > 0)[0]
    
        # Find the next peak after the low point (other than 0 lag).  This bit is
        # not reliable for long signals, due to the desired peak occurring between
        # samples, and other peaks appearing higher.
        i_peak = argmax(corr[start:]) + start
        i_interp = parabolic(corr, i_peak)[0]
    
        return fs / i_interp
        
    @classmethod
    def freq_from_hps(signal, fs):
        """
        Estimate frequency using harmonic product spectrum
    
        Low frequency noise piles up and overwhelms the desired peaks
    
        Doesn't work well if signal doesn't have harmonics
        """
        signal = asarray(signal) + 0.0
    
        N = len(signal)
        signal -= mean(signal)  # Remove DC offset
    
        # Compute Fourier transform of windowed signal
        windowed = signal * kaiser(N, 100)
    
        # Get spectrum
        X = log(abs(rfft(windowed)))
    
        # Remove mean of spectrum (so sum is not increasingly offset
        # only in overlap region)
        X -= mean(X)
    
        # Downsample sum logs of spectra instead of multiplying
        hps = copy(X)
        for h in range(2, 9):  # TODO: choose a smarter upper limit
            dec = decimate(X, h, zero_phase=True)
            hps[:len(dec)] += dec
    
        # Find the peak and interpolate to get a more accurate peak
        i_peak = argmax(hps[:len(dec)])
        i_interp = parabolic(hps, i_peak)[0]
    
        # Convert to equivalent frequency
        return fs * i_interp / N  # Hz


def processAudio(signal=None, targetFreq=None, sr=None, threshold=None):
    window = [];
    sig_len = len(signal);
    # signal = np.float64(signal);
    i = 0;
    for x in signal:
        window[i] = np.float64(x) * (0.54 - 0.46*math.cos(2*math.pi*np.float64(i)/np.float64(sig_len-1)));
        i += 1;
    window = np.asarray(window, dtype=np.float64);
    fftData = np.fft.fft(window);
    fftData_re = np.real(fftData);
    fftData_im = np.imag(fftData);
    targetIndex =  int(np.float64(len(fftData_re)) * np.float64(targetFreq) / np.float64(sr));
    # magnitude = np.sqrt(fftData_re^2 + fftData_im^2);
    magnitude = np.absolute(np.complex128(fftData_re[targetIndex]));
    
    if magnitude > threshold:
        print("Significant data detected {}Hz, magnitude: {%.2f}\n".format(targetFreq, magnitude));
    else:
        print("\rNo significant magnitude at {}Hz\n", targetFreq)


if __name__ == "__main__":
	run();