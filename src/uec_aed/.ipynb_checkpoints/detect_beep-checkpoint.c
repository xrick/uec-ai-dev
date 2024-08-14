#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

double parabolic(double* corr, int index);
int find_first_positive(double* d, int length);
int argmax(double* arr, int length);
void correlate(double* signal, int length, double* corr);
double mean(double* signal, int length);
double freq_from_autocorr(double* signal, int length, double fs);
void find(int* condition, int size, int** res, int* res_size);
void detectAudio(double* signal, int sig_len, int sr, double wanted_freq, double magthreshold, double freqthreshold, double (*freq_func)(double*, int));

/************************FFT implementation***********************/
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

void fft(complex double *x, int N);

// int main() {
//     // Example signal initialization
//     int N = 8; // Size of the signal
//     complex double *signal = malloc(N * sizeof(complex double));

//     // Initialize the signal with some values (example)
//     for (int i = 0; i < N; i++) {
//         signal[i] = i + 0.0 * I; // Replace with actual signal values
//     }

//     fft(signal, N);

//     // Free allocated memory
//     free(signal);
//     return 0;
// }

// FFT implementation (Cooley-Tukey algorithm)
void fft(complex double *x, int N) {
/*
    x: signal
*/
    if (N <= 1) return;

    // Divide
    complex double *even = malloc(N / 2 * sizeof(complex double));
    complex double *odd = malloc(N / 2 * sizeof(complex double));
    for (int i = 0; i < N / 2; i++) {
        even[i] = x[i * 2];
        odd[i] = x[i * 2 + 1];
    }

    // Conquer
    fft(even, N / 2);
    fft(odd, N / 2);

    // Combine
    for (int k = 0; k < N / 2; k++) {
        complex double t = cexp(-2.0 * I * M_PI * k / N) * odd[k];
        x[k] = even[k] + t;
        x[k + N / 2] = even[k] - t;
    }

    // Free allocated memory
    free(even);
    free(odd);
}



/***********C codes for freq_from_autocorr***********************/

double mean(double* signal, int length) {
    double sum = 0.0;
    for (int i = 0; i < length; i++) {
        sum += signal[i];
    }
    return sum / length;
}

void correlate(double* signal, int length, double* corr) {
    for (int lag = 0; lag < length; lag++) {
        corr[lag] = 0.0;
        for (int i = 0; i < length - lag; i++) {
            corr[lag] += signal[i] * signal[i + lag];
        }
    }
}

int find_first_positive(double* d, int length) {
    for (int i = 0; i < length; i++) {
        if (d[i] > 0) {
            return i;
        }
    }
    return -1; // Not found
}

int argmax(double* arr, int length) {
    int max_index = 0;
    for (int i = 1; i < length; i++) {
        if (arr[i] > arr[max_index]) {
            max_index = i;
        }
    }
    return max_index;
}

double parabolic(double* corr, int index) {
    double a = corr[index - 1];
    double b = corr[index];
    double c = corr[index + 1];
    return index + (b - a) / (2 * (b - 2 * a + c));
}

double freq_from_autocorr(double* signal, int length, double fs) {
    double* corr = (double*)malloc(2 * length * sizeof(double));
    double* d = (double*)malloc((length - 1) * sizeof(double));
    
    double dc_offset = mean(signal, length);
    for (int i = 0; i < length; i++) {
        signal[i] -= dc_offset; // Remove DC offset
    }

    correlate(signal, length, corr);
    
    // Throw away the negative lags
    for (int i = 0; i < length; i++) {
        corr[i] = corr[length + i];
    }

    for (int i = 0; i < length - 1; i++) {
        d[i] = corr[i + 1] - corr[i];
    }

    int start = find_first_positive(d, length - 1);
    int i_peak = argmax(corr + start, length - start) + start;
    double i_interp = parabolic(corr, i_peak);
    
    free(corr);
    free(d);
    
    return fs / i_interp;
}
/***********End C codes for freq_from_autocorr***********************/

/******************utility functions********************************/

void find(int* condition, int size, int** res, int* res_size) {
    *res_size = 0;
    *res = (int*)malloc(size * sizeof(int));
    
    for (int i = 0; i < size; i++) {
        if (condition[i]) {
            (*res)[(*res_size)++] = i;
        }
    }
    *res = (int*)realloc(*res, (*res_size) * sizeof(int));
}



//main function of detecting beep and alarm
/*original python codes
def detectAudio(signal=None, sr=None, wanted_freq=None, magthreshold=None, freqthreshold=None, freq_func=None):
    testFreq = freq_func(signal, sr);#freq_estimate.freq_from_fft(signal, sr);
    print(f"testFreq calculated:{testFreq}")
    diff_freq = np.abs(int(np.abs(wanted_freq)-np.abs(testFreq)));
    print(f"diff_freq calculated:{diff_freq}")
    
    sig_len = len(signal);
    window = [0] * sig_len;
    # signal = np.float64(signal);
    i = 0;
    for x in signal:
        window[i] = np.float64(x) * (0.54 - 0.46*math.cos(2*math.pi*np.float64(i)/np.float64(sig_len-1)));
        i += 1;
    window = np.asarray(window, dtype=np.float64);
    fftData = np.fft.fft(window);
    fftData_re = np.real(fftData);
    fftData_im = np.imag(fftData);
    targetIndex =  int(np.float64(len(fftData_re)) * np.float64(testFreq) / np.float64(sr));
    # magnitude = np.sqrt(fftData_re^2 + fftData_im^2);
    magnitude = np.absolute(np.complex128(fftData_re[targetIndex]));

    if diff_freq < freqthreshold:
        if magnitude > magthreshold:
            print("Wanted Frequencey detected:{}Hz and magnitude: {}\n".format(testFreq, magnitude));
        else:
            print("Wanted Frequencey detected:{}Hz but no significant magnitude: {}\n".format(testFreq, magnitude));
            
    else:
        print(f"wanted frequency:{wanted_freq} is not found, found frequency:{testFreq} and magnitude is {magnitude}");
*/
void detectAudio(double* signal, int sig_len, int sr, double wanted_freq, double magthreshold, double freqthreshold, double (*freq_func)(double*, int)) {
    double testFreq = freq_func(signal, sr);
    printf("testFreq calculated: %f\n", testFreq);
    
    double diff_freq = fabs(wanted_freq - testFreq);
    printf("diff_freq calculated: %f\n", diff_freq);
    
    double* window = (double*)malloc(sig_len * sizeof(double));
    if (window == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return;
    }

    for (int i = 0; i < sig_len; i++) {
        window[i] = signal[i] * (0.54 - 0.46 * cos(2 * M_PI * i / (sig_len - 1)));
    }

    double complex* fftData = (double complex*)malloc(sig_len * sizeof(double complex));
    if (fftData == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        free(window);
        return;
    }

    // Perform FFT (this requires a FFT library, e.g., FFTW)
    // fftw_execute_dft_r2c(plan, window, fftData);

    double* fftData_re = (double*)malloc(sig_len * sizeof(double));
    double* fftData_im = (double*)malloc(sig_len * sizeof(double));
    if (fftData_re == NULL || fftData_im == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        free(window);
        free(fftData);
        return;
    }
    for (int i = 0; i < sig_len; i++) {
        fftData_re[i] = creal(fftData[i]);
        fftData_im[i] = cimag(fftData[i]);
    }
    int targetIndex = (int)((double)sig_len * testFreq / sr);
    double magnitude = cabs(fftData_re[targetIndex] + I * fftData_im[targetIndex]);

    if (diff_freq < freqthreshold) {
        if (magnitude > magthreshold) {
            printf("Wanted Frequency detected: %f Hz and magnitude: %f\n", testFreq, magnitude);
        } else {
            printf("Wanted Frequency detected: %f Hz but no significant magnitude: %f\n", testFreq, magnitude);
        }
    } else {
        printf("wanted frequency: %f is not found, found frequency: %f and magnitude is %f\n", wanted_freq, testFreq, magnitude);
    }
    free(window);
    free(fftData);
    free(fftData_re);
    free(fftData_im);
}

