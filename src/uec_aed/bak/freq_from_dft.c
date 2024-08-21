#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

// Function prototypes
void kaiser_window(double* window, int n, double beta);
void freq_from_fft(const double* signal, int n, double fs, double* result);
double parabolic_interpolation(const double* spectrum, int peak_index);

void freq_from_fft(const double* signal, int n, double fs, double* result) {
    // Allocate memory for the FFTW complex array
    fftw_complex *in, *out;
    fftw_plan p;
    
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
    
    if (!in || !out) {
        fprintf(stderr, "Memory allocation failed\n");
        return;
    }

    // Apply Kaiser window to the signal
    double *windowed_signal = (double*) malloc(n * sizeof(double));
    if (!windowed_signal) {
        fprintf(stderr, "Memory allocation failed\n");
        fftw_free(in);
        fftw_free(out);
        return;
    }
    
    kaiser_window(windowed_signal, n, 100.0);
    for (int i = 0; i < n; i++) {
        in[i][0] = signal[i] * windowed_signal[i]; // Real part
        in[i][1] = 0.0; // Imaginary part
    }

    // Create a plan for FFT
    p = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    // Find the peak in the magnitude spectrum
    double *magnitude = (double*) malloc((n / 2 + 1) * sizeof(double));
    if (!magnitude) {
        fprintf(stderr, "Memory allocation failed\n");
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
        free(windowed_signal);
        return;
    }
    
    int peak_index = 0;
    double max_magnitude = 0.0;
    
    for (int i = 0; i < n / 2 + 1; i++) {
        magnitude[i] = sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]);
        if (magnitude[i] > max_magnitude) {
            max_magnitude = magnitude[i];
            peak_index = i;
        }
    }
    
    // Interpolate the peak position
    double i_interp = parabolic_interpolation(magnitude, peak_index);

    // Calculate the frequency
    *result = fs * i_interp / n;

    // Free allocated memory
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    free(windowed_signal);
    free(magnitude);
}

// Kaiser window function
void kaiser_window(double* window, int n, double beta) {
    double I0 = 1.0; // Modified Bessel function of the first kind
    double x;
    double alpha = (n - 1) / 2.0;
    
    for (int i = 0; i < n; i++) {
        x = (i - alpha) / alpha;
        window[i] = 0.0;
        if (fabs(x) <= 1.0) {
            window[i] = I0 * (beta * sqrt(1 - x * x));
        }
    }
}

// Parabolic interpolation
double parabolic_interpolation(const double* spectrum, int peak_index) {
    double a = (peak_index > 0) ? spectrum[peak_index - 1] : spectrum[peak_index];
    double b = spectrum[peak_index];
    double c = (peak_index < n / 2) ? spectrum[peak_index + 1] : spectrum[peak_index];
    
    // Interpolation formula
    return peak_index - 0.5 * (c - a) / (2 * b - c - a);
}

int main() {
    // Example usage
    int n = 1024; // Length of the signal
    double fs = 1000.0; // Sampling frequency
    double signal[1024] = { /* signal values */ };
    double freq;

    freq_from_fft(signal, n, fs, &freq);
    printf("Estimated Frequency: %f Hz\n", freq);

    return 0;
}
