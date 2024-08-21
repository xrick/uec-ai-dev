#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Function prototypes
void autocorrelation(const double* signal, double* corr, int n);
int find_first_valley(const double* diff_corr, int n);
int find_peak(const double* corr, int start, int n);
double parabolic_interpolation(const double* corr, int peak_index);

void freq_from_autocorr(const double* signal, int n, double fs, double* result) {
    // Allocate memory for autocorrelation and its difference
    double* corr = (double*)malloc(2 * n * sizeof(double));
    double* diff_corr = (double*)malloc((2 * n - 1) * sizeof(double));
    if (!corr || !diff_corr) {
        fprintf(stderr, "Memory allocation failed\n");
        return;
    }

    // Calculate autocorrelation
    autocorrelation(signal, corr, n);

    // Compute the difference of the autocorrelation
    for (int i = 1; i < 2 * n - 1; i++) {
        diff_corr[i - 1] = corr[i] - corr[i - 1];
    }

    // Find the first valley in the autocorrelation
    int start = find_first_valley(diff_corr, 2 * n - 2);

    // Find the peak after the valley
    int peak_index = find_peak(corr, start, 2 * n - 1);

    // Interpolate the peak
    double i_interp = parabolic_interpolation(corr, peak_index);

    // Calculate the frequency
    *result = fs / i_interp;

    // Free allocated memory
    free(corr);
    free(diff_corr);
}

// Example function to calculate autocorrelation (using a naive approach)
void autocorrelation(const double* signal, double* corr, int n) {
    for (int i = 0; i < 2 * n - 1; i++) {
        corr[i] = 0.0;
        for (int j = 0; j < n; j++) {
            if (i - j >= 0 && i - j < n) {
                corr[i] += signal[j] * signal[i - j];
            }
        }
    }
}

// Find the first valley in the difference of autocorrelation
int find_first_valley(const double* diff_corr, int n) {
    for (int i = 0; i < n; i++) {
        if (diff_corr[i] > 0) {
            return i + 1; // Return the index of the first valley
        }
    }
    return -1; // Return -1 if no valley is found
}

// Find the peak in the autocorrelation after the start index
int find_peak(const double* corr, int start, int n) {
    int peak_index = start;
    for (int i = start + 1; i < n; i++) {
        if (corr[i] > corr[peak_index]) {
            peak_index = i;
        }
    }
    return peak_index;
}

// Parabolic interpolation to find a more accurate peak
double parabolic_interpolation(const double* corr, int peak_index) {
    double a = (peak_index > 0) ? corr[peak_index - 1] : corr[peak_index];
    double b = corr[peak_index];
    double c = (peak_index < 2 * n - 2) ? corr[peak_index + 1] : corr[peak_index];
    
    // Interpolation formula
    return peak_index - 0.5 * (c - a) / (2 * b - c - a);
}

int main() {
    // Example usage
    int n = 10; // Length of the signal
    double fs = 1000.0; // Sampling frequency
    double signal[10] = { /* signal values */ };
    double freq;

    freq_from_autocorr(signal, n, fs, &freq);
    printf("Estimated Frequency: %f Hz\n", freq);

    return 0;
}
