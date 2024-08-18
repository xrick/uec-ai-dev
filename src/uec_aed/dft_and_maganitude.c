#include <stdio.h>
#include <math.h>
#include <complex.h>

#define N 8 // Number of input samples

// Function to compute the DFT
void compute_dft(double complex X[], double x[], int N) {
    int k, n;
    for (k = 0; k < N; k++) {
        X[k] = 0.0 + 0.0*I; // Initialize to zero
        for (n = 0; n < N; n++) {
            double angle = -2 * M_PI * k * n / N;
            X[k] += x[n] * (cos(angle) + I * sin(angle));
        }
    }
}

// Function to compute the magnitude of the DFT output
void compute_magnitude(double complex X[], double magnitude[], int N) {
    int k;
    for (k = 0; k < N; k++) {
        magnitude[k] = cabs(X[k]); // cabs() returns the magnitude of the complex number
    }
}

int main() {
    double x[N] = {1, 2, 3, 4, 5, 6, 7, 8}; // Input signal
    double complex X[N]; // DFT output
    double magnitude[N]; // Magnitude of DFT output
    
    // Compute DFT
    compute_dft(X, x, N);
    
    // Compute magnitude of DFT output
    compute_magnitude(X, magnitude, N);
    
    // Print results
    printf("DFT Magnitudes:\n");
    for (int k = 0; k < N; k++) {
        printf("Magnitude of X[%d]: %.4f\n", k, magnitude[k]);
    }

    return 0;
}
