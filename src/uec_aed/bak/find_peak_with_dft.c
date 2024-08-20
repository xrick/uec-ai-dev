#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#define WAV_HEADER_SIZE 44
#define PI 3.14159265358979323846

// Function to read a 16-bit sample from the WAV file
int16_t read_sample(FILE *file) {
    int16_t sample;
    fread(&sample, sizeof(int16_t), 1, file);
    return sample;
}

// Function to perform DFT
void dft(double complex *input, double complex *output, int N) {
    for (int k = 0; k < N; k++) {
        output[k] = 0;
        for (int n = 0; n < N; n++) {
            output[k] += input[n] * cexp(-2.0 * PI * I * k * n / N);
        }
    }
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <wav_file>\n", argv[0]);
        return 1;
    }

    FILE *file = fopen(argv[1], "rb");
    if (!file) {
        perror("Error opening file");
        return 1;
    }

    // Read WAV header to get sample rate
    fseek(file, 24, SEEK_SET);
    int32_t sample_rate;
    fread(&sample_rate, sizeof(int32_t), 1, file);

    // Skip to data
    fseek(file, WAV_HEADER_SIZE, SEEK_SET);

    // Read samples
    int N = 8192; // Number of samples to process (power of 2 for efficiency)
    double complex input[N];
    for (int i = 0; i < N; i++) {
        input[i] = read_sample(file);
    }

    fclose(file);

    // Perform DFT
    double complex output[N];
    dft(input, output, N);

    // Find peak frequency
    double max_magnitude = 0;
    int peak_index = 0;
    for (int i = 0; i < N / 2; i++) {
        double magnitude = cabs(output[i]);
        if (magnitude > max_magnitude) {
            max_magnitude = magnitude;
            peak_index = i;
        }
    }

    // Calculate peak frequency
    double peak_frequency = (double)peak_index * sample_rate / N;

    printf("Peak frequency: %.2f Hz\n", peak_frequency);

    return 0;
}