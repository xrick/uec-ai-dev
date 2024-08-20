#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#define WAV_HEADER_SIZE 44

// Function to read a 16-bit sample from the WAV file
int16_t read_sample(FILE *file) {
    int16_t sample;
    fread(&sample, sizeof(int16_t), 1, file);
    return sample;
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

    // Skip WAV header
    fseek(file, WAV_HEADER_SIZE, SEEK_SET);

    int16_t max_sample = 0;
    int16_t min_sample = 0;
    int64_t sample_count = 0;
    int16_t sample;

    while (fread(&sample, sizeof(int16_t), 1, file) == 1) {
        if (sample > max_sample) max_sample = sample;
        if (sample < min_sample) min_sample = sample;
        sample_count++;
    }

    fclose(file);

    double peak_amplitude = fmax(fabs((double)max_sample), fabs((double)min_sample)) / 32768.0;

    printf("Number of samples: %lld\n", sample_count);
    printf("Peak positive sample value: %d\n", max_sample);
    printf("Peak negative sample value: %d\n", min_sample);
    printf("Peak amplitude (normalized): %f\n", peak_amplitude);

    return 0;
}