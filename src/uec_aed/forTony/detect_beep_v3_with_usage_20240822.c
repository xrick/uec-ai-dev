// #pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "beep_clear_clip_01_22K.h" //from beep_clear_clip_01_22K.wav

#define PI          3.1415926535897932384626433832795


void compute_dft(float fftData[], float signal[], int N);
void compute_magnitude(double complex X[], double magnitude[], int N);
void dft_1(float x[], float result[], int num_elems);

int find_first_positive(float* d, int length);
void find(int* condition, int size, int** res, int* res_size);
int argmax(float* arr, int length);
void correlate(float* signal, int length, float* corr);
float mean(float* signal, int length);

///////
void autocorrelation(const double* signal, double* corr, int n);
int find_first_valley(const double* diff_corr, int n);
int find_peak(const double* corr, int start, int n);
double parabolic_interpolation(const double* corr, int peak_index, int n);
double freq_from_autocorr(const double* signal, int n, double fs);
///////
// void freq_from_fft(float* signal, int N, int dft_elements, float fs, float* result);
float freq_from_fft(float* signal, int N, int dft_elements, float fs, float* mag);

int detectAudio(float* signal, int sig_len, float sr, int dft_len, float wanted_freq, float magthreshold, float freqthreshold);


float mean(float* signal, int length) 
{
    float sum = 0.0;
    for (int i = 0; i < length; i++) 
    {
        sum += signal[i];
    }
    return sum / length;
}

void correlate(float* signal, int length, float* corr) {
    for (int lag = 0; lag < length; lag++) {
        corr[lag] = 0.0;
        for (int i = 0; i < length - lag; i++) {
            corr[lag] += signal[i] * signal[i + lag];
        }
    }
}

int find_first_positive(float* d, int length) 
{
    for (int i = 0; i < length; i++) 
    {
        if (d[i] > 0) 
        {
            return i;
        }
    }
    return -1; // Not found
};

int argmax(float* arr, int length)
{
    int max_index = 0;
    for (int i = 1; i < length; i++) 
    {
        if (arr[i] > arr[max_index]) 
        {
            max_index = i;
        }
    }
    return max_index;
};

float parabolic(float* corr, int index) 
{
    float a = corr[index - 1];
    float b = corr[index];
    float c = corr[index + 1];
    return index + (b - a) / (2 * (b - 2 * a + c));
};

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
/*************************** DFT Calculation Section************************/

void compute_dft(float fftData[], float signal[], int N) {
/*
N: length of dft length or can set the legth of signal
*/
    int k, n;
    for (k = 0; k < N; k++) {
        fftData[k] = 0.0 + 0.0*I; // Initialize to zero
        for (n = 0; n < N; n++) {
            double angle = -2 * M_PI * k * n / N;
            fftData[k] += signal[n] * (cos(angle) + I * sin(angle));
        }
    }
}

/*********************** End Of DFT Calculation Section*********************/
/*********************** Maganitude Calculation Section*********************/
void compute_magnitude(double complex X[], double magnitude[], int N) {
    int k;
    for (k = 0; k < N; k++) {
        magnitude[k] = cabs(X[k]); // cabs() returns the magnitude of the complex number
    }
}
/*********************** End Of Maganitude Section*********************/
/************************** Freqency Estimation  ***************************/
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
    printf("find peak index:%d\n", peak_index);
    return peak_index;
}

// Parabolic interpolation to find a more accurate peak
double parabolic_interpolation(const double* corr, int peak_index, int n) {
    double a = (peak_index > 0) ? corr[peak_index - 1] : corr[peak_index];
    double b = corr[peak_index];
    double c = (peak_index < 2 * n - 2) ? corr[peak_index + 1] : corr[peak_index];
    
    // Interpolation formula
    return peak_index - 0.5 * (c - a) / (2 * b - c - a);
}

double freq_from_autocorr(const double* signal, int n, double fs) {
/*
signal: input wav signal
n: length of signal
fs: sampling rate of input signal
*/
    // Allocate memory for autocorrelation and its difference
    double* corr = (double*)malloc(2 * n * sizeof(double));
    double* diff_corr = (double*)malloc((2 * n - 1) * sizeof(double));
    if (!corr || !diff_corr) {
        fprintf(stderr, "Memory allocation failed\n");
        return -1;;
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
    float i_interp = parabolic_interpolation(corr, peak_index, n);

    // Calculate the frequency
    double result = fs / i_interp;

    // Free allocated memory
    free(corr);
    free(diff_corr);
    return result;
}
//////////
// frequency using fft
float freq_from_fft(float* signal, int N, int dft_elements, float fs, float* mag) {
    /*
    Estimate frequency from peak of FFT
    Pros: Accurate, usually even more so than zero crossing counter
    (1000.000004 Hz for 1000 Hz, for instance).  Due to parabolic
    interpolation being a very good fit for windowed log FFT peaks?
    Accuracy also increases with signal length
    Cons: Doesn't find the right value if harmonics are stronger than
    fundamental, which is common.
    signal: input wav;
    N: length of signal;
    dft_elements: dft length we want to compute.
    fs: sampling rate, currently use 20K
    */
    // printf("N:%d,\ndft_elements:%d,\nsampling rate:%f\n",N,dft_elements,fs);
    //set dft_elements for testing
    // dft_elements = 2048;
    /*************************/ 
    int computed_dft_len = (int)(dft_elements/2 + 1);   
    float* windowed = (float*)malloc(N * sizeof(float));
    float* f = (float*)malloc((dft_elements + 1) * sizeof(float));
    // float* log_abs_f = (float*)malloc((N/2 + 1) * sizeof(float));
    // float* abs_f = (float*)malloc((N/2 + 1) * sizeof(float));
    float* log_abs_f = (float*)malloc(computed_dft_len * sizeof(float));
    float* abs_f = (float*)malloc(computed_dft_len * sizeof(float));
    //we don't use kaiser window here.
    //please call cmsis dft function.
    compute_dft(f, signal, dft_elements);
    // Find the peak
    int i_peak = 0;
    float max_val = 0.0;
    for (int i = 0; i < computed_dft_len; i++) {
        abs_f[i] = sqrt(f[i] * f[i]); // Assuming f is complex
        // printf("%d.abs_f[i] is %f\n",i,abs_f[i]);
        if (abs_f[i] > max_val) {
            max_val = abs_f[i];
            i_peak = i;
        }
    }
    *mag = max_val;
    printf("i_peak is %d\n",i_peak);
    // Log and interpolate
    for (int i = 0; i < computed_dft_len; i++) {
        log_abs_f[i] = log(abs_f[i]);
    }
    float i_interp = parabolic(log_abs_f, i_peak);
    // Convert to equivalent frequency
    float result = fs * i_interp / dft_elements; // Hz
    // Free allocated memory
    free(windowed);
    free(f);
    free(log_abs_f);
    free(abs_f);
    return result;
}
/************************** End Of Freqency Estimation  ***************************/
// Note: You need to implement or include the kaiser, rfft, and parabolic functions.

int detectAudio(float* signal, int sig_len, float sr, int dft_len, float wanted_freq, float magthreshold, float freqthreshold)
{
/*
Parameters:
    signal: input wav signal, need to be float.
    sig_len: length of signal.
    sr: sampling, need to be float for float computation.
    dft_len: the points we want to compute in dft algorithm.
    wanted_freq: golden frequency of standard sound we want to detect,
                 for example, beep:3078.918213(estimated by dft)
    magthreshold: the threshold of magnitude.
    freqthreshold:the threshold is to test whether the absolute value of
                  golden sample frequence minus input signal frequence
                  smaller than it.
return:
    1: The input sound matches the wanted_freq and maganitude is also equal or greater than the magthreshold
    98: Wanted Frequency detected, but no significant magnitude.
    99: wanted frequency is not found
*/
    float _mag = 0.0;
    float testFreq = freq_from_fft(signal, sig_len, dft_len, sr, &_mag);

    printf("testFreq calculated: %f\n", testFreq);
    
    float diff_freq = fabs(wanted_freq - testFreq);
    printf("diff_freq calculated: %f\n", diff_freq);
    
    float* window = (float*)malloc(sig_len * sizeof(float));
    if (window == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return -1;
    }

    for (int i = 0; i < sig_len; i++) {
        window[i] = signal[i] * (0.54 - 0.46 * cos(2 * M_PI * i / (sig_len - 1)));
    }

    // float complex* fftData = (float complex*)malloc(sig_len * sizeof(float complex));
    float* fftData = (float*)malloc(sig_len * sizeof(float));
    if (fftData == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        free(window);
        return -1;
    }
    compute_dft(fftData, window, dft_len);
    // int targetIndex = (int)((float)sig_len * testFreq / sr);// get the bin we want to calculate maganitude
    int targetIndex = (int)((float)dft_len * testFreq / sr);
    float magnitude = sqrt(fftData[targetIndex] * fftData[targetIndex]);//cbs(fftData[targetIndex]);
    free(window);
    free(fftData);
    if (diff_freq < freqthreshold) {
        if (magnitude > magthreshold) {
            printf("Wanted Frequency detected: %f Hz and magnitude: %f\n", testFreq, magnitude);
            return 1;
        } else {
            printf("Wanted Frequency detected: %f Hz but no significant magnitude: %f\n", testFreq, magnitude);
            return 98;
        }
    } else {
        printf("wanted frequency: %f is not found, found frequency: %f and magnitude is %f\n", wanted_freq, testFreq, magnitude);
        return 99;
    }
}


int main(void)
{
    /*
        signature
        freq_from_fft(float* signal, int N, int dft_elements, float fs, float* result);
    */
    int wav_len = sizeof(wav_array)/sizeof(wav_array[0]);
    int _dft_len = 2048;
    float _wanted_f = 3078.0;
    float _sr = 20000.0;
    float _magthreshold = 2000;//162590.0;
    float _freqthreshold = 3078.0;
    // detectAudio(float* signal, int sig_len, int sr, int dft_len, float wanted_freq, float magthreshold, float freqthreshold)
    float* wav_sig = (float*)malloc(wav_len * sizeof(float));
    for (int i = 0; i < wav_len; i++) {
        // wav_sig[i] = wav_array[i] * (0.54 - 0.46 * cos(2 * M_PI * i / (wav_len - 1)));
        wav_sig[i] = (float)wav_array[i];
    }
    int res = detectAudio(wav_sig, wav_len, _sr, _dft_len, _wanted_f,_magthreshold,_freqthreshold);
    printf("test result is %d",res);
    return 0;
}




