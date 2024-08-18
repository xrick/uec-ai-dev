#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "beep_clear_clip_01_22K.h"

// #ifndef _RFFT_H 
// #define _RFFT_H
// #define M_PI        3.14159265358979323846
// #define M_SQRT2     1.41421356237309504880
// #define RFFT_LEN  1024
// #define RFFT_ORDER 10
// #endif //end of _RFFT_H
#define PI          3.1415926
#define pi          3.14159265358979


void dft(float x[], float result[], uint32_t num_elems)
void freq_from_fft(float* signal, int dft_elements, float fs, float* result);
float parabolic(float* corr, int index);
int find_first_positive(float* d, int length);
int kaiser( float beta, int M, float *window );
int argmax(float* arr, int length);
void correlate(float* signal, int length, float* corr);
double mean(float* signal, int length);
double freq_from_autocorr(float* signal, int length, float fs);
void find(int* condition, int size, int** res, int* res_size);
void detectAudio(float* signal, int sig_len, int sr, float wanted_freq, float magthreshold, float freqthreshold, float (*freq_func)(float*, int));

double mean(float* signal, int length) 
{
    float sum = 0.0;
    for (int i = 0; i < length; i++) 
    {
        sum += signal[i];
    }
    return sum / length;
}

void correlate(float* signal, int length, double* corr) {
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

int kaiser( float beta, int M, float *window )
{
  // # Docstring adapted from NumPy's kaiser function
  // if _len_guards(M):
  //     return np.ones(M)
  // M, needs_trunc = _extend(M, sym)

  int result = true;

  // n = np.arange(0, M)
  float n[M];
  for ( int i = 0; i < M; i++ )
    n[i] = i;
  // alpha = (M - 1) / 2.0
  float alpha = (M - 1) / 2.0;
  // w = (special.i0(beta * np.sqrt(1 - ((n - alpha) / alpha) ** 2.0)) /
  //      special.i0(beta))
  for ( int i = 0; i < M; i++ ) {
    float p = pow( (n[i] - alpha) / alpha, 2 );
    window[i] = i0( beta * sqrt(1 - p) ) / i0( beta );
  }

  return result;
}



void dft(float x[], float result[], uint32_t num_elems) {
  // See: "modified C code" from https://batchloaf.wordpress.com/2013/12/07/simple-dft-in-c/ 
  // to simplify does not use pre-computed cos/sin(z)
  for(uint32_t k = 0; k < num_elems; k++) {
    float xre[num_elems]; // Real component
    float xim[num_elems]; // Imaginary component
    for(uint64_t n = 0; n < num_elems; n++) {
      float z = (2 * M_PI * k * n) / num_elems;
      xre[n] += x[n] * cos(z);
      xim[n] -= x[n] * sin(z);
    }
    result[k] = (xre[k] * xre[k]) * (xim[k] * xim[k]);
  }
}



void freq_from_fft(float* signal, int dft_elements, float fs, float* result) {
    /*
    Estimate frequency from peak of FFT
    Pros: Accurate, usually even more so than zero crossing counter
    (1000.000004 Hz for 1000 Hz, for instance).  Due to parabolic
    interpolation being a very good fit for windowed log FFT peaks?
    Accuracy also increases with signal length
    Cons: Doesn't find the right value if harmonics are stronger than
    fundamental, which is common.
    */
    
    //set dft_elements for testing
    dft_elements = 2048;
    /*************************/    
    float* windowed = (float*)malloc(N * sizeof(float));
    float* f = (float*)malloc((dft_elements + 1) * sizeof(float));
    float* log_abs_f = (float*)malloc((N/2 + 1) * sizeof(float));
    float* abs_f = (float*)malloc((N/2 + 1) * sizeof(float));
    float* kaiser_window = (float*)malloc(N * sizeof(float));
    float beta = 100.0;
    kaiser(beta, N, kaiser_window);
    // Apply Kaiser window
    for (int n = 0; n < N; n++) {
        windowed[n] = signal[n] * kaiser_window[n];
    }

    // apply dft to windowed signal with dft_elements
    dft(windowed, f, dft_elements);

    // Find the peak
    int i_peak = 0;
    float max_val = 0.0;
    for (int i = 0; i < (dft_elements/2 + 1); i++) {
        abs_f[i] = sqrt(windowed[i] * windowed[i]); // Assuming f is complex
        if (abs_f[i] > max_val) {
            max_val = abs_f[i];
            i_peak = i;
        }
    }
    // Log and interpolate
    for (int i = 0; i < (N/2 + 1); i++) {
        log_abs_f[i] = log(abs_f[i]);
    }
    float i_interp = parabolic(log_abs_f, i_peak);
    // Convert to equivalent frequency
    *result = fs * i_interp / N; // Hz

    // Free allocated memory
    free(windowed);
    // free(f);
    free(log_abs_f);
    free(abs_f);
}

// Note: You need to implement or include the kaiser, rfft, and parabolic functions.


// void detectAudio(double* signal, int sig_len, int sr, double wanted_freq, double magthreshold, double freqthreshold, double (*freq_func)(double*, int)) 
void detectAudio(float* signal, int sig_len, int sr, float wanted_freq, float magthreshold, float freqthreshold)
{
    // double testFreq = freq_func(signal, sr);
    float result=0.0;
    float testFreq = freq_from_fft2(signal, sig_len, sr, &result);
    int _dft_elements = 2048;
    printf("testFreq calculated: %f\n", testFreq);
    
    float diff_freq = fabs(wanted_freq - testFreq);
    printf("diff_freq calculated: %f\n", diff_freq);
    
    floate* window = (float*)malloc(sig_len * sizeof(float));
    if (window == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return;
    }

    for (int i = 0; i < sig_len; i++) {
        window[i] = signal[i] * (0.54 - 0.46 * cos(2 * M_PI * i / (sig_len - 1)));
    }

    // double complex* fftData = (double complex*)malloc(sig_len * sizeof(double complex));
    float* fftData = (float*)malloc(sig_len * sizeof(float));
    if (fftData == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        free(window);
        return;
    }
    dft(window, fftData, _dft_elements);
    // int targetIndex = (int)((double)sig_len * testFreq / sr);// get the bin we want to calculate maganitude
    int targetIndex = (int)((double)_dft_elements * testFreq / sr);
    float magnitude = cbs(fftData(targetIndex));

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
    // free(fftData_re);
    // free(fftData_im);
}


