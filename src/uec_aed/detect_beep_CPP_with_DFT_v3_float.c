// #pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "complex.h"
#include "beep_clear_clip_01_22K.h"

// #ifndef _RFFT_H 
// #define _RFFT_H
#define M_PI        3.14159265358979323846
#define PI          3.1415926535897932384626433832795
// #define M_SQRT2     1.41421356237309504880
// #define RFFT_LEN  1024
// #define RFFT_ORDER 10
// #endif //end of _RFFT_H
// #define PI          3.1415926
// #define pi          3.14159265358979


float chbevl(float x, float array[], int n);
float i0(float x);
void dft_1(float x[], float result[], int num_elems);
complex *dft(complex *h, ulint N);
float parabolic(float* corr, int index);
int find_first_positive(float* d, int length);
int kaiser( float beta, int M, float *window );
int argmax(float* arr, int length);
void correlate(float* signal, int length, float* corr);
float mean(float* signal, int length);
float freq_from_autocorr(float* signal, int length, float fs);
// void freq_from_fft(float* signal, int N, int dft_elements, float fs, float* result);
float freq_from_fft(float* signal, int N, int dft_elements, float fs);
void find(int* condition, int size, int** res, int* res_size);
void detectAudio(float* signal, int sig_len, int sr, int dft_len, float wanted_freq, float magthreshold, float freqthreshold);

static float A[] = {
  -4.41534164647933937950E-18,
  3.33079451882223809783E-17,
  -2.43127984654795469359E-16,
  1.71539128555513303061E-15,
  -1.16853328779934516808E-14,
  7.67618549860493561688E-14,
  -4.85644678311192946090E-13,
  2.95505266312963983461E-12,
  -1.72682629144155570723E-11,
  9.67580903537323691224E-11,
  -5.18979560163526290666E-10,
  2.65982372468238665035E-9,
  -1.30002500998624804212E-8,
  6.04699502254191894932E-8,
  -2.67079385394061173391E-7,
  1.11738753912010371815E-6,
  -4.41673835845875056359E-6,
  1.64484480707288970893E-5,
  -5.75419501008210370398E-5,
  1.88502885095841655729E-4,
  -5.76375574538582365885E-4,
  1.63947561694133579842E-3,
  -4.32430999505057594430E-3,
  1.05464603945949983183E-2,
  -2.37374148058994688156E-2,
  4.93052842396707084878E-2,
  -9.49010970480476444210E-2,
  1.71620901522208775349E-1,
  -3.04682672343198398683E-1,
  6.76795274409476084995E-1
};

/* Chebyshev coefficients for exp(-x) sqrt(x) I0(x)
   in the inverted interval [8,infinity].

   lim(x->inf){ exp(-x) sqrt(x) I0(x) } = 1/sqrt(2pi).
*/
static float B[] = {
  -7.23318048787475395456E-18,
  -4.83050448594418207126E-18,
  4.46562142029675999901E-17,
  3.46122286769746109310E-17,
  -2.82762398051658348494E-16,
  -3.42548561967721913462E-16,
  1.77256013305652638360E-15,
  3.81168066935262242075E-15,
  -9.55484669882830764870E-15,
  -4.15056934728722208663E-14,
  1.54008621752140982691E-14,
  3.85277838274214270114E-13,
  7.18012445138366623367E-13,
  -1.79417853150680611778E-12,
  -1.32158118404477131188E-11,
  -3.14991652796324136454E-11,
  1.18891471078464383424E-11,
  4.94060238822496958910E-10,
  3.39623202570838634515E-9,
  2.26666899049817806459E-8,
  2.04891858946906374183E-7,
  2.89137052083475648297E-6,
  6.88975834691682398426E-5,
  3.36911647825569408990E-3,
  8.04490411014108831608E-1
};

float chbevl(float x, float array[], int n)
{
  float b0, b1, b2, *p;
  int i;
  p = array;
  b0 = *p++;
  b1 = 0.0;
  i = n - 1;
  do {
    b2 = b1;
    b1 = b0;
    b0 = x * b1 - b2 + *p++;
  }
  while (--i);
  return (0.5 * (b0 - b2));
}

float i0(float x)
{
  float y;

  if (x < 0)
    x = -x;
  if (x <= 8.0) {
    y = (x / 2.0) - 2.0;
    return (float)(exp(x) * chbevl(y, A, 30));
  }
  return (float)(exp(x) * chbevl(32.0 / x - 2.0, B, 25) / sqrt(x));
}



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
void dft_1(float x[], float result[], int num_elems) {
  // See: "modified C code" from https://batchloaf.wordpress.com/2013/12/07/simple-dft-in-c/ 
  // to simplify does not use pre-computed cos/sin(z)
  for(int k = 0; k < num_elems; k++) {
    float xre[num_elems]; // Real component
    float xim[num_elems]; // Imaginary component
    for(int n = 0; n < num_elems; n++) {
      float z = (2 * M_PI * k * n) / num_elems;
      xre[n] += x[n] * cos(z);
      xim[n] -= x[n] * sin(z);
    }
    result[k] = (xre[k] * xre[k]) * (xim[k] * xim[k]);
  }
}

/* DFT: h is the input vector of lenth N */
complex *dft_2(complex *h, ulint N){
        
        complex *H, aux;
        register ulint n, k;
                
        H=makecm(N, 1);
        cm_zero(H, N, 1);
        
        for(n=0; n<N; n++){
            for(k=0; k<N; k++){
                aux.Re=cos((double)(2*PI*k*n)/N);
                aux.Im=sin((double)(2*PI*k*n)/N);
                H[n]=cadd(H[n], cmult(h[k], aux));
            }
        }
        
        return H;
}

/*********************** End Of DFT Calculation Section*********************/
int kaiser( float beta, int M, float *window )
{
  // # Docstring adapted from NumPy's kaiser function
  // if _len_guards(M):
  //     return np.ones(M)
  // M, needs_trunc = _extend(M, sym)
  int result = 1;
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



float freq_from_fft(float* signal, int N, int dft_elements, float fs) {
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
    printf("N:%d,\ndft_elements:%d,\nsampling rate:%f\n",N,dft_elements,fs);
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
    float* kaiser_window = (float*)malloc(N * sizeof(float));
    float beta = 100.0;
    kaiser(beta, N, kaiser_window);
    // Apply Kaiser window
    for (int n = 0; n < N; n++) {
        windowed[n] = signal[n] * kaiser_window[n];
    }
    
    // apply dft to windowed signal with dft_elements
    //please call cmsis dft function.
    dft(windowed, f, dft_elements);
    
    // Find the peak
    int i_peak = 0;
    float max_val = 0.0;
    for (int i = 0; i < computed_dft_len; i++) {
        abs_f[i] = sqrt(windowed[i] * windowed[i]); // Assuming f is complex
        printf("%d.abs_f[i] is %f\n",i,abs_f[i]);
        if (abs_f[i] > max_val) {
            max_val = abs_f[i];
            i_peak = i;
        }
    }
    printf("i_peak is %d\n",i_peak);
    // Log and interpolate
    for (int i = 0; i < computed_dft_len; i++) {
        log_abs_f[i] = log(abs_f[i]);
    }
    float i_interp = parabolic(log_abs_f, i_peak);
    // Convert to equivalent frequency
    float result = fs * i_interp / dft_elements; // Hz
    printf("i_interp is %f\n",i_interp);
    printf("result is %f\n",result);
    // Free allocated memory
    free(windowed);
    free(f);
    free(log_abs_f);
    free(abs_f);
    return result;
}

// Note: You need to implement or include the kaiser, rfft, and parabolic functions.


// void detectAudio(float* signal, int sig_len, int sr, float wanted_freq, float magthreshold, float freqthreshold, float (*freq_func)(float*, int)) 
// void detectAudio(float* signal, int sig_len, int sr, int dft_len, float wanted_freq, float magthreshold, float freqthreshold, float (*freq_func)(float*, int));

void detectAudio(float* signal, int sig_len, int sr, int dft_len, float wanted_freq, float magthreshold, float freqthreshold)
{
    // float testFreq = freq_func(signal, sr);
    // float testFreq=0.0;
    float testFreq = freq_from_fft(signal, sig_len, dft_len, sr);
    int _dft_elements = 2048;
    printf("testFreq calculated: %f\n", testFreq);
    
    float diff_freq = fabs(wanted_freq - testFreq);
    printf("diff_freq calculated: %f\n", diff_freq);
    
    float* window = (float*)malloc(sig_len * sizeof(float));
    if (window == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return;
    }

    for (int i = 0; i < sig_len; i++) {
        window[i] = signal[i] * (0.54 - 0.46 * cos(2 * M_PI * i / (sig_len - 1)));
    }

    // float complex* fftData = (float complex*)malloc(sig_len * sizeof(float complex));
    float* fftData = (float*)malloc(sig_len * sizeof(float));
    if (fftData == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        free(window);
        return;
    }
    dft(window, fftData, _dft_elements);
    // int targetIndex = (int)((float)sig_len * testFreq / sr);// get the bin we want to calculate maganitude
    int targetIndex = (int)((float)_dft_elements * testFreq / sr);
    float magnitude = sqrt(fftData[targetIndex] * fftData[targetIndex]);//cbs(fftData[targetIndex]);

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

void changeArrayTest(int ary1[])
{
    ary1[0] = 100;
    ary1[1] = 200;
}
// int ageArray[] = {2, 8, 4, 12};
// printf("ageArray[0] is %d before changeArrayTest\n",ageArray[0]);
// changeArrayTest(ageArray);
// printf("ageArray[0] is %d",ageArray[0]);

int main(void)
{
    /*
        signature
        freq_from_fft(float* signal, int N, int dft_elements, float fs, float* result);
    */
    int wav_len = sizeof(wav_array)/sizeof(wav_array[0]);
    float* float_wav = (float*)malloc(wav_len * sizeof(float));
    //convert int array to float array
    for(int idx=0; idx<wav_len; idx++)
    {
        float_wav[idx] = (float)wav_array[idx];
    }
    int _dft_len = 2048;
    float _sr = 20000.0;
    // float ret_freq = 0.0;
    float ret_freq = freq_from_fft(float_wav, wav_len, _dft_len, _sr);
    printf("ret_freq is %f", ret_freq);
    return 0;
}



