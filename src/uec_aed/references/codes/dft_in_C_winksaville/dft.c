// Compute Descrete Fouier Transform (DFT)
//
// from: https://en.wikipedia.org/wiki/Discrete_Fourier_transform
//
// X[k] = sum(n=0, n=N-1, x[n] * (e^(-((i * 2 * PI) / N) * k * n)))
//      = sum(n=0, n=N-1, x[n] * (cos(2 * PI * k * n) / N) - (i * sin(2 *( PI * k * n) / N)))

#include <stdint.h>
#include <math.h>

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

