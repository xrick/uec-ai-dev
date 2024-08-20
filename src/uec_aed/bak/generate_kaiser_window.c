#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979323846

double bessel_i0(double x) {
    double sum = 1.0;
    double term = 1.0;
    double x_2 = x / 2.0;
    
    for (int i = 1; i <= 50; i++) {
        term *= (x_2 * x_2) / (i * i);
        sum += term;
        if (term < 1e-12 * sum) break;
    }
    
    return sum;
}

void kaiser_window(double *window, int size, double beta) {
    double arg;
    double denom = bessel_i0(beta);
    
    for (int n = 0; n < size; n++) {
        arg = beta * sqrt(1 - pow((2.0 * n / (size - 1) - 1), 2));
        window[n] = bessel_i0(arg) / denom;
    }
}

int main() {
    int size = 51;
    double beta = 5.0;
    double window[size];
    
    kaiser_window(window, size, beta);
    
    for (int i = 0; i < size; i++) {
        printf("window[%d] = %f\n", i, window[i]);
    }
    
    return 0;
}