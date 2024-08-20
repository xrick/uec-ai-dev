#include<stdio.h>
#include<math.h>
#define PI 3.14159265
 
struct DFT_Coefficient {
    double real, img;
} dft_val;
int main(int argc, char **argv) {
    int N = 10;
    float a, b, c;
    int i;
    double function[N];
    int k;
    double cosine[N];
    double sine[N];
 
    printf("Calculation DFT Coefficients\n");
    printf("Enter the coefficient of simple linear function:\n");
    printf("ax + by = c\n");
    scanf("%f", &a);
    scanf("%f", &b);
    scanf("%f", &c);
    for (i = 0; i < N; i++) {
        function[i] = (((a * (double) i) + (b * (double) i)) - c);
        //System.out.print( "  "+function[i] + "  ");
    }
 
    printf("Enter the max K value: ");
    scanf("%d", &k);
 
    for (i = 0; i < N; i++) {
        cosine[i] = cos((2 * i * k * PI) / N);
        sine[i] = sin((2 * i * k * PI) / N);
    }
 
    printf("The coefficients are: ");
 
    for (i = 0; i < N; i++) {
        dft_val.real += function[i] * cosine[i];
        dft_val.img += function[i] * sine[i];
    }
    printf("( %e )-( %e i )", dft_val.real, dft_val.img);
    return 0;
}