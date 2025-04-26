//#include <stdio.h>

#define ITERATIONS 25
#define PI 3.14159265358979323846264338327950288419

unsigned long long fact(int x) {
    unsigned long long fact = 1;

    for (int i = 1; i <= x; i++)
        fact *= i;

    return fact;
}

long int pow(long int x, unsigned long n) {
    long int a = x, p = 1;
    while (n > 0) {
        if ((n & 1) != 0)
            p *= a;
        a *= a;
        n >>= 1;
    }
    return p;
}

double __appr_sin(int x) {
    return (16.02385 * x * (PI - x))/(5.00198 * PI * PI - 4 * x * (PI - x));
}

double __x_sin(double x) {
    double r = 0;
    short int sign = 1;
    const double sqrx = x * x;
    double powv = x;
    for (int i = 1; i < ITERATIONS; i += 2) {
        r += sign * powv / fact(i);
        sign = -sign;
        powv *= sqrx;
    }
    return r;
}

double __x_cos(double x) {
    double r = 0;
    short int sign = 1;
    const double sqrx = x * x;
    double powv = x;
    for (int i = 0; i < ITERATIONS; i += 2) {
        r += sign * (powv / fact(i));
        sign = -sign;
        powv *= sqrx;
    }
    return r;
}

double __x_ex(double x) {
    double r = 1.0;
    double powv = x;
    for (int n = 1; n < ITERATIONS; n++) {
        r += (float)powv / fact(n);
        powv *= x;
    }
    return r;
}

void dft(double *xr, double *xi, int *a, int len) {
    int N = len;

    for (int k = 0; k < N; k++) {
        xr[k] = 0;
        xi[k] = 0;
        for (int n = 0; n < len; n++) {
            xr[k] += (double)a[n] * (double)cos((2 * PI * k * n) / N);
            xi[k] += (double)a[n] * (double)sin((2 * PI * k * n) / N);
        }
    }
}

#define cas(x) (double)(sin(x) + cos(x))
void dht(double *r, int *a, int len) {
    int N = len;

    for (int k = 0; k < N; k++) {
        r[k] = 0;
        for (int n = 0; n < len; n++) {
            r[k] += (double)a[n] * cas(((2 * PI) / N) * n * k);
        }
        //r[k] /= (double)N;
        r[k] /= sqrt(N);
    }
}

void fht(double *r, int *a, int len) {
    int N = len;

    int nN1[len];
    int nN2[len];

    for (int i = 0; i < len; i++) {
        if (i % 2)
            nN1[i] = a[i];
        else
            nN1[i] = 0;
    }

    for (int i = 0; i < len; i++) {
        if (i % 2) {
            nN2[i] = 0;
            continue;
        }
        nN2[i] = a[i];
    }

    double N1[] = {0.0, 0.0, 0.0, 0.0};
    double N2[] = {0.0, 0.0, 0.0, 0.0};

    dht(&N1, nN1, 4);
    dht(&N2, nN2, 4);

    for (int k = 0; k < N; k++) {
        r[k] = N1[k] + (N2[k] * cos((2 * PI) / N) * k) + (N2[N-k] * sin((2 * PI) / N) * k);
    }
}

int main(void) {
    printf("%f\n", __x_sin(1));
    printf("%f\n", __x_cos(1));
    printf("%f\n", __x_ex(2));
    
    int a[] = {1, 4, 9, 16};
    double fht_a[] = {0.0, 0.0, 0.0, 0.0};
    double dht_a[] = {0.0, 0.0, 0.0, 0.0};

    double dft_ar[] = {0.0, 0.0, 0.0, 0.0};
    double dft_ai[] = {0.0, 0.0, 0.0, 0.0};
    
    dht(&dht_a, a, 4);

    for (int k = 0; k < 4; k++) {
        printf("dht - (%f)\n", dht_a[k]);
    }
    printf("\n");

    fht(&fht_a, a, 4);

    for (int k = 0; k < 4; k++) {
        printf("fht - (%f)\n", fht_a[k]);
    }
    printf("\n");

    dft(&dft_ar, &dft_ai, a, 4);

    for (int k = 0; k < 4; k++) {
        printf("dft - (%f) - j(%f)\n", dft_ar[k], dft_ai[k]);
    }
}