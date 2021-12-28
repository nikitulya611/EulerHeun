#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <random>
#include <ctime>
#include <math.h>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long* idum)
{
    int j;
    long k;
    static long idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0) {
        if (-(*idum) < 1) *idum = 1;
        else *idum = -(*idum);
        idum2 = (*idum);
        for (j = NTAB + 7; j >= 0; j--) {
            k = (*idum) / IQ1;
            *idum = IA1 * (*idum - k * IQ1) - k * IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum) / IQ1;
    *idum = IA1 * (*idum - k * IQ1) - k * IR1;
    if (*idum < 0) *idum += IM1;
    k = idum2 / IQ2;
    idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
    if (idum2 < 0) idum2 += IM2;
    j = iy / NDIV;
    iy = iv[j] - idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp = AM * iy) > RNMX) return RNMX;
    else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


float gasdev(long* idum)
{
    float ran2(long* idum);
    static int iset = 0;
    static float gset;
    float fac, rsq, v1, v2;

    if (iset == 0) {
        do {
            v1 = 2.0 * ran2(idum) - 1.0;
            v2 = 2.0 * ran2(idum) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        iset = 1;
        return v2 * fac;
    }
    else {
        iset = 0;
        return gset;
    }
}

using namespace std;

long seed = time(NULL);

//  dy/dt = a - sin(y)
//  a = 0,5; a = 1,5

double foo(double x, double a)
{
    return a - sin(x);
}

void eiler(int n, double a, double h, string filename)
{
    ofstream out;
    out.open(filename);
    if (!out.is_open())
        cerr << "eiler -";

    double y0 = 0, t0 = 0, fi;

    for (int i = 0; i < n; i++)
    {
        t0 += h;
        fi = foo(y0, a);
        y0 = y0 + h * fi;
        
        out << t0 << '\t' << y0 << endl;
    }

    out.close();
}

void eilerWithRNG(int n, double a, double h, double D, string filename)
{
    ofstream out;
    out.open(filename);
    if (!out.is_open())
        cout << "eiler -";

    double y0 = 0, t0 = 0, fi;
    double c = sqrt(2 * D * h);

    for (int i = 0; i < n; i++)
    {
        t0 += h;
        fi = foo(y0, a);
        y0 = y0 + c * gasdev(&seed) + h * fi;
        
        out << t0 << '\t' << y0 << endl;
    }

    out.close();
}

void heunMethod(int n, double a, double h, string filename)
{
    ofstream out;
    out.open(filename);
    if (!out.is_open())
        cout << "heun -";

    double y0 = 0, t0 = 0, fi;

    for (int i = 0; i < n; i++)
    {
        t0 += h;
        fi = foo(y0, a);

        double y_p = y0 + h * fi;
        y0 = y0 + h / 2 * (fi + a - sin(y_p));

        out << t0 << '\t' << y0 << endl;
    }

    out.close();
}

void heunMethodRNG(int n, double a, double h, double D, string filename)
{
    ofstream out;
    out.open(filename);
    if (!out.is_open())
        cout << "heun -";

    double y0 = 0, t0 = 0;
    double c = sqrt(2 * D * h);

    for (int i = 0; i < n; i++)
    {
        t0 += h;

        double tempF = a - sin(y0);
        double y_p = y0 + c * gasdev(&seed) + h * tempF;
        double y = y0 + c * gasdev(&seed) + h / 2 * (tempF + a - sin(y_p));
        y0 = y;

        out << t0 << '\t' << y0 << endl;
    }

    out.close();
}

int main()
{
    eiler(50000, 0.5, 0.01, "d:\\eiler1.dat");
    eiler(50000, 1.5, 0.01, "d:\\eiler2.dat");
    heunMethod(50000, 0.5, 0.01, "d:\\heun1.dat");
    heunMethod(50000, 1.5, 0.01, "d:\\heun2.dat");


    eilerWithRNG(50000, 0.5, 0.01, 0.3, "d:\\eiler1RNG.dat");
    eilerWithRNG(50000, 1.5, 0.01, 0.3, "d:\\eiler2RNG.dat");
    heunMethodRNG(50000, 0.5, 0.01, 0.3, "d:\\heun1RNG.dat");
    heunMethodRNG(50000, 1.5, 0.01, 0.3, "d:\\heun2RNG.dat");
    return 0;
}