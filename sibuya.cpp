/***********************************************************/
/* Sibuya family for r=1 (a=1/2)
 * Method 2
 * For details, see the paper.
 *
 * This version
 * 2 May 2023
 * Ad Ridder
 *
 */
/***********************************************************/

 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <complex.h>

using namespace std; 

#ifndef PI
# define PI	3.14159265358979323846264338327950288
#endif

/********************* auxiliary *****************/
double intrealpow(double x, int n)
/* x^n */
{
    int i;
    double y;

    y = 1.0;
    for (i=0; i<n; i++)
        y *= x;
    return y;
}

complex<double> intcomplexpow(complex<double> z, int n)
/* z^n */
{
    int i;
    complex<double> w;

    w = 1.0;
    for (i=0; i<n; i++)
        w *= z;
    return w;
}


/*********************** the FFT procedures ***************************/

complex<double> cspiegel(complex<double> &z)
/* Input z = a+b*i, output w = b+a*i
*/
{
    return complex<double>(imag(z),real(z));
}

void fft(int n, complex<double> *a)
/* Fast fourier algoritme.
   To use in combination with procedure terug_zetten().
   Input: 
   two-power integer n
   array a of length n of complex numbers.
   Output: 
   same array a with the fourier transformed values
*/
{
    complex<double> aux,u,w,t,tt;
    int i,j,k,ip,ke2,m,*pow2,m1,n1,ke1;
    double aux1;


    m = (int) floor(log(n)/log(2.0) + 0.5);
    pow2 = new int [m+1];
    for(i=1,j=0; j<=m; i*=2,j++)
        pow2[j]=i;

    m1 = m-1;
    n1 = n-1;
    ke1 = pow2[m1];

    j = 0;
    for(i=0; i<n1; i++) {
        if (i < j) {
            aux  = a[i];
            a[i] = a[j];
            a[j] = aux;
        }
        k = ke1;
        ip = m1;
        while (k <= j) {
            j -= k;
            k = pow2[--ip];
        }
        j += k;
    }

    for (i=0; i<n; i+=2) {
        ip = i+1;
        t = a[ip];
        a[ip] = a[i] - t;
        a[i]  += t;
    }

    for (k=2; k<=m; k++) {
        ke2 = pow2[k-2];
        ke1 = 2*ke2;
        u = 1.0;
        aux1 = PI/ke1;
        w = complex<double>(cos(aux1), sin(aux1));

        i = 0;
        while (i < n) {
            ip = i+ke1;
            t  = a[ip];
                a[ip] = a[i] - t;
                a[i]  += t;
                ip += ke2;
                i  += ke2;
                t = conj(cspiegel(a[ip]));
                a[ip] = a[i] + t;
                a[i]  -= t;
                i += ke1+ke2;
        }
        u = w;
        for(j=1; j<ke2; j++) {
            i = j;
            while (i < n) {
                    ip = i+ke1;
                    t = a[ip] * u;
                    a[ip] = a[i] - t;
                    a[i]  += t;
                    ip += ke2;
                    i  += ke2;
                    tt = a[ip]*u;
                    t = conj(cspiegel(tt));
                    a[ip] = a[i] + t;
                    a[i]  -= t;
                    i += ke2+ke1;
            }
            u *= w;
        }
    }
    delete pow2;
}


void terug_zetten(int n, complex<double> *a)
/* Puts the fourier transformd numbers in the proper order
   To use after procedure fft().
*/
{
    int i;
    double sumi = 0.0;
    double sumr = 0.0;

    for(i=0; i<n; i++) {
        sumr += fabs(real(a[i]));
          sumi += fabs(imag(a[i]));
    }
    if (sumi > sumr*1.0e-10) {
        cout << "the imaginary part of the result seems too large\n";
        cout << "error in input is probable, check input please\n";
        cout << "program stops\n";
        return;
    }

    for(i=0; i<n; i++)
    a[i] = complex<double>(real(a[i]), real(a[i]));

    a[0] = complex<double>(real(a[0])/n, imag(a[0]));
    for (i=n-1; i>0; i--)
        a[n-i]=complex<double>(imag(a[i])/n, imag(a[n-i]));
}

/*********************** EDM functions ***********************************/

double funV(double m)
{
    return 2.0 * (1.0 - sqrt(1.0-m));
}

double funA(double m)
{
    return 1.0 - sqrt(1.0-m) - log(funV(m)) + log(m);
}

double funpsi(double m)
{
    double u = 1.0 - sqrt(1.0-m);
    return log(u) - u;
    /* return log(funV(m)) - (1.0 - sqrt(1.0-m)); */
}


double funphi(double m)
{
    /*
    m1 = 1 - m;
    return 0.5*m + 0.3333333*(1.0 - m1*sqrt(m1)); 
    */
    double u = 1.0 - sqrt(1.0-m);
    return u * (2.0 - 1.5*u + 0.33333333*u*u);
}

double funVp(double m, double p)
{
    return p * funV(m/p);
}
double funpsip(double m, double p)
{
    return funpsi(m/p);
}


double funphip(double m, double p)
{
    return p * funphi(m/p);
}

/******************** moments **********************/
void momentsexact(double m, double p)
{
    double mp1 = 1.0 - m / p;
    double v = 2.0 * p * (1.0 - sqrt(mp1));
    double vd1 = 1.0 / sqrt(mp1);
    double vd2 = 0.5 / (mp1 * sqrt(mp1) * p);

    printf("from the variance function\n");
    printf("mean     = %lf\n", m);
    printf("variance = %lf\n", v);
    printf("skewness = %lf\n", vd1 / sqrt(v));
    printf("kurtosis = %lf\n", vd1 * vd1 / v + vd2);
}

void momentsdist(double pr[], int N)
{
    int i;
    double m,s2,s,m3,m4,b1,g2;

    m = 0.0;
    for (i=1; i<N+1; i++) 
        m += i * pr[i];

    s2 = 0.0;
    m3 = 0.0;
    m4 = 0.0;

    for (i=0; i<N+1; i++) {
        s2 += intrealpow(i-m,2) * pr[i];
        m3 += intrealpow(i-m,3) * pr[i];
        m4 += intrealpow(i-m,4) * pr[i];
    }

    s = sqrt(s2);
    b1 = m3 / intrealpow(s,3); 
    g2 = m4 / intrealpow(s2,2) - 3.0; /* klopt */  

    printf("\nfrom the distribution\n");
    printf("mean     = %.6lf\n", m);
    printf("variance = %.6lf\n", s2);
    printf("skewness = %.6lf\n", b1); 
    printf("kurtosis = %.6lf\n", g2);
}

   

/********************* GF **************************/
complex<double> ell(complex<double> z, int n)
/* series upto n-th term */
{
    int i;
    double k,kpow;
    complex<double> w,t,zkf;

    zkf = z; /* z^k / k! */
    w = 2.0 * zkf;

    zkf = zkf * z / 2.0;
    w = w + zkf;

    for (i=3; i<=n; i++) {
        k = 1.0 * i;
        zkf = zkf * z / k;
        kpow = intrealpow(k,i-3);  /* k^{k-3} */
        t = zkf * kpow;
        w = w + 2.0 * t;
    }
    return w;
}


complex<double> GF(complex<double> z, double exth, int n)
{
    complex<double> v = exth * z;
    complex<double> x = ell(v,n);
    return exp(x);
}

complex<double> GFp(complex<double> z, double p, double exth, int n)
{
    complex<double> v = exth * z;
    complex<double> x = ell(v,n);
    return exp(p * x);
}

void eenheidscirkel(double p, double exth, int n, int k, complex<double> *a)
/* Given a the generating function G and integer n,
   this function computes the n-th unit roots G(exp(2*pi*k*i/n))
   and stores in array a.
   The series G is computed upto the k-th term
*/
{
    int j;
    double u;
    complex<double> z,w,x;
    
    z = 1.0;
    w = GFp(z,p,exth,k);
    a[0] = real(w);
    u = 2.0*PI/((double) n);
    for(j=1; j<n/2; j++) {
        x = complex<double>(0.0,j*u);
        z = exp(x);
        w = GFp(z,p,exth,k);        
        a[j] = w;
        a[n-j] = conj(w);
    }
    z = -1.0;
    w = GFp(z,p,exth,k);
    a[n/2] = real(w);
}



void FastFourier(double m, double p, int n, int k, complex<double> *a)
/* This procedure executes the FFT on the generating function 
*/
{
    double exth;
    exth = exp(funpsip(m,p));
    eenheidscirkel(p,exth,n,k,a);
    fft(n,a);
    terug_zetten(n,a);
}

/******************** fit **********************/
double bisection(double m, double v)
/* given mean m and variance v,
   determines the dispersion parameter p
*/
{
    double pa,pb,pc,fa,fb,fc,t;

    t = 0.1;
    pa = m + t;
    fa = funVp(m,pa) - v;
    while (fa<0.0) {
        t *= 0.5;
        pa = m + t;
        fa = funVp(m,pa) - v;
    }
    t = 1.0;
    pb = m + t;    
    fb = funVp(m,pb) - v;
    while (fb>0.0) {
        t += 1.0;
        pb = m + t;
        fb = funVp(m,pb) - v;
    }
    while (pb-pa>1.0e-8) {
        pc = 0.5 * (pa+pb);
        fc = funVp(m,pc) - v;
        if (fc>0.0)
            pa = pc;
        else
            pb = pc;
    }
    return 0.5 * (pa+pb);
}

void distribution(double m, double p, int k, int N)
{
    int i,n;
    double exkap,*pr,sumpr;
    complex<double> *a;     

    n = 256; /* length of the Fourier series */
    a = new complex<double>[n];
    FastFourier(m,p,n,k,a);
    exkap = exp(-funphip(m,p)); 

    sumpr = 0.0;
    pr = new double [N+1];
    for (i=0; i<N+1; i++) {
       pr[i] = real(a[i]) * exkap;
       sumpr += pr[i];
       printf("%3d  %.10lf\n", i, pr[i]); 
    }
    printf("\nsum = %.10lf\n", sumpr);
    momentsexact(m,p);
    momentsdist(pr,N);
} 

/******************** main ********************/
int main(void)
{
    int N,k;
    double m,p,R,theta,v;

    m = 3.4; /* mean */
    v = 6.2; /* variance */
    k = 100; /* highest power in generating function */
    N = 100; /* highest value for the distribution */ 

    p = bisection(m,v);
    printf("dispersion parameter p = %lf\n", p);
    theta = funpsip(m,p);
    R = exp(-theta-1.0);
    printf("radius of convergence: %lf\n", R);
    if (R>1.0) 
        distribution(m,p,k,N);
    else
        printf("something is wrong\n");
    return 0; 
} 

