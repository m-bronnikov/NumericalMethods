#ifndef COMPLEX_H
#define COMPLEX_H

#include <stdio.h>
#include <math.h>

#define EPSILON 0.00000000000001

typedef struct Complex{
    double re;
    double im;
}Complex;

Complex complex_num(){
    Complex ans;
    ans.re = ans.im = 0.0;
    return ans;
}

void print_complex(Complex number){
    printf("%lg", number.re);
    if(fabs(number.im) > EPSILON){
        if(number.im > 0.0){
            printf("+");
        }
        printf("%lg", number.im);
    }
}

double modulo_complex(Complex number){
    return sqrt(number.re*number.re + number.im*number.im);
}

Complex multiply_complex(Complex z1, Complex z2){
    Complex ans;
    ans.re = z1.re*z2.re - z1.im*z2.im;
    ans.im = z1.re*z2.im + z2.re*z1.im;
    return ans;
}

Complex plus_complex(Complex z1, Complex z2){
    Complex ans;
    ans.re = z1.re + z2.im;
    ans.im = z1.im + z2.im;
    return ans;
}

Complex minus_complex(Complex z1, Complex z2){
    Complex ans;
    ans.re = z1.re - z2.im;
    ans.im = z1.im - z2.im;
    return ans;
}

void quadrat_solver_first(double a, double b, double c, Complex* ans1, Complex* ans2){
    double D = b*b - 4*a*c;
    *ans2 = *ans1 = complex_num();
    ans2->re = ans1->re = -b / (2.0 * a); 
    if(D < 0.0){
        ans1->im = sqrt(-D) / (2.0 * a);
        ans2->im = -ans1->im;
    }else{
        double incr = sqrt(D) / (2.0 * a);
        ans1->re += incr;
        ans2->re -= incr;
    }
}

#endif