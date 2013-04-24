#include "complex.h"
#include <cmath>

complex complex::assign(double r, double i) {
    complex res;
    res.re = r;
    res.im = i;
    return res;
}

complex complex::complex_from_polar(double r, double theta_radians) {//creating complex from polar representation
    complex res;
    res.re = r * cos(theta_radians);
    res.im = r * sin(theta_radians);
    return res;
}

double complex::complex_magnitude(complex c) {
    return sqrt(c.re*c.re + c.im*c.im);
}

complex complex::complex_add(complex left, complex right) {//adding two complex numbers
    complex res;
    res.re = left.re + right.re;
    res.im = left.im + right.im;
    return res;
}

complex complex::complex_sub(complex left, complex right) {//substraction
    complex res;
    res.re = left.re - right.re;
    res.im = left.im - right.im;
    return res;
}

complex complex::complex_mult(complex left, complex right) {//multiplication
    complex res;
    res.re = left.re*right.re - left.im*right.im;
    res.im = left.re*right.im + left.im*right.re;
    return res;
}

void complex::operator=(const complex r) {//assigning
    re = r.re;
    im = r.im;
}
