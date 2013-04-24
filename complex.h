#ifndef COMPLEX_H
#define COMPLEX_H

#include <iostream>
#include <cmath>
using namespace std;
//this class implements complex numbers
class complex {
    private:
        double re;
        
        double im;
        
    public:
        complex assign(double r, double i);//assign values to real and imaginary part
        
        complex complex_from_polar(double r, double theta_radians);//form the complex number using r and theta
        
        double  complex_magnitude(complex c);//returns the magnitude of complex number
        
        complex complex_add(complex left, complex right); //adding two complex numbers
        
        complex complex_sub(complex left, complex right);//substracting two complex numbers
        
        complex complex_mult(complex left, complex right);//multiplying two complex numbers
        
        void operator=(const complex right);//assigning
        
        double get_real() {return re;}
};
#endif
