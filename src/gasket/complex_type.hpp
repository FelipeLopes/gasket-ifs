#pragma once

#include <complex>
#include <gmpxx.h>

namespace gasket {

typedef std::complex<double> cx;

template <typename T>
class Complex {
public:
    Complex(): real(T(0)), imag(T(0)) { }
    Complex(T a): real(a), imag(T(0)) { }
    Complex(T a, T b): real(a), imag(b) { }
    T real, imag;
    T norm() const {
        return real*real + imag*imag;
    }
    Complex<T> conj() const {
        return Complex<T>(real, -imag);
    }
    cx toCxDouble() const;
    Complex<double> toComplexDouble() const;
};

template<typename T>
double toDouble(T x);

template <typename T>
Complex<double> Complex<T>::toComplexDouble() const {
    return Complex<double>(toDouble(real), toDouble(imag));
}

template <typename T>
bool operator==(Complex<T> a, Complex<T> b) {
    return a.real == b.real && a.imag == b.imag;
}

template <typename T>
bool operator!=(Complex<T> a, Complex<T> b) {
    return a.real != b.real || a.imag != b.imag;
}

template <typename T>
Complex<T> operator+(Complex<T> a, Complex<T> b) {
    return Complex<T>(a.real+b.real, a.imag+b.imag);
}

template <typename T>
Complex<T> operator-(Complex<T> a) {
    return Complex<T>(-a.real, -a.imag);
}

template <typename T>
Complex<T> operator-(Complex<T> a, Complex<T> b) {
    return Complex<T>(a.real - b.real, a.imag - b.imag);
}

template <typename T>
Complex<T> operator*(T a, Complex<T> b) {
    return Complex<T>(a*b.real,a*b.imag);
}

template <typename T>
Complex<T> operator*(Complex<T> a, Complex<T> b) {
    return Complex<T>(a.real*b.real-a.imag*b.imag, a.real*b.imag+a.imag*b.real);
}

template <typename T>
Complex<T> operator/(Complex<T> a, Complex<T> b) {
    T scale = T(1)/b.norm();
    return scale*(a*b.conj());
}

template <typename T>
T squareRoot(T x);

template <typename T>
double toDouble(T x);

template <typename T>
T abs(T x);

template <typename T>
T abs(T x) {
    return x > 0 ? x : -x;
}

template <typename T>
T expPrec(T x, T prec) {
    T c = 1;
    T d = 0;
    T f = 1;
    c = 2-x + (2*x)/c;
    d = 1/(2-x+(2*x)*d);
    T x2 = x*x;
    T b = 6;
    T nf = c*d*f;
    while (abs<T>(nf-f) > prec) {
        f = nf;
        c = b + x2/c;
        d = 1/(b+x2*d);
        b = b + 4;
        nf = c*d*f;
    }
    return nf;
}

}