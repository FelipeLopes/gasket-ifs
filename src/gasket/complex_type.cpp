#include <gmpxx.h>
#include "complex_type.hpp"

namespace gasket {

template <>
Complex<double>::Complex(): real(0), imag(0) {

}

template <>
Complex<double>::Complex(double a): real(a), imag(0) {

}

template <>
Complex<double>::Complex(double a, double b): real(a), imag(b) {

}

template <>
double Complex<double>::norm() const {
    return real*real + imag*imag;
}

template <>
Complex<double> Complex<double>::conj() const {
    return Complex<double>(real, -imag);
}

template <>
cx Complex<mpq_class>::toCxDouble() const {
    mpf_class a(real);
    mpf_class b(imag);
    return a.get_d() + 1i*b.get_d();
}

template <>
double toDouble<mpq_class>(mpq_class x) {
    mpf_class f(x);
    return f.get_d();
}

template <>
cx Complex<double>::toCxDouble() const {
    return cx(real, imag);
}

template <>
mpz_class squareRoot<mpz_class>(mpz_class s) {
    if (s < 0) {
        throw std::invalid_argument("Square root of negative integer.");
    }
    if (s == 0 || s == 1) {
        return s;
    }
    mpz_class x0 = s / 2;
    mpz_class x1 = (x0 + s / x0) / 2;
    while (x1 < x0) {
        x0 = x1;
        x1 = (x0 + s / x0) / 2;
    }
    if (x0*x0 != s) {
        throw std::invalid_argument("Integer is not a perfect square.");
    }
    return x0;
}

template <>
mpq_class squareRoot<mpq_class>(mpq_class x) {
    x.canonicalize();
    auto n = squareRoot<mpz_class>(x.get_num());
    auto d = squareRoot<mpz_class>(x.get_den());
    return mpq_class(n.get_str()+"/"+d.get_str());
}

template<>
Complex<mpq_class> squareRoot<Complex<mpq_class>>(Complex<mpq_class> z) {
    auto x = z.real;
    auto y = z.imag;
    auto l = squareRoot<mpq_class>(z.norm());
    auto u = squareRoot<mpq_class>((l+x)/2);
    auto v = squareRoot<mpq_class>((l-x)/2);
    if (y < 0) {
        v = -v;
    }
    return Complex<mpq_class>(u,v);
}

template<>
Complex<double> squareRoot<Complex<double>>(Complex<double> z) {
    cx u = sqrt(cx(z.real,z.imag));
    return Complex<double>(u.real(), u.imag());
}

template <>
Complex<double> Complex<double>::toComplexDouble() const {
    return *this;
}

template <>
mpz_class squareRoot<mpz_class>(mpz_class s);

template <>
mpq_class squareRoot<mpq_class>(mpq_class x);

template <>
Complex<mpq_class> squareRoot<Complex<mpq_class>>(Complex<mpq_class> z);

}