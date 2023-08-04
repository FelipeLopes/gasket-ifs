#pragma once

#include "complex_type.hpp"
#include "mobius.hpp"

namespace gasket {

template <typename T>
class Shape {
public:
    Shape(T r1, T r2, Complex<T> f, bool flip = false) {
        if (r2 > r1) {
            throw std::invalid_argument("First radius parameter should be greater than "
                "or equal to the second.");
        } else if (r1 + r2 > 1) {
            throw std::invalid_argument("Radii sum cannot be larger than 1.");
        } else if (r1 < 0 || r2 < 0) {
            throw std::invalid_argument("All radii should be positive.");
        }

        if (f.norm() != 1) {
            throw std::invalid_argument("Phase must have unit norm.");
        }

        T a = -1;
        T b = 1/r1;
        T c = 1/r2;
        T s1 = a + b + c;
        T s2 = a*b + b*c + c*a;

        T e = 2*squareRoot<T>(s2);
        T d = s1 - e;
        if (d <= 0) {
            d =  s1 + e;
        }
        if (d < c) {
            throw std::invalid_argument("Radii given are not the largest circles.");
        }

        T l1 = r1 + r2;
        T l2 = 1 - r1;
        T l3 = 1 - r2;
        T cosx = (l2*l2+l3*l3-l1*l1)/(2*l2*l3);
        T sinx = squareRoot<T>(1-cosx*cosx);

        Complex<T> ii(0,1);
        Complex<T> p1 = f;
        Complex<T> p2 = p1*(Complex<T>(cosx)+Complex<T>(sinx)*ii);

        Complex<T> v1 = (Complex<T>(1)-Complex<T>(r1))*p1;
        Complex<T> v2 = (Complex<T>(1)-Complex<T>(r2))*p2;
        Complex<T> v3 = (r2*v1+r1*v2)/Complex<T>(r1+r2);

        auto m = Mobius<T>::fromPointsToPoints(
            Complex<T>(1),
            Complex<T>(-1),
            Complex<T>(0),
            p1,p2,v3);

        if (flip) {
            m = m.flip();
        }

        tr = Mobius<T>(Complex<T>(1),Complex<T>(0),
            Complex<T>(0,2),Complex<T>(1)).conjugate(m);
        rot = Mobius<T>::fromPointsToPoints(
            Complex<T>(0), Complex<T>(1), Complex<T>(-1),
            Complex<T>(1), Complex<T>(-1), Complex<T>(0)).conjugate(m);

        pa = m.apply(Complex<T>(-1));
        pb = m.apply(Complex<T>(1));
        pc = m.apply(Complex<T>(0));
    }
    std::array<Complex<T>, 3> startingPoints(bool inverseDive) const {
        if (inverseDive) {
            return {pa, pc, pb};
        } else {
            return {pa, pb, pc};
        }
    }
    std::array<Mobius<T>, 3> diveArray(bool inverseDive) const {
        if (inverseDive) {
            return {tr.inverse(), tr.inverse().conjugate(rot),
                tr.inverse().conjugate(rot.inverse())};
        } else {
            return {tr, tr.conjugate(rot), tr.conjugate(rot.inverse())};
        }
    }
    std::array<Mobius<double>, 6> doubleSidedTransforms(T scale, Complex<T> center) const {
        auto s = Mobius<T>::scaling(scale).compose(Mobius<T>::translation(-center));
        std::array<Mobius<double>, 6> ans;
        ans[0] = tr.conjugate(s).toMobiusDouble();
        ans[1] = tr.conjugate(rot).conjugate(s).toMobiusDouble();
        ans[2] = tr.conjugate(rot.inverse()).conjugate(s).toMobiusDouble();
        ans[3] = tr.inverse().conjugate(s).toMobiusDouble();
        ans[4] = tr.inverse().conjugate(rot).conjugate(s).toMobiusDouble();
        ans[5] = tr.inverse().conjugate(rot.inverse()).conjugate(s).toMobiusDouble();
        return ans;
    }
private:
    Complex<T> pa, pb, pc;
    Mobius<T> tr, rot;
};

}