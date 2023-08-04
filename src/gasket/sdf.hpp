#pragma once

#include <array>
#include "complex_type.hpp"

namespace gasket {

template <typename T>
class Sdf {
public:
    Sdf(T a_, T b_, T c_, T d_):
        a(a_),b(b_),c(c_),d(d_) { }

    bool inside(Complex<T> z) {
        T x = z.real;
        T y = z.imag;
        return a*(x*x+y*y) + b*x + c*y + d < 0;
    }

    bool rectInside(Complex<T> p, T w, T h) {
        Complex<T> ii(0,1);

        Complex<T> ra = p - Complex<T>(w/2) - ii*Complex<T>(h/2);
        Complex<T> rb = p + Complex<T>(w/2) - ii*Complex<T>(h/2);
        Complex<T> rc = p - Complex<T>(w/2) + ii*Complex<T>(h/2);
        Complex<T> rd = p + Complex<T>(w/2) + ii*Complex<T>(h/2);

        if (a < 0 && circRectCollision(p, w, h)) {
            return false;
        } else {
            return inside(ra) && inside(rb) && inside(rc) && inside(rd);
        }
    }

    bool circRectCollision(Complex<T> p, T w, T h) {
        T nb = b/a;
        T nc = c/a;
        T nd = d/a;
        double xc = -toDouble(nb)/2;
        double yc = -toDouble(nc)/2;
        double r = sqrt(xc*xc+yc*yc-toDouble(nd));
        double xp = toDouble(p.real);
        double yp = toDouble(p.imag);
        double dx = fabs(xc-xp);
        double dy = fabs(yc-yp);
        if ((dx > w/2 + r) || (dy > h/2 + r)) {
            return false;
        } else if ((dx < w/2) || (dy < h/2)) {
            return true;
        } else {
            return (dx-w/2)*(dx-w/2)+(dy-h/2)*(dy-h/2) < r*r;
        }
    }
    Sdf flip() {
        return Sdf(-a,-b,-c,-d);
    }
    static Sdf fromPoints(Complex<T> p, Complex<T> q, Complex<T> r) {
        T xp = p.real;
        T yp = p.imag;
        T np = xp*xp + yp*yp;
        T xq = q.real;
        T yq = q.imag;
        T nq = xq*xq + yq*yq;
        T xr = r.real;
        T yr = r.imag;
        T nr = xr*xr + yr*yr;

        T a = det33({std::array<T,3>
            {xp, yp, 1},
            {xq, yq, 1},
            {xr, yr, 1}
        });
        T b = -det33({std::array<T,3>
            {np, yp, 1},
            {nq, yq, 1},
            {nr, yr, 1}
        });
        T c = det33({std::array<T,3>
            {np, xp, 1},
            {nq, xq, 1},
            {nr, xr, 1}
        });
        T d = -det33({std::array<T,3>
            {np, xp, yp},
            {nq, xq, yq},
            {nr, xr, yr}
        });
        return Sdf<T>(a, b, c, d);
    }
private:
    T a, b, c, d;
    static T det33(std::array<std::array<T,3>,3> m) {
        return m[0][0]*m[1][1]*m[2][2]+m[0][1]*m[1][2]*m[2][0]+m[0][2]*m[1][0]*m[2][1]-
            m[0][0]*m[1][2]*m[2][1]-m[0][1]*m[1][0]*m[2][2]-m[0][2]*m[1][1]*m[2][0];
    }
};

}