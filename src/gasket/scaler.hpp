#pragma once

#include <vector>
#include "complex_type.hpp"

namespace gasket {

template <typename T>
class Scaler {
public:
    Scaler(T iniLogscale_, T step_, int numSteps_, int precDigits):
        iniLogscale(iniLogscale_), step(step_), numSteps(numSteps_) {

        mpq_class prec(1);
        for (int i=0; i<precDigits; i++) {
            prec = prec / 10;
        }
        base = expPrec<T>(iniLogscale, prec);
        lookup.resize(32-__builtin_clz(numSteps));
        for (int i=0; i<lookup.size(); i++) {
            lookup[i] = expPrec<T>((1<<i)*step, prec);
        }
    }
    T lookupExp(int n) const {
        T ans = base;
        while (n > 0) {
            int bits = (n & -n);
            ans = ans * lookup[__builtin_ctz(bits)];
            n -= bits;
        }
        return ans;
    }
    const T iniLogscale, step;
    const int numSteps;
private:
    T base;
    std::vector<T> lookup;
};

}