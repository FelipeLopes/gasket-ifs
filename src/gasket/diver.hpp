#pragma once

#include "mobius.hpp"
#include <random>

namespace gasket {

template<typename T>
class Diver {
public:
    Diver() { }
    virtual int chooseDive(Mobius<T> accumulator) const = 0;
    virtual int getDepth() const = 0;
    virtual ~Diver() { }
};

}