#include "mobius.hpp"

namespace gasket {

template <>
Mobius<double> Mobius<double>::toMobiusDouble() const {
    return *this;
}

}