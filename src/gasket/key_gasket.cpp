#include <cmath>
#include <gmpxx.h>
#include <memory>
#include "key_gasket.hpp"
#include "mobius.hpp"

namespace gasket {

using std::vector;

KeyGasket::KeyGasket(vector<Mobius<double>> ifsTransforms_, int level_):
    level(level_), ifsTransforms(ifsTransforms_) {

}

int KeyGasket::numTransforms() const {
    return ifsTransforms.size();
}

}