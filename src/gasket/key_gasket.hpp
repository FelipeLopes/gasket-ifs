#pragma once

#include "color_params.hpp"
#include "mobius.hpp"
#include <memory>
#include <vector>

namespace gasket {

class KeyGasket {
public:
    KeyGasket(std::vector<Mobius<double>> ifsTransforms, int level);
    int level = 0;
    int numTransforms() const;
private:
    std::vector<Mobius<double>> ifsTransforms;
};

}