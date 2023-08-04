#pragma once

#include "color_params.hpp"
#include "key_gasket.hpp"
#include <algorithm>
#include <boost/gil.hpp>

namespace gasket {

class Colorer {

public:
    Colorer() { }
    virtual void keyGaskets(const std::map<double, KeyGasket>& keyGaskets) = 0;
    virtual ColorParams color(double logscale, int diveTransform) const = 0;
    virtual ~Colorer() { }
};

}