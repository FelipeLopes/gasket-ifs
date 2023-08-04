#pragma once

#include "colorer.hpp"
#include "complex_type.hpp"
#include "diver.hpp"
#include "key_gasket.hpp"
#include "scaler.hpp"
#include "searcher.hpp"
#include "shape.hpp"

namespace gasket {

template <typename T, typename DiverT, typename ColorerT>
class Zoom {

public:
    class Builder {
    public:
        Builder() {
            static_assert(std::is_base_of<Diver<T>, DiverT>::value,
                "DiverT must implement Diver<T> interface");
            static_assert(std::is_base_of<Colorer, ColorerT>::value,
                "ColorerT must implement Colorer interface");
        }
        Builder& withShape(T r1_, T r2_, Complex<T> f_, bool flip_ = false) {
            initShape = true;
            r1 = r1_;
            r2 = r2_;
            f = f_;
            flip = flip_;
            return *this;
        }
        Builder& withScales(T iniLogscale_, T step_, int numSteps_, int precDigits_ = 10) {
            initScales = true;
            iniLogscale = iniLogscale_;
            step = step_;
            numSteps = numSteps_;
            precDigits = precDigits_;
            return *this;
        }
        Builder& withImageSize(int width_, int height_) {
            initImageSize = true;
            width = width_;
            height = height_;
            return *this;
        }
        Zoom build(DiverT diver, ColorerT colorer) {
            if (!initShape) {
                throw std::invalid_argument("Shape not initialized");
            }
            if (!initScales) {
                throw std::invalid_argument("Scales not initialized");
            }
            if (!initImageSize) {
                throw std::invalid_argument("Aspect ratio not initialized");
            }
            Shape<T> shape(r1, r2, f, flip);
            Scaler<T> scaler(iniLogscale, step, numSteps, precDigits);
            return Zoom(shape, diver, scaler, colorer, width, height);
        }
    private:
        bool initShape = false;
        T r1, r2;
        Complex<T> f;
        bool flip;

        bool initScales = false;
        T iniLogscale, step;
        int numSteps, precDigits;

        bool initImageSize = false;
        int width, height;
    };

private:
    Zoom(const Shape<T>& shape_, DiverT diver_, const Scaler<T>& scaler_, ColorerT colorer_,
        int width_, int height_): shape(shape_), diver(diver_), scaler(scaler_), colorer(colorer_),
        width(width_), height(height_) {

        Mobius<T> acc;
        int k = diver.chooseDive(acc);
        bool inverseDive = (k>=3);
        std::vector<int> diveIndices;
        diveIndices.push_back(k);
        auto pts = shape.startingPoints(inverseDive);
        auto arr = shape.diveArray(inverseDive);
        acc = acc.compose(arr[k%3]);
        std::vector<Mobius<T>> zoomTransforms;
        zoomTransforms.push_back(acc);
        for (int i=0; i<diver.getDepth()-1; i++) {
            int k = diver.chooseDive(acc);
            diveIndices.push_back(k);
            acc = acc.compose(arr[k]);
            zoomTransforms.push_back(acc);
        }
        auto center = (acc.apply(pts[0])+acc.apply(pts[1])+acc.apply(pts[2]))/Complex<T>(3);

        T ar(width, height);
        ar.canonicalize();

        Searcher<T> searcher(shape, scaler, center, inverseDive, zoomTransforms, keyGaskets, ar);

        searcher.start();
        searcher.block();

        colorer.keyGaskets(keyGaskets);

        for (auto g: keyGaskets) {
            diveIndicesMap[g.first] = diveIndices[g.second.level+1];
        }
    }
    const Shape<T>& shape;
    DiverT diver;
    const Scaler<T>& scaler;
    ColorerT colorer;
    int width, height;
    std::map<double, KeyGasket> keyGaskets;
    std::map<double, int> diveIndicesMap;
};

}