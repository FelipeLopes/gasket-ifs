#pragma once

#include "key_gasket.hpp"
#include "scaler.hpp"
#include "sdf.hpp"
#include "shape.hpp"
#include <array>
#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>
#include <map>
#include <memory>

namespace gasket {

template <typename T>
class Searcher {
public:
    Searcher(const Shape<T>& shape_,
        const Scaler<T>& scaler_,
        Complex<T> center_,
        bool inverseDive_,
        const std::vector<Mobius<T>>& zoomTransforms_,
        std::map<double, KeyGasket>& keyGaskets_,
        T ar_, int numThreads_ = 4):
        shape(shape_), center(center_), inverseDive(inverseDive_),
        ar(ar_), numThreads(numThreads_), scaler(scaler_), threadPool(numThreads_),
        zoomTransforms(zoomTransforms_), keyGaskets(keyGaskets_), lastPickedUp(numThreads_-1) {

        pts = shape.startingPoints(inverseDive);
        transforms = shape.diveArray(inverseDive);
    }

    void start(){
        Sdf<T> sdf = Sdf<T>::fromPoints(pts[0], pts[1], pts[2]);
        T iniScale = scaler.lookupExp(0);
        T iniHeight = 2/iniScale;
        T iniWidth = iniHeight*ar;
        if (!sdf.rectInside(center, iniWidth, iniHeight)) {
            std::vector<Mobius<double>> gasketTransforms;
            auto transforms = shape.doubleSidedTransforms(scaler.lookupExp(0), center);
            gasketTransforms.insert(gasketTransforms.end(), transforms.begin(), transforms.end());
            double iniLogscaleDouble = toDouble(scaler.iniLogscale);
            KeyGasket g(gasketTransforms, iniLogscaleDouble);
            keyGaskets.insert(std::pair<double, KeyGasket>(iniLogscaleDouble, g));
        }
        for (int i=0; i<numThreads; i++) {
            boost::asio::post(threadPool, [=] {
                task(i);
            });
        }
    }

    void block() {
        threadPool.join();
    }
private:
    void task(int i) {
        Mobius<T> acc = zoomTransforms[i];
        auto qa = acc.apply(pts[0]);
        auto qb = acc.apply(pts[1]);
        auto qc = acc.apply(pts[2]);
        Sdf<T> sdf = Sdf<T>::fromPoints(qa, qb, qc);
        int scaleVal = searchScale(sdf);

        T logscale = scaler.iniLogscale + scaleVal*scaler.step;
        auto s = Mobius<T>::scaling(scaler.lookupExp(scaleVal))
            .compose(Mobius<T>::translation(-center))
            .compose(acc);

        std::vector<Mobius<double>> gasketTransforms;
        for (int i=0; i<3; i++) {
            gasketTransforms.push_back(transforms[i].conjugate(s).toMobiusDouble());
        }
        double logscaleDouble = toDouble(logscale);
        KeyGasket g(gasketTransforms, i);

        lock.lock();
        auto it = keyGaskets.find(logscaleDouble);
        if (it == keyGaskets.end() || i > it->second.level) {
            keyGaskets.insert(std::pair<double, KeyGasket>(logscaleDouble, g));
        }
        if (scaleVal >= scaler.numSteps) {
            foundEnd = true;
        }
        if (!foundEnd) {
            lastPickedUp++;
            boost::asio::post(threadPool, [=] {
                task(lastPickedUp);
            });
        }
        lock.unlock();
    }

    int searchScale(Sdf<T> sdf) {
        int lb = 0;
        int ub = scaler.numSteps;
        while (ub - lb > 1) {
            int m = (lb + ub) / 2;
            T scale = scaler.lookupExp(m);
            T height = 2/scale;
            T width = height*ar;
            if (sdf.rectInside(center, width, height)) {
                ub = m;
            } else {
                lb = m;
            }
        }
        return ub;
    }

    const Shape<T>& shape;
    std::array<Complex<T>, 3> pts;
    Complex<T> center;
    bool inverseDive;
    std::array<Mobius<T>, 3> transforms;
    T ar;
    int numThreads;
    const Scaler<T>& scaler;
    boost::asio::thread_pool threadPool;
    const std::vector<Mobius<T>>& zoomTransforms;
    std::map<double, KeyGasket>& keyGaskets;
    std::mutex lock;
    int lastPickedUp;
    bool foundEnd = false;
};

}