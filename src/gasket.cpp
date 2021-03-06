#include <iostream>

#include "gasket.h"

using namespace std;

const double EPS = 1e-10;
const double PI = 2*acos(0);

Mobius::Mobius() {
}

Mobius::Mobius(cx a_, cx b_, cx c_, cx d_) {
    cx sdet = sqrt(a_*d_ - b_*c_);
    a = a_/sdet;
    b = b_/sdet;
    c = c_/sdet;
    d = d_/sdet;
}

cx Mobius::apply(cx z) {
    return (a*z+b)/(c*z+d);
}

cx Mobius::applyInf() {
    return a/c;
}

void Mobius::print() {
    printf("(%lf, %lf) (%lf, %lf)\n(%lf, %lf) (%lf, %lf)\n\n",a.real(),a.imag(),b.real(),b.imag(),
            c.real(),c.imag(),d.real(),d.imag());
}

void printCx(cx z) {
    printf("(%lf, %lf)\n",z.real(),z.imag());
}

Mobius Mobius::scaling(cx a) {
    return Mobius(a, 0, 0, 1);
}

Mobius Mobius::translation(cx b) {
    return Mobius(1, b, 0, 1);
}

Mobius Mobius::invertAndReflect() {
    return Mobius(0, 1, 1, 0);
} 

// Takes point p -> 0, q -> 1, r -> infty
Mobius Mobius::fromPoints(cx p, cx q, cx r) {
    return Mobius(q-r,-p*(q-r),(q-p),-r*(q-p));
}

// Inverts the transformation
Mobius Mobius::inverse() {
    return Mobius(d,-b,-c,a);
}

// Composes the transformation with n,
// returning the transformation given by applying n
// before the transformation
Mobius Mobius::compose(Mobius n) {
    return Mobius(a*n.a + b*n.c, a*n.b + b*n.d, c*n.a + d*n.c, c*n.b + d*n.d); 
}

// Returns a transformation that takes points
// p1 -> p2, q1 -> q2, r1 -> r2
Mobius Mobius::fromPointsToPoints(cx p1, cx q1, cx r1, cx p2, cx q2, cx r2) {
    return fromPoints(p2, q2, r2).inverse().compose(fromPoints(p1, q1, r1));
}

Mobius Mobius::conjugate(Mobius s) {
    return s.inverse().compose(*this).compose(s);
}

Mobius Mobius::diagonalize() {
    cx trace = a + d;
    cx l = 0.5*(trace + sqrt(trace*trace-4.0));
    return Mobius(b, l-d, trace-l-a, c);
}

pair<cx,cx> Mobius::fixedPoints() {
    cx s = a - d;
    cx e = sqrt((a-d)*(a-d)+4.0*b*c);
    return make_pair((s+e)/(2.0*c),(s-e)/(2.0*c));
}

shared_ptr<Transform> Mobius::decompose() {
    if (fabs(c) > EPS) {
        auto t = make_shared<Spherical>();

        Mobius preInv = Mobius::translation(d/c);
        cx preOri = conj(preInv.apply(0));
        t->pre.o = Point(preOri);
        t->pre.x = Point(conj(preInv.apply(1))-preOri);
        t->pre.y = Point(conj(preInv.apply(1i))-preOri);

        Mobius postInv = Mobius::translation(a/c).compose(Mobius::scaling(-1.0/(c*c)));
        cx postOri = postInv.apply(0);
        t->post.o = Point(postOri);
        t->post.x = Point(postInv.apply(1)-postOri);
        t->post.y = Point(postInv.apply(1i)-postOri);

        return t;
    } else {
        auto t = make_shared<Linear>();

        cx ori = apply(0);
        t->pre.o = Point(ori);
        t->pre.x = Point(apply(1)-ori);
        t->pre.y = Point(apply(1i)-ori);

        return t;
    }
}

Gasket::Gasket(double fu, double v, double fv) {
    double u = sqrt(1+v*v);
    cx ua = fu*2.0*PI*1i;
    cx va = fv*2.0*PI*1i;
    cx uc = u*exp(ua);
    cx vc = v*exp(va);

    auto tr = Mobius::translation(2);
    auto rot = Mobius::scaling(0.5*(-1+sqrt(3)*1i)).conjugate(Mobius(sqrt(3),-sqrt(3),1,1).inverse());
    auto s = Mobius(0,1i,1,0).compose(Mobius(uc,vc,conj(vc),conj(uc)));

    m.push_back(tr.conjugate(s));
    m.push_back(rot.compose(tr.inverse()).conjugate(s));

    auto root = xmlDoc.NewElement("Flames");
    root->SetAttribute("name", "gasket");
    xmlDoc.InsertFirstChild(root);

    auto flame = xmlDoc.NewElement("flame");
    flame->SetAttribute("name", "Gasket");
    flame->SetAttribute("background", "0 0 0");
    flame->SetAttribute("version", "Apophysis 2.09");
    flame->SetAttribute("size", "600 600");
    flame->SetAttribute("center", "0 0");
    flame->SetAttribute("scale", "144");
    flame->SetAttribute("oversample", "1");
    flame->SetAttribute("filter", "0.2");
    flame->SetAttribute("quality", "1");
    flame->SetAttribute("background", "0 0 0");
    flame->SetAttribute("brightness", "4");
    flame->SetAttribute("gamma", "4");
    flame->SetAttribute("gamma_threshold", "0.04");
    root->InsertEndChild(flame);

    for (auto& mob: m) {
        auto t = mob.decompose();
        t->addTransformXML(xmlDoc, flame);
    }

    auto palette = xmlDoc.NewElement("palette");
    palette->SetAttribute("count", "256");
    palette->SetAttribute("format", "RGB");
    string whiteChars(48, 'F');
    string text = "\n";
    for (int i=0; i<32; i++) {
        text += (whiteChars + "\n");
    }
    palette->SetText(text.c_str());
    flame->InsertEndChild(palette);
}

void Gasket::writeXMLFile(string filename) {
    xmlDoc.SaveFile(filename.c_str());
}

int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);

    if (argc < 7) {
        printf("Usage: %s [px] [py] [qx] [qy] [rx] [ry]\n", argv[0]);
        return 1;
    }

    cx p = cx(stod(argv[1]), stod(argv[2]));
    cx q = cx(stod(argv[3]), stod(argv[4]));
    cx r = cx(stod(argv[5]), stod(argv[6]));

    Gasket g(0.1, 0.3, 0.5);

    g.writeXMLFile("test.flame");

    return 0;
}
