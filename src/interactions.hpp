#pragma once

#include <math.h>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "system.hpp"
#include "vectors.hpp"

using namespace std;

struct AtomParameters {
    double eps;
    double sig;
    double q;
    double mass;
};

struct OverwriteParameters {
    unordered_set<string> atoms;
    double eps;
    double sig;
    bool morse;
    double de;
    double a;
    double re;

};

struct InteractionParameters {
    double qbase;
    double tip_dummy_k, tip_dummy_r0;
    unordered_map<string, AtomParameters> atom_parameters;
    vector<OverwriteParameters> overwrite_parameters;
};

/* Mixing rule functions */
inline double mixsig(double sig1, double sig2) {
    return (sig1 + sig2) / 2;
}

inline double mixeps(double eps1, double eps2) {
    return sqrt(eps1 * eps2);
}

class Interaction {
 public:
    virtual ~Interaction() {};
    virtual Vec3d eval(const System& system) = 0;
 private:
};

class LJInteraction: public Interaction {
 public:
    Vec3d eval(const System& system) override;
 private:
    double es6_;
    double es12_;
};

class MorseInteraction: public Interaction {
 public:
    Vec3d eval(const System& system) override;
 private:
    double de_;
    double a_;
    double re_;
};

class HarmonicInteraction: public Interaction {
 public:
    Vec3d eval(const System& system) override;
 private:
    double k_;
    double r0_;
};

class HarmonicAngleInteraction: public Interaction {
 public:
    Vec3d eval(const System& system) override;
 private:
    double k_;
    double theta0_;
};

class GridInteraction: public Interaction {
 public:
    Vec3d eval(const System& system) override;
 private:
};
