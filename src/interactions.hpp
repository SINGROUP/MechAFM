#pragma once

#include <math.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "vectors.hpp"

using namespace std;

struct AtomParameters {
    double eps;
    double sig;
    double q;
    double mass;
};

struct OverwriteParameters {
    unordered_multiset<string> atoms;
    double eps;
    double sig;
    bool morse;
    double de;
    double a;
    double re;

};

struct PossibleBond {
    unordered_multiset<string> atoms;
    double r0;
};

struct InteractionParameters {
    double qbase;
    double tip_dummy_k, tip_dummy_r0;
    double bond_k, angle_k, substrate_k;
    unordered_map<string, AtomParameters> atom_parameters;
    vector<OverwriteParameters> overwrite_parameters;
    vector<PossibleBond> possible_bonds_;
};

/* Mixing rule functions */
inline double mixsig(double sig1, double sig2) {
    return (sig1 + sig2) / 2;
}

inline double mixeps(double eps1, double eps2) {
    return sqrt(eps1 * eps2);
}

class System;

class Interaction {
 public:
    virtual ~Interaction() {};
    virtual void eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, double& energy) const = 0;
 private:
};

class LJInteraction: public Interaction {
 public:
    LJInteraction():
        atom_i1_(0), atom_i2_(0), es6_(0), es12_(0) {};
    LJInteraction(int atom_i1, int atom_i2, double es6, double es12):
        atom_i1_(atom_i1), atom_i2_(atom_i2), es6_(es6), es12_(es12) {};
    void eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, double& energy) const override;

 private:
    int atom_i1_, atom_i2_;   // Atom indices in the system vectors
    double es6_;
    double es12_;
};

class MorseInteraction: public Interaction {
 public:
    MorseInteraction():
        atom_i1_(0), atom_i2_(0), de_(0), a_(0), re_(0) {};
    MorseInteraction(int atom_i1, int atom_i2, double de, double a, double re):
        atom_i1_(atom_i1), atom_i2_(atom_i2), de_(de), a_(a), re_(re) {};
    void eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, double& energy) const override;

 private:
    int atom_i1_, atom_i2_;
    double de_;
    double a_;
    double re_;
};

class CoulombInteraction: public Interaction {
 public:
    CoulombInteraction():
        atom_i1_(0), atom_i2_(0), qq_(0) {};
    CoulombInteraction(int atom_i1, int atom_i2, double qq):
        atom_i1_(atom_i1), atom_i2_(atom_i2), qq_(qq) {};
    void eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, double& energy) const override;

 private:
    int atom_i1_, atom_i2_;
    double qq_;
};

class Harmonic2DInteraction: public Interaction {
 public:
    Harmonic2DInteraction():
        atom_i1_(0), atom_i2_(0), k_(0), r0_(0) {};
    Harmonic2DInteraction(int atom_i1, int atom_i2, double k, double r0):
        atom_i1_(atom_i1), atom_i2_(atom_i2), k_(k), r0_(r0) {};
    void eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, double& energy) const override;

 private:
    int atom_i1_, atom_i2_;
    double k_;
    double r0_;
};

class SubstrateInteraction: public Interaction {
 public:
    SubstrateInteraction():
        atom_i_(0), k_(0), z0_(0) {};
    SubstrateInteraction(int atom_i, double k, double z0):
        atom_i_(atom_i), k_(k), z0_(z0) {};
    void eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, double& energy) const override;

 private:
    int atom_i_;
    double k_;
    double z0_;
};

class HarmonicInteraction: public Interaction {
 public:
    HarmonicInteraction():
        atom_i1_(0), atom_i2_(0), k_(0), r0_(0) {};
    HarmonicInteraction(int atom_i1, int atom_i2, double k, double r0):
        atom_i1_(atom_i1), atom_i2_(atom_i2), k_(k), r0_(r0) {};
    void eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, double& energy) const override;

 private:
    int atom_i1_, atom_i2_;
    double k_;
    double r0_;
};

class HarmonicAngleInteraction: public Interaction {
 public:
    HarmonicAngleInteraction():
        atom_i1_(0), atom_i2_(0), atom_i3_(0), k_(0), theta0_(0) {};
    HarmonicAngleInteraction(int atom_i1, int atom_i2, int atom_i3, double k, double theta0):
        atom_i1_(atom_i1), atom_i2_(atom_i2), atom_i3_(atom_i3), k_(k), theta0_(theta0) {};
    void eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, double& energy) const override;

 private:
    int atom_i1_, atom_i2_, atom_i3_;
    double k_;
    double theta0_;
};

class GridInteraction: public Interaction {
 public:
    void eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, double& energy) const override;

 private:
};