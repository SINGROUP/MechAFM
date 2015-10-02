#pragma once

#include <cmath>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "force_grid.hpp"
#include "globals.hpp"
#include "vectors.hpp"

class System;

using namespace std;

// Contains the base parameters of an atom
struct AtomParameters {
    double eps;
    double sig;
    double q;
    double mass;
};

// Contains the definition for overridden interactions
struct OverwriteParameters {
    unordered_multiset<string> atoms;
    double eps;
    double sig;
    bool morse;
    double de;
    double a;
    double re;

};

// Contains the definition for possible bond between atoms
struct PossibleBond {
    unordered_multiset<string> atoms;
    double r0;  // The expected bond length
};

// Contains all the parameters for simulation interactions
struct InteractionParameters {
    double qbase;
    double tip_dummy_k, tip_dummy_r0;
    double bond_k, angle_k, dihedral_k;
    double substrate_eps, substrate_sig, substrate_lambda, substrate_rc, substrate_k;
    unordered_map<string, AtomParameters> atom_parameters;
    vector<OverwriteParameters> overwrite_parameters;
    vector<PossibleBond> possible_bonds_;
};

// Mixing rule functions
inline double mixsig(double sig1, double sig2) {
    return (sig1 + sig2) / 2;
}

inline double mixeps(double eps1, double eps2) {
    return sqrt(eps1 * eps2);
}


// Pure virtual interface for all the interactions
class Interaction {
 public:
    virtual ~Interaction() {};
    // Evaluates the forces and energies of the interaction for the given positions
    virtual void eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, vector<double>& energies) const = 0;
    // Return whether the interaction is between the tip and the surface or not
    virtual bool isTipSurface() const = 0;
 private:
};

class LJInteraction: public Interaction {
 public:
    LJInteraction():
        atom_i1_(0), atom_i2_(0), es6_(0), es12_(0) {};
    LJInteraction(int atom_i1, int atom_i2, double es6, double es12):
        atom_i1_(atom_i1), atom_i2_(atom_i2), es6_(es6), es12_(es12) {};
    void eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, vector<double>& energies) const override;
    bool isTipSurface() const override {
        // The interactions are build such that this holds
        return atom_i1_ == 1;
    }

 private:
    int atom_i1_, atom_i2_;  // Atom indices in the state vectors
    // Interaction constants
    double es6_;
    double es12_;
};

class MorseInteraction: public Interaction {
 public:
    MorseInteraction():
        atom_i1_(0), atom_i2_(0), de_(0), a_(0), re_(0) {};
    MorseInteraction(int atom_i1, int atom_i2, double de, double a, double re):
        atom_i1_(atom_i1), atom_i2_(atom_i2), de_(de), a_(a), re_(re) {};
    void eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, vector<double>& energies) const override;
    bool isTipSurface() const override {
        // The interactions are build such that this holds
        return atom_i1_ == 1;
    }

 private:
    int atom_i1_, atom_i2_;  // Atom indices in the state vectors
    // Interaction constants
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
    void eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, vector<double>& energies) const override;
    bool isTipSurface() const override {
        // The interactions are build such that this holds
        return atom_i1_ == 1;
    }

 private:
    int atom_i1_, atom_i2_;  // Atom indices in the state vectors
    // Interaction constants
    double qq_;
};

class GridInteraction: public Interaction {
 public:
    GridInteraction(ForceGrid& fg): force_grid(fg) {};
    void eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, vector<double>& energies) const override;
    bool isTipSurface() const override {
        return true;
    }

 private:
    ForceGrid force_grid; // The force grid containing the samples
};

class TipHarmonicInteraction: public Interaction {
 public:
    TipHarmonicInteraction():
        atom_i1_(0), atom_i2_(0), k_(0), r0_(0) {};
    TipHarmonicInteraction(int atom_i1, int atom_i2, double k, double r0):
        atom_i1_(atom_i1), atom_i2_(atom_i2), k_(k), r0_(r0) {};
    void eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, vector<double>& energies) const override;
    bool isTipSurface() const override {
        return false;
    }

 private:
    int atom_i1_, atom_i2_;  // Atom indices in the state vectors
    // Interaction constants
    double k_;
    double r0_;
};

class XYHarmonicInteraction: public Interaction {
 public:
    XYHarmonicInteraction():
        atom_i_(0), k_(0), p0_(0) {};
    XYHarmonicInteraction(int atom_i, double k, Vec2d p0):
        atom_i_(atom_i), k_(k), p0_(p0) {};
    void eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, vector<double>& energies) const override;
    bool isTipSurface() const override {
        return false;
    }

 private:
    int atom_i_;  // Atom indices in the state vectors
    // Interaction constants
    double k_;
    Vec2d p0_;
};

class SubstrateInteraction: public Interaction {
 public:
    SubstrateInteraction(): 
      atom_i_(0), eps_(0), sig_(0), lambda_(0), rc_(0), multiplier_(0), ulj_(0), ushift_(0), z0_(0) {};
    SubstrateInteraction(int atom_i, double eps, double sig, double rc, double lambda):
        atom_i_(atom_i), eps_(eps), sig_(sig), lambda_(lambda), rc_(rc)
        {
        multiplier_ = 2 * PI * eps_ * pow(sig_/lambda_, 2);
        ulj_ = (PI/(lambda_*lambda_)) * 4 * eps_ * ( pow(sig_/rc_, 12) - pow(sig_/rc_, 6) );
	ushift_ = multiplier_ * ( (12./5) * pow(sig_/rc_, 10) - 3 * pow(sig_/rc_, 4) );
        z0_ = -sig_;
        };
    void eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, vector<double>& energies) const override;
    bool isTipSurface() const override {
        return false;
    }

 private:
    int atom_i_;  // Atom indices in the state vectors
    // Interaction constants
    double eps_;
    double sig_;
    double lambda_;
    double rc_;
    double multiplier_;
    double ulj_;
    double ushift_;
    double z0_;
};

class HarmonicInteraction: public Interaction {
 public:
    HarmonicInteraction():
        atom_i1_(0), atom_i2_(0), k_(0), r0_(0) {};
    HarmonicInteraction(int atom_i1, int atom_i2, double k, double r0):
        atom_i1_(atom_i1), atom_i2_(atom_i2), k_(k), r0_(r0) {};
    void eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, vector<double>& energies) const override;
    bool isTipSurface() const override {
        return false;
    }

 private:
    int atom_i1_, atom_i2_;  // Atom indices in the state vectors
    // Interaction constants
    double k_;
    double r0_;
};

class HarmonicAngleInteraction: public Interaction {
 public:
    HarmonicAngleInteraction():
        shared_i_(0), atom_i1_(0), atom_i2_(0), k_(0), theta0_(0) {};
    HarmonicAngleInteraction(int shared_i, int atom_i1, int atom_i2, double k, double theta0):
        shared_i_(shared_i), atom_i1_(atom_i1), atom_i2_(atom_i2), k_(k), theta0_(theta0) {};
    void eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, vector<double>& energies) const override;
    bool isTipSurface() const override {
        return false;
    }

 private:
    int shared_i_, atom_i1_, atom_i2_;  // Atom indices in the state vectors
    // Interaction constants
    double k_;
    double theta0_;
};

class HarmonicDihedralInteraction: public Interaction {
 public:
    HarmonicDihedralInteraction():
        atom_i1_(0), atom_i2_(0), atom_i3_(0), atom_i4_(0), k_(0), sigma0_(0) {};
    HarmonicDihedralInteraction(int atom_i1, int atom_i2, int atom_i3, int atom_i4, double k, double sigma0):
        atom_i1_(atom_i1), atom_i2_(atom_i2), atom_i3_(atom_i3), atom_i4_(atom_i4), k_(k), sigma0_(sigma0) {};
    void eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, vector<double>& energies) const override;
    bool isTipSurface() const override {
        return false;
    }

 private:
    int atom_i1_, atom_i2_, atom_i3_, atom_i4_;  // Atom indices in the state vectors
    // Interaction constants
    double k_;
    double sigma0_;
};
