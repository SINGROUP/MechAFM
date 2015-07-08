#include "interactions.hpp"

#include <cmath>

#include "globals.hpp"
#include "vectors.hpp"

void LJInteraction::eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, vector<double>& energies) const {
    Vec3d r_vec = positions[atom_i1_] - positions[atom_i2_];
    double r_sqr = r_vec.lensqr();
    double r6 = r_sqr * r_sqr * r_sqr;
    double term_a = es12_ / (r6*r6);
    double term_b = es6_ / r6;
    double e = term_a - term_b;
    Vec3d f = (12*term_a - 6*term_b) / r_sqr * r_vec;
    energies[atom_i1_] += e;
    energies[atom_i2_] += e;
    forces[atom_i1_] += f;
    forces[atom_i2_] -= f;
}

void MorseInteraction::eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, vector<double>& energies) const {
    Vec3d r_vec = positions[atom_i1_] - positions[atom_i2_];
    double r = r_vec.len();
    double d_exp = exp(- a_ * (r - re_));
    double e = de_ * (pow(d_exp, 2) - 2 * d_exp + 1);
    double f_abs = 2 * de_ * a_ * (pow(d_exp, 2) - d_exp);
    Vec3d f = f_abs / r * r_vec;
    energies[atom_i1_] += e;
    energies[atom_i2_] += e;
    forces[atom_i1_] += f;
    forces[atom_i2_] -= f;
}

void CoulombInteraction::eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, vector<double>& energies) const {
    Vec3d r_vec = positions[atom_i1_] - positions[atom_i2_];
    double r = r_vec.len();
    double e = qq_ / r;
    Vec3d f = qq_ / (r*r*r) * r_vec;
    energies[atom_i1_] += e;
    energies[atom_i2_] += e;
    forces[atom_i1_] += f;
    forces[atom_i2_] -= f;
}

void Harmonic2DInteraction::eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, vector<double>& energies) const {
    Vec3d r_vec = positions[atom_i1_] - positions[atom_i2_];
    Vec2d r_2d = r_vec.getXY();
    double r = r_2d.len();
    double dr = r - r0_;
    double e = k_ * dr * dr;
    Vec3d f;
    if (r > TOLERANCE) {
        f = Vec3d((-2 * k_ * dr / r * r_2d), 0);
    }
    energies[atom_i1_] += e;
    energies[atom_i2_] += e;
    forces[atom_i1_] += f;
    forces[atom_i2_] -= f;
}

void SubstrateInteraction::eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, vector<double>& energies) const {
    double dz = positions[atom_i_].z - z0_;
    if (dz < 0) {
        energies[atom_i_] = k_ * pow(dz, 2);
        forces[atom_i_].z += -2 * k_ * dz;
    }
}

void HarmonicInteraction::eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, vector<double>& energies) const {
    Vec3d r_vec = positions[atom_i1_] - positions[atom_i2_];
    double r = r_vec.len();
    double dr = r - r0_;
    double e = k_ * pow(dr, 2);
    Vec3d f = -2 * k_ * dr / r * r_vec;
    energies[atom_i1_] += e;
    energies[atom_i2_] += e;
    forces[atom_i1_] += f;
    forces[atom_i2_] -= f;
}

void HarmonicAngleInteraction::eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, vector<double>& energies) const {
    Vec3d r1_vec = positions[atom_i1_] - positions[shared_i_];
    Vec3d r2_vec = positions[atom_i2_] - positions[shared_i_];
    double r1 = r1_vec.len();
    double r2 = r2_vec.len();
    double cos_t = r1_vec.dot(r2_vec) / (r1 * r2);
    double theta = acos(cos_t);
    double sin_t = sin(theta);
    double f_multiplier = -2 * k_ * (theta - theta0_) / sin_t;
    Vec3d f1 = f_multiplier / r1 * (r2_vec / r2 - r1_vec * cos_t / r1);
    Vec3d f2 = f_multiplier / r2 * (r1_vec / r1 - r2_vec * cos_t / r2);
    double e = k_ * pow((theta - theta0_), 2);
    forces[atom_i1_] += f1;
    forces[atom_i2_] += f2;
    forces[shared_i_] -= f1 + f2;
    energies[atom_i1_] += e;
    energies[atom_i2_] += e;
    energies[shared_i_] += e;
}
