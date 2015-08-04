#include "interactions.hpp"

#include <cmath>

#include "force_grid.hpp"
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

void GridInteraction::eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, vector<double>& energies) const {
    Vec3d tip_force;
    double tip_energy;
    force_grid.interpolate(positions[1], tip_force, tip_energy);
    forces[1] += tip_force;
    energies[1] += tip_energy;
}

void TipHarmonicInteraction::eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, vector<double>& energies) const {
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

void XYHarmonicInteraction::eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, vector<double>& energies) const {
    Vec2d r_2d = positions[atom_i_].getXY() - p0_;
    double r = r_2d.len();
    double e = k_ * r * r;
    Vec3d f;
    if (r > TOLERANCE) {
        f = Vec3d((-2 * k_ * r_2d), 0);
    }
    energies[atom_i_] += e;
    forces[atom_i_] += f;
}

void SubstrateInteraction::eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, vector<double>& energies) const {
    double dz = positions[atom_i_].z - z0_;
    double sig_z = sig_ / dz;
    double sig_z2 = sig_z * sig_z;
    double sig_z4 = sig_z2 * sig_z2;
    double sig_z10 = sig_z4 * sig_z4 * sig_z2;
    forces[atom_i_].z += 4 * multiplier_ / dz * (sig_z10 - sig_z4);
    energies[atom_i_] += multiplier_ * (2/5 * sig_z10 - sig_z4);
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
    if (cos_t > 1) {
        cos_t = 1;
    } else if (cos_t < -1) {
        cos_t = -1;
    }
    double theta = acos(cos_t);
    double sin_t = sin(theta);
    double d_theta = theta - theta0_;
    double f_multiplier = -2 * k_ * d_theta / sin_t;
    Vec3d f1 = -f_multiplier / r1 * (r2_vec / r2 - r1_vec * cos_t / r1);
    Vec3d f2 = -f_multiplier / r2 * (r1_vec / r1 - r2_vec * cos_t / r2);
    double e = k_ * pow(d_theta, 2);
    forces[atom_i1_] += f1;
    forces[atom_i2_] += f2;
    forces[shared_i_] -= f1 + f2;
    energies[atom_i1_] += e;
    energies[atom_i2_] += e;
    energies[shared_i_] += e;
}

void HarmonicDihedralInteraction::eval(const vector<Vec3d>& positions, vector<Vec3d>& forces, vector<double>& energies) const {
    Vec3d r12_vec = positions[atom_i1_] - positions[atom_i2_];
    Vec3d r23_vec = positions[atom_i2_] - positions[atom_i3_];
    Vec3d r43_vec = positions[atom_i4_] - positions[atom_i3_];
    double r23 = r23_vec.len();
    Vec3d m_vec = r12_vec.cross(r23_vec);
    Vec3d n_vec = r43_vec.cross(r23_vec);
    double m = m_vec.len();
    double n = n_vec.len();
    double tan_s = n_vec.dot(r12_vec) * r23 / m_vec.dot(n_vec);
    double sigma = atan(tan_s);
    double d_sigma = sigma - sigma0_;
    double f_multiplier = 2 * k_ * d_sigma;
    Vec3d f1 = -f_multiplier * r23 / (m * m) * m_vec;
    Vec3d f2 = -f_multiplier * (r43_vec.dot(r23_vec) / (n * n * r23) * n_vec
                - ((r23*r23) + r12_vec.dot(r23_vec)) / (m * m * r23) * m_vec);
    Vec3d f3 = -f_multiplier * (r12_vec.dot(r23_vec) / (m * m * r23) * m_vec
                + ((r23*r23) - r43_vec.dot(r23_vec)) / (n * n * r23) * n_vec);
    Vec3d f4 = f_multiplier * r23 / (n * n) * n_vec;
    double e = k_ * pow(d_sigma, 2);
    forces[atom_i1_] += f1;
    forces[atom_i2_] += f2;
    forces[atom_i3_] += f3;
    forces[atom_i4_] += f4;
    energies[atom_i1_] += e;
    energies[atom_i2_] += e;
    energies[atom_i3_] += e;
    energies[atom_i4_] += e;
}
