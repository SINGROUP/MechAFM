#include "simulation.hpp"

#include <cmath>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "interactions.hpp"
#include "messages.hpp"
#include "vectors.hpp"

using namespace std;

bool Simulation::rootProcess() {
    return current_process_ == root_process_;
}

void Simulation::run() {
    double x, y, z;

    switch (options_.minimiser_type) {
        case STEEPEST_DESCENT:
            minimiser = unique_ptr<Minimiser>(new SDMinimiser());
            break;
        case FIRE:
            minimiser = unique_ptr<Minimiser>(new FIREMinimiser());
            break;
        default:
            error(*this, "Invalid minimiser type.");
    }
    minimiser->initialize(*this);

    for (int i = 0; i < n_points_.x; ++i) {
        x = i * options_.dx;
        for (int j = 0; j < n_points_.y; ++j) {
            y = j * options_.dy;
            System min_system = system;  // Take a copy for each z approach
            min_system.setDummyXY(x, y);
            for (int k = 0; k < n_points_.z; ++k) {
                z = k * options_.dz;
                min_system.setDummyZ(z);
                minimiser->minimise(min_system);
            }
        }
    }
}

void Simulation::buildInteractions() {
    buildTipDummyInteractions();
    if (options_.rigidgrid) {
        buildTipGridInteractions();
    } else {
        buildTipSurfaceInteractions();
    }
    if (options_.flexible) {
        buildSurfaceSurfaceInteractions();
        buildSubstrateInteractions();
        buildBondInteractions();
    }
}

bool Simulation::findOverwriteParameters(int atom_i1, int atom_i2, OverwriteParameters op){
    unordered_multiset<string> test_set{system.types_[atom_i1], system.types_[atom_i2]};
    for (const auto& overwrite : interaction_parameters_.overwrite_parameters) {
        if (overwrite.atoms == test_set) {
            op = overwrite;
            return true;
        }
    }
    return false;
}

void Simulation::addLJInteraction(int atom_i1, int atom_i2) {
    OverwriteParameters op;
    if (findOverwriteParameters(atom_i1, atom_i2, op)) {
        if (op.morse) {
            interactions_.emplace_back(new MorseInteraction(atom_i1, atom_i2, op.de, op.a, op.re));
        } else {
            double es6 = op.eps * pow(op.sig, 6);
            double es12 = op.eps * pow(op.sig, 12);
            interactions_.emplace_back(new LJInteraction(atom_i1, atom_i2, es6, es12));
        }
    } else {
        unordered_map<string, AtomParameters> ap = interaction_parameters_.atom_parameters;
        auto atom1_it = ap.find(system.types_[atom_i1]);
        auto atom2_it = ap.find(system.types_[atom_i2]);
        if (atom1_it == ap.end()) {
            error("Parameters for atom type %s not found in parameter file!",
                            system.types_[atom_i1].c_str());
        }
        if (atom2_it == ap.end()) {
            error("Parameters for atom type %s not found in parameter file!",
                            system.types_[atom_i2].c_str());
        }
        double m_eps = mixeps(atom1_it->second.eps, atom2_it->second.eps);
        double m_sig = mixsig(atom1_it->second.sig, atom2_it->second.sig);
        double es6 = m_eps * pow(m_sig, 6);
        double es12 = m_eps * pow(m_sig, 12);
        interactions_.emplace_back(new LJInteraction(atom_i1, atom_i2, es6, es12));
    }
}

void Simulation::addCoulombInteraction(int atom_i1, int atom_i2) {
    unordered_map<string, AtomParameters> ap = interaction_parameters_.atom_parameters;

    double q1, q2;
    if (options_.xyz_charges) {
        q1 = system.charges_[atom_i1];
        q2 = system.charges_[atom_i2];
    } else {
        auto atom1_it = ap.find(system.types_[atom_i1]);
        auto atom2_it = ap.find(system.types_[atom_i2]);
        if (atom1_it == ap.end()) {
            error("Parameters for atom type %s not found in parameter file!",
                            system.types_[atom_i1].c_str());
        }
        if (atom2_it == ap.end()) {
            error("Parameters for atom type %s not found in parameter file!",
                            system.types_[atom_i2].c_str());
        }
        q1 = atom1_it->second.q;
        q2 = atom2_it->second.q;
    }
    if (q1 != 0 && q2 != 0) {
        double qq = interaction_parameters_.qbase * q1 * q2;
        interactions_.emplace_back(new CoulombInteraction(atom_i1, atom_i2, qq));
    }
}

void Simulation::buildTipSurfaceInteractions() {
    for (int i = 2; i < system.n_atoms_; ++i) {
        addLJInteraction(1, i);
        if (options_.coulomb) {
            addCoulombInteraction(1, i);
        }
    }
}

void Simulation::buildTipGridInteractions() {
    warning("Rigid grid not yet implemented! Using regular tip-surface interactions.");
    buildTipSurfaceInteractions();
}

void Simulation::buildTipDummyInteractions() {
    // LJ / Morse
    addLJInteraction(0, 1);

    // Harmonic
    double k = interaction_parameters_.tip_dummy_k;
    double r0 = interaction_parameters_.tip_dummy_r0;
    interactions_.emplace_back(new Harmonic2DInteraction(0, 1, k, r0));

    // Calculate the tip and dummy initial distance
    unordered_map<string, AtomParameters> ap = interaction_parameters_.atom_parameters;
    auto dummy = ap.find(system.types_[0])->second;
    auto tip = ap.find(system.types_[1])->second;
    double d = mixsig(dummy.sig, tip.sig) * SIXTHRT2;
    system.setTipDummyDistance(d);
}

void Simulation::buildSurfaceSurfaceInteractions() {
    unordered_map<string, AtomParameters> ap = interaction_parameters_.atom_parameters;
    for (int i = 2; i < system.n_atoms_; ++i) {
        for (int j = i + 1; j < system.n_atoms_; ++j) {
            addLJInteraction(i, j);
            if (options_.coulomb) {
                addCoulombInteraction(i, j);
            }
        }
    }
}

void Simulation::buildSubstrateInteractions() {
    double k = interaction_parameters_.substrate_k;
    for (int i = 2; i < system.n_atoms_; ++i) {
        interactions_.emplace_back(new SubstrateInteraction(i, k, system.positions_[i].z));
    }
}

void Simulation::buildBondInteractions() {
    vector<pair<int, int>> bonds;
    double bond_k = interaction_parameters_.bond_k;
    for (int i = 2; i < system.n_atoms_; ++i) {
        for (int j = i + 1; j < system.n_atoms_; ++j) {
            unordered_multiset<string> test_set{system.types_[i], system.types_[j]};
            double r0 = 0;
            for (const auto& pos_bond : interaction_parameters_.possible_bonds_) {
                if (pos_bond.atoms == test_set) {
                    r0 = pos_bond.r0;
                }
            }
            double atom_d = (system.positions_[i] - system.positions_[j]).len();
            if (atom_d < 1.1 * r0) {
                bonds.emplace_back(i, j);
                interactions_.emplace_back(new HarmonicInteraction(i, j, bond_k, atom_d));
            }
        }
    }

    double angle_k = interaction_parameters_.angle_k;
    for (int i = 0; i < bonds.size(); ++i) {
        for (int j = i + 1; j < bonds.size(); ++j) {
            auto bond1 = bonds[i];
            auto bond2 = bonds[j];
            int shared_atom = -1;
            int atom1, atom2;
            if (bond1.first == bond2.first) {
                shared_atom = bond1.first;
                atom1 = bond1.second;
                atom2 = bond2.second;
            } else if (bond1.first == bond2.second) {
                shared_atom = bond1.first;
                atom1 = bond1.second;
                atom2 = bond2.first;
            } else if (bond1.second == bond2.first) {
                shared_atom = bond1.second;
                atom1 = bond1.first;
                atom2 = bond2.second;
            } else if (bond1.second == bond2.second) {
                shared_atom = bond1.second;
                atom1 = bond1.first;
                atom2 = bond2.first;
            }
            if (shared_atom != -1) {
                Vec3d d1 = system.positions_[atom1] - system.positions_[shared_atom];
                Vec3d d2 = system.positions_[atom2] - system.positions_[shared_atom];
                double cos_t = d1.normalized().dot(d2.normalized());
                double theta0 = acos(cos_t);
                interactions_.emplace_back(new HarmonicAngleInteraction(
                                    shared_atom, atom1, atom2, angle_k, theta0));
            }
        }
    }
}
