#include "simulation.hpp"

#if MPI_BUILD
    #include <mpi.h>
#endif
#include <cmath>
#include <deque>
#include <memory>
#include <string>
#include <tuple>
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
    const double report_interval = 0.1;
    double next_report = report_interval;
    int processed_points = 0;
    const int total_points = n_points_.x * n_points_.y;

    const unsigned int buffer_size = options_.bufsize * n_points_.z;
    vector<OutputData> output_buffer;
    output_buffer.reserve(buffer_size);

#pragma omp parallel for shared(output_buffer, next_report, processed_points)
    for (int i = 0; i < n_points_.x; ++i) {
        double x = i * options_.dx;
        for (int j = 0; j < n_points_.y; ++j) {
            double y = j * options_.dy;

            int current_point = i * n_points_.y + j;
            processed_points++;
            // Check if this point is handled by this process
            if (current_point % n_processes_ != current_process_) {
                continue;
            } else {
                points_per_process_[current_process_]++;
            }

            System min_system = system;  // Take a copy for each z approach
            vector<OutputData> z_data(n_points_.z);
            if (min_system.interactions_ == nullptr) {
                error("System interactions are not given!");
            }
            min_system.setDummyXY(x, y);
            for (int k = 0; k < n_points_.z; ++k) {
                double z = options_.zhigh - k * options_.dz;
                min_system.setDummyZ(z);
                int n = 0;
                switch (options_.minimiser_type) {
                    case STEEPEST_DESCENT:
                        n = SDMinimisation(min_system, options_);
                        break;
                    case FIRE:
                        n = FIREMinimisation(min_system, options_);
                        break;
                    default:
                        error("Unimplemented minimiser type!");
                }
                // if (options_.flexible && current_point == total_points / 2) {
                if (current_point % 100 == 99) {
                    min_system.makeXYZFile();
                }
                z_data[k] = min_system.getOutput();
                z_data[k].indices = Vec3i(i, j, k);
                z_data[k].minimisation_steps = n;
                n_total_ += n;
            } // z
#pragma omp critical(output)
        {
            output_buffer.insert(output_buffer.end(), z_data.begin(), z_data.end());
            if (output_buffer.size() >= buffer_size) {
                writeOutput(output_buffer);
                output_buffer.clear();
            }

            // Report progress every once in a while
            double current_progress = 1.0f * processed_points / total_points;
            if(rootProcess() && current_progress >= next_report) {
                pretty_print("Finished %4.1f %% of the simulation",
                             100 * next_report);
                next_report += report_interval;
            }
        }
        } // y
    } // x
    // Write the remaining data
    writeOutput(output_buffer);
}

void Simulation::writeOutput(vector<OutputData> output_buffer) {
#if MPI_BUILD
    MPI_Status mpi_status;
    // Send the data to the root process
    if (!rootProcess()) {
        MPI_Send(static_cast<void*>(output_buffer.data()),
                 output_buffer.size() * sizeof(OutputData),
                 MPI_CHAR, root_process_, 0, universe);
    }
#endif
    // Receive the data from the daughter processes and write to file
    if (rootProcess()) {
        vector<OutputData> recieve_buffer(options_.bufsize * n_points_.z);
        for (int i = 0; i < n_processes_; ++i) {
            int data_size = 0;
            // On the main processor we only have to copy the data
            if (i == 0) {
                recieve_buffer = output_buffer;
                data_size = output_buffer.size();
            }
#if MPI_BUILD
            // For all other processors we need to receive the data
            else {
                MPI_Recv(static_cast<void*>(recieve_buffer.data()),
                         recieve_buffer.size() * sizeof(OutputData),
                         MPI_CHAR, i, 0, universe, &mpi_status);
                MPI_Get_count(&mpi_status, MPI_CHAR, &data_size);
                data_size /= sizeof(OutputData);
            }
#endif
            // Write data to file (only the root processor can do this)
            // PLEASE NOTE: DATA IS SENT IN STRIPED FORM, THEY ARE NOT ORDERED!
            for (int bi = 0; bi < data_size; ++bi) {
                int fi = recieve_buffer[bi].indices.z;
                // The file buffer can be a gzip pipe or an ASCII file stream
                fprintf(fstreams_[fi], "%d ", recieve_buffer[bi].indices.z);
                fprintf(fstreams_[fi], "%d ", recieve_buffer[bi].indices.x);
                fprintf(fstreams_[fi], "%d ", recieve_buffer[bi].indices.y);
                fprintf(fstreams_[fi], "%6.3f ", recieve_buffer[bi].position.x);
                fprintf(fstreams_[fi], "%6.3f ", recieve_buffer[bi].position.y);
                fprintf(fstreams_[fi], "%6.3f ", recieve_buffer[bi].position.z);
                fprintf(fstreams_[fi], "%8.4f ", recieve_buffer[bi].tip_force.x);
                fprintf(fstreams_[fi], "%8.4f ", recieve_buffer[bi].tip_force.y);
                fprintf(fstreams_[fi], "%8.4f ", recieve_buffer[bi].tip_force.z);
                fprintf(fstreams_[fi], "%6.3f ", recieve_buffer[bi].r_vec.x);
                fprintf(fstreams_[fi], "%6.3f ", recieve_buffer[bi].r_vec.y);
                fprintf(fstreams_[fi], "%6.3f ", recieve_buffer[bi].r_vec.z);
                fprintf(fstreams_[fi], "%6.3f ", recieve_buffer[bi].r);
                fprintf(fstreams_[fi], "%8.4f ", recieve_buffer[bi].angle);
                fprintf(fstreams_[fi], "%8.4f ", recieve_buffer[bi].tip_energy);
                fprintf(fstreams_[fi], "%d\n", recieve_buffer[bi].minimisation_steps);
            }
        }
    }
}

void Simulation::buildInteractions() {
    // Tip-Dummy distance needs to be calculated first since it's required
    // by some of the interactions.
    calculateTipDummyDistance();
    // Grid interactions have to be build first since it currently
    // clears the interaction list.
    if (options_.rigidgrid) {
        buildTipGridInteractions();
    } else {
        buildTipSurfaceInteractions();
    }
    buildTipDummyInteractions();
    if (options_.flexible) {
        buildSurfaceSurfaceInteractions();
        buildSubstrateInteractions();
    }
    // Give the system a pointer to the interaction list
    system.interactions_ = &interactions_;
}

void Simulation::calculateTipDummyDistance() {
    // Calculate the tip and dummy initial distance
    unordered_map<string, AtomParameters> ap = interaction_parameters_.atom_parameters;
    auto dummy = ap.find(system.types_[0])->second;
    auto tip = ap.find(system.types_[1])->second;
    double d = mixsig(dummy.sig, tip.sig) * SIXTHRT2;
    system.setTipDummyDistance(d);
}

bool Simulation::findOverwriteParameters(int atom_i1, int atom_i2, OverwriteParameters& op){
    unordered_multiset<string> test_set{system.types_[atom_i1], system.types_[atom_i2]};
    for (const auto& overwrite : interaction_parameters_.overwrite_parameters) {
        if (overwrite.atoms == test_set) {
            op = overwrite;
            return true;
        }
    }
    return false;
}

void Simulation::addVDWInteraction(int atom_i1, int atom_i2) {
    OverwriteParameters op;
    // Use overwrite parameters to define the interaction if they exist
    if (findOverwriteParameters(atom_i1, atom_i2, op)) {
        if (op.morse) {
            interactions_.emplace_back(new MorseInteraction(atom_i1, atom_i2, op.de, op.a, op.re));
        } else {
            double es6 = 4 * op.eps * pow(op.sig, 6);
            double es12 = 4 * op.eps * pow(op.sig, 12);
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
        double es6 = 4 * m_eps * pow(m_sig, 6);
        double es12 = 4 * m_eps * pow(m_sig, 12);
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
    // Only add the interaction if charges aren't zero
    if (q1 != 0 && q2 != 0) {
        double qq = interaction_parameters_.qbase * q1 * q2;
        interactions_.emplace_back(new CoulombInteraction(atom_i1, atom_i2, qq));
    }
}

void Simulation::buildTipSurfaceInteractions() {
    for (int i = 2; i < system.n_atoms_; ++i) {
        addVDWInteraction(1, i);
        if (options_.coulomb) {
            addCoulombInteraction(1, i);
        }
    }
}

void Simulation::buildTipGridInteractions() {
    // Check that the interaction list is empty before we begin
    if (!interactions_.empty()) {
        error("Interaction list is not empty when building force grid!.");
    }
    // Build the interactions we want to replace with the grid
    buildTipSurfaceInteractions();

    ForceGrid fg;
    const double min_spacing = 0.01;  // A grid smaller than this (in Angstrom) is silly
    double temp_spacing = options_.dx;
    if (options_.dy < temp_spacing) {
        temp_spacing = options_.dy;
    }
    if (options_.dz < temp_spacing) {
        temp_spacing = options_.dz;
    }
    if (temp_spacing < min_spacing) {
        temp_spacing = min_spacing;
    }
    fg.spacing_ = Vec3d(temp_spacing);
    const int border = ceil(fg.edge_ / min(fg.spacing_.x, min(fg.spacing_.y, fg.spacing_.z)));
    fg.grid_points_.x = floor(options_.box.x / fg.spacing_.x) + 2*border + 1;
    fg.grid_points_.y = floor(options_.box.y / fg.spacing_.y) + 2*border + 1;
    fg.grid_points_.z = floor((options_.zhigh - options_.zlow) / fg.spacing_.z) + 2*border + 1;
    int total_points = fg.grid_points_.x * fg.grid_points_.y * fg.grid_points_.z;
    fg.offset_.x = -border * fg.spacing_.x;
    fg.offset_.y = -border * fg.spacing_.y;
    fg.offset_.z = options_.zhigh - system.getTipDummyDistance()
                   - (fg.grid_points_.z - border - 1) * fg.spacing_.z;

    pretty_print("Computing 3D force grid (%d grid points)", total_points);
    pretty_print("");
    // Store the samples in a single list for easy communication
    // Sample order: f.x, f.y, f.z, e
    vector<double> temp_samples(4*total_points, 0);
#pragma omp parallel for default(none) shared(temp_samples, fg)
    for (int i = 0; i < fg.grid_points_.x; ++i) {
        double x = i * fg.spacing_.x + fg.offset_.x;
        for (int j = 0; j < fg.grid_points_.y; ++j) {
            double y = j * fg.spacing_.y + fg.offset_.y;

            // Check if this point is handled by this process
            int current_point = i * fg.grid_points_.y + j;
            if (current_point % n_processes_ != current_process_) {
                continue;
            }

            for (int k = 0; k < fg.grid_points_.z; ++k) {
                double z = k * fg.spacing_.z + fg.offset_.z;
                vector<Vec3d> positions = system.positions_;
                positions[1] = Vec3d(x, y, z);
                vector<Vec3d> forces(system.n_atoms_);
                vector<double> energies(system.n_atoms_, 0);
                for (const auto& interaction : interactions_) {
                    interaction->eval(positions, forces, energies);
                }
                int index = 4*(i * fg.grid_points_.y * fg.grid_points_.z + j * fg.grid_points_.z + k);
                temp_samples[index] = forces[1].x;
                temp_samples[index + 1] = forces[1].y;
                temp_samples[index + 2] = forces[1].z;
                temp_samples[index + 3] = energies[1];
            } // z
        } // y
    } // x

    // Communicate the data to all processes
#if MPI_BUILD
    MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(temp_samples.data()),
                  4 * total_points, MPI_DOUBLE, MPI_SUM, universe);
#endif

    // Initialise and set the sample vectors
    fg.forces_.assign(total_points, Vec3d(0));
    fg.energies_.assign(total_points, 0);
    for (int i = 0; i < total_points; ++i) {
        fg.forces_[i].x = temp_samples[4*i];
        fg.forces_[i].y = temp_samples[4*i + 1];
        fg.forces_[i].z = temp_samples[4*i + 2];
        fg.energies_[i] = temp_samples[4*i + 3];
    }

    // Replace the interactions with the grid
    interactions_.clear();
    interactions_.emplace_back(new GridInteraction(fg));
}

void Simulation::buildTipDummyInteractions() {
    // LJ / Morse
    addVDWInteraction(0, 1);

    // Harmonic constraint
    double k = interaction_parameters_.tip_dummy_k;
    double r0 = interaction_parameters_.tip_dummy_r0;
    interactions_.emplace_back(new Harmonic2DInteraction(0, 1, k, r0));
}

// Checks which atoms are connected by bonds to other atoms within given distance
vector<unordered_set<int>> getConnectedAtoms(vector<unordered_set<int>> adjacent_atoms, int distance) {
    vector<unordered_set<int>> connected_atoms(adjacent_atoms.size());
    // Perform a BFS for each of the atoms to see what it's connected to
    for (unsigned int i = 0; i < adjacent_atoms.size(); ++i) {
        deque<pair<int, int>> atoms_to_check;
        for (int atom_i : adjacent_atoms[i]) {
            atoms_to_check.emplace_back(atom_i, 1);
        }
        while (!atoms_to_check.empty()) {
            pair<int, int> check = atoms_to_check[0];
            connected_atoms[i].insert(check.first);
            atoms_to_check.pop_front();
            if (check.second < distance) {
                for (int atom_i : adjacent_atoms[check.first]) {
                    if (connected_atoms[i].count(atom_i) == 0) {
                        atoms_to_check.emplace_back(atom_i, check.second + 1);
                    }
                }
            }
        }
    }
    return connected_atoms;
}

void Simulation::buildSurfaceSurfaceInteractions() {
    // Set of atoms each atom is connected with by bonds
    vector<unordered_set<int>> adjacent_atoms(system.n_atoms_);

    // Harmonic bond interactions
    vector<pair<int, int>> bonds;  // List of all the bonds
    double bond_k = interaction_parameters_.bond_k;
    for (int i = 2; i < system.n_atoms_; ++i) {
        for (int j = i + 1; j < system.n_atoms_; ++j) {
            unordered_multiset<string> test_set{system.types_[i], system.types_[j]};
            double r0 = 0;
            // Check if there's a possibly a bond between these atoms
            for (const auto& pos_bond : interaction_parameters_.possible_bonds_) {
                if (pos_bond.atoms == test_set) {
                    r0 = pos_bond.r0;
                }
            }
            double atom_d = (system.positions_[i] - system.positions_[j]).len();
            // Check if the atoms distance is small enough to form a bond
            if (atom_d < 1.1 * r0) {
                bonds.emplace_back(i, j);
                adjacent_atoms[i].insert(j);
                adjacent_atoms[j].insert(i);
                interactions_.emplace_back(new HarmonicInteraction(i, j, bond_k, atom_d));
            }
        }
    }

    // Figure out which atoms are connected to each other with bonds
    auto connected_atoms = getConnectedAtoms(adjacent_atoms, system.n_atoms_);

    // Harmonic angle interactions
    vector<vector<int>> angles;  // List of all the angles
    double angle_k = interaction_parameters_.angle_k;
    for (unsigned int i = 0; i < bonds.size(); ++i) {
        for (unsigned int j = i + 1; j < bonds.size(); ++j) {
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
                angles.push_back(vector<int>{shared_atom, atom1, atom2});
                Vec3d d1 = system.positions_[atom1] - system.positions_[shared_atom];
                Vec3d d2 = system.positions_[atom2] - system.positions_[shared_atom];
                double cos_t = d1.normalized().dot(d2.normalized());
                double theta0 = acos(cos_t);
                interactions_.emplace_back(new HarmonicAngleInteraction(
                                    shared_atom, atom1, atom2, angle_k, theta0));
            }
        }
    }

    // Harmonic dihedral interactions
    double dihedral_k = interaction_parameters_.dihedral_k;
    unordered_set<int> improper_centers;  // Set of improper centers we've allready added
    for (unsigned int i = 0; i < angles.size(); ++i) {
        for (unsigned int j = i + 1; j < angles.size(); ++j) {
            auto angle1 = angles[i];
            auto angle2 = angles[j];
            int atom1, atom2, atom3, atom4;
            atom1 = atom2 = atom3 = atom4 = -1;
            if (angle1[0] == angle2[0] && improper_centers.count(angle1[0]) == 0) {
                improper_centers.insert(angle1[0]);
                atom1 = angle1[0];
                atom2 = angle1[1];
                atom3 = angle1[2];
                if (angle2[1] == angle1[1] || angle2[1] == angle1[2]) {
                    atom4 = angle2[2];
                } else {
                    atom4 = angle2[1];
                }
            } else if (angle1[1] == angle2[0] && angle2[1] == angle1[0]) {
                atom1 = angle1[2];
                atom2 = angle1[0];
                atom3 = angle2[0];
                atom4 = angle2[2];
            } else if (angle1[2] == angle2[0] && angle2[1] == angle1[0]) {
                atom1 = angle1[1];
                atom2 = angle1[0];
                atom3 = angle2[0];
                atom4 = angle2[2];
            } else if (angle1[1] == angle2[0] && angle2[2] == angle1[0]) {
                atom1 = angle1[2];
                atom2 = angle1[0];
                atom3 = angle2[0];
                atom4 = angle2[1];
            } else if (angle1[2] == angle2[0] && angle2[2] == angle1[0]) {
                atom1 = angle1[1];
                atom2 = angle1[0];
                atom3 = angle2[0];
                atom4 = angle2[1];
            }
            if (atom1 != -1) {
                Vec3d d12 = system.positions_[atom1] - system.positions_[atom2];
                Vec3d d23 = system.positions_[atom2] - system.positions_[atom3];
                Vec3d d43 = system.positions_[atom4] - system.positions_[atom3];
                Vec3d m = d12.cross(d23);
                Vec3d n = d43.cross(d23);
                double tan_s = n.dot(d12) * d23.len() / (m.dot(n));
                double sigma0 = atan(tan_s);
                interactions_.emplace_back(new HarmonicDihedralInteraction(
                                atom1, atom2, atom3, atom4, dihedral_k, sigma0));
            }
        }
    }

    // Non-bonded interactions
    unordered_map<string, AtomParameters> ap = interaction_parameters_.atom_parameters;
    for (int i = 2; i < system.n_atoms_; ++i) {
        for (int j = i + 1; j < system.n_atoms_; ++j) {
            // Only add non-bonded interactions if atoms aren't connected by bonds
            if (connected_atoms[i].count(j) == 0) {
                addVDWInteraction(i, j);
                if (options_.coulomb) {
                    addCoulombInteraction(i, j);
                }
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
