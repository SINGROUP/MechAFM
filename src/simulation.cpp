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

//debug only
#include <iostream>

#include "cube_io.hpp"
#include "fft.hpp"
#include "interactions.hpp"
#include "messages.hpp"
#include "matrices.hpp"
#include "vectors.hpp"

using namespace std;

bool Simulation::rootProcess() {
    return current_process_ == root_process_;
}

void Simulation::initialize() {
    // Rotate system if z axis is not perpendicular to the surface
    if (options_.normal == NORMAL_X) {
        system.rotateCoordAxes("YZX");
    } else if (options_.normal == NORMAL_Y) {
        system.rotateCoordAxes("ZXY");
    }
    
    // If pbc for vdw is on, get the unit cell vectors from electrostatic potential file
    // or use the ones given in input file
    if (options_.vdw_pbc) {
        if (options_.use_external_potential) {
            CubeReader cube_file(options_.e_potential_file);
            Vec3i n_voxels = cube_file.getNVoxels();
            vector<Vec3d> voxel_vectors = cube_file.getVoxelVectors();
            vector<Vec3d> cell_vectors;
            cell_vectors.push_back(n_voxels.x * voxel_vectors[0]);
            cell_vectors.push_back(n_voxels.y * voxel_vectors[1]);
            cell_vectors.push_back(n_voxels.z * voxel_vectors[2]);
            system.setUnitCell(cell_vectors);
        }
        else {
            vector<Vec3d> cell_vectors;
            cell_vectors.push_back(options_.cell_a);
            cell_vectors.push_back(options_.cell_b);
            cell_vectors.push_back(options_.cell_c);
            system.setUnitCell(cell_vectors);
        }
    }
    
    system.centerMolecule(options_.center);
    system.setMoleculeZ();
    calculateTipDummyDistance();
    buildInteractions();
}

void Simulation::run() {
    const double report_interval = 0.1;
    double next_report = report_interval;
    int processed_points = 0;
    const int total_points = n_points_.x * n_points_.y;

    const unsigned int buffer_size = options_.bufsize * n_points_.z;
    vector<OutputData> output_buffer;
    output_buffer.reserve(buffer_size);
    pretty_print("Starting simulation");

#pragma omp parallel for schedule(dynamic, 1)
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
            if (options_.rigidgrid)
                min_system.setTipPbc(false);
            min_system.setDummyXY(x, y);
            min_system.setDummyZ(options_.zhigh);
            for (int k = 0; k < n_points_.z; ++k) {
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
                if (options_.flexible && current_point == total_points / 2) {
                    min_system.makeXYZFile(options_.outputfolder);
                }
                z_data[k] = min_system.getOutput();
                z_data[k].indices = Vec3i(i, j, k);
                z_data[k].minimisation_steps = n;
                n_total_ += n;
                min_system.lowerTip(options_.dz);
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
            while (rootProcess() && current_progress >= next_report) {
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

void Simulation::addVDWInteraction(int atom_i1, int atom_i2, Vec3d pbc_shift = Vec3d(0)) {
    OverwriteParameters op;
    // Use overwrite parameters to define the interaction if they exist
    if (findOverwriteParameters(atom_i1, atom_i2, op)) {
        if (op.morse) {
            interactions_.emplace_back(new MorseInteraction(atom_i1, atom_i2, op.de, op.a, op.re, pbc_shift));
        } else {
            double es6 = 4 * op.eps * pow(op.sig, 6);
            double es12 = 4 * op.eps * pow(op.sig, 12);
            interactions_.emplace_back(new LJInteraction(atom_i1, atom_i2, es6, es12, pbc_shift));
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
        interactions_.emplace_back(new LJInteraction(atom_i1, atom_i2, es6, es12, pbc_shift));
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
    if (options_.vdw_pbc) {
        Vec3d pbc_shift;
        Mat3d cell_matrix = system.getUnitCell();
        // This is hard coded to consider atoms across one unit cell boundary in the direction
        // of each unit cell vector perpendicular to the surface normal
        for (int cell_a_shift = -1; cell_a_shift <= 1; cell_a_shift++) {
            for (int cell_b_shift = -1; cell_b_shift <= 1; cell_b_shift++) {
                pbc_shift = cell_a_shift*cell_matrix.getColumn(0) + \
                            cell_b_shift*cell_matrix.getColumn(1);
                for (int i = 2; i < system.n_atoms_; ++i) {
                    addVDWInteraction(1, i, pbc_shift);
                }
            }
        }
    }
    else {
        for (int i = 2; i < system.n_atoms_; ++i) {
            addVDWInteraction(1, i);
            if (options_.coulomb) {
                addCoulombInteraction(1, i);
            } 
        }
    }
    
    // Interaction of tip atom with an external electrostatic potential
    if (options_.use_external_potential) {
        pretty_print("Calculating energy and force on grid from external electrostatic potential.");
        DataGrid<double> electrostatic_potential;
        CubeReader cube_file(options_.e_potential_file);
        if (rootProcess()) {
            cube_file.storeToDataGrid(electrostatic_potential);
        } else {
            const Vec3i& n_grid = cube_file.getNVoxels();
            const vector<Vec3d> voxel_vectors = cube_file.getVoxelVectors();
            const Vec3d& origin = cube_file.getOrigin();
            electrostatic_potential.initValues(n_grid.x, n_grid.y, n_grid.z, 0.0);
            electrostatic_potential.setBasis(voxel_vectors);
            electrostatic_potential.setOrigin(origin);
        }
        
#if MPI_BUILD
        // Broadcast electrostatic potential values from root process to others
        vector<double> tmp_values;
        electrostatic_potential.swapValues(tmp_values);
        MPI_Bcast(tmp_values.data(), tmp_values.size(), MPI_DOUBLE, root_process_, universe);
        electrostatic_potential.swapValues(tmp_values);
#endif
        
        // Units for potential are in Hartree units in the case of CP2k cube files
        // Hartree potential is defined for negatively charge electrons -> multiply by -1
        if (options_.units == U_EV)
            electrostatic_potential.scaleValues(-g_hartree_to_eV);
        else if (options_.units == U_KJ)
            electrostatic_potential.scaleValues(-g_hartree_to_kJ);
        else if (options_.units == U_KCAL)
            electrostatic_potential.scaleValues(-g_hartree_to_kcal);
        else
            error("unit conversion of Hartree potential to given units is not implemented");
        
        // Rotate potential grid if z coordinate is not perpendicular to the surface
        if (options_.normal == NORMAL_X) {
            electrostatic_potential.rotateCoordAxes("ZYX");
        } else if (options_.normal == NORMAL_Y) {
            electrostatic_potential.rotateCoordAxes("ZXY");
        }
        
        // Offset the potential data by the same amount as the atomic system
        electrostatic_potential.setOrigin(electrostatic_potential.getOrigin() + system.getOffset());
        
        if (DEBUG_MODE) {
            cout << endl <<"========== Debug ==========" << endl << endl;
            
            cout << "Hello from process " << current_process_ << endl;
            cout << "Number of grid points:" << endl;
            cout << electrostatic_potential.getNGrid().x << " " << electrostatic_potential.getNGrid().y << " " << electrostatic_potential.getNGrid().z << endl;
            cout << "Basis matrix:" << endl;
            electrostatic_potential.getBasis().print();
            cout << "Three first values of electrostatic potential:" << endl;
            cout << electrostatic_potential.at(0,0,0) << " " << electrostatic_potential.at(1,0,0) << " " << electrostatic_potential.at(2,0,0) << endl << endl;

            // Test FFT
            DataGrid<dcomplex> e_potential_kspace;
            chrono::high_resolution_clock::time_point start, end;
            chrono::duration<double> time_interval;
            
            start = chrono::high_resolution_clock::now();
            fft_data_grid(electrostatic_potential, e_potential_kspace);
            end = chrono::high_resolution_clock::now();
            time_interval = end - start;
            cout << "FFT took " << time_interval.count() << " s" << endl;
            
            cout << "Hello from process " << current_process_ << " after FFT" << endl;
            cout << "Number of grid points:" << endl;
            cout << e_potential_kspace.getNGrid().x << " " << e_potential_kspace.getNGrid().y << " " << e_potential_kspace.getNGrid().z << endl;
            cout << "Basis matrix:" << endl;
            e_potential_kspace.getBasis().print();
            cout << "Three first values of electrostatic potential:" << endl;
            cout << e_potential_kspace.at(0,0,0) << ", " << e_potential_kspace.at(1,0,0) << ", " << e_potential_kspace.at(2,0,0) << endl << endl;
            
            start = chrono::high_resolution_clock::now();
            ffti_data_grid(e_potential_kspace, electrostatic_potential);
            end = chrono::high_resolution_clock::now();
            time_interval = end - start;
            cout << "Inverse FFT took " << time_interval.count() << " s" << endl;
            
            cout << "Hello from process " << current_process_ << " after FFT and inverse FFT" << endl;
            cout << "Number of grid points:" << endl;
            cout << electrostatic_potential.getNGrid().x << " " << electrostatic_potential.getNGrid().y << " " << electrostatic_potential.getNGrid().z << endl;
            cout << "Basis matrix:" << endl;
            electrostatic_potential.getBasis().print();
            cout << "Three first values of electrostatic potential:" << endl;
            cout << electrostatic_potential.at(0,0,0) << " " << electrostatic_potential.at(1,0,0) << " " << electrostatic_potential.at(2,0,0) << endl << endl;
            
            cout << "========== End debug ==========" << endl << endl;
        }
        
        // Create the interaction between the tip atom and the electrostatic potential
        interactions_.emplace_back(new ElectrostaticPotentialInteraction(electrostatic_potential, system.charges_[1], g_tip_gaussian_width));
        pretty_print("Done!");
    }
}

void Simulation::buildTipGridInteractions() {
    // Check that the interaction list is empty before we begin
    if (!interactions_.empty()) {
        error("Interaction list is not empty when building force grid!.");
    }
    // Build the interactions we want to replace with the grid
    buildTipSurfaceInteractions();

    Vec3d spacing = Vec3d(options_.dx, options_.dy, options_.dz) / 2;
    Vec3i border;
    border.x = ceil(g_force_grid_margin / spacing.x);
    border.y = ceil(g_force_grid_margin / spacing.y);
    border.z = ceil(g_force_grid_margin / spacing.z);
    Vec3i n_grid;
    n_grid.x = floor(options_.area.x / spacing.x) + 2*border.x + 1;
    n_grid.y = floor(options_.area.y / spacing.y) + 2*border.y + 1;
    n_grid.z = floor((options_.zhigh - options_.zlow) / spacing.z) + 2*border.z + 1;
    int total_points = n_grid.x * n_grid.y * n_grid.z;
    Vec3d offset;
    offset.x = -border.x * spacing.x;
    offset.y = -border.y * spacing.y;
    offset.z = options_.zlow - system.getTipDummyDistance() - border.z * spacing.z;

    pretty_print("Computing 3D force grid: %d, %d, %d (%d grid points)",
        n_grid.x, n_grid.y, n_grid.z, total_points);
    // Initialize temporary sample vectors
    vector<Vec3d> forces;
    vector<double> energies;
    forces.assign(total_points, Vec3d(0));
    energies.assign(total_points, 0);
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < n_grid.x; ++i) {
        double x = i * spacing.x + offset.x;
        for (int j = 0; j < n_grid.y; ++j) {
            double y = j * spacing.y + offset.y;

            // Check if this point is handled by this process
            int current_point = i * n_grid.y + j;
            if (current_point % n_processes_ != current_process_) {
                continue;
            }
            
            System temp_system = system; // Create a copy of system for each (x, y) point
            temp_system.setTipDummyDistance(0);
            temp_system.setDummyXY(x, y);

            for (int k = 0; k < n_grid.z; ++k) {
                double z = k * spacing.z + offset.z;
                temp_system.setDummyZ(z);
                fill(temp_system.forces_.begin(), temp_system.forces_.end(), Vec3d(0));
                fill(temp_system.energies_.begin(), temp_system.energies_.end(), 0);
                for (const auto& interaction : interactions_) {
                    interaction->eval(temp_system.positions_, temp_system.forces_, temp_system.energies_);
                }
                int index = i * n_grid.y * n_grid.z + j * n_grid.z + k;
                forces[index] = temp_system.forces_[1];
                energies[index] = temp_system.energies_[1];
            } // z
        } // y
    } // x

    // Communicate the data to all processes
#if MPI_BUILD
    MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(forces_.data()),
                  total_points * sizeof(Vec3d), MPI_CHAR, MPI_SUM, universe);
    MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(energies_.data()),
                  total_points, MPI_DOUBLE, MPI_SUM, universe);
#endif

    // Initialize ForceGrid object
    ForceGrid fg;
    fg.setNGrid(n_grid);
    fg.setSpacing(spacing);
    fg.setOffset(offset);
    fg.swapForceValues(forces);
    fg.swapEnergyValues(energies);

    // Replace the interactions with the grid
    interactions_.clear();
    interactions_.emplace_back(new GridInteraction(fg));
    pretty_print("Done!");
}

void Simulation::buildTipDummyInteractions() {
    // LJ / Morse
    addVDWInteraction(0, 1);
    
    // Coulomb
    if (options_.tip_dummy_coulomb)
        addCoulombInteraction(0, 1);

    // Harmonic constraint
    double k = interaction_parameters_.tip_dummy_k;
    double r0 = interaction_parameters_.tip_dummy_r0;
    interactions_.emplace_back(new TipHarmonicInteraction(0, 1, k, r0));
}

void Simulation::buildSubstrateInteractions() {
    double eps = interaction_parameters_.substrate_eps;
    double sig = interaction_parameters_.substrate_sig;
    double lambda = interaction_parameters_.substrate_lambda;
    double rc = interaction_parameters_.substrate_rc;
    double k = interaction_parameters_.substrate_k;
    for (int i = 2; i < system.n_atoms_; ++i) {
      interactions_.emplace_back(new SubstrateInteraction(i, eps, sig, rc, lambda));
        if (system.fixed_[i] == 2) {
            interactions_.emplace_back(new XYHarmonicInteraction(i, k,
                                            system.positions_[i].getXY()));
        }
    }
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
