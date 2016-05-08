#pragma once

#include <memory>
#include <vector>

#include "data_grid.hpp"
#include "matrices.hpp"
#include "vectors.hpp"

using namespace std;

struct InputOptions;
class Interaction;

class ForceGrid {
 public:
    ForceGrid();
    ~ForceGrid() {};
    
    void setPeriodic(bool is_periodic) { is_periodic_ = is_periodic; };
    void setNGrid(const Vec3i& n_grid);
    void setBasis(const vector<Vec3d>& basis_vectors);
    void setBasis(const Mat3d& basis_matrix);
    void setSpacing(const Vec3d& spacing);
    void setOffset(const Vec3d& offset);
    
    void swapForceValues(vector<Vec3d>& forces_to_swap);
    void swapForceValues(DataGrid<Vec3d>& forces_to_swap);
    void swapEnergyValues(vector<double>& energies_to_swap);
    void swapEnergyValues(DataGrid<double>& energies_to_swap);
    
    // Calculates the interpolated force and energy at the given position
    void interpolate(const Vec3d& position, Vec3d& force, double& energy) const;

 private:
    // Returns the grid point matching for the given position
    Vec3i getGridPoint(const Vec3d& pos) const;
    // Returns the matching array index for the given grid point
    int getGridPointIndex(const Vec3i& grid_point) const;
    Vec3i pbcGridPoint(const Vec3i& grid_point) const;
    
    bool is_periodic_; // Determines whether the force grid has periodic boundary conditions
    bool is_orthogonal_basis_; // Determines whether the basis vectors of force grid are orthogonal
    Vec3i n_grid_;  // The number of grid points along each basis vector
    Mat3d basis_;   // The basis in which each point of the force grid is represented.
                                    // If is_orthogonal_coord_ == true, this is a 3x3 diagonal matrix
                                    // and the diagonal elements define the spacing between grid points.
    Vec3d offset_;  // The real position of grid point (0, 0, 0)
    vector<Vec3d> forces_;  // List of force samples
    vector<double> energies_;  // List of energy samples
};
