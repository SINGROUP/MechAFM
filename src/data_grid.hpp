/*
 * data_grid.hpp
 * 
 * Generic container class for any data that is represented on a 3D grid.
 * 
 */

#pragma once

#include <string>
#include <vector>

#include "matrices.hpp"
#include "vectors.hpp"

using namespace std;

/** \brief Generic container class for any data that is represented on a 3D grid.
 * 
 * Works with types defined at the bottom of the cpp file. Add more if you need.
 * 
 */

template<typename T>
class DataGrid {

public:
    DataGrid(): is_orthogonal_basis_(false), n_grid_(Vec3i(0)), basis_(Mat3d(0)), origin_(Vec3d(0)) {};
    DataGrid(int nx, int ny, int nz, const T& value);
    ~DataGrid() {};
    
    void initValues(int nx, int ny, int nz, const T& value);
    void scaleValues(double scaling_factor);
    void swapValues(vector<T>& values_to_swap);
    void rotateCoordAxes(const string& new_coord_sequence);
    
    // Returns reference to value at given grid indices
    T& at(int ix, int iy, int iz) { return values_[index(ix, iy, iz)]; };
    // at() with periodic boundary conditions
    T& atPBC(int ix, int iy, int iz);
    // Returns reference to value at given internal storage index
    T& at(int ind) { return values_[ind]; };
    // Returns constant reference to value at given internal storage index
    const T& at(int ind) const { return values_[ind]; };
    // Returns position at given grid indices
    Vec3d positionAt(int ix, int iy, int iz) const;
    
    const Vec3i& getNGrid() const { return n_grid_; };
    const Mat3d& getBasis() const { return basis_; };
    Vec3d getSpacing() const;
    const Vec3d& getOrigin() const { return origin_; };
    
    void setNGrid(int nx, int ny, int nz);
    void setNGrid(const Vec3i& n_grid);
    void setBasis(const vector<Vec3d>& basis_vectors);
    void setBasis(const Mat3d& basis_matrix);
    void setSpacing(double dx, double dy, double dz);
    void setSpacing(const Vec3d& spacing);
    void setOrigin(double x, double y, double z);
    void setOrigin(const Vec3d& origin);
    
private:
    inline int index(int ix, int iy, int iz) { return ix*n_grid_.y*n_grid_.z + iy*n_grid_.z + iz; };

    bool is_orthogonal_basis_; // Determines whether the basis vectors of data grid are orthogonal
    Vec3i n_grid_; // The number of grid points along each basis vector
    Mat3d basis_; // The basis in which each point of the force grid is represented.
                                    // If is_orthogonal_coord_ == true, this is a 3x3 diagonal matrix
                                    // and the diagonal elements define the spacing between grid points.
    Vec3d origin_; // The real position of grid point (0, 0, 0)
    vector<T> values_; // The container for the actual data on the grid points
    
};
