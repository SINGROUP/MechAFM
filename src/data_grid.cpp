#include "data_grid.hpp"

#include <complex>
#include <stdexcept>

using namespace std;


template<typename T>
DataGrid<T>::DataGrid(int nx, int ny, int nz, const T& value) {
    is_orthogonal_basis_ = false;
    basis_ = Mat3d(0);
    origin_ = Vec3d(0);
    initValues(nx, ny, nz, value);
}


template<typename T>
void DataGrid<T>::initValues(int nx, int ny, int nz, const T& value) {
    n_grid_.x = nx;
    n_grid_.y = ny;
    n_grid_.z = nz;
    values_.assign(nx*ny*nz, value);
}


template<typename T>
void DataGrid<T>::scaleValues(double scaling_factor) {
    for (auto it = values_.begin(); it != values_.end(); it++)
        *it *= scaling_factor;
}


template<typename T>
void DataGrid<T>::swapValues(vector<T>& values_to_swap) {
    values_.swap(values_to_swap);
}


template<typename T>
void DataGrid<T>::rotateCoordAxes(const string& new_coord_sequence) {
    Vec3i new_n_grid;
    Mat3d new_basis;
    Vec3d new_origin;
    int new_index;
    vector<T> new_values;
    new_values.resize(values_.size());
    
    if (new_coord_sequence == "ZXY") {
        new_n_grid.x = n_grid_.z;
        new_n_grid.y = n_grid_.x;
        new_n_grid.z = n_grid_.y;
        Mat3d swap_matrix = Mat3d(0);
        swap_matrix.at(2, 0) = 1;
        swap_matrix.at(0, 1) = 1;
        swap_matrix.at(1, 2) = 1;
        new_basis = swap_matrix.transpose().multiply(basis_.multiply(swap_matrix));
        new_origin.x = origin_.z;
        new_origin.y = origin_.x;
        new_origin.z = origin_.y;
        
        for (int ix = 0; ix < new_n_grid.x; ix++) {
            for (int iy = 0; iy < new_n_grid.y; iy++) {
                for (int iz = 0; iz < new_n_grid.z; iz++) {
                    new_index = ix*new_n_grid.y*new_n_grid.z + iy*new_n_grid.z + iz;
                    new_values[new_index] = values_[index(iy, iz, ix)];
                }
            }
        }
    }
    else if (new_coord_sequence == "YZX") {
        new_n_grid.x = n_grid_.y;
        new_n_grid.y = n_grid_.z;
        new_n_grid.z = n_grid_.x;
        Mat3d swap_matrix = Mat3d(0);
        swap_matrix.at(1, 0) = 1;
        swap_matrix.at(2, 1) = 1;
        swap_matrix.at(0, 2) = 1;
        new_basis = swap_matrix.transpose().multiply(basis_.multiply(swap_matrix));
        new_origin.x = origin_.y;
        new_origin.y = origin_.z;
        new_origin.z = origin_.x;
        
        for (int ix = 0; ix < new_n_grid.x; ix++) {
            for (int iy = 0; iy < new_n_grid.y; iy++) {
                for (int iz = 0; iz < new_n_grid.z; iz++) {
                    new_index = ix*new_n_grid.y*new_n_grid.z + iy*new_n_grid.z + iz;
                    new_values[new_index] = values_[index(iz, ix, iy)];
                }
            }
        }
    }
    else {
        throw "The coordinate sequence for rotateCoordAxes must be either ZXY or YZX.";
    }
    
    n_grid_ = new_n_grid;
    basis_ = new_basis;
    origin_ = new_origin;
    values_.swap(new_values);
}


template<typename T>
T& DataGrid<T>::atPBC(int ix, int iy, int iz) {
    int n;
    
    // Check if ix is outside the grid
    if ((ix < 0) || (ix >= n_grid_.x)) {
        n = floor(double(ix)/n_grid_.x);
        ix -= n*n_grid_.x;
    }
    
    // Check if iy is outside the grid
    if ((iy < 0) || (iy >= n_grid_.y)) {
        n = floor(double(iy)/n_grid_.y);
        iy -= n*n_grid_.y;
    }
    
    // Check if iz is outside the grid
    if ((iz < 0) || (iz >= n_grid_.z)) {
        n = floor(double(iz)/n_grid_.z);
        iz -= n*n_grid_.z;
    }
    
    return at(ix, iy, iz);
}


template<typename T>
Vec3d DataGrid<T>::positionAt(int ix, int iy, int iz) const {
    Vec3d position;
    if (is_orthogonal_basis_) {
        position.x = origin_.x + ix*basis_.at(0, 0);
        position.y = origin_.y + iy*basis_.at(1, 1); 
        position.z = origin_.z + iz*basis_.at(2, 2);
    } else {
        Vec3d grid_point = Vec3d(ix, iy, iz);
        position = origin_ + basis_.multiply(grid_point);
    }
    return position;
}


template<typename T>
Vec3d DataGrid<T>::getSpacing() const {
    Vec3d spacing;
    
    if (is_orthogonal_basis_) {
        spacing.x = basis_.at(0, 0);
        spacing.y = basis_.at(1, 1);
        spacing.z = basis_.at(2, 2);
    } else {
        throw runtime_error("The basis vectors are not orthogonal and thus the grid spacing is undefined.");
    }
    
    return spacing;
}


template<typename T>
void DataGrid<T>::setNGrid(int nx, int ny, int nz) {
    n_grid_.x = nx;
    n_grid_.y = ny;
    n_grid_.z = nz;
}


template<typename T>
void DataGrid<T>::setNGrid(const Vec3i& n_grid) {
    n_grid_.x = n_grid.x;
    n_grid_.y = n_grid.y;
    n_grid_.z = n_grid.z;
}


template<typename T>
void DataGrid<T>::setBasis(const vector<Vec3d>& basis_vectors) {
    for (int i = 0; i < 3; i++) {
        basis_.at(0, i) = basis_vectors[i].x;
        basis_.at(1, i) = basis_vectors[i].y;
        basis_.at(2, i) = basis_vectors[i].z;
    }
    
    is_orthogonal_basis_ = basis_.isDiagonal();
}


template<typename T>
void DataGrid<T>::setBasis(const Mat3d& basis_matrix) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            basis_.at(i, j) = basis_matrix.at(i, j);
        }
    }
    
    is_orthogonal_basis_ = basis_.isDiagonal();
}


template<typename T>
void DataGrid<T>::setSpacing(double dx, double dy, double dz) {
    is_orthogonal_basis_ = true;
    basis_.at(0, 0) = dx;
    basis_.at(1, 1) = dy;
    basis_.at(2, 2) = dz;
}


template<typename T>
void DataGrid<T>::setSpacing(const Vec3d& spacing) {
    is_orthogonal_basis_ = true;
    basis_.at(0, 0) = spacing.x;
    basis_.at(1, 1) = spacing.y;
    basis_.at(2, 2) = spacing.z;
}


template<typename T>
void DataGrid<T>::setOrigin(double x, double y, double z) {
    origin_.x = x;
    origin_.y = y;
    origin_.z = z;
}


template<typename T>
void DataGrid<T>::setOrigin(const Vec3d& origin) {
    origin_.x = origin.x;
    origin_.y = origin.y;
    origin_.z = origin.z;
}


// add extra templates here if needed
template class DataGrid<double>;
template class DataGrid<int>;
template class DataGrid<complex<double>>;
template class DataGrid<Vec3d>;
template class DataGrid<Vec3i>;
