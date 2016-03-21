#include "data_grid.hpp"

#include <complex>

using namespace std;


template<typename T>
DataGrid<T>::DataGrid(int nx, int ny, int nz, const T& value) {
    spacing_ = Vec3d(0);
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
    Vec3d new_spacing;
    Vec3d new_origin;
    int new_index;
    vector<T> new_values;
    new_values.resize(values_.size());
    
    if (new_coord_sequence == "ZXY") {
        new_n_grid.x = n_grid_.z;
        new_n_grid.y = n_grid_.x;
        new_n_grid.z = n_grid_.y;
        new_spacing.x = spacing_.z;
        new_spacing.y = spacing_.x;
        new_spacing.z = spacing_.y;
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
        new_spacing.x = spacing_.y;
        new_spacing.y = spacing_.z;
        new_spacing.z = spacing_.x;
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
    spacing_ = new_spacing;
    origin_ = new_origin;
    values_.swap(new_values);
}


template<typename T>
Vec3d DataGrid<T>::positionAt(int ix, int iy, int iz) const {
    Vec3d position;
    position.x = origin_.x + ix*spacing_.x;
    position.y = origin_.y + iy*spacing_.y; 
    position.z = origin_.z + iz*spacing_.z;
    return position;
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
void DataGrid<T>::setSpacing(double dx, double dy, double dz) {
    spacing_.x = dx;
    spacing_.y = dy;
    spacing_.z = dz;
}


template<typename T>
void DataGrid<T>::setSpacing(const Vec3d& spacing) {
    spacing_.x = spacing.x;
    spacing_.y = spacing.y;
    spacing_.z = spacing.z;
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
