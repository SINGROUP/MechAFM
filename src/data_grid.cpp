#include "data_grid.hpp"

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


template<typename T>
void DataGrid<T>::swapValues(vector<T>& values_to_swap) {
    values_.swap(values_to_swap);
}


// add extra templates here if needed
template class DataGrid<double>;
template class DataGrid<int>;
template class DataGrid<Vec3d>;
template class DataGrid<Vec3i>;
