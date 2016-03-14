/*
 * data_grid.hpp
 * 
 * Copyright 2016 Juha Ritala <juha.ritala@aalto.fi>
 * 
 * Generic container class for any data that is represented on a 3D grid.
 * 
 */

#pragma once

#include <vector>

#include "vectors.hpp"

using namespace std;

/** \brief Generic container class for any data that is represented on a 3D grid.
 * 
 */

template<typename T>
class DataGrid {

public:
    DataGrid(): n_grid_(Vec3i(0)), spacing_(Vec3d(0)), origin_(Vec3d(0)) {};
    DataGrid(int nx, int ny, int nz, const T& value);
    ~DataGrid() {};
    
    T& at(int ix, int iy, int iz) { return values_[index(ix, iy, iz)]; };
    void initValues(int nx, int ny, int nz, const T& value);
    
    const Vec3i& getNGrid() const { return n_grid_; };
    const Vec3d& getSpacing() const { return spacing_; };
    const Vec3d& getOrigin() const { return origin_; };
    
    void setNGrid(int nx, int ny, int nz);
    void setNGrid(const Vec3i& n_grid);
    void setSpacing(double dx, double dy, double dz);
    void setSpacing(const Vec3d& spacing);
    void setOrigin(double x, double y, double z);
    void setOrigin(const Vec3d& origin);
    void swapValues(vector<T>& values_to_swap);
    
private:
    inline int index(int ix, int iy, int iz) { return ix*n_grid_.y*n_grid_.z + iy*n_grid_.z + iz; };

    Vec3i n_grid_;
    Vec3d spacing_;
    Vec3d origin_;
    vector<T> values_;
    
};
