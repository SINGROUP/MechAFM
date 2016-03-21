/*
 * cube_io.hpp
 * 
 * Copyright 2016 Juha Ritala <juha.ritala@aalto.fi>
 * 
 * CubeReader class represents a read access to the volumetric data of a cube file
 * and contains the metadata of the file and the atom types and positions.
 * 
 */

#pragma once

#include <fstream>
#include <vector>

#include "data_grid.hpp"
#include "vectors.hpp"

const double bohr_to_angst = 0.52917721092;

using namespace std;

/** \brief A cube file parser.
 *  
 * Represents a read access to the volumetric data of a cube file
 * and contains the metadata of the file and the atom types and positions.
 */

class CubeReader {
public:
    CubeReader(const string& filepath);
    ~CubeReader() { cube_file_.close(); };
    
    const Vec3i& getNVoxels() const { return n_voxels_; };
    //TODO: return voxel_vectors in either Bohr or Angstrom
    const vector<Vec3d>& getVoxelVectors() const { return voxel_vectors_; };
    Vec3d getVoxelSpacing() const;
    const Vec3d& getOrigin() const { return origin_; };
    int getNAtoms() const { return n_atoms_; };
    
    // allocates a vector and returns it after filling it with volumetric data from the file
    vector<double> readVolumetricData();
    // stores the volumetric data from the file to the preallocated vector given as a reference
    void readVolumetricData(vector<double>& volumetric_data);
    // stores all contents to a DataGrid object (only for orthogonal voxels)
    void storeToDataGrid(DataGrid<double>& data_grid, const Vec3d& offset = Vec3d(0.0));

private:
    bool isVoxelsOrthogonal() const;
    bool readMetadata();

    ifstream cube_file_;
    int volumetric_data_pos_;
    string comment_lines_;
    Vec3i n_voxels_;
    vector<Vec3d> voxel_vectors_;
    Vec3d origin_;
    int n_atoms_;
    vector<int> atom_numbers_;
    vector<Vec3d> atom_positions_;
};
