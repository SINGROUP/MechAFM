#include "cube_io.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace std;


CubeReader::CubeReader(const string& filepath) {
    voxel_vectors_.assign(3, Vec3d(0));
    
    // open file access and read metadata
    cube_file_.exceptions(ifstream::failbit | ifstream::badbit);
    cube_file_.open(filepath);
    bool success = readMetadata();

    if (not success) {
        cube_file_.close();
        throw runtime_error("Could not read metadata from the cube file. Check the format of the file.");
    }
}


//TODO: let user to choose whether to return spacing in units of Bohr or Angstrom
Vec3d CubeReader::getVoxelSpacing() const {
    Vec3d voxel_spacing;
    if (isVoxelsOrthogonal()) {
        voxel_spacing.x = voxel_vectors_[0].x;
        voxel_spacing.y = voxel_vectors_[1].y;
        voxel_spacing.z = voxel_vectors_[2].z;
    } else {
        throw runtime_error("The voxel vectors are not orthogonal and thus the spacing between "
            "voxels is undefined.");
    }
    return voxel_spacing;
}


vector<double> CubeReader::readVolumetricData() {
    vector<double> volumetric_data;
    volumetric_data.assign(n_voxels_.x * n_voxels_.y * n_voxels_.z, 0.0);
    try {
        cube_file_.seekg(volumetric_data_pos_);
        for (int i = 0; i < n_voxels_.x * n_voxels_.y * n_voxels_.z; i++) {
            cube_file_ >> volumetric_data[i];
        }
    }
    catch (ifstream::failure) {
        throw runtime_error("Could not read the volumetric data from the cube file. "
            "Check the format of the file.");
    }
    return volumetric_data;
}


void CubeReader::readVolumetricData(vector<double>& volumetric_data) {
    try {
        cube_file_.seekg(volumetric_data_pos_);
        for (int i = 0; i < n_voxels_.x * n_voxels_.y * n_voxels_.z; i++) {
            cube_file_ >> volumetric_data[i];
        }
    }
    catch (ifstream::failure) {
        throw runtime_error("Could not read the volumetric data from the cube file. "
            "Check the format of the file.");
    }
}


void CubeReader::storeToDataGrid(DataGrid<double>& data_grid, const Vec3d& offset) {
    vector<double> volumetric_data = readVolumetricData();
    data_grid.setNGrid(n_voxels_);
    data_grid.setBasis(voxel_vectors_);
    data_grid.setOrigin(origin_ + offset);
    data_grid.swapValues(volumetric_data);
}


bool CubeReader::isVoxelsOrthogonal() const {
    return ( abs(voxel_vectors_[0].y) + abs(voxel_vectors_[0].z) + \
        abs(voxel_vectors_[1].x) + abs(voxel_vectors_[1].z) + \
        abs(voxel_vectors_[2].x) + abs(voxel_vectors_[2].y) ) \
        < \
        1.0e-10 * (abs(voxel_vectors_[0].x) + abs(voxel_vectors_[1].y) + abs(voxel_vectors_[2].z));
}


bool CubeReader::readMetadata() {
    string line;
    double temp;
    bool is_bohr;
    
    try {
        // read the two comment lines
        getline(cube_file_, line);
        comment_lines_ = line + '\n';
        getline(cube_file_, line);
        comment_lines_ += line;
        
        // read number of atoms and grid metadata
        cube_file_ >> n_atoms_ >> origin_.x >> origin_.y >> origin_.z;
        cube_file_ >> n_voxels_.x >> voxel_vectors_[0].x >> voxel_vectors_[0].y >> voxel_vectors_[0].z;
        cube_file_ >> n_voxels_.y >> voxel_vectors_[1].x >> voxel_vectors_[1].y >> voxel_vectors_[1].z;
        cube_file_ >> n_voxels_.z >> voxel_vectors_[2].x >> voxel_vectors_[2].y >> voxel_vectors_[2].z;
        
        // length in Angstroms
        if (n_voxels_.x < 0) {
            is_bohr = false;
            n_voxels_ *= -1;
        }
        // length in Bohr
        else {
            is_bohr = true;
            for (int i = 0; i < 3; i++)
                voxel_vectors_[i] *= bohr_to_angst;
            origin_ *= bohr_to_angst;
        }
        
        // reserve memory for atomic data
        atom_numbers_.assign(n_atoms_, 0);
        atom_positions_.assign(n_atoms_, Vec3d(0));
        
        // read atom types and positions
        for (int ia = 0; ia < n_atoms_; ia++) {
            cube_file_ >> atom_numbers_[ia] >> temp >> atom_positions_[ia].x >> atom_positions_[ia].y >> atom_positions_[ia].z;
            if (is_bohr)
                atom_positions_[ia] *= bohr_to_angst;
        }
        
        // record the position of file at which the volumetric data begins
        volumetric_data_pos_ = cube_file_.tellg();
    }
    catch (ifstream::failure) {
        return false;
    }
    
    return true;
}
