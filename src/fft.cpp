#include "fft.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

#include "matrices.hpp"
#include "vectors.hpp"


void fft_data_grid(const DataGrid<double>& data_grid_in, DataGrid<dcomplex>& data_grid_out) {
    const Vec3i& n_grid = data_grid_in.getNGrid();
    Mat3d basis = data_grid_in.getBasis();
    int n_grid_total = n_grid.x * n_grid.y * n_grid.z;
    
    vector<kiss_fft_cpx> realspace_in, kspace_out;
    realspace_in.assign(n_grid_total, {0.0, 0.0});
    kspace_out.assign(n_grid_total, {0.0, 0.0});
    
    // Copy values from data_grid_in to realspace_in (from double to kiss_fft_cpx type)
    for (int i = 0; i < n_grid_total; i++) {
        realspace_in[i].r = data_grid_in.at(i);
    }
    
    // Do the FFT
    int dims[] = {n_grid.x, n_grid.y, n_grid.z};
    int n_dims = 3;
    int is_inverse = 0;
    kiss_fftnd_cfg st = kiss_fftnd_alloc(dims, n_dims, is_inverse, nullptr, nullptr);
    kiss_fftnd(st, realspace_in.data(), kspace_out.data());
    free(st);
    
    // Copy the values from kspace_out to data_grid_out (from kiss_fft_cpx to complex<double> type)
    data_grid_out.initValues(n_grid.x, n_grid.y, n_grid.z, dcomplex(0.0, 0.0));
    for (int i = 0; i < n_grid_total; i++) {
        data_grid_out.at(i).real(kspace_out[i].r);
        data_grid_out.at(i).imag(kspace_out[i].i);
    }
    
    // Set the basis in k-space
    Mat3d k_basis;
    Mat3d n_grid_scaling = Mat3d(0);
    n_grid_scaling.at(0, 0) = n_grid.x;
    n_grid_scaling.at(1, 1) = n_grid.y;
    n_grid_scaling.at(2, 2) = n_grid.z;
    k_basis = basis.multiply(n_grid_scaling).inverse().transpose();
    data_grid_out.setBasis(k_basis);
}


void ffti_data_grid(const DataGrid<dcomplex>& data_grid_in, DataGrid<double>& data_grid_out) {
    const Vec3i& n_grid = data_grid_in.getNGrid();
    Mat3d k_basis = data_grid_in.getBasis();
    int n_grid_total = n_grid.x * n_grid.y * n_grid.z;
    
    vector<kiss_fft_cpx> kspace_in, realspace_out;
    kspace_in.assign(n_grid_total, {0.0, 0.0});
    realspace_out.assign(n_grid_total, {0.0, 0.0});
    
    // Copy values from data_grid_in to kspace_in (from complex<double> to kiss_fft_cpx type)
    for (int i = 0; i < n_grid_total; i++) {
        kspace_in[i].r = data_grid_in.at(i).real();
        kspace_in[i].i = data_grid_in.at(i).imag();
    }
    
    // Do the inverse FFT
    int dims[] = {n_grid.x, n_grid.y, n_grid.z};
    int n_dims = 3;
    int is_inverse = 1;
    kiss_fftnd_cfg st = kiss_fftnd_alloc(dims, n_dims, is_inverse, nullptr, nullptr);
    kiss_fftnd(st, kspace_in.data(), realspace_out.data());
    free(st);
    
    // Copy the values from kspace_out to data_grid_out (from kiss_fft_cpx to double type)
    // Normalize at the same time with 1/(nx*ny*nz)
    data_grid_out.initValues(n_grid.x, n_grid.y, n_grid.z, 0.0);
    for (int i = 0; i < n_grid_total; i++) {
        data_grid_out.at(i) = realspace_out[i].r / n_grid_total;
        if (realspace_out[i].i / n_grid_total > 1.0e-8) {
            char buf[100];
            sprintf(buf, "Imaginary part of inverse FFT exceeded threshold 1.0e-8! Value: %e", realspace_out[i].i);
            throw runtime_error(buf);
        }
    }
    
    // Set the basis in real space
    Mat3d basis;
    Mat3d n_grid_scaling = Mat3d(0);
    n_grid_scaling.at(0, 0) = n_grid.x;
    n_grid_scaling.at(1, 1) = n_grid.y;
    n_grid_scaling.at(2, 2) = n_grid.z;
    basis = k_basis.multiply(n_grid_scaling).inverse().transpose();
    data_grid_out.setBasis(basis);
}
