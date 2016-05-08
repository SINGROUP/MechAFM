/*
 * fft.hpp
 * 
 * A wrapper for Kiss FFT library.
 * 
 */

#pragma once

#include <complex>

#include "data_grid.hpp"
#include "kiss_fftnd.h"

using namespace std;

typedef complex<double> dcomplex;


/** \brief Does FFT on data_grid_in and stores the result to data_grid_out
 * 
 *  data_grid_out is allocated by the function, so you don't need to do that beforehand.
 *  Automatically changes grid basis from real space to k-space.
 */
void fft_data_grid(const DataGrid<double>& data_grid_in, DataGrid<dcomplex>& data_grid_out);


/** \brief Does inverse FFT on data_grid_in and stores the result to data_grid_out
 * 
 *  data_grid_out is allocated by the function, so you don't need to do that beforehand.
 *  The result is scaled using factor 1/total_number_of_grid_points. Automatically
 *  changes grid basis from k-space to real space.
 */
void ffti_data_grid(const DataGrid<dcomplex>& data_grid_in, DataGrid<double>& data_grid_out);
