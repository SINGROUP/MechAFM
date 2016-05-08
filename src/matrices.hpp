/*
 * matrices.hpp
 * 
 * Contains template for a 3x3 matrix (Matrix3) and implements common operations
 * for the matrix itself and between Matrix3 and Vector3 templates.
 * Defines type Mat3d (double type) for easy template access.
 * 
 */

#pragma once

#include <cmath>
#include <iomanip>
#include <iostream>

#include "vectors.hpp"

/** \brief A template for a 3x3 matrix
 * 
 * Defines template for a 3x3 matrix and implements common operations for
 * for the matrix itself and between Matrix3 and Vector3 templates.
 * 
 */

template<typename T>
class Matrix3 {

public:
    Matrix3() {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                values_[i][j] = 0;
            }
        }
    }
    
    Matrix3(const T& val) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                values_[i][j] = val;
            }
        }
    }
    
    ~Matrix3() {};
    
    
    // Element access
    T& at(int i, int j) { return values_[i][j]; };
    const T& at(int i, int j) const { return values_[i][j]; };
    
    // Column access
    Vec3d getColumn(int j) const {
        Vec3d column_vec;
        column_vec.x = values_[0][j];
        column_vec.y = values_[1][j];
        column_vec.z = values_[2][j];
        return column_vec;
    }
    
    // Print elements
    void print() const {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                std::cout << std::setw(10) << std::setprecision(4) << values_[i][j];
            }
            std::cout << std::endl;
        }
    }
    
    // Check if the matrix is diagonal
    bool isDiagonal() {
        T diag_sum = 0;
        T nondiag_sum = 0;
        
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (i == j)
                    diag_sum += std::abs(values_[i][i]);
                else
                    nondiag_sum += std::abs(values_[i][j]);
            }
        }
        
        return nondiag_sum < 1.0e-10*diag_sum;
    }
    
    // Inverse (https://en.wikipedia.org/wiki/Invertible_matrix#Inversion_of_3.C3.973_matrices)
    Matrix3 inverse() const {
        Matrix3<T> inv_matrix;
        
        T inv_det = 1/determinant();
        
        T a = values_[0][0];
        T b = values_[0][1];
        T c = values_[0][2];
        T d = values_[1][0];
        T e = values_[1][1];
        T f = values_[1][2];
        T g = values_[2][0];
        T h = values_[2][1];
        T i = values_[2][2];
        
        inv_matrix.at(0, 0) = inv_det*(e*i - f*h);
        inv_matrix.at(0, 1) = -inv_det*(b*i - c*h);
        inv_matrix.at(0, 2) = inv_det*(b*f - c*e);
        inv_matrix.at(1, 0) = -inv_det*(d*i - f*g);
        inv_matrix.at(1, 1) = inv_det*(a*i - c*g);
        inv_matrix.at(1, 2) = -inv_det*(a*f - c*d);
        inv_matrix.at(2, 0) = inv_det*(d*h - e*g);
        inv_matrix.at(2, 1) = -inv_det*(a*h - b*g);
        inv_matrix.at(2, 2) = inv_det*(a*e - b*d);
        
        return inv_matrix;
    }
    
    // Transpose
    Matrix3 transpose() const {
        Matrix3<T> t_matrix;
        
        t_matrix.at(0, 0) = values_[0][0];
        t_matrix.at(0, 1) = values_[1][0];
        t_matrix.at(0, 2) = values_[2][0];
        t_matrix.at(1, 0) = values_[0][1];
        t_matrix.at(1, 1) = values_[1][1];
        t_matrix.at(1, 2) = values_[2][1];
        t_matrix.at(2, 0) = values_[0][2];
        t_matrix.at(2, 1) = values_[1][2];
        t_matrix.at(2, 2) = values_[2][2];
        
        return t_matrix;
    }
    
    // Determinant
    T determinant() const {
        T det = 0;
        
        det = values_[0][0] * (values_[1][1]*values_[2][2] - values_[1][2]*values_[2][1]);
        det += values_[0][1] * (values_[1][2]*values_[2][0] - values_[1][0]*values_[2][2]);
        det += values_[0][2] * (values_[1][0]*values_[2][1] - values_[1][1]*values_[2][0]);
        
        return det;
    }
    
    // Matrix-vector multiplication
    Vector3<T> multiply(Vector3<T> rhs_vec) const {
        Vector3<T> result_vec;
        
        result_vec.x = values_[0][0]*rhs_vec.x + values_[0][1]*rhs_vec.y + values_[0][2]*rhs_vec.z;
        result_vec.y = values_[1][0]*rhs_vec.x + values_[1][1]*rhs_vec.y + values_[1][2]*rhs_vec.z;
        result_vec.z = values_[2][0]*rhs_vec.x + values_[2][1]*rhs_vec.y + values_[2][2]*rhs_vec.z;
        
        return result_vec;
    }
    
    // Matrix-matrix multiplication
    Matrix3 multiply(Matrix3 rhs_mat) const {
        Matrix3<T> result_mat;
        
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    result_mat.at(i, j) += values_[i][k]*rhs_mat.at(k, j);
                }
            }
        }
        
        return result_mat;
    }


private:
    T values_[3][3];
    
};


typedef Matrix3<double> Mat3d;
