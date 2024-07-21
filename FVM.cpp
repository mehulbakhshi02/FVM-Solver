#include "FVM.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "Eigen/Dense"
using namespace std;

FVM::FVM(double error_req)
    : error_req(error_req) {
    initialize();
}

void FVM::initialize() {
    phi = Eigen::MatrixXd::Zero(ny, nx);
    phi_old = Eigen::MatrixXd::Zero(ny, nx);
    A_P = Eigen::MatrixXd::Zero(ny, nx);
    A_E = Eigen::MatrixXd::Constant(ny, nx, k*dy/dx);
    A_W = Eigen::MatrixXd::Constant(ny, nx, k*dy/dx);
    A_N = Eigen::MatrixXd::Constant(ny, nx, k*dx/dy);
    A_S = Eigen::MatrixXd::Constant(ny, nx, k*dx/dy);
    b = Eigen::MatrixXd::Zero(ny, nx);

    A_E.col(nx-1).setZero();
    A_W.col(0).setZero();
    A_N.row(ny-1).setZero();
    A_S.row(0).setZero();
    A_P = A_E+A_W+A_S+A_N;

}

void FVM::solve() {
    // Implementation of B.C.
    int row, col;

   // West boundary
   col = 0;
   for (row = 0; row < ny; row++)
   {
       b(row,col) += dir_w*(phi_w*2*k*dy/dx)+neu_w*(k*dphi_w*dy);
        A_P(row,col) += dir_w*(2*dy/dx*k);        
   }
   // East boundary
   col = nx-1;
   for (row = 0; row < ny; row++)
   {
       b(row,col) += dir_e*(phi_e*2*k*dy/dx)+neu_e*(k*dphi_e*dy);
        A_P(row,col) += dir_e*(2*dy/dx*k);
   }

   // North boundary
   row = ny-1;
   for (col = 0; col < nx; col++)
   {
       b(row,col) += dir_n*(phi_n*2*k*dx/dy)+neu_n*(k*dphi_n*dx);
        A_P(row,col) += dir_n*(2*dx/dy*k);
   }

   // South boundary
   row = 0;
   for (col = 0; col < nx; col++)
   {
       b(row,col) += dir_s*(phi_s*2*k*dx/dy)+neu_s*(k*dphi_s*dx);
       A_P(row,col) += dir_s*(2*dx/dy*k);
   }

    double error_mag = 1.0;
    int iterations = 0;
    double S, E, W, N;

    while (error_mag > error_req) {
        for (row=0; row<ny; row++){
            for (col=0; col<nx; col++){
               S=0.0; N=0.0; E = 0.0; W = 0.0;

                // South face cells
               if(row!=0)
                   S = A_S(row, col) * phi(row-1, col);

               // North face cells
               if(row!=ny-1)
                   N = A_N(row, col) * phi(row+1, col);

               // West face cells
               if(col!=0)
                   W = A_W(row, col) * phi(row, col-1);

               // East face cells
               if(col!=nx-1)
                   E = A_E(row, col) * phi(row, col+1);

                phi(row, col) = (E + W + N + S + b(row, col)) / A_P(row, col);    
                    
            }
        }
        error_mag = compute_error();
        phi_old = phi;
        iterations++;
    }

    std::cout << "Iterations: " << iterations << std::endl;
}

double FVM::compute_error() {
    double error_mag = 0.0;
    Eigen::MatrixXd phi_error;
    phi_error = Eigen::MatrixXd::Zero(ny, nx);
    for (int row = 0; row<ny; row++) {
        for (int col=0; col<nx; col++){
                phi_error(row,col)= abs(phi(row, col)-phi_old(row, col));
       	        error_mag=error_mag+phi_error(row,col);             //Error_mag is summation of all errors across all points       
            }
    }
    return error_mag;
}

void FVM::write_csv(string filename) {
    ofstream myFile(filename);
    for (int row = 0; row<ny; row++) {
        for (int col=0; col<nx; col++)
        {
            myFile << row << "," << col << "," << phi(row,col) << endl;
        }
    }
    myFile.close();
}
