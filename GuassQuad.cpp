#include <iostream>
#include <cmath>
#include <functional>
#include "FVM.h"
#include "Eigen/Dense"

using namespace std;

FVM::FVM(double error_req)
    : error_req(error_req) {
    initialize(), verify();
}

// Function that we take in place of phi
double FVM::phifunction(double x, double y){
    return M_PI*(u*cos(M_PI*x)*cos(M_PI*y)-v*sin(M_PI*x)*sin(M_PI*y)+2*M_PI*sin(M_PI*x)*cos(M_PI*y));
}

// 3 point Guass Quad
double FVM::gaussQuad(function<double(double, double)> func, double llim_x, double ulim_x, double llim_y, double ulim_y){
    double weight[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};
    double nodes[3] = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)};

    //transformation
    double a1_x = (llim_x+ulim_x)/2;
    double a2_x = (ulim_x-llim_x)/2;
    double a1_y = (llim_y+ulim_y)/2;
    double a2_y = (ulim_y-llim_y)/2;

    double integral = 0.0;

    //x is j and y is i
    for (int i=0; i < 3; i++){
        for (int j=0; j < 3; j++){
            integral += weight[i]*weight[j]*(func(a1_x+a2_x*nodes[j], a1_y+a2_y*nodes[i]));
        }
    }
    integral = integral*a2_x*a2_y;

    return integral;

}

void FVM::initialize() {
    phi = Eigen::MatrixXd::Zero(ny, nx);
    phi_old = Eigen::MatrixXd::Zero(ny, nx);
    A_P = Eigen::MatrixXd::Constant(ny, nx, rho*(u*dy-u*dy+v*dx-v*dx));
    A_E = Eigen::MatrixXd::Constant(ny, nx, k*dy/dx-(rho*u*dy)/2);
    A_W = Eigen::MatrixXd::Constant(ny, nx, k*dy/dx+(rho*u*dy)/2);
    A_N = Eigen::MatrixXd::Constant(ny, nx, k*dx/dy-(rho*v*dx)/2);
    A_S = Eigen::MatrixXd::Constant(ny, nx, k*dx/dy+(rho*v*dx)/2);
    b = Eigen::MatrixXd::Zero(ny, nx);

    A_E.col(nx-1).setZero();
    A_W.col(0).setZero();
    A_N.row(ny-1).setZero();
    A_S.row(0).setZero();
    A_P = A_E+A_W+A_S+A_N;
}

void FVM::verify(){
    int row, col;
    
    for (row=0; row<ny; row++){
            for (col=0; col<nx; col++){
                double llim_x = col*dx;
                double ulim_x = (col+1)*dx;
                double llim_y = row*dy;
                double ulim_y = (row+1)*dy;
                double result = gaussQuad(bind(&FVM::phifunction, this, placeholders::_1, placeholders::_2), llim_x, ulim_x, llim_y, ulim_y);
                b(row, col) = result;
            }
    }
}