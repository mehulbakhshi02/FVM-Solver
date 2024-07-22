#ifndef FVM_H
#define FVM_H

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "Eigen/Dense"
using namespace std;


//#define HeatDiff1D
//#define HeatDiff2D
#define HeatConv2D

#if defined HeatDiff1D

class FVM {
public:
    FVM(double error_req);
    void solve();
    void write_csv(string filename);

private:
//-------------------------------------------------User Defined Parameters--------------------------------
    int nx=100;
    double xL=1;
    double k=1.0;
    double Su=0.0;
    double Sp=0.0;

    #define Dirichlet_E
    double phi_e=1.0;
    #define Dirichlet_W
    double phi_w=0.0;
    #define Neumann_N
    double dphi_n=0.0;
    #define Neumann_S
    double dphi_s=0.0;

//----------------------------------------------------------------------------------------------------------

    int ny=1;
    double yL=0.1;
    double dx = xL/nx;
    double dy = yL/ny;
    double error_req;
    Eigen::MatrixXd phi;
    Eigen::MatrixXd phi_old;
    Eigen::MatrixXd A_P;
    Eigen::MatrixXd A_E;
    Eigen::MatrixXd A_W;
    Eigen::MatrixXd A_N;
    Eigen::MatrixXd A_S;
    Eigen::MatrixXd b;

    void initialize();
    double compute_error();
    };

#endif

#if defined HeatDiff2D

class FVM {
public:
    FVM(double error_req);
    void solve();
    void write_csv(string filename);

private:
//-------------------------------------------------User Defined Parameters--------------------------------
    int nx=100;
    int ny=100;
    double xL=1.0;
    double yL=1.0;
    double k=1.0;
    double Su=1.0;
    double Sp=1.0;
    double rho=1.0;
    double u=0.0;
    double v=0.0;

    //boundary conditions
    #define Dirichlet_E
    double phi_e=1.0;
    #define Dirichlet_W
    double phi_w=0.0;
    #define Dirichlet_N
    double phi_n=0.0;
    #define Dirichlet_S
    double phi_s=0.0;

    //Convection boundary conditions
    int in_w=0;
    int out_w=0;
    int in_e=0;
    int out_e=0;
    int in_n=0;
    int out_n=0;
    int in_s=0;
    int out_s=0;

//-------------------------------------------------------------------------------------------------------

    double error_req;
    double dx = xL/nx;
    double dy = yL/ny;
    Eigen::MatrixXd phi;
    Eigen::MatrixXd phi_old;
    Eigen::MatrixXd A_P;
    Eigen::MatrixXd A_E;
    Eigen::MatrixXd A_W;
    Eigen::MatrixXd A_N;
    Eigen::MatrixXd A_S;
    Eigen::MatrixXd b;

    void initialize();
    double compute_error();
    };

#endif

#if defined HeatConv2D

class FVM {
public:
    FVM(double error_req);
    void solve();
    void write_csv(string filename);

private:
//-------------------------------------------------User Defined Parameters--------------------------------
    int nx=100;
    int ny=100;
    double xL=1.0;
    double yL=1.0;
    double k=1.0;
    double Su=0.0;
    double Sp=0.0;
    double rho=1.0;
    double u=10.0;
    double v=10.0;

    //Diffusion boundary conditions
    #define Dirichlet_E
    double phi_e=0.0;
    #define Dirichlet_W
    double phi_w=1.0;
    #define Dirichlet_N
    double phi_n=1.0;
    #define Dirichlet_S
    double phi_s=0.0;

    //Convection boundary conditions
    int in_w=1;
    int out_w=0;
    int in_e=0;
    int out_e=1;
    int in_n=1;
    int out_n=0;
    int in_s=0;
    int out_s=1;

//-------------------------------------------------------------------------------------------------------

    double error_req;
    double dx = xL/nx;
    double dy = yL/ny;
    Eigen::MatrixXd phi;
    Eigen::MatrixXd phi_old;
    Eigen::MatrixXd A_P;
    Eigen::MatrixXd A_E;
    Eigen::MatrixXd A_W;
    Eigen::MatrixXd A_N;
    Eigen::MatrixXd A_S;
    Eigen::MatrixXd b;

    void initialize();
    double compute_error();
    };

#endif

//Dirichlet and Neumann BC
#if defined Dirichlet_E
int dir_e = 1;
int neu_e = 0;
double dphi_e = 0.0;
#endif

#if defined Dirichlet_W
int dir_w = 1;
int neu_w = 0;
double dphi_w = 0.0;
#endif

#if defined Dirichlet_N
int dir_n = 1;
int neu_n = 0;
double dphi_n = 0.0;
#endif

#if defined Dirichlet_S
int dir_s = 1;
int neu_s = 0;
double dphi_s = 0.0;
#endif

#if defined Neumann_E
int dir_e = 0;
int neu_e = 1;
double phi_e = 0.0;
#endif

#if defined Neumann_W
int dir_w = 0;
int neu_w = 1;
double phi_w = 0.0;
#endif

#if defined Neumann_N
int dir_n = 0;
int neu_n = 1;
double phi_n = 0.0;
#endif

#if defined Neumann_S
int dir_s = 0;
int neu_s = 1;
double phi_s = 0.0;
#endif

#endif