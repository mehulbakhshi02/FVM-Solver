#include "FVM.h"
#include "FVM.cpp"

int main() {
    double error_req = 1e-2;

    FVM fvm(error_req);
    fvm.solve();
    fvm.write_csv("2DHD.csv");
    
    //Mehul git test

    return 0;
}
