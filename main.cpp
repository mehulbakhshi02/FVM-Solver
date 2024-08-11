#include "FVM.h"
#include "FVM.cpp"

int main() {
    double error_req = 1e-6;

    FVM fvm(error_req);
    fvm.solve();
    fvm.write_csv("2DHD.csv");

    return 0;
}
