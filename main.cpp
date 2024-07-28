#include "FVM.h"
#include "FVM.cpp"

int main() {
    double error_req = 1e-4;

    FVM fvm(error_req);
    fvm.solve();
    fvm.write_csv("1DHD.csv");

    return 0;
}
