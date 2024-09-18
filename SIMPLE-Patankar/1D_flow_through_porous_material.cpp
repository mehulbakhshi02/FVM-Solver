#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <functional>
using namespace std;

int main () {

//Given parameters
double AB = 5.0;
double AC = 4.0;
double CB = 0.25;
double CC = 0.2;
double p1 = 200.00;
double p3 = 38;
double deltax = 2.0;

//under-relaxation, tolerance
double tolerance = 1.0e-6;
double alphaP = 0.8;
double alphaU = 0.9;
int imax = 200;

//intial guess values (Step 1)
double uB=15.0;
double uC=15.0;
double p2=120.0;

for (int i=0; i<imax; i++) {
	//Coefficient calculation
	double aB = CB*abs(uB)*deltax;
	double aC = CC*abs(uC)*deltax;

	//Step 2
	uB = ((p1-p2)+(1-alphaU)/alphaU*aB*uB)*(alphaU/aB);
	uC = ((p2-p3)+(1-alphaU)/alphaU*aC*uC)*(alphaU/aC);

	//Step 4
	double p2prime = (AB*uB-AC*uC)/(AC/aC+AB/aB);

	//Step 5
	double uBprime = -p2prime/aB;
	double uCprime = p2prime/aC;
	uB = uB+uBprime;
	uC = uC+uCprime;
	p2 = p2 + alphaP*p2prime;

	//Step 6
	double b = AB*uB - AC*uC;

	//Step 8
	double u_residual = abs(uB*aB/alphaU - ((p1-p2)+(1-alphaU)/alphaU*aB*uB)) + abs(uC*aC/alphaU - ((p2-p3)+(1-alphaU)/alphaU*aC*uC));
	u_residual = u_residual/(abs(uB*aB/alphaU + uC*aC/alphaU));

	double c_residual = abs(AB*uB - AC*uC);
	c_residual = c_residual/(0.5*(abs(AB*uB)+abs(AC*uC)));

	if (i==0){
		cout << "Iterations " << "uB " << "uC " << "p2 " << endl;
	}
	cout << i << " " << uB << " " << uC << " " << p2 << endl;

	if (u_residual+c_residual < tolerance){
		cout << "converged solutions is uB: " << uB << " uC: " << uC << " p2: " << p2 << endl;
		break;
	}

}
}


