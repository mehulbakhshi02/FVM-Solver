#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <functional>
using namespace std;

int main () {

//Given parameters
double rho = 1.0;
double AA=3.0;
double AB=1.0;
double p1 = 28.0;
double p3 =0.0;

//under-relaxation, tolerance
double tolerance = 1.0e-6;
double alphaP = 0.8;
double alphaU = 0.1;
int imax = 200;

//intial guess values (Step 1)
double uA=5.0/3.0;
double uB=5.0;
double p2=25.0;

for (int i=0; i<imax; i++) {
	//Coefficient calculation
	double FA = rho*uA*AA;
	double FB = rho*uB*AB;

	//Step 2
	uA = (AA*(p1-p2)+(1-alphaU)/alphaU*FA*uA)*alphaU/FA;
	uB = (FA*uA+AB*(p2-p3)+(1-alphaU)/alphaU*FB*uB)*alphaU/FB;

	//Step 4
	double dA = AA/FA;
	double dB = AB/FB;
	double p2prime = ((uA*AA-uB*AB)/(dA*AA+dB*AB));

	//Step 5
	double uAprime = -AA*p2prime/FA;
	double uBprime = AB*p2prime/FB;
	uA = uA+uAprime;
	uB = uB+uBprime;
	p2 = p2 + alphaP*p2prime;

	//Step 6
	double b = rho*uB*AB-rho*uA*AA;

	//Step 8
	double u_residual = abs(FA*uA/alphaU - (AA*(p1-p2)+(1-alphaU)/alphaU*FA*uA)) + abs(FB*uB/alphaU - (FA*uA+AB*(p2-p3)+(1-alphaU)/alphaU*FB*uB));
	//u_residual = u_residual/(abs(FA*uA/alphaU + FB*uB/alphaU));

	double c_residual = abs(b);
	//c_residual = c_residual/(0.5*(abs(uA*AA)+abs(uB*AB)));

	if (i==0){
		cout << "Iterations " << "uA " << "uB " << "p2 " << endl;
	}
	cout << i << " " << uA << " " << uB << " " << p2 << endl;

	if (u_residual+c_residual < tolerance){
		cout << "converged solutions is uA: " << uA << " uB: " << uB << " p2: " << p2 << endl;
		break;
	}

}
}


