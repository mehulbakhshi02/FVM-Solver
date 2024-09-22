#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <functional>
using namespace std;

int main () {

//Given parameters
double deltax = 2.0;
double A1 = 6.0;
double A2 = 4.0;
double A3 = 2.0;
double u1 = 10.0;
double u3 = 30.0;
double c1 = 10.0;
double c2 = 10.0;
double c3 = 10.0;
double cB = 10.0;
double cC = 10.0;

//tolerance
double tolerance = 1.0e-6;
double alphaP = 0.8;
double alphaU = 0.8;
int imax = 200;

//intial guess values (only for pressure as linear momentum equations) (Step 1)
double uB = 15.0;
double uC = 15.0;
double pB = 120.0;
double pC = 120.0;
double p1 = 120.0;
double p2 = 120.0;
double p3 = 120.0;

double pCprime = 0.0;

for (int i=0; i<imax; i++) {
	//Coefficient calculation
	double aB = cB*abs(uB)*deltax/alphaU;
    double aC = cC*abs(uC)*deltax/alphaU;
    double a1 = c1*abs(u1)*deltax/(2.0*alphaU);
    double a3 = c3*abs(u3)*deltax/(2.0*alphaU);

    double bB = (1.0-alphaU)*aB*uB;
    double bC = (1.0-alphaU)*aC*uC;
    double b1 = (1.0-alphaU)*a1*u1;
    double b3 = (1.0-alphaU)*a3*u3;

	//Step 2 (momentum eqn solving only at cell centers)
    uB = ((p1-p2)+bB)/aB;
    uC = ((p2-p3)+bC)/aC;

    // Hat velocity calculation
    double uBhat = bB/aB;
    double uChat = bC/aC;
    double u1hat = b1/a1;
    double u3hat = b3/a3;
    double u2hat = (uBhat+uChat)/2.0;

    double dB = 1.0/aB;
    double dC = 1.0/aC;
    double d1 = 1.0/a1;
    double d3 = 1.0/a3;
    double d2 = (dB+dC)/2.0;

    //Step 3
    double u2 = u2hat + d2*(pB-pC);

	//Step 4
	double pBprime = (u1*A1-u2*A2)/(d2*A2);

	//Step 5
    // correctors
	double p1prime = pBprime;
	double p3prime = pCprime;
    double p2prime = (pBprime+pCprime)/2;

    // face velocity correction
    double u2prime = d2*(pBprime - pCprime);
	u2 = u2 + u2prime;

    // cell pressure corrections
    pC = pC + alphaP*pCprime;
    pB = pB + alphaP*pBprime;

    // cell velocity corrections
    uB = uB + dB*(p1prime-p2prime);
    uC = uC + dC*(p2prime-p3prime);

    // face pressure corrections
    p1 = pB + (u1 - u1hat)/d1;
    p3 = pC - (u3 - u3hat)/d3;
    p2 = (pB+pC)/2.0;

	//Step 8
	double u_residual = abs(uB*aB-(p1-p2)-bB) + abs(uC*aC-(p2-p3)-bC);
	u_residual = u_residual/(abs(uB*aB) + abs(uC*aC));

	double c_residual = abs(u1*A1 - u2*A2)+abs(u2*A2-u3*A3);
	c_residual = c_residual/(0.5*(abs(A1*u1)+abs(A2*u2)));

	if (i==0){
		cout << "Iterations " << "uB " << "uC " << "pB " << "pC " << "p1 " << "p3 " << endl;
	}
	cout << i << " " << uB << " " << uC << " " << pB << " " << pC << " " << p1 << " " << p3 << endl;

	if (u_residual+c_residual < tolerance){
		cout << "converged solutions is uB: " << uB << " uC: " << uC << " pB: " << pB << " pC: " << pC << " p1: " << p1 << " p3: " << p3 << endl;
		break;
	}

}
}
