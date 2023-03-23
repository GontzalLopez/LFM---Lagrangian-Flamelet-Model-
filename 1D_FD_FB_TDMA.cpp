/*
	This code uses the Crank-Nicolson scheme. 
	It constructs a tridiagonal matrix to solve the system of equations at each time step, with the boundary conditions 
	imposed separately. The output is written to a file `output.dat

	The erfc_inv function uses the bisection algorithm to approximate the inverse complementary error function. 

*/

#include <iostream>
#include <cmath>
#include <limits>
#include <fstream>

using namespace std;

#include <erfc_bisection.H>

// Diffusion coefficient function 
double D(double x, double chiSt, const double chiInfSt) 
{
    double compErfc = erfc_inv(2*x);
    double chiInf = exp(-2*pow(compErfc,2));
    double Fst = chiInf/(chiInfSt);

    return 0.5*chiSt*Fst;
}

// Compute chi profile
double chiProf(double x) 
{
    double InvErfc = erfc_inv(2*x);
    cout<<"InvErfc for "<< x <<" ="<< InvErfc << endl;
    double a = 50;

    return (a/3.1416)*exp(-2*pow(InvErfc,2));
}

// Parameters
const double L = 1.0;   // Length of the domain
const double T = 0.01;   // Final time
const int N = 200;      // Number of grid points
const double alpha = 0.0; // Reaction rate
const double Ns = 3.0; // Number of species

const double H2_L = 2.34261e-02;  // Left boundary condition
const double H2_R = 0.0;  // Right boundary condition
const double O2_L = 0.0;  // Left boundary condition
const double O2_R = 0.237;  // Right boundary condition
const double N2_L = 0.0;  // Left boundary condition
const double N2_R = 1-O2_R;  // Right boundary condition

const double chi_st=0.5; //Stoichiometric scalar dissipation
const double Zst = 0.4789;


int main() 
{
    // Initialize variables
    double dx = L/(N-1);
    double dt = 1e-5;
    double r = dt/(2*dx*dx);
    double t = 0.0;

    const double compErfc0 = erfc_inv(2*Zst);
    const double chiInfSt = exp(-2*pow(compErfc0,2));
    
    // Create arrays
    double u[N][3];
    double Dv[N];
    double chiProfile[N];

    double A[N-2], B[N-2], C[N-2], F[N-2];
    
    // Initialize u and Dv
    for (int i=0; i<N; i++) 
    {
        u[i][0] = 0.0;
        u[i][1] = 0.0;
        u[i][2] = 0.0;
        Dv[i] = D(i*dx,chi_st,chiInfSt);
	chiProfile[i] = chiProf(i*dx);
    }

    //Boundary conditions

    //H2	
    u[0][0] = H2_L;
    u[N-1][0] = H2_R;

    //O2	
    u[0][1] = O2_L;
    u[N-1][1] = O2_R;

    //N2	
    u[0][2] = N2_L;
    u[N-1][2] = N2_R;
    
    // Open output file
    ofstream outfile;
    ofstream outfile2;
    outfile.open("output.dat");
    outfile2.open("chiProfile.dat");
    
    // Time loop
    while (t < T) 
    {
        // Construct the tridiagonal matrix

    //	cout << "Solving for time step "<<t<<endl;
     //   cout << "... "<<endl;
	for(int k=0;k<Ns; k++)
	{
		for (int i=0; i<N-2; i++) 
		{

		    double Dm = (Dv[i]+Dv[i+1])/2.0;
		    double Dp = (Dv[i+1]+Dv[i+2])/2.0;

		    dt = std::min(dt,pow(dx,2)/(2*Dp));

		    A[i] = -r*Dm;
		    B[i] = 1.0 + 2.0*r*Dm + dt*alpha;
		    C[i] = -r*Dp;
		    F[i] = u[i+1][k] + r*Dm*u[i][k] + r*Dp*u[i+2][k] + dt*alpha*u[i+1][k];
		}
		
		// Solve the tridiagonal system
		for (int i=1; i<(N-1); i++) 
		{
		    double m = A[i-1]/B[i-1];
		    B[i] -= m*C[i-1];
		    F[i] -= m*F[i-1];
		}

		u[N-2][k] = F[N-3]/B[N-3];

		for (int i=N-3; i>=1; i--) 
		{
		    u[i][k] = (F[i-1] - C[i-1]*u[i+1][k])/B[i-1];
		}
	}
        t += dt;

    }

    // Update time and output
    for (int i=0; i<N; i++) 
    {
	outfile << i*dx << "	" << t << "	" << u[i][0] << "	" << u[i][1] << "	" << u[i][2] <<endl;
	outfile2 << i*dx << "	" << chiProfile[i] << endl;
    }

    // Close output file
    outfile.close();
    
    return 0;
}

