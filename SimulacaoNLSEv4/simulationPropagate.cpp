
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <memory.h>
#include <fftw3.h>
#include <time.h>
#include <netcdfcpp.h>
//#include <clapack.h>
using namespace std;

#include "library.h"
#include "simulationVariables.h"

//Prototype
void ApplyTransverseLaplacian(double step);


/////////////////////////////////////////////////////////////
// Helper functions for nonlinear contribution calculation
//

inline double drho(double current_rho,double Emag_square)
{
	return sigmaK*(rho_neutral-current_rho)*pow(Emag_square,K)+
				                     sigma/Ui*current_rho*Emag_square;
}


inline double complex ionization(double current_rho,double Emag_square)
{
	return -1*I*M_PI/(lambdaZero*n*n*rhoC)*current_rho
			                      -sigma/2*current_rho
								  -Ui*sigmaK/2*pow(Emag_square,K-1)*(rho_neutral-current_rho);
}

inline double cmag_square(double complex number)
{
	double real = creal(number),imag = cimag(number);
	return real*real+imag*imag;
}
	
/////////////////////////////////////////////////////////////////////////////
//
//  Propagate
// 
//  Calculates E using the first order split step method with step 'step'
//	Uses E as initial condition, stores fft in fftE
void 	Propagate(double step)
{
	double 		omega,k,Emag_square;
	double 		tmpRho1,tmpRho2,tmpRho3,tmpRho4;
	double 		complex 	bZ;
	int 		i,j;

	k = 2*M_PI/lambdaZero;
	
//1) Calculate phi(f,z0)=F[phi(t,z0)]
// execute fftw forward: fftPhi = F[phi]
	/*fftw_execute(forward);

//2) Apply dispersive propagation z0->z0+h


	for(i = 0;i < NPOINTS_T;i++)
	{
		if(i >= NPOINTS_T/2)
			omega = 2*M_PI*(i-NPOINTS_T)/NPOINTS_T/DELTAT;
		else omega = 2*M_PI*i/NPOINTS_T/DELTAT;
		fftPhi[i] = fftPhi[i]*cexp(-GVD/2.*omega*omega*step/2.*I); 
	}		

//3) Calculate phi(t,z0+h/2)=F-1[phi(f,z0+h/2)]
// execute fftw backward: phi = F[fftPhi]
	fftw_execute(backward); */

//4) Apply nonlinear propagation z0->z0+h
// Apply nonlinear term with correction
// Also solve for rho(t) at that point and time (RungeKutta4)

	for(j = 0;j < NPOINTS_R;j++)
	{
		i = 0;
		
		rho[0+j*NPOINTS_T] = 0;
		
		Emag_square = cmag_square(E[i+j*NPOINTS_T]);
		bZ = -absorptionCalc[j]-I*omegaZero/c*n2*Emag_square+ionization(rho[0+j*NPOINTS_T],Emag_square);
		E[0+j*NPOINTS_T] = E[0+j*NPOINTS_T]*cexp(bZ*step); 
		
		for(i = 1;i < NPOINTS_T;i++)
		{
			Emag_square = cmag_square(E[i+j*NPOINTS_T]);
			
			//Calc rho(t)
			tmpRho1 = drho(rho[i+j*NPOINTS_T]                   ,Emag_square);
			tmpRho2 = drho(rho[i+j*NPOINTS_T]+0.5*DELTAT*tmpRho1,Emag_square); //TODO important!
			tmpRho3 = drho(rho[i+j*NPOINTS_T]+0.5*DELTAT*tmpRho2,Emag_square); //faltam updates a Emag_square
			tmpRho4 = drho(rho[i+j*NPOINTS_T]+    DELTAT*tmpRho3,Emag_square);
				
			rho[i+j*NPOINTS_T] = rho[i-1+j*NPOINTS_T] + 
				       DELTAT/6.*(tmpRho1+2*tmpRho2+2*tmpRho3+tmpRho4);
					   			
			//calc ionization nonlinear terms (without the E, because we are doing an exp.)
			
			bZ = -absorptionCalc[j]-I*omegaZero/c*n2*Emag_square+ionization(rho[i+j*NPOINTS_T],Emag_square);
			E[i+j*NPOINTS_T] = E[i+j*NPOINTS_T]*cexp(bZ*step); 
		}	
	}

//1) Calculate phi(f,z0)=F[phi(t,z0)]
// execute fftw forward: fftPhi = F[phi]
	for(i = 0;i < NPOINTS_R;i++)
		fftw_execute(forward[i]);

//2) Apply dispersive propagation z0->z0+h

	for(j = 0;j < NPOINTS_R;j++)
		for(i = 0;i < NPOINTS_T;i++)
		{
			if(i >= NPOINTS_T/2)
				omega = 2*M_PI*(i-NPOINTS_T)/NPOINTS_T/DELTAT;
			else omega = 2*M_PI*i/NPOINTS_T/DELTAT;
			fftE[i+j*NPOINTS_T] = fftE[i+j*NPOINTS_T]*cexp(-GVD/2.*omega*omega*step*I); 
		}	
		
	ApplyTransverseLaplacian(step);

	for(i = 0;i < NPOINTS_R;i++)
		fftw_execute(backward[i]); 

//4.7.1) Renormalize (fftw multiplies data per sqrt(N))
	for(i = 0;i < NPOINTS_T*NPOINTS_R;i++)
		E[i] = E[i]/NPOINTS_T;

	//result = propagated phi and calculated fftPhi
}


//http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
void Tridiagonal()
{
#define n NPOINTS_R
	int i;
	triRight[0] /= triDiag[0];				
	triIn[0] /= triDiag[0];				
	for(i = 1; i < n; i++)
	{
		double complex id = (triDiag[i] - triRight[i-1] * triLeft[i]);	
		triRight[i] /= id;				
		triIn[i] = (triIn[i] - triIn[i-1] * triLeft[i])/id;
	}
 	triOut[n - 1] = triIn[n - 1];
	for(i = n - 2; i >= 0; i--)
		triOut[i] = triIn[i] - triRight[i] * triOut[i + 1];
}

void CalcTriIn(int i)
{
	int j;
	for(j = 1;j < NPOINTS_R-1;j++)
		triIn[j] = triLeft[j]*fftE[i+(j-1)*NPOINTS_T]+
					triDiag[j]*fftE[i+j*NPOINTS_T]+
					triRight[j]*fftE[i+(j+1)*NPOINTS_T];
	triIn[0] = triLeft[0]*fftE[i+(NPOINTS_R-1)*NPOINTS_T]+
					triDiag[0]*fftE[i+0*NPOINTS_T]+
					triRight[0]*fftE[i+(0+1)*NPOINTS_T];
	triIn[NPOINTS_R-1] = triLeft[NPOINTS_R-1]*fftE[i+(NPOINTS_R-2)*NPOINTS_T]+
					     triDiag[NPOINTS_R-1]*fftE[i+(NPOINTS_R-1)*NPOINTS_T]+
					     triRight[NPOINTS_R-1]*fftE[i+0*NPOINTS_T];
}



void ApplyTransverseLaplacian(double step)
{
	int i,j; //i e a frequencia, j e o raio
	double omega,invJ;
	double complex D;
	
	
	for(i = 0;i < NPOINTS_T;i++)
		{
			if(i >= NPOINTS_T/2)
				omega = 2.*M_PI*(i-NPOINTS_T)/NPOINTS_T/DELTAT;
			else omega = 2.*M_PI*i/NPOINTS_T/DELTAT;
			
			D = step/2.*I*lambdaZero/(2.*M_PI)/2.*(1+omega/omegaZero)/DELTAR/DELTAR; 
			
			triLeft[0]  = 0;
			triDiag[0]  = 1.0-D*4.;
			triRight[0] = 4.0*D;
			
			for(j = 1;j < NPOINTS_R;j++)
			{
				invJ = 1./(double)j;
								
				triLeft[j]  = D*(1.-0.5*invJ);
				triDiag[j]  = 1.0-D*2.;
				triRight[j] = D*(1.+0.5*invJ);
			}
			triRight[NPOINTS_R-1] = 0;
								
			CalcTriIn(i);
			
			triLeft[0]  = 0;
			triDiag[0]  = 1.0+D*4.;
			triRight[0] = -4.0*D;
			
			for(j = 1;j < NPOINTS_R;j++)
			{
				invJ = 1./(double)j;
								
				triLeft[j]  = -1*D*(1.-0.5*invJ);
				triDiag[j]  = 1.0+D*2.;
				triRight[j] = -1*D*(1.+0.5*invJ);
			}
			triRight[NPOINTS_R-1] = 0;
			
			Tridiagonal();
			
			for(j = 0;j < NPOINTS_R;j++)
			{
				fftE[i+j*NPOINTS_T] = triOut[j];
			}
		}	
}
