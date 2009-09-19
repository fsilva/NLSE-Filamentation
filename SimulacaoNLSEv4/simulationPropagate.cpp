
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
#include <netcdf/netcdfcpp.h>
//#include <clapack.h>
using namespace std;

#include "library.h"
#include "simulationVariables.h"

//Prototype
void ApplyTransverseLaplacian(double step);


/////////////////////////////////////////////////////////////
// Helper functions for nonlinear contribution calculation
//
//
inline double multiply_sigmaK(double Emag_square,int n)
{
	return sigmaK*pow(Emag_square,n);
       // if we just did sigmaK*pow(Emag_square,n) we would be multiplying a small number by a bi number
       // this way my multiply two numbers closer to one
       double x = Emag_square * 1e-15;
       double sigmaK_higher = sigmaK*pow(1e15,n);
       return sigmaK_higher*pow(x,n);       
}

inline double drho(double current_rho,double Emag_square)
{
	return 0;
	return (rho_neutral-current_rho)*multiply_sigmaK(Emag_square,n)+
	                     sigma/Ui*current_rho*Emag_square;
}


inline double complex ionization(double current_rho,double Emag_square)
{
	return -sigma/2*current_rho//-0*1.*I*M_PI/(lambdaZero*n*n*rhoC)*current_rho
			  -Ui/2.*multiply_sigmaK(Emag_square,K-1)*(rho_neutral-current_rho);
}

inline double cmag_square(double complex number)
{
	double real = creal(number),imag = cimag(number);
	return real*real+imag*imag;
}

inline double complex polarization_terms(int j,double current_rho, double complex E)
{
	return E*omegaZero/c*n2*cmag_square(E)*I;
//        return E*(-absorptionCalc[j]+I*omegaZero/c*n2*cmag_square(E)+ionization(current_rho,cmag_square(E)));
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
	double complex	tmpE1,tmpE2,tmpE3,tmpE4;
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
//#pragma omp parallel for 
	for(j = 0;j < NPOINTS_R;j++)
	{
		for(i = 0;i < NPOINTS_T;i++)
		{
		   
		    /*tmpE1  = polarization_terms(0,0,E[i+j*NPOINTS_T]);
		    tmpE2  = polarization_terms(0,0,E[i+j*NPOINTS_T]+tmpE1*step*0.5);
		    tmpE3  = polarization_terms(0,0,E[i+j*NPOINTS_T]+tmpE2*step*0.5);
		    tmpE4  = polarization_terms(0,0,E[i+j*NPOINTS_T]+tmpE3*step);
		    E[i+j*NPOINTS_T]  += step/6.*(tmpE1+2*tmpE2+2*tmpE3+tmpE4);*/
		E[i+j*NPOINTS_T] *= cexp(step*omegaZero/c*n2*cmag_square(E[i+j*NPOINTS_T])*I);	
		}
		/*i = 0;
		
		//rho[0+j*NPOINTS_T] = 0;
	
		//tmpRho1 = drho(rho[0+j*NPOINTS_T]  ,cmag_square(E[i+j*NPOINTS_T]));
		tmpE1   = polarization_terms(j,rho[i-1+j*NPOINTS_T],E[i+j*NPOINTS_T]);

                //Second step
		//tmpRho2 = drho(rho[0+j*NPOINTS_T]+DELTAT*tmpRho1,cmag_square(E[i+j*NPOINTS_T]+tmpE1*DELTAT));
		tmpE2   = polarization_terms(j,rho[0+j*NPOINTS_T]+DELTAT*tmpRho1,E[i+j*NPOINTS_T]+tmpE1*step); 
		//Final
	        E[i+j*NPOINTS_T] = E[i+j*NPOINTS_T] + step/2.*(tmpE1+tmpE2);	

		for(i = 1;i < NPOINTS_T;i++)
		{		
   			//Runge Kutta 4 (for rho and E at the same time)
			//
			//First step
		//	tmpRho1 = drho(rho[i-1+j*NPOINTS_T]  ,cmag_square(E[i+j*NPOINTS_T]));
			tmpE1   = polarization_terms(j,rho[i-1+j*NPOINTS_T],E[i+j*NPOINTS_T]);

                        //Second step
		//	tmpRho2 = drho(rho[i-1+j*NPOINTS_T]+DELTAT*0.5*tmpRho1,cmag_square(E[i+j*NPOINTS_T]+tmpE1*DELTAT*0.5));
			tmpE2   = polarization_terms(j,rho[i-1+j*NPOINTS_T]+0.5*DELTAT*tmpRho1,E[i+j*NPOINTS_T]+tmpE1*step*0.5); 

                        //Third step
		//	tmpRho3 = drho(rho[i-1+j*NPOINTS_T]+0.5*DELTAT*tmpRho2,cmag_square(E[i+j*NPOINTS_T]+tmpE2*DELTAT*0.5));
			tmpE3   = polarization_terms(j,rho[i-1+j*NPOINTS_T]+0.5*DELTAT*tmpRho2,E[i+j*NPOINTS_T]+tmpE2*step*0.5); 

                        //Fourth step
		//	tmpRho4 = drho(rho[i-1+j*NPOINTS_T]+DELTAT*tmpRho3,cmag_square(E[i+j*NPOINTS_T]+tmpE3*DELTAT));
			tmpE4   = polarization_terms(j,rho[i-1+j*NPOINTS_T]+DELTAT*tmpRho3,E[i+j*NPOINTS_T]+tmpE3*step); 
			
			//Final
		  //      rho[i+j*NPOINTS_T] = rho[i-1+j*NPOINTS_T] + DELTAT/6.*(tmpRho1+2*tmpRho2+2*tmpRho3+tmpRho4);	

		        E[i+j*NPOINTS_T] = E[i+j*NPOINTS_T] + step/6.*(tmpE1+2*tmpE2+2*tmpE3+tmpE4);

		}	*/
	}

//1) Calculate phi(f,z0)=F[phi(t,z0)]
// execute fftw forward: fftPhi = F[phi]

//#pragma omp parallel for 
	for(i = 0;i < NPOINTS_R;i++)
		fftw_execute(forward[i]);

//2) Apply dispersive propagation z0->z0+h

//#pragma omp parallel for 
	for(j = 0;j < NPOINTS_R;j++)
		for(i = 0;i < NPOINTS_T;i++)
		{
			if(i >= NPOINTS_T/2)
				omega = 2*M_PI*(i-NPOINTS_T)/NPOINTS_T/DELTAT;
			else omega = 2*M_PI*i/NPOINTS_T/DELTAT;
			fftE[i+j*NPOINTS_T] = fftE[i+j*NPOINTS_T]*cexp(-GVD/2.*omega*omega*step*I); 
		}
		
	ApplyTransverseLaplacian(step);

//#pragma omp parallel for 
	for(i = 0;i < NPOINTS_R;i++)
		fftw_execute(backward[i]); 

//4.7.1) Renormalize (fftw multiplies data per sqrt(N))

#pragma omp parallel for 
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
	
//#pragma omp parallel for	
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
