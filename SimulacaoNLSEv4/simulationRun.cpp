
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
using namespace std;

#include "library.h"
#include "simulationVariables.h"

void	OpenDirectory();
void	SetInitialConditions();
void	StartSimulation();
void	Cleanup();


void RunSimulation()
{
	startTime = time (0l);
	
	OpenDirectory(); //2)
	
	SetInitialConditions(); //3)

	StartSimulation(); //4)

	Cleanup(); //5)
}

//1) was opening the configuration file, in main()















void	OpenDirectory()
{
//2) 
//   Open Directory with name based on time and date and copy config file there.
//	Save name in simDir
	string  command;
	time_t rawtime;

  	time ( &rawtime );
  	timeinfo = localtime ( &rawtime );

	std::stringstream out;
	out << setfill('0');
	out << "Results_" << setw(2) << timeinfo->tm_year-100 << "_" 
			  << setw(2) << timeinfo->tm_mon+1 << "_" 
			  << setw(2) << timeinfo->tm_mday  << "___" 
			  << setw(2) << timeinfo->tm_hour << "_" 
			  << setw(2) << timeinfo->tm_min << "_" 
			  << setw(2) << timeinfo->tm_sec; 

	simDir = out.str();
	
	command = "mkdir " + simDir;
	int ret = system(command.c_str());
	command = "cp " + configFileName +"   " + simDir +"/";
	ret = system(command.c_str());
}
















void	SetInitialConditions()
{
	int i,j;
//3) Initialization - FFTS and E intial conditions
//   Initialize FFTs  //TODO TODO TODO change E arrays to fftw_mallocatted to match SIMD instructions
//   Set initial conditions
//	Save them with netCDF to initial.nc
	cout << "Initializing FFTW......." << flush;
	for(i = 0; i < NPOINTS_R;i++)
	{
		forward[i]  = fftw_plan_dft_1d(NPOINTS_T, &E[i*NPOINTS_T], &fftE[i*NPOINTS_T], FFTW_FORWARD, FFTW_PATIENT);
		backward[i] = fftw_plan_dft_1d(NPOINTS_T, &fftE[i*NPOINTS_T], &E[i*NPOINTS_T], FFTW_BACKWARD, FFTW_PATIENT);
	}
	cout << "Done\n";


//  E initial conditions:

	if(Energy > 1e-15)  //If energy != 0, set right Pmax
	{
		//TODO!
	}

        double w_l,R_l,z_l,z0,w0;
        double complex invQ1,invQ2;
        double r,sigmaR,sigmaT,t;

        sigmaR = w_l = spot/2;
        R_l = -curvature;
        //calculate z_l,z0,w0
	z_l = 1/R_l;
        invQ1 = 1./R_l-1j*lambdaZero/M_PI/(w_l*w_l);
       
        invQ2 = invQ1/(1+z_l*invQ1);

        w0 = sqrt(-1./cimag(invQ2)/M_PI*lambdaZero);

        z0 = M_PI*w0*w0/lambdaZero;
        //cout << "w0=" << w0 << "   z0=" << z0;
	

	
	//sigmaR = spot/2; //spot is the diameter at 1/e^2
	sigmaT = pulseT/1.17741;
	for(i = 0; i < NPOINTS_R;i++)
	{
		r = (double)(i)*(double)DELTAR;
		
		for(j = 0; j < NPOINTS_T;j++)
		{
			t = (double)(j-NPOINTS_T/2)*(double)DELTAT;
			
			initialE[i*NPOINTS_T+j] = sqrt(2*Pmax/M_PI/sigmaR/sigmaR)*cexp( //I*t*t*1000e-30+ //TODO:remove this or add to config file
											-r*r/w_l/w_l//sigmaR/sigmaR
											-t*t/sigmaT/sigmaT-I*2.*M_PI/lambdaZero*r*r/2/R_l-I*atan(z_l/z0)); 
		}
	}

	//copy initial conditions into E
	memcpy(E,initialE,sizeof(double complex)*NPOINTS_T*NPOINTS_R);
	
// Precalculate absorption coeficients = f(r)
//    to allow 'fake' absorptive boundary layers
	double r_min, r_max,alpha;
	r_max = NPOINTS_R*DELTAR;
	r_min = (1.0-boundaryRatio)*r_max;
	for(i = 0; i < NPOINTS_R;i++)
	{
		r = i*DELTAR;
		if(r < r_min)
			alpha = absorption;
		else 
			alpha = absorption+boundaryI0*exp(boundarySlope*((r_max-r_min)/(r_max-r)-1))-1;
		absorptionCalc[i] = alpha;
		//debug : cout << i << " r=" << r << " r_min=" << r_min <<" alpha="<<alpha << "   " << (r_max-r_min)/(r_max-r) << "\n";
	}
}


























void StartSimulation()
{
//4) 
//	Start simulation - open simulation.nc and set dims and data
//
	NcVar *dataSpectrum=0,*dataImpulse=0,*dataZStep=0,*dataError=0,*dataRho;

	string filename = simDir + "/simulation.nc";
	NcFile dataFile(filename.c_str(), NcFile::Replace);
	int 	i,j,outputIndex,last_time;
	double 	distance,zStep,uC,uF,EReal,EImag,error,lastDistance,zNextStep;

	if(!dataFile.is_valid())
	{
		cout << "Error: Could not open simulation.nc.";
		return;
	}
	else
	{	
		NcDim* rDim = dataFile.add_dim("r", NPOINTS_R);
		NcDim* tDim = dataFile.add_dim("time", NPOINTS_T*2);
                NcDim* tDim_rho = dataFile.add_dim("time_rho", NPOINTS_T);
		NcDim* fDim = dataFile.add_dim("omega", NPOINTS_T*2);
		//NcDim* zDim = dataFile.add_dim("z",int(zDistance/zOutputStep));
                NcDim* zDim = dataFile.add_dim("z");  

		dataSpectrum = dataFile.add_var("Spectrum", ncDouble, zDim, rDim, fDim);		
		dataImpulse  = dataFile.add_var("Impulse", ncDouble, zDim, rDim, tDim);

		dataRho  = dataFile.add_var("ElectronDensity", ncDouble, zDim, rDim, tDim_rho);

		dataZStep  = dataFile.add_var("zStep", ncDouble, zDim);
		dataError  = dataFile.add_var("Error", ncDouble, zDim);

		dataFile.add_att("version",VERSION);

		dataFile.add_att("date", asctime(timeinfo));
		dataFile.add_att("deltaT",(double)DELTAT);
		dataFile.add_att("deltaR",(double)DELTAR);
		dataFile.add_att("Energy",Energy);
		dataFile.add_att("Pmax",Pmax);
		dataFile.add_att("pulseT",pulseT);
		dataFile.add_att("curvature",curvature);
		dataFile.add_att("timescale",TIMESCALE);
		dataFile.add_att("distancescale",DISTANCESCALE);
		dataFile.add_att("zDistance",zDistance);
		dataFile.add_att("GVD",GVD);
		dataFile.add_att("absorption",absorption);
		dataFile.add_att("lambda0",lambdaZero);
		dataFile.add_att("desiredError",desiredError);
		dataFile.add_att("spot",spot);
		dataFile.add_att("boundaryRatio",boundaryRatio);
		dataFile.add_att("boundarySlope",boundarySlope);
		dataFile.add_att("boundaryI0",boundaryI0);
		dataFile.add_att("n2",n2);
		dataFile.add_att("n",n);
		dataFile.add_att("sigma",sigma);
		dataFile.add_att("Ui",Ui);
		dataFile.add_att("sigmaK",sigmaK);
		dataFile.add_att("K",K);
		dataFile.add_att("rho_neutral",rho_neutral);
		dataFile.add_att("tau_r",tau_r);
		dataFile.add_att("tau_c",tau_c);
		dataFile.add_att("p",p);
	}

	distance = 0;
	lastDistance = 0;
	outputIndex = 0;
	zStep = zOutputStep;
	zNextStep = 0;

//4.0) Output first data point
//
	//dataZStep->set_cur(0.);
	dataZStep->put_rec(&zStep,0);	
	//dataError->set_cur(0.);
	dataError->put_rec(&error,0);	

	/*for(i = 0;i < NPOINTS_R;i++)
	{
		dataImpulse->set_cur(0,i,0);
		dataImpulse->put(reinterpret_cast<double*>(E+NPOINTS_T*i),NPOINTS_T*2,1,1);
	
		fftw_execute(forward[i]);
		
		dataSpectrum->set_cur(0,i,0);
		dataSpectrum->put(reinterpret_cast<double*>(fftE+NPOINTS_T*i),NPOINTS_T*2,1,1);	
		
		dataRho->set_cur(0,i,0);
		dataRho->put(reinterpret_cast<double*>(rho+NPOINTS_T*i),NPOINTS_T,1,1);	
	}*/
        dataImpulse->put_rec(reinterpret_cast<double*>(E),0);
	for(i = 0;i < NPOINTS_R;i++)
	     fftw_execute(forward[i]);
	dataSpectrum->put_rec(reinterpret_cast<double*>(fftE),0);	
	dataRho->put_rec(reinterpret_cast<double*>(rho),0);	
	
//	outputIndex++;

	float percentage = (float)(distance)/(float)(zDistance);
	cout << setw(5) << percentage*100  << " %      z = " << setw(10) << distance << " m " << endl;
	last_time = time(0l);

	for(distance = 0;distance < zDistance;)
	{ cout << "         Current Step Size = " << zStep << "     z=" << distance << endl;
//4.1)
//    Adaptive Step Size trial propagations

//4.1.1) Do 1x zStep propagation
		memcpy(E,initialE,sizeof(double complex)*NPOINTS_T*NPOINTS_R);
		Propagate(zStep);
		memcpy(h2E,E,sizeof(double complex)*NPOINTS_T*NPOINTS_R);

//4.1.2) Do 2x zStep/2 propagation
		memcpy(E,initialE,sizeof(double complex)*NPOINTS_T*NPOINTS_R);
		Propagate(zStep*0.5);
		Propagate(zStep*0.5);

//4.3) 
//     Calculate relative local error and change step accordingly
		
//4.3.1) Calculate time integrals  //TODO FIX! make adaptive step size care about all R
		uF = 0;
		uC = 0;
		for(i = 0;i < NPOINTS_T;i++)
		{
			for(j = 0;j < NPOINTS_R;j++)
			{
				EReal = creal(E[i+j*NPOINTS_T]);
				EImag = cimag(E[i+j*NPOINTS_T]);
				uF += EReal*EReal+EImag*EImag; 
			}

			for(j = 0;j < NPOINTS_R;j++)
			{
				EReal = creal(h2E[i+j*NPOINTS_T]);
				EImag = cimag(h2E[i+j*NPOINTS_T]);
				uC += EReal*EReal+EImag*EImag; 
			}
		}
		uF = sqrt(uF);
		uC = sqrt(uC);
//4.3.2) Calculate Error, decide what to do

//		cout << "11zStep =" << zStep << "      outputIndex=" << outputIndex << "      distance=" << distance << "      lastDistance=" << lastDistance<<endl;
//		cout << setw(15) << uF << " " << setw(15) <<uC << endl;

		error =  fabs(uF-uC)/uF;

		if(fabs(zStep) < zOutputStep*1e-20)
		{	
			cout << "zStep converged to Zero. desiredError too low. Aborting." << endl;
			return;
		}

		if(isnan(error) || error > 2*desiredError)
		{
			zStep *= 0.5;
			//cout << "error>2*desirederror -> zStep= " << zStep << endl;
			continue;
		}

// Step Accepted
//4.4) Save data, output to console, increase distance
		distance += zStep;

		if(distance-lastDistance+1e-7 > zOutputStep && outputIndex < int(zDistance/zOutputStep))
		{
			lastDistance = distance;

			/*dataZStep->set_cur(outputIndex);
			dataZStep->put(&zStep,1);	
			dataError->set_cur(outputIndex);
			dataError->put(&error,1);	
			

			for(i = 0;i < NPOINTS_R;i++)
			{
				dataImpulse->set_cur(0,i,outputIndex);
				dataImpulse->put(reinterpret_cast<double*>(&E[i*NPOINTS_T]),NPOINTS_T*2,1,1);
				
				fftw_execute(forward[i]);
		
				dataSpectrum->set_cur(0,i,outputIndex);
				dataSpectrum->put(reinterpret_cast<double*>(fftE+NPOINTS_T*i),NPOINTS_T*2,1,1);	
				
				dataRho->set_cur(0,i,outputIndex);
				dataRho->put(reinterpret_cast<double*>(rho+NPOINTS_T*i),NPOINTS_T,1,1);	
			}*/
	                dataZStep->put_rec(&zStep,outputIndex);	
	                dataError->put_rec(&error,outputIndex);
                        dataImpulse->put_rec(reinterpret_cast<double*>(E),outputIndex);
                        for(i = 0;i < NPOINTS_R;i++)
	                    fftw_execute(forward[i]);
	                dataSpectrum->put_rec(reinterpret_cast<double*>(fftE),outputIndex);	
	                dataRho->put_rec(reinterpret_cast<double*>(rho),outputIndex);	

//			dataFile.sync();
			outputIndex++;
			//cout << outputIndex <<endl;

			float percentage = (float)(distance)/(float)(zDistance),dt = time(0l)-startTime;
			float remaining = dt/percentage-dt;
			last_time = time(0l);
			
			cout << setw(5) << percentage*100  << " %      z = " << setw(10) << distance << " m " << "   error = " << error;
		       if(remaining < 100)
		       	     cout  << "     Remaining: " << setw(4) << remaining << " s" << endl;
		       else if(remaining < 3600)
		       	     cout  << "     Remaining: " << setw(4) << remaining/60. << " min" << endl;
		       else if(remaining < 24*3600)
		       	     cout  << "     Remaining: " << setw(4) << remaining/3600. << " h" << endl;
		       else
		       	     cout  << "     Remaining: " << setw(4) << remaining/3600./24 << " days"  << endl;
		       
		}else if(time(0l)-last_time > 5) //output to console in case > 5 seconds have passed
		{
			last_time = time(0l);
			
			float percentage = (float)(distance)/(float)(zDistance),dt = time(0l)-startTime;
			float remaining = dt/percentage-dt;

			cout << setw(5) << percentage*100  << " %      z = " << setw(10) << distance << " m " << "   error = " << error;
		       if(remaining < 100)
		       	     cout  << "     Remaining: " << setw(4) << remaining << " s" << endl;
		       else if(remaining < 3600)
		       	     cout  << "     Remaining: " << setw(4) << remaining/60. << " min" << endl;
		       else if(remaining < 24*3600)
		       	     cout  << "     Remaining: " << setw(4) << remaining/3600. << " h" << endl;
		       else
		       	     cout  << "     Remaining: " << setw(4) << remaining/3600./24 << " days"  << endl;
				
		}

		if(zStep > zOutputStep)
			zStep = zOutputStep;

		if(distance+zStep > lastDistance+zOutputStep && outputIndex < int(zDistance/zOutputStep))
		{
			zStep = zOutputStep-(distance-lastDistance); 
			if(zStep < 0)
			{
				cout << "negative zStep - makes no sense. Aborting" << (distance-lastDistance)<< endl; 
				return;
			}	
		}

		if(error > desiredError)
		{
			zStep /= pow(2.,1./3.);
			//cout << "error>desiredError   -> zStep= " << zStep << endl;
		}else if(error < 0.5*desiredError && zStep < zOutputStep/pow(2.,1./3.))
		{
			zStep *= pow(2.,1./3.);
			//cout << "error<desiredError/2 -> zStep= " << zStep << endl;
		}
//cout << zStep << endl;
		//if(zStep < zOutputStep*1e-3)
			//zStep = zOutputStep*1e-3;	

//		zStep = 100;

		//cout << "22zStep =" << zStep  << "      outputIndex=" << outputIndex << "      distance=" << distance << "      lastDistance=" << lastDistance<<endl;

		memcpy(initialE,E,sizeof(double complex)*NPOINTS_T*NPOINTS_R);
	}

//Save last set of points, if not saved in the last cycle
	if(outputIndex < int(zDistance/zOutputStep))
	{
		/*dataZStep->set_cur(outputIndex);
		dataZStep->put(&zStep,1);	
		dataError->set_cur(outputIndex);
		dataError->put(&error,1);
		for(i = 0;i < NPOINTS_R;i++)
		{
			dataImpulse->set_cur(0,i,outputIndex);
			dataImpulse->put(reinterpret_cast<double*>(E+NPOINTS_T*i),NPOINTS_T*2,1,1);
	
			fftw_execute(forward[i]);
		
			dataSpectrum->set_cur(0,i,outputIndex);
			dataSpectrum->put(reinterpret_cast<double*>(fftE+NPOINTS_T*i),NPOINTS_T*2,1,1);	
			
			dataRho->set_cur(0,i,outputIndex);
			dataRho->put(reinterpret_cast<double*>(rho+NPOINTS_T*i),NPOINTS_T,1,1);	
		}*/
	        dataZStep->put_rec(&zStep,outputIndex);	
	        dataError->put_rec(&error,outputIndex);
                dataImpulse->put_rec(reinterpret_cast<double*>(E),outputIndex);
                for(i = 0;i < NPOINTS_R;i++)
	            fftw_execute(forward[i]);
	        dataSpectrum->put_rec(reinterpret_cast<double*>(fftE),outputIndex);	
	        dataRho->put_rec(reinterpret_cast<double*>(rho),outputIndex);	
	}
	dataFile.add_att("finished","yes");

}





























void	Cleanup()
{
	int i;
//5) 
//    Cleanup
	for(i = 0; i < NPOINTS_R;i++)
	{
		fftw_destroy_plan(forward[i]);
		fftw_destroy_plan(backward[i]);
	}

}
