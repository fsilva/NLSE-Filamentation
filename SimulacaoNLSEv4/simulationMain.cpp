
#define __MAINCPP__

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


///////////////////////////////////////////////////
//
//
//		Main
//
//
int main(int argc, char *argv[])
{
	cout << "***********************************************\n";
	cout << "    Non linear Envelope Equation Solver\n";
	cout << "***********************************************\n";


// Open Configuration file. Should be argv[1].
// if no argument assume 'default.conf'

	ConfigurationParser configParser;
	
//Process configuration file.
	
	if(argc == 1)
		configFileName = "default.conf";
	else
		configFileName = argv[1];
	cout << configFileName.c_str() << endl;

	if(configParser.Process(configFileName.c_str()))
	{	
		cout << "ERROR: Cannot open " << configFileName << "\nExiting.\n";
		return 0;
	}

	int  currentSim,numSims;

	numSims = configParser.GetNumberofConfigurations();

	for(currentSim = 0;currentSim < numSims;currentSim++)
	{

		configParser.SetCurrentConfiguration(currentSim);

	//Input parameter values
		try
		{
			zDistance	= configParser.GetNumber("zDistance");
			zOutputStep	= configParser.GetNumber("zOutputStep");
			Energy 		= configParser.GetNumber("Energy");
			Pmax 		= configParser.GetNumber("Pmax");
			pulseT   	= configParser.GetNumber("pulseT");
			spot    	= configParser.GetNumber("spot");
			curvature       = configParser.GetNumber("curvature");
			GVD 		= configParser.GetNumber("GVD");
			absorption 	= configParser.GetNumber("absorption");
			n2 		    = configParser.GetNumber("n2");
			lambdaZero 	= configParser.GetNumber("lambdaZero");
			desiredError 	= configParser.GetNumber("desiredError");
			boundaryRatio	= configParser.GetNumber("boundaryRatio");
			boundarySlope	= configParser.GetNumber("boundarySlope");
			boundaryI0		= configParser.GetNumber("boundaryI0");
			n               = configParser.GetNumber("n");
			sigma           = configParser.GetNumber("sigma");
			Ui              = configParser.GetNumber("Ui");
			sigmaK          = configParser.GetNumber("sigmaK");
			K               = configParser.GetNumber("K");
			rho_neutral     = configParser.GetNumber("rho_neutral");
			sigmaK          = configParser.GetNumber("sigmaK");
			tau_r           = configParser.GetNumber("tau_r");
			tau_c           = configParser.GetNumber("tau_c");
			p 				= configParser.GetNumber("p");
			
			omegaZero = 2*M_PI*lightVelocity/lambdaZero;
			rhoC	  = epslon0*m_electron*omegaZero*omegaZero/q_electron/q_electron;
			
			n2     *= p;
			sigma  *= p;
			sigmaK *= p; //TODO: ver se isto faz algum sentido
			
		}catch(string error) //some parameter was not found
		{
			cout << error << endl;
			cout << "Exiting.\n";
			return 0;
		}

		// Run Simulation with the current parameters
		RunSimulation();

		cout << "Simulation Done - " << currentSim+1 << "/" << numSims << "  -   " <<  time(0l)-startTime << " Seconds Elapsed"<<endl;
		sleep(1);
	}

//END. Success.
	cout << "Program Ended Successfully   -   " <<  time(0l)-startTime << " Seconds Elapsed"<<endl;

	return 0;
}

