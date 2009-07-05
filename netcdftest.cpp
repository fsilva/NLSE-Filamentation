
#include <iostream>
#include <complex.h>
#include <string>
#include <netcdfcpp.h>
using namespace std;

#define NPOINTS_T 	(1024)
#define NPOINTS_R 	(1024)
#define NZSTEPS          2

double 		E[NPOINTS_T*NPOINTS_R];          //
double 		fftE[NPOINTS_T*NPOINTS_R];       // envelopes' transform


int main()
{
        NcVar *dataSpectrum=0,*dataImpulse=0;

	string filename = "./simulation.nc";
	NcFile dataFile(filename.c_str(), NcFile::Replace);//,0l,0,NcFile::Netcdf4Classic);
	int 	i,j;
	//double 	distance,zStep,uC,uF,EReal,EImag,error,lastDistance,zNextStep;

	if(!dataFile.is_valid())
	{
		cout << "Error: Could not open simulation.nc.";
		return 0;
	}
	else
	{	
		NcDim* rDim = dataFile.add_dim("r", NPOINTS_R);
		NcDim* tDim = dataFile.add_dim("time", NPOINTS_T);
		NcDim* fDim = dataFile.add_dim("omega", NPOINTS_T);
		

                //old way
		/*
                NcDim* zDim = dataFile.add_dim("z",NZSTEPS); 
                dataSpectrum = dataFile.add_var("Spectrum", ncDouble, fDim, rDim, zDim);		
		dataImpulse  = dataFile.add_var("Impulse", ncDouble, tDim, rDim, zDim);*/
                //new way
                NcDim* zDim = dataFile.add_dim("z"); 
	//	dataSpectrum = dataFile.add_var("Spectrum", ncDouble, zDim, fDim, rDim);		
		dataImpulse  = dataFile.add_var("Impulse", ncDouble, zDim, tDim, rDim);

         }

         for(j = 0;j < NZSTEPS;j++)
{
         cout << j << " ";
         for(int o = 0;o < NPOINTS_R*NPOINTS_T;o++)
         {
             E[o] = fftE[o] = j;
         }
         cout <<j<<"   "<< E[2];

         //old way
         /*for(i = 0;i < NPOINTS_R;i++)
			{
				dataImpulse->set_cur(0,i,j);
				dataImpulse->put(reinterpret_cast<double*>(&E[i*NPOINTS_T]),NPOINTS_T,1,1);
				dataSpectrum->set_cur(0,i,j);
				dataSpectrum->put(reinterpret_cast<double*>(fftE+NPOINTS_T*i),NPOINTS_T,1,1);	
			}*/
         
         //new way
         if(!dataImpulse->put_rec(E,j))
              cout << "error" << endl;
	 //if(!dataSpectrum->put(fftE,j))	
         //     cout << "error" << endl;
	 cout <<"written " << j << endl;
         dataFile.sync();
}
       return 0;
}