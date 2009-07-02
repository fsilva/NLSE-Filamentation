
#ifndef LIBRARY__H
#define LIBRARY__H

#include <string>
#include <vector>
#include <complex.h>
using namespace std;

//init array (with 'size') to gaussian centered at point 'center' with 'hwhm', 'max' is the value at the peak
void	GaussianSlow(double complex* array,int size,double complex max,int center,double complex hwhm);

class ConfigurationParser
{
// fixed fields
	vector<string> 		Names,Values;
// list fields
	vector<string> 		listNames;
	vector<vector<double> > listValues;
// iteration fields
	vector<string> 		iterNames;
	vector<double> 		MaxValues;
	vector<double> 		MinValues;
	vector<double> 		StepValues;
	vector<double> 		NumValues;
// Variables
	int			CurrentConfiguration,
				NumberOfConfigurations;
public:	
	ConfigurationParser();
	
	int 		Process(const char* fileName);

	int		GetNumberofConfigurations();
	void		SetCurrentConfiguration(int number);
	
	string 		GetString(string name);		
	double 		GetNumber(string name);		
};


#endif
