
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <math.h>
#include <complex.h>
using namespace std;

#include "library.h"



//init array (with 'size') to gaussian centered at point 'center' with 'fwhm',
//            'max' is the value at the peak
void	GaussianSlow(double complex* array,int size,double complex max,int center,double complex fwhm)
{
	double complex x;
	fwhm /= 2.354820;
	for(int i = 0;i < size;i++)
	{
		x = double(i-center)/fwhm;
		array[i] = max*cexp(-1/2.*x*x);
		x = double(i-center/2.0)/fwhm;
	}
}



//*************************************************************************
//*************************************************************************
//*************************************************************************
//*************************************************************************
//*************************************************************************

ConfigurationParser::ConfigurationParser()
{
	Names.clear();
	Values.clear();
	listNames.clear();
	listValues.clear();
	iterNames.clear();
	MaxValues.clear();
	MinValues.clear();
	StepValues.clear();
	CurrentConfiguration = 0;
	NumberOfConfigurations = 0;
}
	
int 		ConfigurationParser::Process(const char* fileName)
{
	fstream file;

	file.open(fileName,fstream::in);

	if(file.fail()) //error occurred
		return 1;

	NumberOfConfigurations = 1;

	string str,name,value,min,max,step,listVal;
	unsigned int i,j,k;
	double  vmin,vmax,vstep;

	while(!file.eof())
	{
		getline(file,str);
	//	cout << str <<endl;

//remove whitespace,tabs
		i =  str.find_first_not_of(" 	");
	//	cout <<"test2"<<endl;
		if(i == string::npos || i > str.length())
		{
			continue; //Characters not found, ignore line
		}

		if(str[i] == '#' || str[i]=='\n') //it is a comment or empty line
		        continue;	
//find = or whitespace	
		j = str.find_first_of(" =	",i);
		if(j == string::npos)
			continue; //Characters not found, ignore line

//remove whitespace and =
		k =  str.find_first_not_of(" 	=",j);
		if(k == string::npos)
			continue; //Characters not found, ignore line
                

		name = str.substr(i,j-i);
		value = str.substr(k,str.length()-k);
	
		if(value[0] == '[') //iteration parameter
		{ 
			i = 1;
			j = value.find_first_of(",");
			min = value.substr(i,j-i);
			
			i = j+1;
			j = value.find_first_of(",",i);
			max = value.substr(i,j-i);

			i = j+1;
			j = value.find_first_of("]",i);
			step = value.substr(i,j-i);

			vmin = atof(min.c_str());
			vmax = atof(max.c_str());
			vstep = atof(step.c_str());

			iterNames.push_back(name);
			MinValues.push_back(vmin);
			MaxValues.push_back(vmax);
			StepValues.push_back(vstep);
			NumValues.push_back(1+floor((vmax-vmin)/vstep));
			
			NumberOfConfigurations *= 1+floor((vmax-vmin)/vstep);

		}else if(value[0] == '(') //list parameter
		{
			bool end = false;
			i = 1;
			k = 0;
			listValues.push_back(vector<double>());
			do{
				j = value.find_first_of(",",i);
				if(j > value.length())
		       	        {
					j = value.find_first_of(")",i);	
					if(j == i || j > value.length())
					{
						end = true;
						continue;
					}
				}
				listVal = value.substr(i,j-i);
				listValues.back().push_back(atof(listVal.c_str()));
				k++;
				i = j+1;

			}while(!end);

			listNames.push_back(name);
			if(k > 0)
				NumberOfConfigurations *= k;

		}else //fixed parameter
		{
			Names.push_back(name);
			Values.push_back(value);
		}
	}
	file.close();
	return 0; //OK
}

int		ConfigurationParser::GetNumberofConfigurations()
{
	return NumberOfConfigurations;
}

void		ConfigurationParser::SetCurrentConfiguration(int number)
{
	if(number > NumberOfConfigurations || number < 0)
		cout << "Error in ConfigurationParser::SetCurrentConfiguration - number has invalid value" << endl;
	else	
		CurrentConfiguration = number;
}
	
string 		ConfigurationParser::GetString(string name)
{
	unsigned int i;
	for(i = 0;i < Names.size();i++)
	{
		if(Names[i] == name)
			return Values[i];
	}
	throw string("Error: No configuration entry for parameter ") + name;
}

double 		ConfigurationParser::GetNumber(string name)
{
	unsigned int i;
	for(i = 0;i < Names.size();i++)
	{
		if(Names[i] == name)
			return atof(Values[i].c_str());
	}

	int index = CurrentConfiguration;


	for(i = 0;i < iterNames.size();i++)
	{
		//cout << " index : "<<index << endl;
		if(iterNames[i] == name)
			return MinValues[i] + StepValues[i]*int(index%int(NumValues[i]));

		index /= NumValues[i];
		//cout << "i: " << i << " index : "<<index << "    " << int(NumValues[i]) <<  endl;
	}

	for(i = 0;i < listNames.size();i++)
	{
		if(listNames[i] == name)
			return listValues[i][index%listValues[i].size()];
		if(listValues[i].size() > 0)
			index /= listValues[i].size();
		//cout << "i: " << i << " index : "<<index << endl;
	}

	throw string("Error: No configuration entry for parameter ") + name;
}

