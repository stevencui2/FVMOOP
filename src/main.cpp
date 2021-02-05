#include <iostream>
#include "forAllOperations.h"
#include "Fields.h"

int main(int argc, char **argv){
	std::cout<<"Hello CFD"<<std::endl;
	int N=7,M=7;double length=1.0;

	Grid mygrid_(N,M,length);

	std::cout<<"mygrid dx is :"<<mygrid_.retdx()<<std::endl;
	for(unsigned int i=0;i<mygrid_.X.size();i++)
		std::cout<<mygrid_.X[i]<<',';
	std::cout<<std::endl;
	for(unsigned int i=0;i<mygrid_.XC.size();i++)
		std::cout<<mygrid_.XC[i]<<',';
	std::cout<<std::endl;
	for(unsigned int i=0;i<mygrid_.FX.size();i++)
		std::cout<<mygrid_.FX[i]<<',';

	return 0;
}
