#include "FiniteMatrix.h"
#include <iomanip>

FiniteMatrix::FiniteMatrix():
value(0.0), aevalue(0.0), awvalue(0.0), asvalue(0.0), anvalue(0.0),svalue(0.0)
{

}

FiniteMatrix::~FiniteMatrix(){

}
void FiniteMatrix::print2dMat(finiteMat&)
{
  for(unsigned int i=0;i<vec.size();i++){
		for(unsigned int j=0;j<vec[i].size();j++){
			std::cout<<sdt::setprecision(3)<<vec[i][j].value<<',';
		}
		cout<<endl;
	}
}
