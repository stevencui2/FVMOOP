#include <iostream>
#include "forAllOperations.h"
#include "Fields.h"
#include "FiniteMatrix.h"

int main(int argc, char **argv){
	std::cout<<"Hello CFD"<<std::endl;
	#include "createGridFieldsSol.h"
	#include "createFiniteMatrices.h"


	fieldsOper.initializeInternalFields(P,0.2);
	fieldsOper.linearextrapolateCondition(P,mygrid_.FX,mygrid_.FY,north);
	fieldsOper.linearextrapolateCondition(P,mygrid_.FX,mygrid_.FY,south);
	fieldsOper.linearextrapolateCondition(P,mygrid_.FX,mygrid_.FY,east);
	fieldsOper.linearextrapolateCondition(P,mygrid_.FX,mygrid_.FY,west);
	//fieldsOper.print2dmat(P);

FiniteMatrix::finiteMat AE2(AE);
 forAll(AE){
	 AE[i][j].value=2.0;
	 AE2[i][j].value=3.0;
 }

const double two=2.0;
 FiniteMatrix::finiteMat AE3(two*AE);

 finiteobj.print2dMat(AE3);
	//std::cout<<"mygrid dx is :"<<mygrid_.retdx()<<std::endl;
	//for(unsigned int i=0;i<mygrid_.X.size();i++)
	//	std::cout<<mygrid_.X[i]<<',';
	//std::cout<<std::endl;
	//for(unsigned int i=0;i<mygrid_.XC.size();i++)
	//	std::cout<<mygrid_.XC[i]<<',';
	//std::cout<<std::endl;
	//for(unsigned int i=0;i<mygrid_.FX.size();i++)
	//	std::cout<<mygrid_.FX[i]<<',';

	return 0;
}
