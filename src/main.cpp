#include <iostream>
#include "forAllOperations.h"
#include "Fields.h"

int main(int argc, char **argv){
	std::cout<<"Hello CFD"<<std::endl;
	#include "createGridFieldsSol.h"

	fieldsOper.initializeInternalFields(P,0.2);
	fieldsOper.linearextrapolateCondition(P,mygrid_.FX,mygrid_.FY,north);
	fieldsOper.linearextrapolateCondition(P,mygrid_.FX,mygrid_.FY,south);
	fieldsOper.linearextrapolateCondition(P,mygrid_.FX,mygrid_.FY,east);
	fieldsOper.linearextrapolateCondition(P,mygrid_.FX,mygrid_.FY,west);
	fieldsOper.print2dmat(P);

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
