#ifndef FIELDS_H
#define FIELDS_H
#include "forAllOperations.h"
#include "Grid.h"
#include "Solution.h"
#include <string>

class Fields{

	public: Fields();
		Fields(int&, int&);
	virtual ~Fields();

	typedef vector<Fields> vec1dfields;
	typedef vector<vec1dfields> vectorfields;

	void getGridInfoPassed(Fields::vectorfields&, Grid&, Solution&);
	void setVectorFieldGridFeatures();
	void initializeFields(Fields::vectorfields&,double);
	void initializeInternalFields(Fields::vectorfields&,double);

	void print2dmat(Fields::vectorfields&);
	void copyInternalField(Fields::vectorfields&, Fields::vectorfields&);
	//boundary conditions (1) -which modifieds the value of Fields (not matrices)
	void inletboundaryCondition(Fields::vectorfields&, string&, double);
	void linearextrapolateCondition(Fields::vectorfields&, vector<double>&, vector<double>&, string&);
	
	

	double value;
	int NI,NJ,NIM,NJM;
	double X,XC,Y,YC,FXE,FXP,FYN,FYP,DXPtoE,DYPtoN;
	double Se,Sn,visc,density,volume;

	
	

};

#endif
