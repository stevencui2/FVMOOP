#ifndef FIELDS_H
#define FIELDS_H
#include "forAllOperations.h"
#include "Grid.h"
#include "Solution.h"
#include "string"

class Fields{

	public: Fields();
		Fields(int&, int&);
	virtual ~Fields();

	typedef vector<Fields> vec1dfields;
	typedef vector<vector<vec1dfields>> vectorfields;

	double value;
	int NI,NJ,NIM,NJM;
	double X,XC,Y,YC,FXE,FXP,FYN,FYP,DXPtoE,DYPtoN;
	double Se,Sn,visc,density,volume;

	
	
	void getGridInfoPassed(Fields::vectorfields&, Grid&, Solution&);
	void setVectorFieldGridFeatures();
	void copyInternalField(Fields::vectorfields&, Fields::vectorfields&);
	void initalizeFields(Fields::vectorfields&,double&);
	//boundary conditions (1) -which modifieds the value of Fields (not matrices)
	void inletboundaryCondition(Fields::vectorfields&, string&, double&);
	void linearextrapolateCondition(Fields::vectorfields&);

};

#endif
