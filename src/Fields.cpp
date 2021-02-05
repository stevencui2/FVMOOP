#include "Fields.h"

Fields::Fields():
value(0.0)
{
}

Fields::Fields(int& NI_ , int& NJ_):
value(0.0),NI(NI_),NJ(NJ_),NIM(NI-1),NJM(NJ-1)
{
}

Fields::~Fields(){
}
void Fields::getGridInfoPassed(Fields::vectorfields& f, Grid& grid_, Solution& sol_){

	forAll(f){
		f[i][j].X=grid_.X[i];
		f[i][j].Y=grid_.Y[j];
		f[i][j].XC=grid_.XC[i];
		f[i][j].YC=grid_.YC[j];
		f[i][j].FXE=grid_.FX[i];
		f[i][j].FYN=grid_.FY[j];
		f[i][j].FXP=1.0-grid_.FX[i];
		f[i][j].FYP=1.0-grid_.FY[j];
		f[i][j].visc=sol_.visc;
		f[i][j].density=sol_.density;
	}

	forAllInternal(f){
	
	f[i][j].DXPtoE=grid_.XC[i+1]-grid_.XC[i];
	f[i][j].DYPtoN=grid_.YC[j+1]-grid_.YC[j];
	f[i][j].Se=grid_.Y[j]-grid_.Y[j-1];
	f[i][j].Sn=grid_.X[i]-grid_.X[i-1];
	f[i][i].volume=f[i][j].Se*f[i][j].Sn;
	}



}
void Fields::setVectorFieldGridFeatures(){
}
void Fields::copyInternalField(Fields::vectorfields&, Fields::vectorfields&){
}
void Fields::initalizeFields(Fields::vectorfields&,double&){
}
//boundary conditions (1) -which modifieds the value of Fields (not matrices)
void Fields::inletboundaryCondition(Fields::vectorfields&, string&, double&){
}
void Fields::linearextrapolateCondition(Fields::vectorfields&){
}

