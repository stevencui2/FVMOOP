#include "Grid.h"

Grid::Grid(int& N_,int& M_,double& length_):
	N(N_),M(M_),NI(N_+2),NJ(M_+2),NIM(NI-1),NJM(NJ-1),length(length_), 
	X(N_+2,0.0),Y(M_+2,0.0),XC(N_+2,0.0),YC(M_+2,0.0),FX(N_+2,0.0),FY(M_+2,0.0)
{
	dx=length/N;
	dy=length/M;
	setX(X);
	setY(Y);
	setXC(X,XC);
	setYC(Y,YC);
	setFX(X,XC,FX);
	setFY(Y,YC,FY);
}

Grid::~Grid(){
}
void Grid::setX(vector<double>& vecX){
	for(unsigned int i=1;i<vecX.size();i++)
		vecX[i]=vecX[i-1]+dx;

	vecX[vecX.size()-1]=vecX[vecX.size()-2];
}
void Grid::setY(vector<double>&vecY){
	for(unsigned int i=1;i<vecY.size();i++)
		vecY[i]=vecY[i-1]+dy;

	vecY[vecY.size()-1]=vecY[vecY.size()-2];
}
void Grid::setXC(vector<double>&vecX,vector<double>&vecXC){
	for(unsigned int i=1;i<vecXC.size();i++)
		vecXC[i]=(vecX[i]+vecX[i-1])*0.5;
}
void Grid::setYC(vector<double>&vecY,vector<double>&vecYC){
	for(unsigned int i=1;i<vecYC.size();i++)
		vecYC[i]=(vecY[i]+vecY[i-1])*0.5;
}
void Grid::setFX(vector<double>&vecX,vector<double>&vecXC,vector<double>&vecFX){
	for(unsigned int i=0;i<vecFX.size();i++)
		vecFX[i]=(vecX[i]-vecXC[i])/(vecXC[i+1]-vecXC[i]);
}
void Grid::setFY(vector<double>&vecY,vector<double>&vecYC,vector<double>&vecFY){
for(unsigned int i=0;i<vecFY.size();i++)
		vecFY[i]=(vecY[i]-vecYC[i])/(vecYC[i+1]-vecYC[i]);
}
