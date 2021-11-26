#ifndef FINITEMATRIX_H
#define FINITEMATRIX_H

#include <vector>
#include <iostream>
#include "forAllOperations.h"
#include "Fields.h"
using namespace std;

class FiniteMatrix
{
public:
	FiniteMatrix();
	~FiniteMatrix();

	typedef vector<vector<FiniteMatrix>> finiteMat;

	double value, aevalue, awvalue, asvalue, anvalue, svalue;

	FiniteMatrix::finiteMat interpolatedFieldEast(FiniteMatrix::finiteMat &, Grid &);
	FiniteMatrix::finiteMat interpolatedFieldNorth(FiniteMatrix::finiteMat &, Grid &);
	Fields::vectorfields correctFaceVelocityEast(Fields::vectorfields &, Fields::vectorfields &, Fields::vectorfields &,FiniteMatrix::finiteMat&, Grid&);
	Fields::vectorfields correctFaceVelocityNorth(Fields::vectorfields &, Fields::vectorfields &, Fields::vectorfields &,FiniteMatrix::finiteMat&,Grid&);

	void correctEastMassFluxes(Fields::vectorfields&, Fields::vectorfields&, FiniteMatrix::finiteMat&);
	void correctNorthMassFluxes(Fields::vectorfields&, Fields::vectorfields&, FiniteMatrix::finiteMat&);
	void print2dMat(finiteMat &);

	void print2dsource(finiteMat &);

	friend FiniteMatrix::finiteMat operator+(const FiniteMatrix::finiteMat &, const FiniteMatrix::finiteMat &);
	friend FiniteMatrix::finiteMat operator-(const FiniteMatrix::finiteMat &, const FiniteMatrix::finiteMat &);
	friend FiniteMatrix::finiteMat operator&&(const FiniteMatrix::finiteMat &, const FiniteMatrix::finiteMat &);
	friend FiniteMatrix::finiteMat operator*(const double, const FiniteMatrix::finiteMat &);
};

#endif
