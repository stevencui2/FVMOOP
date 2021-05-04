#ifndef FINITEMATRIX_H
#define FINITEMATRXI_H

#include <vector>
#include <iostream>
#include "forAllOperations.h"
using namespace std;


class FiniteMatrix {
	public: FiniteMatrix();
		~FiniteMatrix();

typedef vector<vector<FiniteMatrix>> finiteMat;

double value, aevalue, awvalue, asvalue, anvalue,svalue;

void print2dMat(finiteMat&);
void print2dsource(finiteMat&);

friend FiniteMatrix::finiteMat operator+(const FiniteMatrix::finiteMat&, const FiniteMatrix::finiteMat&);
friend FiniteMatrix::finiteMat operator-(const FiniteMatrix::finiteMat&, const FiniteMatrix::finiteMat&);
friend FiniteMatrix::finiteMat operator&&(const FiniteMatrix::finiteMat&, const FiniteMatrix::finiteMat&);
friend FiniteMatrix::finiteMat operator*(const double, const FiniteMatrix::finiteMat&);

};

#endif
