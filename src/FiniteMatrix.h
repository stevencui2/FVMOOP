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

};

#endif
