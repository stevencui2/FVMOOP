#ifndef EQUATION_H
#define EQUATION_H

#include <string>
#include <vector>
#include "Grid.h"
#include "Fields.h"
#include "FiniteMatrix.h"
#include <cmath>

class Equation
{
public:
  Equation(const FiniteMatrix::finiteMat &);
  virtual ~Equation();

  typedef vector<FiniteMatrix> svector1d;
  typedef FiniteMatrix::finiteMat svector;

  void relax(Fields::vectorfields &);
  void resetEqn();
  void noWallShearXBoundaryCondition(Fields::vectorfields &);
  void noWallShearYBoundaryCondition(Fields::vectorfields &);
  void assembleEquation();

  Fields::vectorfields solve(Fields::vectorfields &, FiniteMatrix::finiteMat &, Solution &, int&);

  double value;
  double Residual;
  double RSM;
  double RESOR;

  double URF, rURF;
  string EqnName;
  double SOR;

  FiniteMatrix::finiteMat APinitial;
  FiniteMatrix::finiteMat AP;
  FiniteMatrix::finiteMat AE;
  FiniteMatrix::finiteMat AW;
  FiniteMatrix::finiteMat AN;
  FiniteMatrix::finiteMat AS;
  FiniteMatrix::finiteMat rAP;
  FiniteMatrix::finiteMat APU;
  FiniteMatrix::finiteMat sourceInitial;
  FiniteMatrix::finiteMat sourceB;
  FiniteMatrix::finiteMat sourceRelaxed;
  FiniteMatrix::finiteMat sourceFinial;

private:
  int NI, NJ, NIM, NJM, Literation;
  svector UE, UN, LW, LS, LPR, RES;
};
#endif
