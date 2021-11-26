#include "FiniteMatrix.h"
#include <iomanip>

FiniteMatrix::FiniteMatrix() : value(0.0), aevalue(0.0), awvalue(0.0), asvalue(0.0), anvalue(0.0), svalue(0.0)
{
}

FiniteMatrix::~FiniteMatrix()
{
}
FiniteMatrix::finiteMat operator+(const FiniteMatrix::finiteMat &lhs, const FiniteMatrix::finiteMat &rhs)
{
  FiniteMatrix::finiteMat results(lhs);
  forAllInternal(lhs)
  {
    results[i][j].value += rhs[i][j].value;
    results[i][j].aevalue += rhs[i][j].aevalue;
    results[i][j].awvalue += rhs[i][j].awvalue;
    results[i][j].anvalue += rhs[i][j].anvalue;
    results[i][j].asvalue += rhs[i][j].asvalue;
    results[i][j].svalue += rhs[i][j].svalue;
  }
  return results;
}
FiniteMatrix::finiteMat operator-(const FiniteMatrix::finiteMat &lhs, const FiniteMatrix::finiteMat &rhs)
{
  FiniteMatrix::finiteMat results(lhs);
  forAllInternal(lhs)
  {
    results[i][j].value -= rhs[i][j].value;
    results[i][j].aevalue -= rhs[i][j].aevalue;
    results[i][j].awvalue -= rhs[i][j].awvalue;
    results[i][j].anvalue -= rhs[i][j].anvalue;
    results[i][j].asvalue -= rhs[i][j].asvalue;
    results[i][j].svalue -= rhs[i][j].svalue;
  }
  return results;
}
FiniteMatrix::finiteMat operator&&(const FiniteMatrix::finiteMat &lhs, const FiniteMatrix::finiteMat &rhs)
{
  FiniteMatrix::finiteMat results(lhs);
  forAllInternal(lhs)
  {
    results[i][j].value *= rhs[i][j].value;
    results[i][j].aevalue *= rhs[i][j].aevalue;
    results[i][j].awvalue *= rhs[i][j].awvalue;
    results[i][j].anvalue *= rhs[i][j].anvalue;
    results[i][j].asvalue *= rhs[i][j].asvalue;
    results[i][j].svalue *= rhs[i][j].svalue;
  }
  return results;
}

FiniteMatrix::finiteMat operator*(const double dbvalue, const FiniteMatrix::finiteMat &rhs)
{
  FiniteMatrix::finiteMat results(rhs);
  forAllInternal(results)
  {
    results[i][j].value *= dbvalue;
  }
  return results;
}

void FiniteMatrix::print2dMat(finiteMat &vec)
{
  for (unsigned int i = 0; i < vec.size(); i++)
  {
    for (unsigned int j = 0; j < vec[i].size(); j++)
    {
      std::cout << std::setprecision(3) << vec[i][j].value << ',';
    }
    cout << endl;
  }
}

void FiniteMatrix::print2dsource(finiteMat &vec)
{
  for (unsigned int i = 0; i < vec.size(); i++)
  {
    for (unsigned int j = 0; j < vec[i].size(); j++)
    {
      std::cout << std::setprecision(3) << vec[i][j].svalue << ',';
    }
    cout << endl;
  }
}

FiniteMatrix::finiteMat FiniteMatrix::interpolatedFieldEast(FiniteMatrix::finiteMat &vec, Grid &myGrid_)
{
  FiniteMatrix::finiteMat temp(vec.size(), vector<FiniteMatrix>(vec[0].size()));
  forAllInternalUCVs(vec)
  {
    double FXE = myGrid_.FX[i];
    double FXP = 1.0 - FXE;

    temp[i][j].value = (vec[i + 1][j].value * FXE) + (vec[i][j].value * FXP);
  }
  return temp;
}
FiniteMatrix::finiteMat FiniteMatrix::interpolatedFieldNorth(FiniteMatrix::finiteMat &vec, Grid &myGrid_)
{
  FiniteMatrix::finiteMat temp(vec.size(), vector<FiniteMatrix>(vec[0].size()));
  forAllInternalVCVs(vec)
  {
    double FYN = myGrid_.FY[j];
    double FYP = 1.0 - FYN;

    temp[i][j].value = (vec[i][j + 1].value * FYN) + (vec[i][j].value * FYP);
  }
  return temp;
}
Fields::vectorfields FiniteMatrix::correctFaceVelocityEast(Fields::vectorfields &interpolatedCellFaceVelocity,
                                                           Fields::vectorfields &cellFacePressureGradient,
                                                           Fields::vectorfields &DPXField,
                                                           FiniteMatrix::finiteMat &APinterpolated,
                                                           Grid &myGrid_)
{
  Fields::vectorfields temp(interpolatedCellFaceVelocity.size(), vector<Fields>(interpolatedCellFaceVelocity[0].size()));
  forAllInternalUCVs(temp)
  {
     double DXPE = (myGrid_.XC[i + 1] - myGrid_.XC[i]);
    double sArea = myGrid_.Y[j] - myGrid_.Y[j - 1];
    double volume = DXPE * sArea;

    temp[i][j].value= interpolatedCellFaceVelocity[i][j].value - (APinterpolated[i][j].value * volume * (cellFacePressureGradient[i][j].value - DPXField[i][j].value));
  }
  return temp;
}
Fields::vectorfields FiniteMatrix::correctFaceVelocityNorth(Fields::vectorfields &interpolatedCellFaceVelocity,
                                                            Fields::vectorfields &cellFacePressureGradient,
                                                            Fields::vectorfields &DPYField,
                                                            FiniteMatrix::finiteMat &APinterpolated,
                                                            Grid &myGrid_)
{
  Fields::vectorfields temp(interpolatedCellFaceVelocity.size(), vector<Fields>(interpolatedCellFaceVelocity[0].size()));
  forAllInternalVCVs(temp)
  {
    double DYPN = (myGrid_.YC[j + 1] - myGrid_.XC[j]);
    double sArea = myGrid_.X[i] - myGrid_.X[i - 1];
    double volume = DYPN * sArea;

    temp[i][j].value = interpolatedCellFaceVelocity[i][j].value - (APinterpolated[i][j].value * volume * (cellFacePressureGradient[i][j].value - DPYField[i][j].value));
  }
  return temp;
}
void FiniteMatrix::correctEastMassFluxes(Fields::vectorfields &massFE, Fields::vectorfields &PressureCorr, FiniteMatrix::finiteMat &AEmat)
{
  forAllInternalUCVs(massFE)
  {
    massFE[i][j].value = massFE[i][j].value + (AEmat[i][j].value * (PressureCorr[i + 1][j].value - PressureCorr[i][j].value));
  }
}
void FiniteMatrix::correctNorthMassFluxes(Fields::vectorfields &massFN, Fields::vectorfields &PressureCorr, FiniteMatrix::finiteMat &ANmat)
{
  forAllInternalVCVs(massFN)
  {
    massFN[i][j].value = massFN[i][j].value + (ANmat[i][j].value * (PressureCorr[i][j+1].value - PressureCorr[i][j].value));
  }
}