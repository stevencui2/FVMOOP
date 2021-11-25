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
