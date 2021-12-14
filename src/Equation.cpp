#include "Equation.h"

Equation::Equation(const FiniteMatrix::finiteMat &Fmatrix) : APinitial(Fmatrix),
                                                             URF(0.8),
                                                             NI(APinitial.size()),
                                                             NJ(APinitial[0].size()),
                                                             NIM(NI - 1),
                                                             NJM(NJ - 1),
                                                             AE(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             AW(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             AN(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             AS(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             AP(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             APU(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             rAP(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             sourceInitial(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             sourceB(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             sourceFinial(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             sourceRelaxed(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             UE(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             UN(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             LW(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             LS(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             LPR(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             RES(Fmatrix.size(), vector<FiniteMatrix>(Fmatrix[0].size())),
                                                             value(0.0), Residual(0.0), RSM(0.0), Literation(5), RESOR(4),
                                                             EqnName("U-momentum"), SOR(0.2)
{
  forAllInternal(AP)
  {
    AP[i][j].value = 0.0;
    AE[i][j].value = 0.0;
    AW[i][j].value = 0.0;
    AN[i][j].value = 0.0;
    AS[i][j].value = 0.0;
    sourceInitial[i][j].value = 0.0;
  }

  forAllInternal(AP)
  {
    AE[i][j].value = APinitial[i][j].aevalue;
    AW[i][j].value = APinitial[i][j].awvalue;
    AS[i][j].value = APinitial[i][j].asvalue;
    AN[i][j].value = APinitial[i][j].anvalue;
    sourceInitial[i][j].value = APinitial[i][j].svalue;
  }
}

Equation::~Equation()
{
}

void Equation::noWallShearXBoundaryCondition(Fields::vectorfields &vec)
{
  // south boundary
  for (int i = 1; i < vec.size() - 1; i++)
  {
    int j = 1;
    double diffcoeff = vec[i][j].visc * (vec[i][j].X - vec[i - 1][j].X) / (vec[i][j].YC - vec[i][j - 1].YC);
    APU[i][j].value += diffcoeff;
    sourceB[i][j].value = diffcoeff * vec[i][j - 1].value;
  }
  // north boundary
  for (int i = 1; i < vec.size() - 1; i++)
  {
    int j = vec[i].size() - 2;
    double diffcoeff = vec[i][j].visc * (vec[i][j].X - vec[i - 1][j].X) / (vec[i][j + 1].YC - vec[i][j].YC);
    APU[i][j].value += diffcoeff;
    sourceB[i][j].value = diffcoeff * vec[i][j + 1].value;
  }
}
void Equation::noWallShearYBoundaryCondition(Fields::vectorfields &vec)
{
  // west boundary
  for (int j = 1; j < vec[0].size() - 1; j++)
  {
    int i = 1;
    double diffcoeff = vec[i][j].visc * (vec[i][j].Y - vec[i][j - 1].Y) / (vec[i][j].XC - vec[i - 1][j].XC);
    APU[i][j].value += diffcoeff;
    sourceB[i][j].value = diffcoeff * vec[i - 1][j].value;
  }
  // east boundary
  for (int j = 1; j < vec[0].size() - 1; j++)
  {
    int i = vec.size() - 2;
    double diffcoeff = vec[i][j].visc * (vec[i][j].Y - vec[i][j - 1].Y) / (vec[i + 1][j].XC - vec[i][j].XC);
    APU[i][j].value += diffcoeff;
    sourceB[i][j].value = diffcoeff * vec[i + 1][j].value;
  }
}

void Equation::assembleEquation()
{
  forAllInternal(AP)
  {
    AP[i][j].value = -AE[i][j].value - AW[i][j].value - AS[i][j].value - AN[i][j].value + APU[i][j].value;
    sourceFinial[i][j].value = sourceInitial[i][j].value + sourceB[i][j].value;
    // source inital  is the source matrix without boundary source and without under relaxation
  }
}
void Equation::relax(Fields::vectorfields &vec)
{
  forAllInternal(AP)
  {
    AP[i][j].value = (1.0 / URF) * AP[i][j].value;
    sourceRelaxed[i][j].value = sourceFinial[i][j].value + (1.0 - URF) * (AP[i][j].value * vec[i][j].value);
    rAP[i][j].value = 1.0 / AP[i][j].value;
  }
}
// strongly implicit solver (SIP) will be replace by author's parallel AMG solver soon
Fields::vectorfields Equation::solve(Fields::vectorfields &phi, FiniteMatrix::finiteMat &source1, Solution &sol_, int &iterations){
  Fields::vectorfields phitemp(phi.size(), vector<Fields>(phi[0].size()));

  //..COEFFICIENTS OF UPPER AND LOWER TRIANGULAR MATRICES
  for (unsigned int i = 1; i < phi.size() - 1; i++)
  {
    for (unsigned int j = 1; j < phi[i].size() - 1; j++)
    {
      LW[i][j].value = AW[i][j].value / (1.0 + (sol_.Alfa * UN[i - 1][j].value));
      LS[i][j].value = AS[i][j].value / (1.0 + (sol_.Alfa * UE[i][j - 1].value));
      double P1 = sol_.Alfa * LW[i][j].value * UN[i - 1][j].value;
      double P2 = sol_.Alfa * LS[i][j].value * UE[i][j - 1].value;
      LPR[i][j].value = 1.0 / (AP[i][j].value + P1 + P2 - LW[i][j].value * UE[i - 1][j].value - LS[i][j].value * UN[i][j - 1].value);
      UN[i][j].value = (AN[i][j].value - P1) * LPR[i][j].value;
      UE[i][j].value = (AE[i][j].value - P2) * LPR[i][j].value;
    }
  }

  //.....CALCULATE RESIDUAL AND OVERWRITE IT BY INTERMEDIATE VECTOR
  for (int L = 0; L < iterations; L++)
  {
    Residual = 0.0;

    for (unsigned int i = 1; i < phi.size() - 1; i++)
    {
      for (unsigned int j = 1; j < phi[i].size() - 1; j++)
      {
        RES[i][j].value = source1[i][j].value - (AN[i][j].value * phi[i][j + 1].value) - (AS[i][j].value * phi[i][j - 1].value) -
                          (AE[i][j].value * phi[i + 1][j].value) - (AW[i][j].value * phi[i - 1][j].value) - (AP[i][j].value * phi[i][j].value);

        Residual += abs(RES[i][j].value);
        RES[i][j].value = (RES[i][j].value - (LS[i][j].value * RES[i][j - 1].value) - (LW[i][j].value * RES[i - 1][j].value)) * LPR[i][j].value;
        //    cout << std::setprecision(3)  << RES[i][j].value << ' ';
      }
    }

    double small = 1e-20;
    if (L == 0)
    {
      RESOR = Residual;
    }
    RSM = Residual / (RESOR + small);

    cout << EqnName << " Inner It: " << L << " and Residual --> " << Residual << " RSM " << RSM << endl;

    // Back Subsitution and Correction
    for (unsigned int i = phi.size() - 2; i >= 1; --i)
    {
      for (unsigned int j = phi[i].size() - 2; j >= 1; --j)
      {

        RES[i][j].value = RES[i][j].value - (UN[i][j].value * RES[i][j + 1].value) - (UE[i][j].value * RES[i + 1][j].value);
        phi[i][j].value = phi[i][j].value + RES[i][j].value;
      }
    }
  }

  forAll(phitemp)
  {
    phitemp[i][j].value = phi[i][j].value;
  }
  return phitemp;
}