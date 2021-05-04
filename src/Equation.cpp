#include "Equation.h"

Equation::Equation(const FiniteMatrix::finiteMat& Fmatrix):
APinitial(Fmatrix),
URF(0.8),
NI(APinitial.size()),
NJ(APinitial[0].size()),
NIM(NI-1),
NJM(NJ-1),
AE(Fmatrix.size(),vector<FiniteMatrix>(Fmatrix[0].size())),
AW(Fmatrix.size(),vector<FiniteMatrix>(Fmatrix[0].size())),
AN(Fmatrix.size(),vector<FiniteMatrix>(Fmatrix[0].size())),
AS(Fmatrix.size(),vector<FiniteMatrix>(Fmatrix[0].size())),
AP(Fmatrix.size(),vector<FiniteMatrix>(Fmatrix[0].size())),
sourceInitial(Fmatrix.size(),vector<FiniteMatrix>(Fmatrix[0].size())),
sourceFinial(Fmatrix.size(),vector<FiniteMatrix>(Fmatrix[0].size())),
sourceRelaxed(Fmatrix.size(),vector<FiniteMatrix>(Fmatrix[0].size())),
UE(Fmatrix.size(),vector<FiniteMatrix>(Fmatrix[0].size())),
UN(Fmatrix.size(),vector<FiniteMatrix>(Fmatrix[0].size())),
LW(Fmatrix.size(),vector<FiniteMatrix>(Fmatrix[0].size())),
LS(Fmatrix.size(),vector<FiniteMatrix>(Fmatrix[0].size())),
LPR(Fmatrix.size(),vector<FiniteMatrix>(Fmatrix[0].size())),
RES(Fmatrix.size(),vector<FiniteMatrix>(Fmatrix[0].size())),
value(0.0), Residual(0.0),RSM(0.0),Literation(5),RESOR(4),
EqnName("U-momentum"), SOR(0.2)
{
  forAllInternal(AP){
    AP[i][j].value=0.0;
    AE[i][j].value=0.0;
    AW[i][j].value=0.0;
    AN[i][j].value=0.0;
    AS[i][j].value=0.0;
    sourceInitial[i][j].value=0.0;
  }

  forAllInternal(AP){
    AE[i][j].value=APinitial[i][j].aevalue;
    AW[i][j].value=APinitial[i][j].awvalue;
    AS[i][j].value=APinitial[i][j].asvalue;
    AN[i][j].value=APinitial[i][j].anvalue;
    sourceInitial[i][j].value=APinitial[i][j].svalue;
  }
}

Equation::~Equation(){

}


//strongly implicit solver (SIP)
