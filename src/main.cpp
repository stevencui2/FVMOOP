#include <iostream>
#include "forAllOperations.h"
#include "Fields.h"
#include "Equation.h"
#include "finiteVolumeOperations.h"

int main(int argc, char **argv)
{
	std::cout << "Hello CFD" << std::endl;
#include "createGridFieldsSol.h"
#include "createFiniteMatrices.h"

	// U velocity inlet boundaryCondition
	fieldsOper.inletboundaryCondition(U, north, 1.0);
	// cout << " initial values before entering the iterations " << endl;
	// fieldsOper.print2dmat(U);
	int iter = 0;
	int maxit = 1;
	for (iter = 0; iter < maxit; iter++)
	{
		// pressure boundary conditions should be updated every iteration

		fieldsOper.linearextrapolateCondition(P, mygrid_.FX, mygrid_.FY, north);
		fieldsOper.linearextrapolateCondition(P, mygrid_.FX, mygrid_.FY, south);
		fieldsOper.linearextrapolateCondition(P, mygrid_.FX, mygrid_.FY, east);
		fieldsOper.linearextrapolateCondition(P, mygrid_.FX, mygrid_.FY, west);

		// convection -blend factor
		double blendFactor = 1.0;

		int eastdir = 1, northdir = 2;
		// equation takes 1 argu, cosnt finite matrix
		Equation *UEqn = new Equation(
			fvm::diffusiveTerm(U) + fvm::convectionTerm(U, massFluxE, massFluxN, blendFactor) + fvm::pressureGrad(P, eastdir));

		UEqn->noWallShearXBoundaryCondition(U);
		UEqn->assembleEquation();
		UEqn->relax(U);

		FiniteMatrix::finiteMat SU(UEqn->sourceRelaxed);
		int Uiter = 1;
		Fields::vectorfields Utemp = UEqn->solve(U, SU, sol, Uiter);
		fieldsOper.copyInternalField(Utemp, U);
		cout << " U velocity updated after solving  " << endl;
		fieldsOper.print2dmat(U);
		RAPU = UEqn->rAP;

		// V momentum Eqn
		Equation *VEqn = new Equation(
			fvm::diffusiveTerm(V) + fvm::convectionTerm(V, massFluxE, massFluxN, blendFactor) + fvm::pressureGrad(P, northdir) // DPY
		);

		VEqn->noWallShearYBoundaryCondition(V);
		VEqn->assembleEquation();
		VEqn->relax(V);

		FiniteMatrix::finiteMat SV(VEqn->sourceRelaxed);
		int Viter = 1;
		Fields::vectorfields Vtemp = VEqn->solve(V, SV, sol, Viter);
		fieldsOper.copyInternalField(Vtemp, V);
		cout << " V velocity updated after solving  " << endl;
		fieldsOper.print2dmat(V);
		RAPV = VEqn->rAP;
		// pressure -poisson RAPU and RAPV
	}

	return 0;
}