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
		VEqn->EqnName = "V-momentum";
		Fields::vectorfields Vtemp = VEqn->solve(V, SV, sol, Viter);
		fieldsOper.copyInternalField(Vtemp, V);
		RAPV = VEqn->rAP;


		
		// 1 - DPX,DPY
		// 2 - field, U, V, or P, interpolate East North
		// 3 - Compute Mass FLuxes(MassFE, massFN);
		// 4 - interpolate RAPU and RAPV
		// 5 - correct mass fluxs
		// grad P cell center
		fieldsOper.computeCellCenterPressureGrad(P, DPX, DPY);
		//grad p at cell face contributed for sum of face  gradPf_avg
		Fields::vectorfields gradPeAvg(fieldsOper.interpolatedFieldEast(DPX, mygrid_));
		Fields::vectorfields gradPnAvg(fieldsOper.interpolatedFieldNorth(DPY, mygrid_));
		// velocity at cell face 
		Fields::vectorfields UVelPtoE(fieldsOper.interpolatedFieldEast(U, mygrid_));
		Fields::vectorfields VVelPtoN(fieldsOper.interpolatedFieldNorth(V, mygrid_));

		FiniteMatrix::finiteMat RAPUPtoE(finiteobj.interpolatedFieldEast(RAPU, mygrid_));
		FiniteMatrix::finiteMat RAPVPtoN(finiteobj.interpolatedFieldNorth(RAPV, mygrid_));

		// cell face pressure gradient contributed from PP and PF gradPf
		Fields::vectorfields gradPe(fieldsOper.cellFaceGradientEast(P, mygrid_));
		Fields::vectorfields gradPn(fieldsOper.cellFaceGradientNorth(P, mygrid_));
		// correct face velocities
		Fields::vectorfields UE(finiteobj.correctFaceVelocityEast(UVelPtoE, gradPe, gradPeAvg, RAPUPtoE, mygrid_));
		Fields::vectorfields VN(finiteobj.correctFaceVelocityNorth(VVelPtoN, gradPn, gradPnAvg, RAPVPtoN, mygrid_));
		fieldsOper.computeEastMassFluxes(massFluxE, UE);
		fieldsOper.computeNorthMassFluxes(massFluxN, VN);

		Equation *PEqn = new Equation(
			fvm::HTerm(RAPUPtoE, RAPVPtoN, U, V) + fvm::divPhi(massFluxE, massFluxN));
		PEqn->assembleEquation();
		double zero = 0.0;
		forAllInternal(PP)
		{
			PP[i][j].value = 0.0;
		}
		FiniteMatrix::finiteMat SP(PEqn->sourceFinial);
		int Piter = 6;
		PEqn->EqnName = "P-Eqn";
		PEqn->URF = 0.2;
		FiniteMatrix::finiteMat PAE(PEqn->AE);
		FiniteMatrix::finiteMat PAN(PEqn->AN);
		PP = PEqn->solve(PP, SP, sol, Piter);

		fieldsOper.linearextrapolateCondition(PP, mygrid_.FX, mygrid_.FY, north);
		fieldsOper.linearextrapolateCondition(PP, mygrid_.FX, mygrid_.FY, south);
		fieldsOper.linearextrapolateCondition(PP, mygrid_.FX, mygrid_.FY, east);
		fieldsOper.linearextrapolateCondition(PP, mygrid_.FX, mygrid_.FY, west);

		// reference location to obtain pressure
		int IPR = 4;
		int JPR = 4;
		double PPO = PP[IPR][JPR].value;
		cout << "pressure PP at ref location 4,4 = " << PPO << endl;
		finiteobj.correctEastMassFluxes(massFluxE, PP, PAE);
		finiteobj.correctNorthMassFluxes(massFluxN, PP, PAN);

		// finial correct U V and pressure
		forAllInternal(U)
		{
			double DX = mygrid_.X[i] - mygrid_.X[i - 1];
			double DY = mygrid_.Y[j] - mygrid_.Y[j - 1];

			double pressureEastFace = (PP[i + 1][j].value * PP[i][j].FXE) + (PP[i][j].value * PP[i][j].FXP);
			double pressureWestFace = (PP[i][j].value * PP[i - 1][j].FXE) + (PP[i - 1][j].value * PP[i - 1][j].FXP);
			double pressureNorthFace = (PP[i][j + 1].value * PP[i][j].FYN) + (PP[i][j].value * PP[i][j].FYP);
			double pressureSouthFace = (PP[i][j].value * PP[i][j - 1].FYN) + (PP[i][j - 1].value * PP[i][j - 1].FYP);

			double URF_P = 0.2;

			// U V P
			U[i][j].value = U[i][j].value - (pressureEastFace - pressureWestFace) * DY * RAPU[i][j].value;
			V[i][j].value = V[i][j].value - (pressureNorthFace - pressureSouthFace) * DX * RAPV[i][j].value;
			P[i][j].value = P[i][j].value + URF_P * (PP[i][j].value - PPO);
		}

		// pressure -poisson RAPU and RAPV
	}

	return 0;
}