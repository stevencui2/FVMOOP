#ifndef FINITEVOLUMEOPERATIONS_H
#define FINITEVOLUMEOPERATIONS_H

#include "Grid.h"
#include "Fields.h"
#include "Solution.h"
#include "forAllOperations.h"
#include "FiniteMatrix.h"

namespace fvm
{
  // Momentum Eqns
  // Diffusion(surface integrate), conveciton(surface integrate), Presure source(volume integrate), Time derivative(volume integrate)
  FiniteMatrix::finiteMat diffusiveTerm(Fields::vectorfields &vec)
  {
    FiniteMatrix::finiteMat APtemp(vec.size(), vector<FiniteMatrix>(vec[0].size()));

    // towards east setVectorFieldGridFeatures
    forAllInternalUCVs(vec)
    {
      APtemp[i][j].aevalue = -(vec[i][j].visc * vec[i][j].Se) / (vec[i][j].DXPtoE);
      APtemp[i + 1][j].awvalue = -(vec[i][j].visc * vec[i][j].Se) / (vec[i][j].DXPtoE);
    }
    // towards north size setVectorFieldGridFeatures
    forAllInternalVCVs(vec)
    {
      APtemp[i][j].anvalue = -(vec[i][j].visc * vec[i][j].Sn) / (vec[i][j].DYPtoN);
      APtemp[i][j + 1].asvalue = -(vec[i][j].visc * vec[i][j].Sn) / (vec[i][j].DYPtoN);
    }
    return APtemp;
  }
  FiniteMatrix::finiteMat convectionTerm(Fields::vectorfields &vec, Fields::vectorfields &massFluxEast, Fields::vectorfields &massFluxNorth, double &blendFac)
  {

    FiniteMatrix::finiteMat APtemp(vec.size(), vector<FiniteMatrix>(vec[0].size()));
    // towards to East Side
    forAllInternalUCVs(vec)
    {
      double cnvFluxNodeE = min(massFluxEast[i][j].value, 0.0);
      double cnvFluxNodeP = max(massFluxEast[i][j].value, 0.0);

      APtemp[i][j].aevalue = cnvFluxNodeE;
      APtemp[i + 1][j].awvalue = -cnvFluxNodeP;

      double fluxUDSapprox = cnvFluxNodeP * vec[i][j].value + cnvFluxNodeE * vec[i + 1][j].value;
      double fluxCDSapprox = massFluxEast[i][j].value * (vec[i + 1][j].value * vec[i][j].FXE + vec[i][j].value * vec[i][j].FXP);

      double resultvalue = blendFac * (fluxUDSapprox - fluxCDSapprox);

      // source term
      APtemp[i][j].svalue = APtemp[i][j].svalue + resultvalue;
      APtemp[i + 1][j].svalue = APtemp[i][j].svalue - resultvalue;
    }
    // towards to North Side
    forAllInternalVCVs(vec)
    {
      double cnvFluxNodeN = min(massFluxNorth[i][j].value, 0.0);
      double cnvFluxNodeP = max(massFluxNorth[i][j].value, 0.0);

      APtemp[i][j].anvalue = cnvFluxNodeN;
      APtemp[i][j + 1].asvalue = -cnvFluxNodeP;

      double fluxUDSapprox = cnvFluxNodeP * vec[i][j].value + cnvFluxNodeN * vec[i][j + 1].value;
      double fluxCDSapprox = massFluxNorth[i][j].value * (vec[i][j + 1].value * vec[i][j].FYN + vec[i][j].value * vec[i][j].FYP);

      double resultvalue = blendFac * (fluxUDSapprox - fluxCDSapprox);

      APtemp[i][j].svalue = APtemp[i][j].svalue + resultvalue;
      APtemp[i][j + 1].svalue = APtemp[i][j].svalue - resultvalue;
    }
    return APtemp;
  } // end convection term

  FiniteMatrix::finiteMat pressureGrad(Fields::vectorfields &vec, int &direction)
  {
    FiniteMatrix::finiteMat APtemp(vec.size(), vector<FiniteMatrix>(vec[0].size()));
    forAllInternal(vec)
    {
      double DX = vec[i][j].X - vec[i - 1][j].X;
      double DY = vec[i][j].Y - vec[i][j - 1].Y;

      // calculate interpolated face pressure
      double pressureEastFace = (vec[i + 1][j].value * vec[i][j].FXE) + (vec[i][j].value * vec[i][j].FXP);
      double pressureWestFace = (vec[i][j].value * vec[i - 1][j].FXE) + (vec[i - 1][j].value * vec[i - 1][j].FXP);
      double pressureNorthFace = (vec[i][j + 1].value * vec[i][j].FYN) + (vec[i][j].value * vec[i][j].FYP);
      double pressureSouthFace = (vec[i][j].value * vec[i][j - 1].FYN) + (vec[i][j - 1].value * vec[i][j - 1].FYP);

      double pressureEastGrad = (pressureEastFace - pressureWestFace) / DX;
      double pressureNorthGrad = (pressureNorthFace - pressureSouthFace) / DY;
      if (direction == 1) // U momentum Eqn
      {
        APtemp[i][j].svalue = -pressureEastGrad * vec[i][j].volume;
      }
      else if (direction == 2) // V momentum Eqn
      {
        APtemp[i][j].svalue = -pressureNorthGrad * vec[i][j].volume;
      }
    }
    return APtemp;
  } // end pressureGrad

  // pressure poisson equation
  FiniteMatrix::finiteMat HTerm(const FiniteMatrix::finiteMat &RAPU, const FiniteMatrix::finiteMat &RAPV, Fields::vectorfields &vec, Fields::vectorfields vec2)
  {
    FiniteMatrix::finiteMat APtemp(vec.size(), vector<FiniteMatrix>(vec[0].size()));
    // towards to east side
    forAllInternalUCVs(APtemp)
    {
      APtemp[i][j].aevalue = -vec[i][j].density * vec[i][j].Se * RAPU[i][j].value * vec2[i][j].Se;
      APtemp[i + 1][j].awvalue = -vec[i][j].density * vec[i][j].Se * RAPU[i][j].value * vec2[i][j].Se;
    }
    // towards to north side
    forAllInternalVCVs(APtemp)
    {
      APtemp[i][j].anvalue = -vec[i][j].density * vec[i][j].Sn * RAPV[i][j].value * vec2[i][j].Sn;
      APtemp[i][j + 1].asvalue = -vec[i][j].density * vec[i][j].Sn * RAPV[i][j].value * vec2[i][j].Sn;
    }
    return APtemp;
  } // end HTerm
  FiniteMatrix::finiteMat divPhi(const Fields::vectorfields &Feast, const Fields::vectorfields &Fnorth)
  {
    FiniteMatrix::finiteMat APtemp(Feast.size(), vector<FiniteMatrix>(Feast[0].size()));
    forAllInternalUCVs(APtemp)
    {
      APtemp[i][j].value = Feast[i - 1][j].value - Feast[i][j].value + Fnorth[i][j - 1].value - Fnorth[i][j].value;
    }
    return APtemp;
  }
} // end namespace fvm

#endif
