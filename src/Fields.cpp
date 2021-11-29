#include "Fields.h"

Fields::Fields() : value(0.0)
{
}

Fields::Fields(int &NI_, int &NJ_) : value(0.0), NI(NI_), NJ(NJ_), NIM(NI - 1), NJM(NJ - 1)
{
}

Fields::~Fields()
{
}
void Fields::getGridInfoPassed(Fields::vectorfields &f, Grid &grid_, Solution &sol_)
{

	forAll(f)
	{
		f[i][j].X = grid_.X[i];
		f[i][j].Y = grid_.Y[j];
		f[i][j].XC = grid_.XC[i];
		f[i][j].YC = grid_.YC[j];
		f[i][j].FXE = grid_.FX[i];
		f[i][j].FYN = grid_.FY[j];
		f[i][j].FXP = 1.0 - grid_.FX[i];
		f[i][j].FYP = 1.0 - grid_.FY[j];
		f[i][j].visc = sol_.visc;
		f[i][j].density = sol_.density;
	}

	forAllInternal(f)
	{

		f[i][j].DXPtoE = grid_.XC[i + 1] - grid_.XC[i];
		f[i][j].DYPtoN = grid_.YC[j + 1] - grid_.YC[j];
		f[i][j].Se = grid_.Y[j] - grid_.Y[j - 1];
		f[i][j].Sn = grid_.X[i] - grid_.X[i - 1];
		f[i][j].volume = f[i][j].Se * f[i][j].Sn;
	}
}
void Fields::setVectorFieldGridFeatures()
{
}
void Fields::initializeFields(Fields::vectorfields &vec, double val)
{
	forAll(vec)
	{
		vec[i][j].value = val;
	}
}

void Fields::initializeInternalFields(Fields::vectorfields &vec, double val)
{
	forAllInternal(vec)
	{
		vec[i][j].value = val;
	}
}
void Fields::print2dmat(Fields::vectorfields &vec)
{
	for (unsigned int i = 0; i < vec.size(); i++)
	{
		for (unsigned int j = 0; j < vec[i].size(); j++)
		{
			cout << vec[i][j].value << ',';
		}
		cout << endl;
	}
}
void Fields::copyInternalField(Fields::vectorfields &from, Fields::vectorfields &to)
{
	forAllInternal(from)
	{
		to[i][j].value = from[i][j].value;
	}
}
// boundary conditions (1) -which modifieds the value of Fields (not matrices)
void Fields::inletboundaryCondition(Fields::vectorfields &vec, string &wallname, double bvalue)
{
	if (wallname == "East")
	{
		forEastBoundary(vec)
		{
			vec[i][j].value = bvalue;
		}
	}
	if (wallname == "West")
	{
		forWestBoundary(vec)
		{
			vec[i][j].value = bvalue;
		}
	}
	if (wallname == "North")
	{
		forNorthBoundary(vec)
		{
			vec[i][j].value = bvalue;
		}
	}
	if (wallname == "South")
	{
		forSouthBoundary(vec)
		{
			vec[i][j].value = bvalue;
		}
	}
}
void Fields::linearextrapolateCondition(Fields::vectorfields &vec, vector<double> &FXvec, vector<double> &FYvec, string &wallname)
{
	if (wallname == "East")
	{
		for (unsigned int j = 1; j < vec[0].size() - 1; j++)
		{
			unsigned int i = vec[0].size() - 1;
			vec[i][j].value = vec[i - 1][j].value + (vec[i - 1][j].value - vec[i - 2][j].value) * FXvec[i - 1];
		}
	}
	if (wallname == "West")
	{
		for (unsigned int j = 1; j < vec[0].size() - 1; j++)
		{
			unsigned int i = 0;
			vec[i][j].value = vec[i + 1][j].value + (vec[i + 1][j].value - vec[i + 2][j].value) * FXvec[i + 1];
		}
	}
	if (wallname == "North")
	{
		for (unsigned int i = 1; i < vec[0].size() - 1; i++)
		{
			unsigned int j = vec[0].size() - 1;
			vec[i][j].value = vec[i][j - 1].value + (vec[i][j - 1].value - vec[i][j - 2].value) * FYvec[j - 1];
		}
	}
	if (wallname == "South")
	{
		for (unsigned int i = 1; i < vec[0].size() - 1; i++)
		{
			unsigned int j = 0;
			vec[i][j].value = vec[i][j + 1].value + (vec[i][j + 1].value - vec[i][j + 2].value) * FYvec[j + 1];
		}
	}
}
void Fields::computeCellCenterPressureGrad(Fields::vectorfields &vec, Fields::vectorfields &dvecx, Fields::vectorfields &dvecy)
{
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
		dvecx[i][j].value = pressureEastGrad;
		dvecy[i][j].value = pressureNorthGrad;
	}
}
Fields::vectorfields Fields::interpolatedFieldEast(Fields::vectorfields &vec, Grid &myGrid_)
{
	Fields::vectorfields temp(vec.size(), vector<Fields>(vec[0].size()));
	forAllInternalUCVs(vec)
	{
		double FXE = myGrid_.FX[i];
		double FXP = 1.0 - FXE;

		temp[i][j].value = (vec[i + 1][j].value * FXE) + (vec[i][j].value * FXP);
	}
	return temp;
}
Fields::vectorfields Fields::interpolatedFieldNorth(Fields::vectorfields &vec, Grid &myGrid_)
{
	Fields::vectorfields temp(vec.size(), vector<Fields>(vec[0].size()));
	forAllInternalVCVs(vec)
	{
		double FYN = myGrid_.FY[j];
		double FYP = 1.0 - FYN;

		temp[i][j].value = (vec[i][j + 1].value * FYN) + (vec[i][j].value * FYP);
	}
	return temp;
}
Fields::vectorfields Fields::cellFaceGradientEast(Fields::vectorfields &vec, Grid &myGrid_)
{
	Fields::vectorfields temp(vec.size(), vector<Fields>(vec[0].size()));
	forAllInternalUCVs(vec)
	{
		double DXPE = (myGrid_.XC[i + 1] - myGrid_.XC[i]);
		temp[i][j].value = (vec[i + 1][j].value - vec[i][j].value) / DXPE;
	}
	return temp;
}

Fields::vectorfields Fields::cellFaceGradientNorth(Fields::vectorfields &vec, Grid &myGrid_)
{
	Fields::vectorfields temp(vec.size(), vector<Fields>(vec[0].size()));
	forAllInternalVCVs(vec)
	{
		double DYPN = (myGrid_.YC[j + 1] - myGrid_.YC[j]);
		temp[i][j].value = (vec[i][j + 1].value - vec[i][j].value) / DYPN;
	}
	return temp;
}

void Fields::computeEastMassFluxes(Fields::vectorfields &vec, Fields::vectorfields &corrU)
{
	forAllInternalUCVs(vec)
	{
		double sArea = vec[i][j].Se;
		double density = vec[i][j].density;
		vec[i][j].value = sArea * density * corrU[i][j].value;
	}
}

void Fields::computeNorthMassFluxes(Fields::vectorfields &vec, Fields::vectorfields &corrV)
{
	forAllInternalVCVs(vec)
	{
		double sArea = vec[i][j].Sn;
		double density = vec[i][j].density;
		vec[i][j].value = sArea * density * corrV[i][j].value;
	}
}