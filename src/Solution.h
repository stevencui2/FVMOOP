#ifndef SOLUTION_H
#define SOLUTION_H

class Solution
{
public:
	Solution();
	virtual ~Solution();

	double dt;
	double R;
	double visc;
	double density;
	int nsteps, maxit;
	int Imonitor, Jmonitor;
	double URFUVel, URFVVel, URFPressure;
	double URFU, URFV, URFP;
	double Alfa; // strongly Implicit Procedure
};

#endif
