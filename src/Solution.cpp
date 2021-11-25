#include "Solution.h"

Solution::Solution() : dt(0.001), R(1.0), visc(1e-03), density(1.0),
                       nsteps(1), maxit(2), Imonitor(1), Jmonitor(1), URFUVel(0.8), URFVVel(0.8), URFPressure(0.2),
                       Alfa(0.92)
{
}

Solution::~Solution()
{
}
