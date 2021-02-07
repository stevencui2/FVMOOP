// FVMCPP.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "Grid.h"
int main()
{
    std::cout << "Hello CFD World!\n";

    int NX = 7, NY = 7;
    double length = 1.0;

    Grid mygrid(NX, NY, length);

    std::cout << "Grid spaceing dx " << mygrid.retdx() << std::endl;
    
    for (int i = 0; i < mygrid.X.size(); i++)
        std::cout << mygrid.X[i] <<' ';
    std::cout << std::endl;
    for (int i = 0; i < mygrid.XC.size(); i++)
        std::cout << mygrid.XC[i] << ' ' ;
    std::cout << std::endl;
    for (int i = 0; i < mygrid.XC.size(); i++)
        std::cout << mygrid.FX[i] << ' ';
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
