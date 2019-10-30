#ifndef SOLVER_H
#define SOLVER_H

#include "data2d.h"
#include <iostream>
#include <fstream>
class Solver
{
    Matrix initialization(int type, double k1, double k2);
    Matrix convertColumnToField(Matrix& X);
    Matrix convertFieldToColumn(Matrix& X);

public:
    int M; //x
    int N; //y

    double hx, hy;

    bool logging;

    Solver(int M, int N, double hx, double hy, bool flag);

    Matrix Analytic();
    Matrix CraigMethod(int initializationType, std::ofstream &s, int BoundaryConditionsType, double k1, double k2);
    Matrix SymmetrizedConjugateGradients(int initializationType, std::ofstream &fout11, int BoundaryConditionsType, double k1, double k2);
};

#endif // SOLVER_H
