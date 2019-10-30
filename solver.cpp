#include "solver.h"

//#define _USE_MATH_DEFINES
//#include <cmath>
#define M_PI 3.14159265358979323846

#include <chrono>

Matrix Solver::initialization(int type, double k1, double k2)
{
    int s = 0;

    int n = N * M;
    Matrix A_flag(M, N);
    Matrix A(n, n);


    double alpha = k2 / k1 * hx / 2.0 / hy / hy;
    double beta = k2 / k1 * hy / 2.0 / hx / hx;

    double a = 1.0 / hx / hx;
    double b = 1.0 / hy / hy;


    for (int i = 0; i < M; i++)
        for (int j = 0; j < N; j++)
        {
            s = i * M + j;
            if (A_flag.type(i, j) == ANGLE)
            {
                A.setValue(s, i * M + j, 1.0);
            }
            else if (A_flag.type(i, j) == LEFT)
            {
                switch (type) {
                case 1:
                    // second order of approximation
                    A.setValue(s, 0 * M + j, -1.0 / hx);
                    A.setValue(s, 1 * M + j - 1,  alpha);
                    A.setValue(s, 1 * M + j, -2.0 * alpha + 1.0 / hx);
                    A.setValue(s, 1 * M + j + 1,  alpha);
                    break;
                case 2:
                    // first order of approximation
                    A.setValue(s, 1 * M + j, 1.0 / hx);
                    A.setValue(s, 0 * M + j, -1.0 / hx);
                    break;
                default:
                    // non-monotonic approximation
                    A.setValue(s, 0 * M + j, -3.0 / 2.0 / hx);
                    A.setValue(s, 1 * M + j, 2.0 / hx);
                    A.setValue(s, 2 * M + j, -1.0 / 2.0 / hx);
                    break;
                }
            }
            else if (A_flag.type(i, j) == RIGHT)
            {
                switch (type) {
                case 1:
                    // second order of approximation
                    A.setValue(s, (N - 1) * M + j, 1.0 / hx);
                    A.setValue(s, (N - 2) * M + j - 1, -alpha);
                    A.setValue(s, (N - 2) * M + j, 2.0 *alpha - 1.0 / hx);
                    A.setValue(s, (N - 2) * M + j + 1, -alpha);
                    break;
                case 2:
                    // first order of approximation
                    A.setValue(s, (N - 1) * M + j, 1.0 / hx);
                    A.setValue(s, (N - 2) * M + j, -1.0 / hx);
                    break;
                default:
                    // non-monotonic approximation
                    A.setValue(s, (N - 1) * M + j, 3.0 / 2.0 / hx);
                    A.setValue(s, (N - 2) * M + j, -2.0 / hx);
                    A.setValue(s, (N - 3) * M + j, 1.0 / 2.0 / hx);
                    break;
                }
            }
            else if (A_flag.type(i, j) == BOTTOM)
            {
                switch (type) {
                case 1:
                    A.setValue(s, i * M + 0, -1.0 / hy);
                    A.setValue(s, (i - 1) * M + 1, beta);
                    A.setValue(s, i * M + 1, - 2.0 *beta + 1.0 / hy);
                    A.setValue(s, (i + 1) * M + 1, beta);
                    break;
                case 2:
                    A.setValue(s, i * M + 1, 1.0 / hy);
                    A.setValue(s, i * M + 0, -1.0 / hy);
                    break;
                default:
                    // non-monotonic approximation
                    A.setValue(s, i * M + 0, -3.0 / 2.0 / hy);
                    A.setValue(s, i * M + 1, 2.0 / hy);
                    A.setValue(s, i * M + 2, -1.0 / 2.0 / hy);
                    break;
                }
            }
            else if (A_flag.type(i, j) == TOP)
            {
                switch (type) {
                case 1:
                    A.setValue(s, i * M + N - 1, 1.0 / hy);
                    A.setValue(s, (i - 1) * M + N - 2, -beta);
                    A.setValue(s, i * M + N - 2,2.0 *beta - 1.0 / hy);
                    A.setValue(s, (i + 1) * M + N - 2, -beta);
                    break;
                case 2:
                    A.setValue(s, i * M + N - 1, 1.0 / hy);
                    A.setValue(s, i * M + N - 2, -1.0 / hy);
                    break;
                default:
                    // non-monotonic approximation
                    A.setValue(s, i * M + N - 1, 3.0 / 2.0 / hy);
                    A.setValue(s, i * M + N - 2, -2.0 / hy);
                    A.setValue(s, i * M + N - 3, 1.0 / 2.0 / hy);
                    break;
                }
            }
            else if (A_flag.type(i, j) == CENTER)
            {
                A.setValue(s, i * M + j - 1, k2*b);
                A.setValue(s, i * M + j + 1, k2*b);
                A.setValue(s, i * M + j, -2.0 * k1 * a - 2.0 * k2 * b);
                A.setValue(s, (i - 1) * M + j,  k1*a);
                A.setValue(s, (i + 1) * M + j, k1*a);
            }
        }
    return A;
}

Matrix Solver::convertColumnToField(Matrix &X)
{
    Matrix res(M, N);
    for(int i = 0; i < M; i++)
    {
        for(int j = 0; j < N; j++)
        {
            //std::cout << X.getValue(res.getPointIndex(i, j), 0);
            res.setValue(i, j, X.getValue(res.getPointIndex(i, j), 0));
        }
    }
    return res;
}

Matrix Solver::convertFieldToColumn(Matrix &X)
{
    Matrix res(M*N, 1);
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            res.setValue(X.getPointIndex(i, j), 0, X.getValue(i, j));
        }
    }
    return res;
}

Solver::Solver(int M, int N, double hx, double hy, bool logging = false)
{
    this->M = M;
    this->N = N;
    this->hx = hx;
    this->hy = hy;
    this->logging = logging;
}

Matrix Solver::Analytic()
{
    Matrix x_0(M, N);

    for (int i = 0; i < M; i++)
        for (int j = 0; j < N; j++)
            //x_0.setValue(i, j, cos(2 * M_PI*i*hx)*exp(2 * M_PI*j*hy));
            //x_0.setValue(i, j, 1.0);
            //x_0.setValue(i, j, (i*hx - 1.0 / 2.0)*(i*hx - 1.0 / 2.0) - (j*hy - 1.0 / 2.0)*(j*hy - 1.0 / 2.0));
            x_0.setValue(i, j, i*i*i*hx*hx*hx - 3.0*i*hx*j*j*hy*hy);

    return x_0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


Matrix Solver::SymmetrizedConjugateGradients(int initializationType, std::ofstream &fout01, int BoundaryConditionsType, double k1, double k2)
{
    Matrix A(M*N, M*N);
    Matrix b(M*N, 1);
    Matrix b2(M*N, 1);

    A = initialization(initializationType, k1, k2);

    Matrix x_0(M, N);

    Matrix Ax_k(M*N, 1);
    Matrix x_k(M*N, 1);

    Matrix R_km1(M*N, 1);
    Matrix R_k(M*N, 1);
    Matrix R_1(M*N, 1);

    Matrix q_k(M*N, 1);
    Matrix q_km1(M*N, 1);

    Matrix p_k(M*N, 1);

    Matrix A_T(M*N, M*N);

    Matrix x_temp(M*N, 1);
    Matrix zij(M, N);

    Matrix stop(M*N, 1);

    double z_l = 0.0, z_l2 = 0.0;

    int s = 0;




    if (BoundaryConditionsType == 1)
    {
        for (int i = -7; i < 7; i++)
        {
            b.setValue(x_0.getPointIndex(i + M / 2, 0), 0, k2*50.0);//низ
        }

        for (int i = -7; i < 7; i++)
        {
            b.setValue(x_0.getPointIndex(i + M / 2, N - 1), 0, k1*50.0);//верх
        }

        ///////////////////////////
        s = x_0.getPointIndex(0, N / 2);
        for (int i = 0; i < M * N; i++)
            A.setValue(s, i, 0.0);
        A.setValue(s, s, 1.0);
        b.setValue(s, 0, -10.0);
    }
    else if (BoundaryConditionsType == 5)
    {
        for (int i = -3; i < 4; i++)
        {
            b.setValue(x_0.getPointIndex(i + M / 2, 0), 0, k2*100.0);//низ
        }

        for (int i = -7; i < 7; i++)
        {
            b.setValue(x_0.getPointIndex(i + M / 2, N - 1), 0, k1*50.0);//верх
        }

        ///////////////////////////
        s = x_0.getPointIndex(0, N / 2);
        for (int i = 0; i < M * N; i++)
            A.setValue(s, i, 0.0);
        A.setValue(s, s, 1.0);
        b.setValue(s, 0, 1.0);
    }
    else if (BoundaryConditionsType == 6)
    {
        for (int i = -7; i < 7; i++)
        {
            b.setValue(x_0.getPointIndex(i + M / 2, 0), 0, k2*50.0);//низ
        }

        for (int i = -3; i < 4; i++)
        {
            b.setValue(x_0.getPointIndex(i + M / 2, N - 1), 0, k1*100.0);//верх
        }

        ///////////////////////////
        s = x_0.getPointIndex(0, N / 2);
        for (int i = 0; i < M * N; i++)
            A.setValue(s, i, 0.0);
        A.setValue(s, s, 1.0);
        b.setValue(s, 0, 1.0);
    }
    else if (BoundaryConditionsType == 2)
    {
        for (int j = 1; j < N - 1; j++)
        {

            //x^3 - 3xy^2
            b.setValue(x_0.getPointIndex(0, j), 0, -3.0*hy*hy*j*j);//лево
            b.setValue(x_0.getPointIndex(M - 1, j), 0, 3.0 - 3.0*hy*hy*j*j);//право
        }
        for (int i = 1; i < M - 1; i++)
        {
            //x^3 - 3xy^2
            b.setValue(x_0.getPointIndex(i, N - 1), 0, -6.0*hx*i);//верх
            b.setValue(x_0.getPointIndex(i, 0), 0, 0.0);//низ
        }

        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                x_0.setValue(i, j, i*i*i*hx*hx*hx - 3.0 * i*hx * j*j*hy*hy);
        //////////////////////////////////////////////////////////////////////////////////////////
        //x^3-3xy^2

        s = x_0.getPointIndex(0, N / 2);
        for (int i = 0; i < M * N; i++)
            A.setValue(s, i, 0.0);
        A.setValue(s, s, 1.0);
        b.setValue(s, 0, 0.0);

    }
    else if (BoundaryConditionsType == 3)
    {
        for (int j = 1; j < N - 1; j++)
        {
            //(x-1/2)^2-(y-1/2)^2
            b.setValue(x_0.getPointIndex(0, j), 0, -1.0);//лево
            b.setValue(x_0.getPointIndex(M - 1, j), 0, 1.0);//право
        }

        for (int i = 1; i < M - 1; i++)
        {
            //(x-1/2)^2-(y-1/2)^2
            b.setValue(x_0.getPointIndex(i, N - 1), 0, -1.0);//верх
            b.setValue(x_0.getPointIndex(i, 0), 0, 1.0);//низ
        }

        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                x_0.setValue(i, j, (i*hx - 1.0 / 2.0)*(i*hx - 1.0 / 2.0) - (j*hy - 1.0 / 2.0)*(j*hy - 1.0 / 2.0));

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //(x-1/2)^2-(y-1/2)^2
        s = x_0.getPointIndex(0, N / 2);
        for (int i = 0; i < M * N; i++)
            A.setValue(s, i, 0.0);
        A.setValue(s, s, 1.0);
        b.setValue(s, 0, 1.0 / 4.0 - (N / 2 * hy - 1.0 / 2.0)*(N / 2 * hy - 1.0 / 2.0));

    }
    else if (BoundaryConditionsType == 4)
    {

        for (int i = 1; i < M - 1; i++)
        {
            //cos*exp
            b.setValue(x_0.getPointIndex(i, N - 1), 0, 2 * M_PI*cos(2 * M_PI*i*hx)*exp(2 * M_PI));//верх
            b.setValue(x_0.getPointIndex(i, 0), 0, 2 * M_PI*cos(2 * M_PI*i*hx));//низ
        }

        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                x_0.setValue(i, j, cos(2 * M_PI*i*hx)*exp(2 * M_PI*j*hy));

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        //cos*exp
        s = x_0.getPointIndex(M / 2 - 2, 0);
        for (int i = 0; i < M * N; i++)
            A.setValue(s, i, 0.0);
        A.setValue(s, s, 1.0);
        b.setValue(s, 0, cos(2 * M_PI * (M / 2 - 2) * hx));
        s = x_0.getPointIndex(M / 2 - 1, 0);
        for (int i = 0; i < M * N; i++)
            A.setValue(s, i, 0.0);
        A.setValue(s, s, 1.0);
        b.setValue(s, 0, cos(2 * M_PI * (M / 2 - 1) * hx));
        s = x_0.getPointIndex(M / 2, 0);
        for (int i = 0; i < M * N; i++)
            A.setValue(s, i, 0.0);
        A.setValue(s, s, 1.0);
        b.setValue(s, 0, cos(2 * M_PI * (M / 2) * hx));
        s = x_0.getPointIndex(M / 2 + 1, 0);
        for (int i = 0; i < M * N; i++)
            A.setValue(s, i, 0.0);
        A.setValue(s, s, 1.0);
        b.setValue(s, 0, cos(2 * M_PI * (M / 2 + 1) * hx));
        s = x_0.getPointIndex(M / 2 + 2, 0);
        for (int i = 0; i < M * N; i++)
            A.setValue(s, i, 0.0);
        A.setValue(s, s, 1.0);
        b.setValue(s, 0, cos(2 * M_PI * (M / 2 + 2) * hx));

    }

//	x_0.printDataToFile("Mya0.txt");

    int flag = 0;

    double alpha, beta;

    //1. Prestep

    trans(A, A_T);

    A.fillNonZeroElementsPos(0.00000001);
    A_T.fillNonZeroElementsPos(0.00000001);



    multiplication(A, x_k, Ax_k);
    summAlpha(Ax_k, b, -1.0, R_km1);
    multiplication(A_T, R_km1, R_k);

    R_1.copyFromData(R_k);

    std::ofstream fout("11.txt");

    while ((R_k.normL2() / R_1.normL2()) > 0.00001)
    {
        if (flag % 100000 == 0)
        {
            std::cout << "FLAG = " << flag << " Error = " << R_k.normL2() / R_1.normL2() << std::endl;
        }

        multiplication(A, x_k, stop);
        summAlpha(stop, b, -1.0, stop);

        if (BoundaryConditionsType == 2 || BoundaryConditionsType == 3 || BoundaryConditionsType == 4)
        {

            if (dotProduct(stop, stop) < 0.000001)
                break;

            x_temp = convertColumnToField(x_k);
            summAlpha(x_temp, x_0, -1.0, zij);
            z_l = zij.max();
            z_l2 = zij.normL2();
            z_l2 *= hx;

            fout << flag << " " << z_l2 << std::endl;
        }
        flag++;

        if(logging)
            x_k.printData();

        alpha = 1.0 / dotProduct(R_k, R_k);
        summAlpha(p_k, R_k, alpha, p_k);

        multiplication(A, p_k, q_km1);
        multiplication(A_T, q_km1, q_k);

        beta = 1.0 / dotProduct(q_k, p_k);
        summAlpha(x_k, p_k, -beta, x_k);

        summAlpha(R_k, q_k, -beta, R_k);

        q_km1.copyFromData(q_k);
    }

    //fout.close();

    Matrix z_n(M*N, 1);
    Matrix A_n(M*N, M*N);
    Matrix x_n(M*N, 1);
    Matrix b_n(M*N, 1);

    A_n.copyFromData(A);
    sqr(A_n);
    x_n.copyFromData(x_k);
    sqr(x_n);
    b_n.copyFromData(b);
    sqr(b_n);

    A_n.fillNonZeroElementsPos(0.00000001);

    multiplication(A_n, x_n, z_n);

    double z = 0;

    for (int i = 0; i < M*N; i++)
        z += z_n.getValue(i, 1) + b_n.getValue(i, 1);

    double r_1 = 1.0 / dotProduct(R_1, R_1);

    long double eps = 1.0 / 10000000000000000;

    double err = eps * eps * z * r_1;

    //fout.open("21.txt");

    while (err <= 1)
    {
        //fout << flag << ';' << err << std::endl;

        if (flag % 100000 == 0)
        {
            std::cout << "FLAG = " << flag << " Error = " << err << std::endl;
        }

        multiplication(A, x_k, stop);
        summAlpha(stop, b, -1.0, stop);

        if (BoundaryConditionsType == 2 || BoundaryConditionsType == 3 || BoundaryConditionsType == 4)
        {

            if (dotProduct(stop, stop) < 0.000001)
                break;

            x_temp = convertColumnToField(x_k);
            summAlpha(x_temp, x_0, -1.0, zij);
            z_l = zij.max();
            z_l2 = zij.normL2();
            z_l2 *= hx;

            fout << flag << " " << z_l2 << std::endl;
        }

        flag++;


        if (logging)
            x_k.printData();

        alpha = 1.0 / dotProduct(R_k, R_k);
        summAlpha(p_k, R_k, alpha, p_k);
        multiplication(A, p_k, q_km1);
        multiplication(A_T, q_km1, q_k);
        beta = 1.0 / dotProduct(q_k, p_k);

        summAlpha(x_k, p_k, -beta, x_k);
        summAlpha(R_k, q_k, -beta, R_k);
        q_km1.copyFromData(q_k);

        r_1 += 1.0 / dotProduct(R_k, R_k);

        err = eps * eps * z * r_1;
    }

    fout.close();

    // ЗАБИВАЕМ УГЛЫ В x_k
    x_k.setValue(x_0.getPointIndex(0, 0), 0, (x_k.getValue(x_0.getPointIndex(0, 1), 0) + x_k.getValue(x_0.getPointIndex(1, 0), 0)) / 2.0 ); // верх лево
    x_k.setValue(x_0.getPointIndex(0, N-1), 0, (x_k.getValue(x_0.getPointIndex(0, N-2), 0) + x_k.getValue(x_0.getPointIndex(1, N-1), 0)) / 2.0); // верх право
    x_k.setValue(x_0.getPointIndex(M-1, 0), 0, (x_k.getValue(x_0.getPointIndex(M-2, 0), 0) + x_k.getValue(x_0.getPointIndex(M-1, 1), 0)) / 2.0); // низ лево
    x_k.setValue(x_0.getPointIndex(M-1, N-1), 0, (x_k.getValue(x_0.getPointIndex(M-2, N-1), 0) + x_k.getValue(x_0.getPointIndex(M-1, N-2), 0)) / 2.0); // низ право


    //convertColumnToField(x_k).printDataToFile("Mya1.txt");

    if (BoundaryConditionsType == 2 || BoundaryConditionsType == 3 || BoundaryConditionsType == 4)
    {
        x_temp = convertColumnToField(x_k);
        summAlpha(x_temp, x_0, -1.0, zij);
        z_l = zij.max();
        z_l2 = zij.normL2();
        z_l2 *= hx;
        std::cout << "Okonchatelino " << z_l << ' ' << z_l2 << std::endl;
        fout01 << N << " " << flag << " " << z_l2 << std::endl;
    }

    std::cout << "Flag = " << flag << std::endl;

    return x_k;
}

Matrix Solver::CraigMethod(int initializationType, std::ofstream &fout11, int BoundaryConditionsType, double k1, double k2)
{
    Matrix A(M*N, M*N);
    Matrix b(M*N, 1);
    Matrix b2(M*N, 1);

    A = initialization(initializationType, k1, k2);

    Matrix x_0(M, N);

    Matrix Ax_k(M*N, 1);
    Matrix x_k(M*N, 1);

    Matrix R_km1(M*N, 1);
    Matrix R_k(M*N, 1);
    Matrix R_1(M*N, 1);

    Matrix q_k(M*N, 1);
    Matrix q_km1(M*N, 1);

    Matrix p_k(M*N, 1);

    Matrix A_T(M*N, M*N);

    Matrix x_temp(M*N, 1);
    Matrix zij(M, N);

    Matrix stop(M*N, 1);

    double z_l = 0.0, z_l2 = 0.0;

    int flag = 0;

    double alpha, beta;
    int s ;



    if (BoundaryConditionsType == 5)
    {
        for (int i = -7; i < 7; i++)
        {
            b.setValue(x_0.getPointIndex(i + M / 2, 0), 0, k1*50.0);//низ
        }

        for (int i = -3; i < 4; i++)
        {
            b.setValue(x_0.getPointIndex(M - 1, i + N / 2), 0, k2*100.0);//верх
        }

        ///////////////////////////
        s = x_0.getPointIndex(0, N / 2);
        for (int i = 0; i < M * N; i++)
            A.setValue(s, i, 0.0);
        A.setValue(s, s, 1.0);
        b.setValue(s, 0, 1.0);
    }
    else if (BoundaryConditionsType == 1)
    {
        for (int i = -7; i < 7; i++)
        {
            b.setValue(x_0.getPointIndex(0, i + N / 2), 0, k2*100.0);//низ
        }

        for (int i = -14; i < 14; i++)
        {
            b.setValue(x_0.getPointIndex(i + N / 2, M - 1), 0, k1*50.0);//право
        }


        ///////////////////////////
        s = x_0.getPointIndex(M / 2, 0);
        for (int i = 0; i < M * N; i++)
            A.setValue(s, i, 0.0);
        A.setValue(s, s, 1.0);
        b.setValue(s, 0, 1.0);
    }
    else if (BoundaryConditionsType == 6)
    {
        for (int i = -7; i < 7; i++)
        {
            b.setValue(x_0.getPointIndex(i + M / 2, N-1), 0, k2*50.0);//низ
        }

        for (int i = -3; i < 4; i++)
        {
            b.setValue(x_0.getPointIndex(0, i + N / 2), 0, k1*50.0);//верх
        }

        ///////////////////////////
        s = x_0.getPointIndex(0, N / 2);
        for (int i = 0; i < M * N; i++)
            A.setValue(s, i, 0.0);
        A.setValue(s, s, 1.0);
        b.setValue(s, 0, 1.0);
    }
    else if (BoundaryConditionsType == 2)
    {
        for (int j = 1; j < N - 1; j++)
        {

            //x^3 - 3xy^2
            b.setValue(x_0.getPointIndex(0, j), 0, -3.0*hy*hy*j*j);//лево
            b.setValue(x_0.getPointIndex(M - 1, j), 0, 3.0 - 3.0*hy*hy*j*j);//право
        }
        for (int i = 1; i < M - 1; i++)
        {
            //x^3 - 3xy^2
            b.setValue(x_0.getPointIndex(i, N - 1), 0, -6.0*hx*i);//верх
            b.setValue(x_0.getPointIndex(i, 0), 0, 0.0);//низ
        }

        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                x_0.setValue(i, j, i*i*i*hx*hx*hx - 3.0 * i*hx * j*j*hy*hy);
        //////////////////////////////////////////////////////////////////////////////////////////
        //x^3-3xy^2

        s = x_0.getPointIndex(0, N / 2);
        for (int i = 0; i < M * N; i++)
            A.setValue(s, i, 0.0);
        A.setValue(s, s, 1.0);
        b.setValue(s, 0, 0.0);

    }
    else if (BoundaryConditionsType == 3)
    {
        for (int j = 1; j < N - 1; j++)
        {
            //(x-1/2)^2-(y-1/2)^2
            b.setValue(x_0.getPointIndex(0, j), 0, -1.0);//лево
            b.setValue(x_0.getPointIndex(M - 1, j), 0, 1.0);//право
        }

        for (int i = 1; i < M - 1; i++)
        {
            //(x-1/2)^2-(y-1/2)^2
            b.setValue(x_0.getPointIndex(i, N - 1), 0, -1.0);//верх
            b.setValue(x_0.getPointIndex(i, 0), 0, 1.0);//низ
        }

        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                x_0.setValue(i, j, (i*hx - 1.0 / 2.0)*(i*hx - 1.0 / 2.0) - (j*hy - 1.0 / 2.0)*(j*hy - 1.0 / 2.0));

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //(x-1/2)^2-(y-1/2)^2
        s = x_0.getPointIndex(0, N / 2);
        for (int i = 0; i < M * N; i++)
            A.setValue(s, i, 0.0);
        A.setValue(s, s, 1.0);
        b.setValue(s, 0, 1.0 / 4.0 - (N / 2 * hy - 1.0 / 2.0)*(N / 2 * hy - 1.0 / 2.0));

    }
    else if (BoundaryConditionsType == 4)
    {

        for (int i = 1; i < M - 1; i++)
        {
            //cos*exp
            b.setValue(x_0.getPointIndex(i, N - 1), 0, 2 * M_PI*cos(2 * M_PI*i*hx)*exp(2 * M_PI));//верх
            b.setValue(x_0.getPointIndex(i, 0), 0, 2 * M_PI*cos(2 * M_PI*i*hx));//низ
        }

        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                x_0.setValue(i, j, cos(2 * M_PI*i*hx)*exp(2 * M_PI*j*hy));

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        //cos*exp
        s = x_0.getPointIndex(M / 2 - 2, 0);
        for (int i = 0; i < M * N; i++)
            A.setValue(s, i, 0.0);
        A.setValue(s, s, 1.0);
        b.setValue(s, 0, cos(2 * M_PI * (M / 2 - 2) * hx));
        s = x_0.getPointIndex(M / 2 - 1, 0);
        for (int i = 0; i < M * N; i++)
            A.setValue(s, i, 0.0);
        A.setValue(s, s, 1.0);
        b.setValue(s, 0, cos(2 * M_PI * (M / 2 - 1) * hx));
        s = x_0.getPointIndex(M / 2, 0);
        for (int i = 0; i < M * N; i++)
            A.setValue(s, i, 0.0);
        A.setValue(s, s, 1.0);
        b.setValue(s, 0, cos(2 * M_PI * (M / 2) * hx));
        s = x_0.getPointIndex(M / 2 + 1, 0);
        for (int i = 0; i < M * N; i++)
            A.setValue(s, i, 0.0);
        A.setValue(s, s, 1.0);
        b.setValue(s, 0, cos(2 * M_PI * (M / 2 + 1) * hx));
        s = x_0.getPointIndex(M / 2 + 2, 0);
        for (int i = 0; i < M * N; i++)
            A.setValue(s, i, 0.0);
        A.setValue(s, s, 1.0);
        b.setValue(s, 0, cos(2 * M_PI * (M / 2 + 2) * hx));

    }

    //x_0.printDataToFile("Mya0.txt");

    //1. Prestep
    trans(A, A_T);

    A.fillNonZeroElementsPos(0.00000001);
    A_T.fillNonZeroElementsPos(0.00000001);

    // x_k при к = 0 любое, пусть x_1 = 0 || x_1 = b
    /*
    for (int i = 0; i < 9; i++)
        x_k.setValue(i, 0, b.getValue(i, 0));
    */

    multiplication(A, x_k, Ax_k);
    summAlpha(Ax_k, b, -1.0, R_k);
    R_1.copyFromData(R_k);

    std::ofstream fout2("12.txt");

    while ((R_k.normL2() / R_1.normL2()) > 0.00001)
    {

        if (flag % 100000 == 0)
        {
            std::cout << "FLAG = " << flag << " Error = " << R_k.normL2() / R_1.normL2() << std::endl;
        }

        multiplication(A, x_k, stop);
        summAlpha(stop, b, -1.0, stop);

        if (dotProduct(stop, stop) < 0.000001)
            break;

        if (BoundaryConditionsType == 2 || BoundaryConditionsType == 3 || BoundaryConditionsType == 4)
        {

            x_temp = convertColumnToField(x_k);
            summAlpha(x_temp, x_0, -1.0, zij);
            //z_l = zij.max();
            z_l2 = zij.normL2();
            z_l2 *= hx;

            fout2 << flag << " " << z_l2 << std::endl;
        }

        flag++;

        alpha = 1.0 / dotProduct(R_k, R_k);
        summAlpha(p_k, R_k, alpha, p_k);

        multiplication(A_T, p_k, q_k);

        beta = 1.0 / dotProduct(q_k, q_k);
        summAlpha(x_k, q_k, -beta, x_k);

        multiplication(A, x_k, Ax_k);
        summAlpha(Ax_k, b, -1.0, R_k);
    }

    //fout2.close();

    Matrix z_n(M*N, 1);
    Matrix A_n(M*N, M*N);
    Matrix x_n(M*N, 1);
    Matrix b_n(M*N, 1);

    A_n.copyFromData(A);
    sqr(A_n);
    x_n.copyFromData(x_k);
    sqr(x_n);
    b_n.copyFromData(b);
    sqr(b_n);

    A_n.fillNonZeroElementsPos(0.00000001);

    multiplication(A_n, x_n, z_n);

    double z;
    z = dotProduct(z_n, b_n);

    long double eps = 1.0 / 10000000000000000;

    double r_1 = 1.0 / dotProduct(R_1, R_1);

    double err = eps * eps * z * r_1;

    //fout2.open("22.txt");

    while (err <= 1)
    {
        if (flag % 100000 == 0)
        {
            std::cout << "FLAG = " << flag << " Error = " << err << std::endl;
        }

        multiplication(A, x_k, stop);
        summAlpha(stop, b, -1.0, stop);

        if (dotProduct(stop, stop) < 0.000001)
            break;

        if (BoundaryConditionsType == 2 || BoundaryConditionsType == 3 || BoundaryConditionsType == 4)
        {

            x_temp = convertColumnToField(x_k);
            summAlpha(x_temp, x_0, -1.0, zij);
            //z_l = zij.max();
            z_l2 = zij.normL2();
            z_l2 *= hx;

            fout2 << flag << " " << z_l2 << std::endl;
        }

        flag++;

        if (logging)
            x_k.printData();

        alpha = 1.0 / dotProduct(R_k, R_k);
        summAlpha(p_k, R_k, alpha, p_k);

        multiplication(A_T, p_k, q_k);

        beta = 1.0 / dotProduct(q_k, q_k);
        summAlpha(x_k, q_k, -beta, x_k);

        multiplication(A, x_k, Ax_k);
        summAlpha(Ax_k, b, -1.0, R_k);

        r_1 += 1.0 / dotProduct(R_k, R_k);

        err = eps * eps * z * r_1;
    }

    fout2.close();

    // ЗАБИВАЕМ УГЛЫ В x_k
    x_k.setValue(x_0.getPointIndex(0, 0), 0, (x_k.getValue(x_0.getPointIndex(0, 1), 0) + x_k.getValue(x_0.getPointIndex(1, 0), 0)) / 2.0); // верх лево
    x_k.setValue(x_0.getPointIndex(0, N - 1), 0, (x_k.getValue(x_0.getPointIndex(0, N - 2), 0) + x_k.getValue(x_0.getPointIndex(1, N - 1), 0)) / 2.0); // верх право
    x_k.setValue(x_0.getPointIndex(M - 1, 0), 0, (x_k.getValue(x_0.getPointIndex(M - 2, 0), 0) + x_k.getValue(x_0.getPointIndex(M - 1, 1), 0)) / 2.0); // низ лево
    x_k.setValue(x_0.getPointIndex(M - 1, N - 1), 0, (x_k.getValue(x_0.getPointIndex(M - 2, N - 1), 0) + x_k.getValue(x_0.getPointIndex(M - 1, N - 2), 0)) / 2.0); // низ право


    convertColumnToField(x_k).printDataToFile(fout11);
    //convertColumnToField(x_k).printDataToFile("Mya2.txt");

    if (BoundaryConditionsType == 2 || BoundaryConditionsType == 3 || BoundaryConditionsType == 4)
    {
        x_temp = convertColumnToField(x_k);
        summAlpha(x_temp, x_0, -1.0, zij);
        z_l = zij.max();
        z_l2 = zij.normL2();
        z_l2 *= hx;
        std::cout << "Okonchatelino " << z_l << ' ' << z_l2 << std::endl;
        fout11 << N << " " << flag << " " << z_l2 << std::endl;
    }

    std::cout << "Flag = " << flag << std::endl;

    return x_k;
}

