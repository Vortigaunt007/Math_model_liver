#include "data2d.h"
#include <iostream>
#include <fstream>
#include <iomanip>

Matrix::Matrix()
{
    setSize(10, 10);
}

Matrix::Matrix(int m, int n)
{
    setSize(m, n);
}

void Matrix::setSize(int m, int n)
{
    field.resize(m);
    for(int i = 0; i < m; i++)
        field[i].resize(n);
}

void Matrix::setValue(int x, int y, double value)
{
    field[x][y] = value;
}

int Matrix::getXSize()
{
    return field.size();
}

int Matrix::getYSize()
{
    return field[0].size();
}

double Matrix::getValue(int x, int y)
{
    return field[x][y];
}

void Matrix::printData()
{
    for(int i = 0; i < field.size(); i++)
    {
        for(int j = 0; j < field[i].size(); j++)
            std::cout << std::setw(10) << std::setprecision(5) << field[i][j] << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

//void Data2D::printDataToFile(std::string s)
void Matrix::printDataToFile(std::ofstream &fout)
{
    //std::ofstream fout(s);
    for(int i = 0; i < field.size(); i++)
    {
        for(int j = 0; j < field[i].size(); j++)
            fout << field[i][j] << "; ";
        fout << std::endl;
    }
    fout << std::endl;
    fout.close();
}

void Matrix::copyFromData(Matrix &D)
{
    if(!checkSizes(*this, D)) return;

    //#pragma omp parallel for schedule(dynamic) num_threads(4)
    for(int i = 0; i < getXSize(); i++)
        for(int j = 0; j < getYSize(); j++)
            setValue(i, j, D.getValue(i, j));
}

double Matrix::normL2()
{
    return sqrt(dotProduct(*this, *this));
}

double Matrix::max()
{
    double c = 0.0;
    for (int i = 0; i < this->getXSize(); i++)
        for (int j = 0; j < this->getYSize(); j++)
            if (abs(this->getValue(i, j)) > c)
                c = abs(this->getValue(i, j));
    return c;

}

void Matrix::reverse_str(int k, int l)
{
    for (int j = 0; j < field[0].size(); j++)
    {
        double temp = field[k][j];
        field[k][j] = field[l][j];
        field[l][j] = temp;
    }
}

void Matrix::multCoef_str(int i, double coef)
{
    for (int j = 0; j < field[0].size(); j++)
        field[i][j] = field[i][j] * coef;
}

typePoint Matrix::type(int i, int j)
{
    // ПРОВЕРИТЬ N и M + ПРОВЕРИТЬ i и j
    if ((i == 0 && j == 0) ||
        (i == 0 && j == getXSize() - 1) ||
        (i == getYSize() - 1 && j == 0) ||
        (i == getYSize() - 1 && j == getXSize() - 1))
        return ANGLE;
    if (i == 0) return LEFT;
    if (i == getYSize() - 1) return RIGHT;
    if (j == 0) return BOTTOM;
    if (j == getXSize() - 1) return TOP;
    return CENTER;
}


int Matrix::getPointIndex(int i, int j)
{
    return i * getYSize() + j;
}

void Matrix::fillNonZeroElementsPos(double eps)
{
    nonZeroElementsPos.resize(getXSize());
    for(int i = 0; i < getXSize(); i++) {
        for(int j = 0; j < getYSize(); j++)
            if(fabs(getValue(i, j)) > eps)
                nonZeroElementsPos[i].push_back(j);
        //std::cout << "Found " << nonZeroElementsPos[i].size() << "elements in " << i << "'s row" << std::endl;
    }

}

bool checkSizes(Matrix &D1, Matrix &D2)
{
    if(D1.getXSize() != D2.getXSize() || D1.getYSize() != D2.getYSize())
        return false;
    return true;
}

double dotProduct(Matrix &D1, Matrix &D2)
{
    if(!checkSizes(D1, D2)) {
        std::cout << "Cannot calculate dot product, sizes aren't equale!" << std::endl;
        return 0.0;
    }
    double result = 0.0;
    for(int i = 0; i < D1.getXSize(); i++)
        for(int j = 0; j < D1.getYSize(); j++)
            result += D1.getValue(i, j)*D2.getValue(i, j);
    return result;
}

void summAlpha(Matrix &D1, Matrix &D2, double alpha, Matrix &out)
{
    if(!checkSizes(D1, D2)) {
        std::cout << "Cannot calculate summ, sizes aren't equale!" << std::endl;
        return;
    }
    //#pragma omp parallel for schedule(dynamic) num_threads(4)
    for(int i = 0; i < D1.getXSize(); i++)
        for(int j = 0; j < D1.getYSize(); j++)
            out.setValue(i, j, D1.getValue(i, j) + alpha * D2.getValue(i, j));
}

void summAlpha_str(Matrix &D1, int i, int k, double alpha)
{
    for (int j = 0; j < D1.getYSize(); j++)
        D1.setValue(i, j, D1.getValue(i, j) + alpha * D1.getValue(k, j));
}

//matrix transposition
void trans(Matrix &D1, Matrix &D2)
{
    for (int i = 0; i < D1.getXSize(); i++) {
        for (int j = 0; j < D1.getYSize(); j++) {
            D2.setValue(i, j, D1.getValue(j, i));
        }
    }
}


void sqr(Matrix &D1)
{
    double c = 0.0;
    for(int i = 0; i < D1.getXSize(); i++)
        for (int j = 0; j < D1.getYSize(); j++)
        {
            c = D1.getValue(i, j) * D1.getValue(i, j);
            D1.setValue(i, j, c);

        }
}

//matrix multiplication i,j,k -> c[i,j] += a[i,k] * b[k,j]
void	multiplication(Matrix &D1, Matrix &D2, Matrix &D3)
{
    for (int i = 0; i < D1.getXSize(); i++)
        //#pragma omp parallel for schedule(dynamic) num_threads(4)
        for (int j = 0; j < D2.getYSize(); j++)
        {
            double c = 0.0;
            for (int k : D1.nonZeroElementsPos[i]) {
            //for (int k = 0; k < D1.getYSize(); k++) {
                c += D1.getValue(i, k) * D2.getValue(k, j);
            }

            D3.setValue(i, j, c);
        }
}

Matrix Gauss(Matrix &D1)
{
    int n = D1.getXSize();

    Matrix I(n, n);

    double coef = 0.0;

    for (int i = 0; i < n; i++)
        I.setValue(i, i, 1.0);

    for (int i = 0; i < n; i++)
    {
        //std::cout << "Gauss " << i << std::endl;
        if ((D1.getValue(i, i) == 0.0) && (i == n))
        {
            std::cout << "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFU" << std::endl;
            return I;
        }
        else if ((D1.getValue(i, i) == 0.0) && (i < n))
        {
            for (int j = i + 1; j < n; j++)
            {
                if (abs(D1.getValue(j, i)) > 0)
                {
                    D1.reverse_str(i, j);
                    I.reverse_str(i, j);
                    break;
                }
            }
        }

        for (int j = i + 1; j < n; j++)
        {

            coef = D1.getValue(j, i) / D1.getValue(i, i);
            if (D1.getValue(i, i) == 0.0) continue;
            summAlpha_str(D1, j, i, -coef);
            //D1.printData();
            summAlpha_str(I, j, i, -coef);
        }
    }

    for (int i = 0; i < n; i++)
    {
        coef = D1.getValue(i, i);
        //std::cout << coef << std::endl; //с.зн.
        D1.multCoef_str(i, 1.0 / coef);
        //D1.printData();
        I.multCoef_str(i, 1.0 / coef);
    }

    for (int i = n - 2; i >= 0; i--)
        for (int j = i; j >= 0; j--)
        {
            coef = D1.getValue(j, i + 1);
            summAlpha_str(D1, j, i + 1, -coef);
            //D1.printData();
            summAlpha_str(I, j, i + 1, -coef);
        }

    return I;
}
