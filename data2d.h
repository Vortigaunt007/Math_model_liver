#ifndef VECTOR2D_H
#define VECTOR2D_H
#include <omp.h>
#include <vector>
#include <math.h>
#include <string>

enum typePoint { ANGLE, LEFT, RIGHT, BOTTOM, TOP, CENTER };

class Matrix
{
    std::vector< std::vector <double> > field;
public:
    Matrix();
    Matrix(int m, int n);


    //Setters and getters
    void setSize(int m, int n);
    int getXSize();
    int getYSize();

    void setValue(int x, int y, double value);
    double getValue(int x, int y);

    void copyFromData(Matrix &D);

    //Debugging
    void printData();

    void printDataToFile(std::ofstream &fout1);
    //void printDataToFile(std::string s);

    //Some useful functions
    double normL2();
    double max();
    void reverse_str(int i, int j);
    void multCoef_str(int i, double coef);

    typePoint type(int i, int j);
    int getPointIndex(int i, int j);

    //
    std::vector< std::vector<int> > nonZeroElementsPos;
    void fillNonZeroElementsPos(double eps);


};

bool checkSizes(Matrix& D1, Matrix& D2);

// operations
double dotProduct(Matrix& D1, Matrix& D2);
void summAlpha(Matrix& D1, Matrix& D2, double alpha, Matrix& out);
void summAlpha_str(Matrix &D1, int i, int k, double alpha);
void trans(Matrix& D1, Matrix& D2);
void sqr(Matrix& D1);
void multiplication(Matrix &D1, Matrix &D2, Matrix &D3);


//for testing
Matrix Gauss(Matrix& A);

#endif // VECTOR2D_H
