#include <include/matrixdet.h>
#include "test_data.h"
#include <assert.h>

using namespace EASYLINALG;

void testMatrixDet()
{
    // 3*3
    Matrix<double, 3, 3> x3;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            x3[i][j] = in33_0[i][j];
        }
    }
    double det = SymmDetByEigenValues(x3);
    std::cout << "det 3*3: " << det << std::endl;
    assert(fabs(det - 19.99999) < 0.00001);

    // 4*4
    Matrix<double, 4, 4> x41;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            x41[i][j] = in44_1[i][j];
        }
    }
    det = SymmDetByEigenValues(x41);
    std::cout << "det 4*4_1: " << det << std::endl;
    assert(fabs(det - 0.0) < 0.00001);

    Matrix<double, 4, 4> x42;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            x42[i][j] = in44_5[i][j];
        }
    }
    det = SymmDetByEigenValues(x42);
    std::cout << "det 4*4_2: " << det << std::endl;
    assert(fabs(det - -0.005100) < 0.00001);

    // 8*8
    Matrix<double, 8, 8> x8;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            x8[i][j] = in88_0[i][j];
        }
    }
    det = SymmDetByEigenValues(x8);
    std::cout << "det 8*8: " << det << std::endl;
    assert(fabs(det - 0.0) < 0.00001);
}

int main()
{
    testMatrixDet();
}