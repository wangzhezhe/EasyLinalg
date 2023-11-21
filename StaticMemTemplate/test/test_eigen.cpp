#include <include/basic.h>
#include <include/eigen.h>
#include "test_data.h"
#include <assert.h>

using namespace EASYLINALG;

void testInit()
{
    // test vector
    Vec<int, 10> v;
    // the initial array may not be all zero
    v[0] = 1;
    std::cout << "index 0 is " << v[0] << std::endl;
    v.Show();
    v.InitZero();
    v.Show();

    // test matrix
    Matrix<int, 3, 3> m;
    m[0][0] = 1;
    std::cout << "index 0 0 is " << m[0][0] << std::endl;
    m.InitEye();
    m.Show();
    m.Transpose();

    Matrix<int, 3, 3> m1;
    m1.InitEye();
    Matrix<int, 3, 3> m2;
    m2.InitEye();
    // the size can be determined automatically
    auto m3 = MMMultiply(m1, m2);
    m3.Show();
}

template <uint Num>
void testMMVPV()
{
    Matrix<int, Num, Num> A;
    Vec<int, Num> U;
    Vec<int, Num> M;

    for (uint i = 0; i < Num; i++)
    {
        for (uint j = 0; j < Num; j++)
        {
            A[i][j] = i + 1;
        }
        U[i] = 1;
        M[i] = i + 1;
    }

    auto aum = MMVPV(A, U, M);
    for (uint i = 0; i < Num; i++)
    {
        assert(aum[i] == static_cast<int>((i + 1) * (Num + 1)));
    }
}

void testAssignment()
{
    Vec<int, 10> v1;
    for (int i = 0; i < 10; i++)
    {
        v1[i] = i;
    }
    Vec<int, 10> v2;
    v2 = v1;
    for (int i = 0; i < 10; i++)
    {
        assert(v1[i] == v2[i]);
    }

    Matrix<float, 10, 10> m1;
    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 10; j++)
        {
            m1[i][j] = i * j * 0.5 + 1;
        }
    }

    auto m2 = m1;
    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 10; j++)
        {
            assert(m2[i][j] == m1[i][j]);
        }
    }
}

void testBasicOperations()
{
    testAssignment();
    testMMVPV<3>();
    testMMVPV<4>();
    testMMVPV<8>();
}

template <uint Num>
void testQRInner()
{

    printf("---basic_qr_test\n");
    // mat_t R, Q;
    Matrix<double, Num, Num> x;
    Matrix<double, Num, Num> Q;
    Matrix<double, Num, Num> R;

    for (uint i = 0; i < Num; i++)
    {
        for (uint j = 0; j < Num; j++)
        {
            if (Num == 3)
            {
                // x[i][j] = in33_0[i][j];
                x[i][j] = in33_1[i][j];
            }
            else if (Num == 4)
            {
                x[i][j] = in44_0[i][j];
            }
            else if (Num == 8)
            {
                x[i][j] = in88_0[i][j];
            }
            else
            {
                puts("unsupported num value");
                assert(false);
            }
        }
    }

    puts("original matrix");
    x.Show();

    Householder(x, Q, R);

    puts("Q");
    Q.Show();
    puts("R");
    R.Show();

    // R is upper triangular matrix
    // how to decide the value of epsilon?
    bool ifupper = R.IsUpperTriangular(0.00001);
    assert(ifupper == true);

    // Q is check orthoganal matrix
    auto Qt = Q;
    Qt.Transpose();
    auto QtQ = MMMultiply(Qt, Q);
    puts("QtQ is");
    QtQ.Show();

    // the QtQ is supposed to be identical matrix
    for (uint i = 0; i < Num; i++)
    {
        for (uint j = 0; j < Num; j++)
        {
            if (i == j)
            {
                assert(fabs(QtQ[i][j] - 1) < 0.00001);
            }
            else
            {
                assert(fabs(QtQ[i][j] - 0) < 0.00001);
            }
        }
    }

    // checking results of Q*R
    auto m = MMMultiply(Q, R);

    puts("Q * R");
    m.Show();
    assert(m.IsEqual(x) == true);
}

void testQR()
{
    testQRInner<3>();
    testQRInner<4>();
    testQRInner<8>();
}

template <uint Num, void EigenFunc(const Matrix<double, Num, Num> &, double, int, Vec<double, Num> &)>
void testEVNQRInner()
{
    Matrix<double, Num, Num> m1;
    if (Num == 3)
    {
        for (uint i = 0; i < Num; i++)
        {
            for (uint j = 0; j < Num; j++)
            {
                m1[i][j] = in33_1[i][j];
            }
        }
        Vec<double, Num> eigenValues;
        m1.Show();
        EigenFunc(m1, 0.0001, 20, eigenValues);

        for (int i = 0; i < 3; i++)
        {
            printf("%f ", eigenValues[i]);
        }

        printf("\n");
        // the assert only works for debug case
        assert(fabs(eigenValues[0] - 2.8) < 0.00001);
        assert(fabs(eigenValues[1] - 0.0) < 0.00001);
        assert(fabs(eigenValues[2] - 0.0) < 0.00001);
    }
    if (Num == 4)
    {
        for (uint i = 0; i < Num; i++)
        {
            for (uint j = 0; j < Num; j++)
            {
                m1[i][j] = in44_1[i][j];
            }
        }
        Vec<double, Num> eigenValues;
        m1.Show();
        EigenFunc(m1, 0.0001, 20, eigenValues);

        for (int i = 0; i < 4; i++)
        {
            printf("%f ", eigenValues[i]);
        }

        printf("\n");

        // sort array from large to small
        eigenValues.Sort();
        // the assert only works for debug case
        assert(fabs(eigenValues[0] - 6.0) < 0.00001);
        assert(fabs(eigenValues[1] - 0.0) < 0.00001);
        assert(fabs(eigenValues[2] - 0.0) < 0.00001);
        assert(fabs(eigenValues[3] - 0.0) < 0.00001);
    }
    if (Num == 8)
    {
        for (uint i = 0; i < Num; i++)
        {
            for (uint j = 0; j < Num; j++)
            {
                m1[i][j] = in88_0[i][j];
            }
        }
        Vec<double, Num> eigenValues;
        m1.Show();
        EigenFunc(m1, 0.0001, 50, eigenValues);

        eigenValues.Sort();

        std::cout << "--show eigen values" << std::endl;
        for (uint i = 0; i < Num; i++)
        {
            printf("%f ", eigenValues[i]);
        }

        printf("\n");
        // the assert only works for debug case
        assert(fabs(eigenValues[0] - 9.678793) < 0.001);
        assert(fabs(eigenValues[1] - 7.1734) < 0.00001);
        assert(fabs(eigenValues[2] - 6.16789) < 0.00001);
        assert(fabs(eigenValues[3] - 5.0000) < 0.001);
        assert(fabs(eigenValues[4] - 4.0000) < 0.00001);
        assert(fabs(eigenValues[5] - 2.8321) < 0.00001);
        assert(fabs(eigenValues[6] - 1.82659) < 0.00001);
        assert(fabs(eigenValues[7] - (-0.678903)) < 0.00001);
    }
}

void testEigenValuesNaiveQR()
{
    testEVNQRInner<3, SymmEigenValues>();
    testEVNQRInner<4, SymmEigenValues>();
    testEVNQRInner<8, SymmEigenValues>();
}

void testEigenValuesShift()
{
    std::cout << "--testing SymmEigenValuesShift---" << std::endl;
    testEVNQRInner<3, SymmEigenValuesShift>();
    testEVNQRInner<4, SymmEigenValuesShift>();
    // the value here converge faster and more accurate
    testEVNQRInner<8, SymmEigenValuesShift>();
}

void testED3by3()
{
    printf("---test_eigen_vectors_decomposition 3by3\n");
    Matrix<double, 3, 3> x;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            x[i][j] = in33_3[i][j];
        }
    }

    Matrix<double, 3, 3> A = SymmEigenDecomposition(x, 0.0001, 50);
    Matrix<double, 3, 3> ATrans = A;
    ATrans.Transpose();
    Matrix<double, 3, 3> rst = MMMultiply(A, ATrans);
    puts("A");
    A.Show();
    puts("rst");
    rst.Show();
    assert(x.IsEqual(rst) == true);
}

void testED4by4()
{
    printf("---test_eigen_vectors_decomposition 4by4\n");
    Matrix<double, 4, 4> x1;
    Matrix<double, 4, 4> x2;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            x1[i][j] = in44_3[i][j];
            x2[i][j] = in44_4[i][j];
        }
    }

    Matrix<double, 4, 4> A = SymmEigenDecomposition(x1, 0.00001, 20);
    Matrix<double, 4, 4> ATrans = A;
    ATrans.Transpose();
    Matrix<double, 4, 4> rst = MMMultiply(A, ATrans);
    puts("A");
    A.Show();
    puts("rst");
    rst.Show();
    assert(x1.IsEqual(rst) == true);

    A = SymmEigenDecomposition(x2, 0.00001, 20);
    ATrans = A;
    ATrans.Transpose();
    rst = MMMultiply(A, ATrans);
    puts("A");
    A.Show();
    puts("rst");
    rst.Show();
    assert(x2.IsEqual(rst) == true);
}

void testED8by8()
{
    printf("---test_eigen_vectors_decomposition 8by8\n");
    constexpr uint size = 8;
    Matrix<double, size, size> x1;
    for (uint i = 0; i < size; i++)
    {
        for (uint j = 0; j < size; j++)
        {
            x1[i][j] = in88_1[i][j];
        }
    }

    Matrix<double, size, size> A = SymmEigenDecomposition(x1, 0.00001, 30);
    Matrix<double, size, size> ATrans = A;
    ATrans.Transpose();
    Matrix<double, size, size> rst = MMMultiply(A, ATrans);
    puts("A");
    A.Show();
    puts("rst");
    rst.Show();

    // TODO,
    // current precision is a little bit low
    // how to improve this precision further?
    assert(x1.IsEqual(rst, 0.001) == true);
}

void testEigenDecomposition()
{
    testED3by3();
    testED4by4();
    testED8by8();
}

void testDGEMV()
{
    Matrix<double, 4, 4> A;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            A[i][j] = i * 4 + j;
        }
    }

    A.Show();

    Vec<double, 4> y;
    Vec<double, 4> x;
    for (int i = 0; i < 4; i++)
    {
        x[i] = 1.0;
        y[i] = 1.0;
    }

    Vec<double, 4> result = DGEMV(0.5, A, x, 0.5, y);
    result.Show();
}

void testInverse8by8()
{
    puts("---testInverse8by8---");
    constexpr uint size = 8;
    Matrix<double, size, size> x;
    Matrix<double, size, size> x_inv;
    for (uint i = 0; i < size; i++)
    {
        for (uint j = 0; j < size; j++)
        {
            x[i][j] = in88_0[i][j];
        }
    }

    bool invertOk = SymmInvertMatrixInner(x, x_inv);

    assert(invertOk == true);

    puts("x");
    x.Show();
    puts("x_inv");
    x_inv.Show();

    auto result1 = MMMultiply(x, x_inv);

    puts("x*x_inv");
    result1.Show();

    auto result2 = MMMultiply(x_inv, x);

    puts("x_inv*x");
    result2.Show();

    for (uint i = 0; i < size; i++)
    {
        for (uint j = 0; j < size; j++)
        {
            if (i == j)
            {
                assert(fabs(result1[i][j] - 1.0) < 0.00001);
                assert(fabs(result2[i][j] - 1.0) < 0.00001);
            }
            else
            {
                assert(fabs(result1[i][j] - 0.0) < 0.00001);
                assert(fabs(result2[i][j] - 0.0) < 0.00001);
            }
        }
    }
}

void testInverse()
{
    testInverse8by8();
}

void testEigenDecompositionComplicated()
{
    printf("---test_eigen_vectors_decomposition testEigenDecompositionComplicated 8by8\n");
    constexpr uint size = 8;
    Matrix<double, size, size> x1;
    for (uint i = 0; i < size; i++)
    {
        for (uint j = 0; j < size; j++)
        {
            x1[i][j] = in88_2[i][j];
        }
    }

    Matrix<double, size, size> A = SymmEigenDecomposition(x1, 0.00001, 1000);
    Matrix<double, size, size> ATrans = A;
    ATrans.Transpose();
    Matrix<double, size, size> rst = MMMultiply(A, ATrans);
    puts("A");
    A.Show();
    puts("rst");
    rst.Show();

    // TODO,
    // current precision is a little bit low
    // how to improve this precision further?
    assert(x1.IsEqual(rst, 0.001) == true);

    Matrix<double, size, size> x3;
    for (uint i = 0; i < size; i++)
    {
        for (uint j = 0; j < size; j++)
        {
            x3[i][j] = in88_3[i][j];
        }
    }

    Matrix<double, size, size> A3 = SymmEigenDecomposition(x3, 0.00001, 100);
    Matrix<double, size, size> ATrans3 = A3;
    ATrans3.Transpose();
    Matrix<double, size, size> rst3 = MMMultiply(A3, ATrans3);
    puts("A3");
    A3.Show();
    puts("rst3");
    rst3.Show();

    // TODO,
    // current precision is a little bit low
    // how to improve this precision further?
    assert(x3.IsEqual(rst3, 0.001) == true);
}

void testEigenDecompositionComplicatedShift()
{
    std::cout << "---testEigenDecompositionComplicatedShift---" << std::endl;
    constexpr uint size = 8;
    Matrix<double, size, size> x3;
    for (uint i = 0; i < size; i++)
    {
        for (uint j = 0; j < size; j++)
        {
            x3[i][j] = in88_3[i][j];
        }
    }
    // for this A3 matrix, we need to improve iteration to 1000 to get accurate eigen value
    // the tolerance is around 0.0005
    // there is one value at the corner, it is hard to converge to the value under the 0.0005
    // after 500 iteration step
    Matrix<double, size, size> A3 = SymmEigenDecomposition(x3, 0.0005, 1000, true);
    Matrix<double, size, size> ATrans3 = A3;
    ATrans3.Transpose();
    Matrix<double, size, size> rst3 = MMMultiply(A3, ATrans3);
    puts("A3");
    A3.Show();
    puts("rst3");
    rst3.Show();

    // TODO,
    // current precision is a little bit low
    // how to improve this precision further?
    assert(x3.IsEqual(rst3, 0.001) == true);
}

void testEigenGetVector()
{
    std::cout << "---testEigenGetVector---" << std::endl;
    // input eigen value and associated matrix
    constexpr uint size = 4;
    Matrix<double, size, size> x;
    for (uint i = 0; i < size; i++)
    {
        for (uint j = 0; j < size; j++)
        {
            x[i][j] = in44_1[i][j];
        }
    }

    // test eigen
    //    LIAG_FUNC_MACRO Vec<T, Size> ComputeEigenVectors(const Matrix<T, Size, Size> &A, const T &eigenValue)
    Vec<double, size> eigenVector;
    double eigenValue = 6.0;
    eigenVector = ComputeEigenVectors(x, eigenValue, 1000);

    // eigenVector.Show();
    // 0.182574 0.547723 0.365148 0.730297
    //  get eigen vector
    assert(fabs(eigenVector[0] - 0.182574) < 0.00001);
    assert(fabs(eigenVector[1] - 0.547723) < 0.00001);
    assert(fabs(eigenVector[2] - 0.365148) < 0.00001);
    assert(fabs(eigenVector[3] - 0.730297) < 0.00001);
}

int main()
{
    testInverse();
    testQR();
    testEigenValuesNaiveQR();
    testEigenValuesShift();
    testEigenDecomposition();
    testEigenGetVector();

    testEigenDecompositionComplicated();
    testEigenDecompositionComplicatedShift();
}