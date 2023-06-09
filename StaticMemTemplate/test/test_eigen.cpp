#include <include/basic.h>
#include <include/eigen.h>
#include <assert.h>

double in33_0[3][3] = {
    {1, -1, 4},
    {1, 4, -2},
    {1, 4, 2}};

double in33_1[3][3] = {
    {0.20, 0.60, 0.40},
    {0.60, 1.80, 1.20},
    {0.40, 1.20, 0.80}};

double in33_2[3][3] = {
    {1.0, 3.0, 7.0},
    {3.0, 2.0, 6.0},
    {7.0, 6.0, 5.0}};

double in33_3[3][3] = {
    {0.0, 0.0, 0.0},
    {0.0, 20.0, 60.0},
    {0.0, 60.0, 180.0}};

double in44_0[4][4] = {
    {1, -1, 4, 1},
    {1, 4, -2, 1},
    {1, 4, 2, 1},
    {1, -1, 0, 1}};

double in44_1[4][4] = {
    {0.20, 0.60, 0.40, 0.80},
    {0.60, 1.80, 1.20, 2.40},
    {0.40, 1.20, 0.80, 1.60},
    {0.80, 2.40, 1.60, 3.20}};

double in44_2[4][4] = {
    {1.0, 3.0, 7.0, 8.0},
    {3.0, 2.0, 6.0, 7.0},
    {7.0, 6.0, 5.0, 6.0},
    {8.0, 7.0, 6.0, 5.0}};

double in44_3[4][4] = {
    {0.0, 0.0, 0.0, 0.0},
    {0.0, 20.0, 60.0, 40.0},
    {0.0, 60.0, 180.0, 120.0},
    {0.0, 40.0, 120.0, 80.0}};

double in44_4[4][4] = {
    {0.0024507, 0.0026104, 0.0019078, 0.0025373},
    {0.0026104, 0.0028755, 0.0020373, 0.0027783},
    {0.0019078, 0.0020373, 0.0015548, 0.0020305},
    {0.0025373, 0.0027783, 0.0020305, 0.0027494}};

double in88_0[8][8] = {
    {1, 1, 0, 1, 0, 1, 0, 1},
    {1, 2, 1, 0, 1, 0, 1, 0},
    {0, 1, 3, 1, 0, 1, 0, 1},
    {1, 0, 1, 4, 1, 0, 1, 0},
    {0, 1, 0, 1, 5, 1, 0, 1},
    {1, 0, 1, 0, 1, 6, 1, 0},
    {0, 1, 0, 1, 0, 1, 7, 1},
    {1, 0, 1, 0, 1, 0, 1, 8}};

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

    for (int i = 0; i < Num; i++)
    {
        for (int j = 0; j < Num; j++)
        {
            A[i][j] = i + 1;
        }
        U[i] = 1;
        M[i] = i + 1;
    }

    auto aum = MMVPV(A, U, M);
    for (int i = 0; i < Num; i++)
    {
        assert(aum[i] == (i + 1) * (Num + 1));
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

    for (int i = 0; i < Num; i++)
    {
        for (int j = 0; j < Num; j++)
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
    for (int i = 0; i < Num; i++)
    {
        for (int j = 0; j < Num; j++)
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

template <uint Num>
void testEVNQRInner()
{
    Matrix<double, Num, Num> m1;
    if (Num == 3)
    {
        for (int i = 0; i < Num; i++)
        {
            for (int j = 0; j < Num; j++)
            {
                m1[i][j] = in33_1[i][j];
            }
        }
        Vec<double, Num> eigenValues;
        m1.Show();
        SymmEigenValues(m1, 0.0001, 20, eigenValues);

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
        for (int i = 0; i < Num; i++)
        {
            for (int j = 0; j < Num; j++)
            {
                m1[i][j] = in44_1[i][j];
            }
        }
        Vec<double, Num> eigenValues;
        m1.Show();
        SymmEigenValues(m1, 0.0001, 20, eigenValues);

        for (int i = 0; i < 4; i++)
        {
            printf("%f ", eigenValues[i]);
        }

        printf("\n");
        // the assert only works for debug case
        assert(fabs(eigenValues[0] - 6) < 0.00001);
        assert(fabs(eigenValues[1] - 0.0) < 0.00001);
        assert(fabs(eigenValues[2] - 0.0) < 0.00001);
        assert(fabs(eigenValues[3] - 0.0) < 0.00001);
    }
    if (Num == 8)
    {
        for (int i = 0; i < Num; i++)
        {
            for (int j = 0; j < Num; j++)
            {
                m1[i][j] = in88_0[i][j];
            }
        }
        Vec<double, Num> eigenValues;
        m1.Show();
        SymmEigenValues(m1, 0.0001, 20, eigenValues);

        for (int i = 0; i < Num; i++)
        {
            printf("%f ", eigenValues[i]);
        }

        printf("\n");
        // the assert only works for debug case
        assert(fabs(eigenValues[0] - 9.6789) < 0.00001);
        assert(fabs(eigenValues[1] - 7.1717) < 0.00001);
        assert(fabs(eigenValues[2] - 6.169579) < 0.00001);
        assert(fabs(eigenValues[3] - 4.999833) < 0.00001);
        assert(fabs(eigenValues[4] - 4.000179) < 0.00001);
        assert(fabs(eigenValues[5] - 2.832101) < 0.00001);
        assert(fabs(eigenValues[6] - 1.826598) < 0.00001);
        assert(fabs(eigenValues[7] - (-0.678903)) < 0.00001);
    }
}

void testEigenValuesNaiveQR()
{
    testEVNQRInner<3>();
    testEVNQRInner<4>();
    testEVNQRInner<8>();
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

    Matrix<double, 3, 3> A = SymmEigenDecomposition(x, 0.00001, 20);
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

void testEigenDecomposition()
{
    testED3by3();
    testED4by4();
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

int main()
{
    testInit();
    testBasicOperations();
    testQR();

    testEigenValuesNaiveQR();
    testEigenDecomposition();

    testDGEMV();
}