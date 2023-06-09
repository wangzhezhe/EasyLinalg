//eigen value and eigen vectors

#ifndef EASY_LIALG_EIGEN
#define EASY_LIALG_EIGEN

#include "./basic.h"
#include "./blas.h"

template <typename T,
          uint NumRow,
          uint NumCol>
void Householder(
    Matrix<T, NumRow, NumCol> &A,
    Matrix<T, NumRow, NumCol> &Q,
    Matrix<T, NumRow, NumCol> &R)
{
    // current code only works for the symetry matrix
    assert(NumRow == NumCol);
    const int msize = NumRow;
    // Matrix<T, NumRow, NumCol> MList[msize];
    Matrix<T, NumRow, NumCol> MTemp;

    // operator both for const and non-const
    // MList[0] = A;
    MTemp = A;

    R = A;
    double ANorm;

    Vec<T, msize> Ve;
    Ve.InitZero();
    Vec<T, msize> Va;
    Va.InitZero();

    for (uint k = 1; k < msize; k++)
    {
        auto HInit = R.GetMinor(k - 1);
        // printf("hhd iter%d\n", k);
        // HInit.Show();

        Va = HInit.GetColum(k - 1);

        // puts("a colum");
        // Va.Show();

        ANorm = Va.GetNorm();

        // printf("a norm %f\n", ANorm);

        if (A[k - 1][k - 1] > 0)
            ANorm = -ANorm;

        for (uint i = 0; i < msize; i++)
            Ve[i] = (i == (k - 1)) ? 1 : 0;

        // puts("ve is");
        // Ve.Show();
        // puts("va is");
        // Va.Show();

        // Ve = ANorm*Ve+Va
        Ve = AXPY(ANorm, Ve, Va);
        // puts("after vmadd v vec");
        // Ve.Show();

        // TODO, only do this when norm is not zero
        if (Ve.GetNorm() > 0.00001)
        {
            Ve.Div(Ve.GetNorm());
        }

        // puts("after vdiv v vec");
        // Ve.Show();

        // I-v*v
        // MList[k] = IMINUSVVT(Ve);
        MTemp = IMINUSVVT(Ve);
        // puts("mat_list[k] H is");
        // MTemp.Show();

        // update Q
        if (k == 1)
        {
            // Q = MList[k];
            Q = MTemp;
        }
        else
        {
            // Q = MatrixMultiply(Q, MList[k]);
            Q = MMMultiply(Q, MTemp);
        }

        // update R
        // puts("before updating R");
        // puts("MList[k]");
        // MTemp.Show();

        // puts("R");
        // R.Show();

        // R = MatrixMultiply(MList[k], R);
        R = MMMultiply(MTemp, R);

        // puts("update R is");
        // R.Show();
    }
}

template <typename T>
bool SymmInvertMatrixInner(const Matrix<T, 3, 3> &m, Matrix<T, 3, 3> &inv_m)
{
    inv_m[0][0] = m[1][1] * m[2][2] - m[1][2] * m[2][1];

    inv_m[1][0] = m[1][2] * m[2][0] - m[1][0] * m[2][2];

    inv_m[2][0] = m[1][0] * m[2][1] - m[1][1] * m[2][0];

    inv_m[0][1] = m[0][2] * m[2][1] - m[0][1] * m[2][2];

    inv_m[1][1] = m[0][0] * m[2][2] - m[0][2] * m[2][0];

    inv_m[2][1] = m[0][1] * m[2][0] - m[0][0] * m[2][1];

    inv_m[0][2] = m[0][1] * m[1][2] - m[0][2] * m[1][1];

    inv_m[1][2] = m[0][2] * m[1][0] - m[0][0] * m[1][2];

    inv_m[2][2] = m[0][0] * m[1][1] - m[0][1] * m[1][0];

    double det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

    if (det == 0)
    {
        printf("singular matrix\n");
        // matrix_show(m);
        return false;
    }

    det = 1.0 / det;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            inv_m[i][j] = inv_m[i][j] * det;
        }
    }

    return true;
}

template <typename T>
bool SymmInvertMatrixInner(const Matrix<T, 4, 4> &m, Matrix<T, 4, 4> &inv_m)
{
    inv_m[0][0] = m[1][1] * m[2][2] * m[3][3] -
                  m[1][1] * m[2][3] * m[3][2] -
                  m[2][1] * m[1][2] * m[3][3] +
                  m[2][1] * m[1][3] * m[3][2] +
                  m[3][1] * m[1][2] * m[2][3] -
                  m[3][1] * m[1][3] * m[2][2];

    inv_m[1][0] = -m[1][0] * m[2][2] * m[3][3] +
                  m[1][0] * m[2][3] * m[3][2] +
                  m[2][0] * m[1][2] * m[3][3] -
                  m[2][0] * m[1][3] * m[3][2] -
                  m[3][0] * m[1][2] * m[2][3] +
                  m[3][0] * m[1][3] * m[2][2];

    inv_m[2][0] = m[1][0] * m[2][1] * m[3][3] -
                  m[1][0] * m[2][3] * m[3][1] -
                  m[2][0] * m[1][1] * m[3][3] +
                  m[2][0] * m[1][3] * m[3][1] +
                  m[3][0] * m[1][1] * m[2][3] -
                  m[3][0] * m[1][3] * m[2][1];

    inv_m[3][0] = -m[1][0] * m[2][1] * m[3][2] +
                  m[1][0] * m[2][2] * m[3][1] +
                  m[2][0] * m[1][1] * m[3][2] -
                  m[2][0] * m[1][2] * m[3][1] -
                  m[3][0] * m[1][1] * m[2][2] +
                  m[3][0] * m[1][2] * m[2][1];

    inv_m[0][1] = -m[0][1] * m[2][2] * m[3][3] +
                  m[0][1] * m[2][3] * m[3][2] +
                  m[2][1] * m[0][2] * m[3][3] -
                  m[2][1] * m[0][3] * m[3][2] -
                  m[3][1] * m[0][2] * m[2][3] +
                  m[3][1] * m[0][3] * m[2][2];

    inv_m[1][1] = m[0][0] * m[2][2] * m[3][3] -
                  m[0][0] * m[2][3] * m[3][2] -
                  m[2][0] * m[0][2] * m[3][3] +
                  m[2][0] * m[0][3] * m[3][2] +
                  m[3][0] * m[0][2] * m[2][3] -
                  m[3][0] * m[0][3] * m[2][2];

    inv_m[2][1] = -m[0][0] * m[2][1] * m[3][3] +
                  m[0][0] * m[2][3] * m[3][1] +
                  m[2][0] * m[0][1] * m[3][3] -
                  m[2][0] * m[0][3] * m[3][1] -
                  m[3][0] * m[0][1] * m[2][3] +
                  m[3][0] * m[0][3] * m[2][1];

    inv_m[3][1] = m[0][0] * m[2][1] * m[3][2] -
                  m[0][0] * m[2][2] * m[3][1] -
                  m[2][0] * m[0][1] * m[3][2] +
                  m[2][0] * m[0][2] * m[3][1] +
                  m[3][0] * m[0][1] * m[2][2] -
                  m[3][0] * m[0][2] * m[2][1];

    inv_m[0][2] = m[0][1] * m[1][2] * m[3][3] -
                  m[0][1] * m[1][3] * m[3][2] -
                  m[1][1] * m[0][2] * m[3][3] +
                  m[1][1] * m[0][3] * m[3][2] +
                  m[3][1] * m[0][2] * m[1][3] -
                  m[3][1] * m[0][3] * m[1][2];

    inv_m[1][2] = -m[0][0] * m[1][2] * m[3][3] +
                  m[0][0] * m[1][3] * m[3][2] +
                  m[1][0] * m[0][2] * m[3][3] -
                  m[1][0] * m[0][3] * m[3][2] -
                  m[3][0] * m[0][2] * m[1][3] +
                  m[3][0] * m[0][3] * m[1][2];

    inv_m[2][2] = m[0][0] * m[1][1] * m[3][3] -
                  m[0][0] * m[1][3] * m[3][1] -
                  m[1][0] * m[0][1] * m[3][3] +
                  m[1][0] * m[0][3] * m[3][1] +
                  m[3][0] * m[0][1] * m[1][3] -
                  m[3][0] * m[0][3] * m[1][1];

    inv_m[3][2] = -m[0][0] * m[1][1] * m[3][2] +
                  m[0][0] * m[1][2] * m[3][1] +
                  m[1][0] * m[0][1] * m[3][2] -
                  m[1][0] * m[0][2] * m[3][1] -
                  m[3][0] * m[0][1] * m[1][2] +
                  m[3][0] * m[0][2] * m[1][1];

    inv_m[0][3] = -m[0][1] * m[1][2] * m[2][3] +
                  m[0][1] * m[1][3] * m[2][2] +
                  m[1][1] * m[0][2] * m[2][3] -
                  m[1][1] * m[0][3] * m[2][2] -
                  m[2][1] * m[0][2] * m[1][3] +
                  m[2][1] * m[0][3] * m[1][2];

    inv_m[1][3] = m[0][0] * m[1][2] * m[2][3] -
                  m[0][0] * m[1][3] * m[2][2] -
                  m[1][0] * m[0][2] * m[2][3] +
                  m[1][0] * m[0][3] * m[2][2] +
                  m[2][0] * m[0][2] * m[1][3] -
                  m[2][0] * m[0][3] * m[1][2];

    inv_m[2][3] = -m[0][0] * m[1][1] * m[2][3] +
                  m[0][0] * m[1][3] * m[2][1] +
                  m[1][0] * m[0][1] * m[2][3] -
                  m[1][0] * m[0][3] * m[2][1] -
                  m[2][0] * m[0][1] * m[1][3] +
                  m[2][0] * m[0][3] * m[1][1];

    inv_m[3][3] = m[0][0] * m[1][1] * m[2][2] -
                  m[0][0] * m[1][2] * m[2][1] -
                  m[1][0] * m[0][1] * m[2][2] +
                  m[1][0] * m[0][2] * m[2][1] +
                  m[2][0] * m[0][1] * m[1][2] -
                  m[2][0] * m[0][2] * m[1][1];

    double det = m[0][0] * inv_m[0][0] +
                 m[0][1] * inv_m[1][0] +
                 m[0][2] * inv_m[2][0] +
                 m[0][3] * inv_m[3][0];

    if (det == 0)
    {
        printf("singular matrix for 4by4\n");
        // matrix_show(m);
        return false;
    }

    det = 1.0 / det;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            inv_m[i][j] = inv_m[i][j] * det;
        }
    }

    return true;
}

template <typename T, uint Size>
bool SymmInvertMatrixInnerGeneral(const Matrix<T, Size, Size> &m, Matrix<T, Size, Size> &inv_m)
{
    puts("TODO");
    return true;
}

template <typename T,
          uint Size>
bool SymmInvertMatrix(const Matrix<T, Size, Size> &A, Matrix<T, Size, Size> &AInv)
{

    if (Size == 3)
    {
        SymmInvertMatrixInner(A, AInv);
        return true;
    }
    else if (Size == 4)
    {
        SymmInvertMatrixInner(A, AInv);
        return true;
    }
    else
    {
        // TODO adding size = 2
        SymmInvertMatrixInnerGeneral(A, AInv);
    }
    return false;
}

template <typename T, uint Size>
void SymmEigenValues(const Matrix<T, Size, Size> &A, double tol, int maxIter, Vec<T, Size> &eigenValues)
{
    // refer to https://www.andreinc.net/2021/01/25/computing-eigenvalues-and-eigenvectors-using-qr-decomposition
    // mat_t ak = *x;
    Matrix<T, Size, Size> Ak = A;

    // mat_t qq = matrix_new_eye(x->m, x->m);
    Matrix<T, Size, Size> QQ;
    QQ.InitEye();

    Matrix<T, Size, Size> R, Q;
    for (int i = 0; i < maxIter; ++i)
    {

        Householder(Ak, Q, R);
        // printf("iter i %d \nR\n", i);
        // R.Show();
        // puts("Q");
        // Q.Show();

        Ak = MMMultiply(R, Q);
        QQ = MMMultiply(QQ, Q);

        if (Ak.IsUpperTriangular(tol))
        {
            break;
        }
    }

    for (int i = 0; i < Size; i++)
    {
        if (fabs(Ak[i][i] - 0) < tol)
        {
            eigenValues[i] = 0.0;
        }
        else
        {
            eigenValues[i] = Ak[i][i];
        }
    }
}

// using inverse itertaion to compute the eigen vectors
// https://en.wikipedia.org/wiki/Inverse_iteration
template <typename T, uint Size>
Matrix<T, Size, Size> SymmEigenVectors(const Matrix<T, Size, Size> &A, const Vec<T, Size> &eigenValues, int maxIter)
{

    Matrix<T, Size, Size> eigenVectors;
    Matrix<T, Size, Size> A_minus_lambda_i;
    Matrix<T, Size, Size> A_minus_lambda_inv;
    // go through each eigen values
    for (uint j = 0; j < Size; j++)
    {
        if (fabs(eigenValues[j] - 0.0) < 0.00001)
        {
            // the associated eigen vector is all zero
            continue;
        }

        double lambda = eigenValues[j] + 0.000001;
        for (uint ii = 0; ii < Size; ii++)
        {
            for (uint jj = 0; jj < Size; jj++)
            {
                if (ii == jj)
                {
                    A_minus_lambda_i[ii][jj] = A[ii][jj] - lambda;
                }
                else
                {
                    A_minus_lambda_i[ii][jj] = A[ii][jj];
                }
            }
        }

        // puts("a-lambda*I");
        // A_minus_lambda_i.Show();

        bool invertok = SymmInvertMatrix(A_minus_lambda_i, A_minus_lambda_inv);
        assert(invertok == true);

        // puts("a-lambda*I inv");
        // A_minus_lambda_inv.Show();

        Vec<T, Size> bPrev;
        Vec<T, Size> bCurr;
        // init as 1
        for (int i = 0; i < Size; i++)
        {
            bPrev[i] = 1.0;
        }
        uint iterNum = 0;
        while (iterNum < maxIter)
        {
            bCurr = MMVultiply(A_minus_lambda_inv, bPrev);
            // puts("check mmv");
            // bCurr.Show();

            bCurr.Normalize();
            // puts("check norm");
            // bCurr.Show();

            if (bCurr.Equal(bPrev))
            {
                // results convert
                break;
            }
            bPrev = bCurr;
            iterNum++;
        }

        // assign to the matrix
        for (uint i = 0; i < Size; i++)
        {
            eigenVectors[i][j] = bCurr[i];
        }
    }
    return eigenVectors;
}

template <typename T, uint Size>
Matrix<T, Size, Size> SymmEigenDecomposition(const Matrix<T, Size, Size> &A, double tol, int maxIter)
{

    // solve eigen values
    Vec<T, Size> eigenValues;
    eigenValues.InitZero();
    SymmEigenValues(A, tol, 20, eigenValues);

    eigenValues.Show();

    // solve eigen vectors
    Matrix<T, Size, Size> eigenVactors;
    eigenVactors = SymmEigenVectors(A, eigenValues, maxIter);

    // create diagonal matrix
    Matrix<T, Size, Size> diag;
    for (uint i = 0; i < Size; i++)
    {
        diag[i][i] = sqrt(eigenValues[i]);
    }

    Matrix<T, Size, Size> decompA;
    decompA = MMMultiply(eigenVactors, diag);
    return decompA;
}

#endif