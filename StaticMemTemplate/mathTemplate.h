#ifndef MATHTEMPLATE
#define MATHTEMPLATE

#include <iostream>
#include <assert.h>
#include <cmath>
#include <type_traits>

template <typename T, uint Size>
class Vec
{
public:
    static const int NUM_COMPONENTS = Size;
    T Components[NUM_COMPONENTS];

    // getting element by [] operator
    // return a reference, so the value can be updated further
    T &operator[](uint index)
    {
        assert(index >= 0);
        assert(index < this->NUM_COMPONENTS);
        return this->Components[index];
    }

    // if not adding this, there are some error such as
    /*
    /Users/zw1/Documents/cworkspace/src/EasyLinalg/StaticMemTemplate/mathTemplate.h:141:42: error: no viable overloaded operator[] for type 'const Vec<Vec<double, NUM_COLUMNS>, NUM_ROWS>'
                this->Components[i][j] = t.Components[i][j];
                                         ^~~~~~~~~~~~ ~
/Users/zw1/Documents/cworkspace/src/EasyLinalg/StaticMemTemplate/mathTemplate.h:343:15: note: in instantiation of member function 'Matrix<double, 3, 3>::operator=' requested here
            Q = MatrixMultiply(Q, MList[k]);
              ^
/Users/zw1/Documents/cworkspace/src/EasyLinalg/StaticMemTemplate/test.cpp:163:5: note: in instantiation of function template specialization 'Householder<double, 3U, 3U>' requested here
    Householder(x,Q,R);
    ^
/Users/zw1/Documents/cworkspace/src/EasyLinalg/StaticMemTemplate/test.cpp:187:5: note: in instantiation of function template specialization 'testQRInner<3U>' requested here
    testQRInner<3>();
    ^
/Users/zw1/Documents/cworkspace/src/EasyLinalg/StaticMemTemplate/mathTemplate.h:18:8: note: candidate function not viable: 'this' argument has type 'const Vec<Vec<double, NUM_COLUMNS>, NUM_ROWS>', but method is not marked const
    T &operator[](const int index)

    The failure happens because Non-const functions can only be called by non-const objects.
    However, when a function is declared as const, it can be called on any type of object.

    still curious about results here, some cases, the object that call the [] is a const object
    */
    const T &operator[](uint index) const
    {
        assert(index >= 0);
        assert(index < this->NUM_COMPONENTS);
        return this->Components[index];
    }

    Vec &operator=(const Vec &t)
    {
        assert(this->NUM_COMPONENTS == t.NUM_COMPONENTS);
        for (uint i = 0; i < this->NUM_COMPONENTS; i++)
        {
            this->Components[i] = t.Components[i];
        }
        return *this;
    }

    Vec &operator=(Vec &t)
    {
        assert(this->NUM_COMPONENTS == t.NUM_COMPONENTS);
        for (uint i = 0; i < this->NUM_COMPONENTS; i++)
        {
            this->Components[i] = t.Components[i];
        }
        return *this;
    }

    // this is only used for single vector
    // the vector can be in inserted form
    // so we do not set the init value of vector as zero
    void InitZero()
    {
        for (int i = 0; i < this->NUM_COMPONENTS; i++)
        {
            this->Components[i] = 0;
        }
    }

    void Show()
    {
        for (int i = 0; i < this->NUM_COMPONENTS; i++)
        {
            printf("%lf ", 1.0 * this->Components[i]);
        }
        printf("\n");
    }

    // ||x||
    double GetNorm()
    {
        double sum = 0;
        for (uint i = 0; i < this->NUM_COMPONENTS; i++)
            sum = sum + this->Components[i] * this->Components[i];
        return sqrt(sum);
    }

    // y = y/d
    void Div(double d)
    {
        // if the d is a value close to zero
        // just set results to 0
        assert((d - 0.0) > 0.00001);
        for (uint i = 0; i < Size; i++)
        {
            this->Components[i] = static_cast<T>(1.0 * this->Components[i] / d);
        }
        return;
    }
};

template <typename T, uint RowSize, uint ColSize>
class Matrix
{
public:
    static constexpr int NUM_ROWS = RowSize;
    static constexpr int NUM_COLUMNS = ColSize;
    Vec<Vec<T, NUM_COLUMNS>, NUM_ROWS> Components;

    Matrix()
    {
        for (uint i = 0; i < this->NUM_ROWS; i++)
        {
            for (uint j = 0; j < this->NUM_COLUMNS; j++)
            {
                this->Components[i][j] = 0;
            }
        }
    };

    bool IsUpperTriangular(double tol)
    {
        // For now, only treat square matricies.
        if (this->NUM_COLUMNS != this->NUM_ROWS)
        {
            return false;
        }

        for (uint i = 0; i < this->NUM_COLUMNS; i++)
        {
            for (uint j = 0; j < i; j++)
            {
                if (fabs(this->Components[i][j]) > tol)
                {
                    return false;
                }
            }
        }
        return true;
    }

    bool IsEqual(const Matrix &t)
    {
        if (this->NUM_ROWS != t.NUM_ROWS || this->NUM_COLUMNS != t.NUM_COLUMNS)
        {
            return false;
        }
        for (uint i = 0; i < this->NUM_ROWS; i++)
        {
            for (uint j = 0; j < this->NUM_COLUMNS; j++)
            {
                if (fabs(this->Components[i][j] - t.Components[i][j]) > 0.0001)
                {
                    return false;
                }
            }
        }
        return true;
    }

    Vec<T, NUM_COLUMNS> &operator[](uint rowIndex)
    {
        assert(rowIndex >= 0);
        assert(rowIndex < NUM_ROWS);
        return this->Components[rowIndex];
    }

    const Vec<T, NUM_COLUMNS> &operator[](uint rowIndex) const
    {
        assert(rowIndex >= 0);
        assert(rowIndex < NUM_ROWS);
        return this->Components[rowIndex];
    }

    Matrix &operator=(Matrix &t)
    {
        assert(this->NUM_ROWS == t.NUM_ROWS);
        assert(this->NUM_COLUMNS == t.NUM_COLUMNS);

        for (uint i = 0; i < this->NUM_ROWS; i++)
        {
            for (uint j = 0; j < this->NUM_COLUMNS; j++)
            {
                this->Components[i][j] = t.Components[i][j];
            }
        }
        return *this;
    }

    Matrix &operator=(const Matrix &t)
    {
        assert(this->NUM_ROWS == t.NUM_ROWS);
        assert(this->NUM_COLUMNS == t.NUM_COLUMNS);

        for (uint i = 0; i < this->NUM_ROWS; i++)
        {
            for (uint j = 0; j < this->NUM_COLUMNS; j++)
            {
                this->Components[i][j] = t.Components[i][j];
            }
        }
        return *this;
    }

    void InitEye()
    {
        assert(this->NUM_ROWS == this->NUM_COLUMNS);
        for (int i = 0; i < this->NUM_ROWS; i++)
        {
            for (int j = 0; j < this->NUM_COLUMNS; j++)
            {
                if (i == j)
                {
                    this->Components[i][j] = 1;
                }
                else
                {
                    this->Components[i][j] = 0;
                }
            }
        }
    }

    void Show()
    {
        for (int i = 0; i < this->NUM_ROWS; i++)
        {
            for (int j = 0; j < this->NUM_COLUMNS; j++)
            {
                printf(" %8.7f", 1.0 * this->Components[i][j]);
            }
            printf("\n");
        }
        printf("\n");
    }

    void Transpose()
    {
        for (int i = 0; i < this->NUM_ROWS; i++)
        {
            for (int j = 0; j < i; j++)
            {
                T t = this->Components[i][j];
                this->Components[i][j] = this->Components[j][i];
                this->Components[j][i] = t;
            }
        }
    }

    // only works for the symmetric matrix
    Matrix GetMinor(uint d)
    {
        assert(this->NUM_ROWS == this->NUM_COLUMNS);
        Matrix minor;
        for (uint i = 0; i < d; i++)
        {
            minor.Components[i][i] = 1;
        }

        // start form the dth position
        // make it same with the original data

        for (uint i = d; i < NUM_ROWS; i++)
        {
            for (uint j = d; j < NUM_COLUMNS; j++)
            {
                minor.Components[i][j] = this->Components[i][j];
            }
        }

        return minor;
    }

    // take c-th column of matrix, put in vector
    Vec<T, NUM_COLUMNS> GetColum(uint c)
    {
        Vec<T, NUM_COLUMNS> vec;
        for (uint i = 0; i < this->NUM_ROWS; i++)
            vec[i] = this->Components[i][c];
        return vec;
    }
};

// m = I - 2* v v^T
template <typename T, uint Size>
Matrix<T, Size, Size> IMINUSVVT(Vec<T, Size> &v)
{
    // this member become the static const member?
    Matrix<T, Size, Size> Mat;
    for (uint i = 0; i < Size; i++)
        for (uint j = 0; j < Size; j++)
            Mat[i][j] = -2 * v[i] * v[j];
    // diagnal elements add 1
    for (uint i = 0; i < Size; i++)
        Mat[i][i] = Mat[i][i] + 1;

    return Mat;
}

template <typename T,
          uint NumRow,
          uint NumCol,
          uint NumInternal>
Matrix<T, NumRow, NumCol> MatrixMultiply(
    Matrix<T, NumRow, NumInternal> &leftFactor,
    Matrix<T, NumInternal, NumCol> &rightFactor)
{
    Matrix<T, NumRow, NumCol> result;

    for (uint i = 0; i < NumRow; i++)
    {
        for (uint j = 0; j < NumCol; j++)
        {
            // init specific position of z as zero
            result[i][j] = 0;
            for (uint k = 0; k < NumInternal; k++)
            {
                result[i][j] = result[i][j] + leftFactor[i][k] * rightFactor[k][j];
            }
        }
    }
    return result;
}

// matrix multiply vector plus vector
// A*U+M
template <typename T,
          uint NumRow,
          uint NumCol>
Vec<T, NumRow> MMVPV(
    Matrix<T, NumRow, NumCol> &A,
    Vec<T, NumCol> &U,
    Vec<T, NumRow> &M)
{
    Vec<T, NumRow> AUM;
    AUM.InitZero();
    // Be careful about it, init the vector to zero
    for (int i = 0; i < NumRow; i++)
    {
        // adding vector into the matrix
        for (int j = 0; j < NumCol; j++)
        {
            AUM[i] = AUM[i] + (A[i][j] * U[j]);
        }
        AUM[i] = AUM[i] + M[i];
    }

    return AUM;
}

// results = a*x + b
template <typename T, uint Size>
Vec<T, Size> AXPY(double a, Vec<T, Size> &x, Vec<T, Size> &b)
{
    Vec<T, Size> result;
    for (uint i = 0; i < Size; i++)
    {
        result[i] = a * x[i] + b[i];
    }
    return result;
}

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
            Q = MatrixMultiply(Q, MTemp);
        }

        // update R
        // puts("before updating R");
        // puts("MList[k]");
        // MTemp.Show();

        // puts("R");
        // R.Show();

        // R = MatrixMultiply(MList[k], R);
        R = MatrixMultiply(MTemp, R);

        // puts("update R is");
        // R.Show();
    }
}

template <typename T, uint Size>
void SymmEigenvalues(const Matrix<T, Size, Size> &A, double tol, int maxIter, Vec<T, Size> &eigenValues)
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

        Ak = MatrixMultiply(R, Q);
        QQ = MatrixMultiply(QQ, Q);

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

#endif