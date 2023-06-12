#ifndef EASY_LIALG_BASIC
#define EASY_LIALG_BASIC

#include <iostream>
#include <assert.h>
#include <cmath>
#include <type_traits>

//basic defination of template based matrix and vector

template <typename T, uint Size>
class Vec
{
public:
    static const int NUM_COMPONENTS = Size;
    T Components[NUM_COMPONENTS];

    Vec(){}

    Vec(T init){
        for(uint i=0;i<this->NUM_COMPONENTS;i++){
            this->Components[i]=init;
        }
    }

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

    bool Equal(const Vec &t)
    {
        if (this->NUM_COMPONENTS != t.NUM_COMPONENTS)
        {
            return false;
        }

        for (uint i = 0; i < this->NUM_COMPONENTS; i++)
        {
            if (fabs(this->Components[i] - t.Components[i]) > 0.0001)
            {
                return false;
            }
        }
        return true;
    }

    void Normalize()
    {
        double sum = 0;
        for (uint i = 0; i < this->NUM_COMPONENTS; i++)
        {
            sum += this->Components[i] * this->Components[i];
        }
        assert(fabs(sum - 0.0) > 0.000001);
        for (uint i = 0; i < this->NUM_COMPONENTS; i++)
        {
            this->Components[i] = this->Components[i] / sqrt(sum);
        }
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

    bool IsUpperTriangular(double tol) const
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

    bool IsEqual(const Matrix &t, double tol=0.00001)
    {
        if (this->NUM_ROWS != t.NUM_ROWS || this->NUM_COLUMNS != t.NUM_COLUMNS)
        {
            return false;
        }
        for (uint i = 0; i < this->NUM_ROWS; i++)
        {
            for (uint j = 0; j < this->NUM_COLUMNS; j++)
            {
                if (fabs(this->Components[i][j] - t.Components[i][j]) > tol)
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





#endif