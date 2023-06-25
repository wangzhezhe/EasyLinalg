#include <include/basic.h>
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

    // init vector by fixed value
    Vec<double, 10> v_init(1.0);
    std::cout << "v init is" << std::endl;
    v_init.Show();

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
}

void testOperator()
{
    Matrix<int, 3, 3> m1;
    m1.InitEye();
    Matrix<int, 3, 3> m2;
    m2.InitEye();
    // operator +
    m1 = m1 + m2;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (i == j)
            {
                assert(m1[i][j] == 2);
            }
            else
            {
                assert(m1[i][j] == 0);
            }
        }
    }
    m1.Show();

    Matrix<int, 3, 3> m3;
    m3.InitEye();
    //operator -
    m2 = m2 - m3;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            assert(m2[i][j] == 0);
        }
    }
    m2.Show();

    return;
}
int main()
{
    testInit();
    testOperator();
}