#include <include/basic.h>
#include <include/kde.h>
#include <assert.h>

using namespace EASYLINALG;

double testArray[6] = {23.40953488, 20.06452322, 14.63870771, 21.61275585, 29.82232875, 11.23601927};

void testBandwidthScott()
{
    Vec<int, 10> v;
    for (int i = 0; i < 10; i++)
    {
        v[i] = i + 1;
    }

    double bwScott = BandWidthScott(v);
    assert(fabs(bwScott - 4.90830) < 0.00001);
}

void testKernel()
{
    double x0 = 0.0;
    double g0 = GaussianKernel(x0);
    double x1 = 0.5;
    double g1 = GaussianKernel(x1);
    double x2 = 1.0;
    double g2 = GaussianKernel(x2);

    assert(fabs(g0 - 0.398942) < 0.00001);
    assert(fabs(g1 - 0.352065) < 0.00001);
    assert(fabs(g2 - 0.241970) < 0.00001);
}

void testKDE()
{

    Vec<double, 6> v(0);
    for (int i = 0; i < 6; i++)
    {
        v[i] = testArray[i];
    }

    double x0 = 10.0;
    double g0 = KDE1D(v, x0);
    double x1 = 15.0;
    double g1 = KDE1D(v, x1);
    double x2 = 20.0;
    double g2 = KDE1D(v, x2);

    assert(fabs(g0 - 0.02197104) < 0.00001);
    assert(fabs(g1 - 0.02659957) < 0.00001);
    assert(fabs(g2 - 0.02839158) < 0.00001);
}

void testKDECDF1D()
{
    Vec<double, 6> v(0);
    for (int i = 0; i < 6; i++)
    {
        v[i] = testArray[i];
    }

    double cdfValue = KDECDF1D(v, 30);

    std::cout << "cdfValue " << cdfValue << std::endl;
}

int main()
{
    testBandwidthScott();
    testKernel();
    testKDE();
    testKDECDF1D();
}