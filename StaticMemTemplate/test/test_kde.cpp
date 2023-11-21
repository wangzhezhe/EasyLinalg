#include <include/basic.h>
#include <include/kde.h>
#include <assert.h>

using namespace EASYLINALG;

double testArray[6] = {23.40953488, 20.06452322, 14.63870771, 21.61275585, 29.82232875, 11.23601927};
double testArray2[80] = {0.38051, 5.37015, 5.39432, 0.42799, 5.37231, 5.38274, 0.35108, 0.40848, 5.35569, 5.38413, 0.37381, 0.3363, 5.39313, 5.391, 0.36787, 0.35847, 0.36968, 0.39352, 0.32975, 5.35745, 0.36587, 0.3432, 5.37259, 0.35032, 5.37496, 5.36477, 5.33913, 0.36177, 5.33206, 5.34491, 0.38014, 5.41331, 0.33376, 0.36947, 5.37308, 5.36478, 5.37552, 0.35477, 0.37243, 0.35559, 0.37561, 0.37756, 5.36344, 5.37984, 0.38726, 5.32108, 5.35691, 5.3622, 0.3725, 0.40043, 5.3801, 0.3643, 5.37245, 0.38601, 0.36525, 5.39771, 0.36984, 5.35069, 5.37312, 5.32356, 0.35683, 5.38209, 5.37106, 5.40225, 5.35648, 0.36374, 0.38218, 0.35696, 0.37202, 0.39191, 5.3624, 0.41786, 0.38808, 0.36177, 5.34164, 0.36941, 0.34972, 5.32496, 5.38209, 0.36764};

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

void testKDECDF1D_2()
{
    std::cout << "test testKDECDF1D_2" << std::endl;
    Vec<double, 80> v(0);
    for (int i = 0; i < 80; i++)
    {
        v[i] = testArray2[i];
    }

    for (int i = 0; i < 80; i++)
    {
        //this input value [0, n] here should match with the actual range in the test data set
        double cdfValue = KDE1D(v, (5*i)/80.0);
        std::cout << cdfValue << ",";
    }
    std::cout << std::endl;

}

int main()
{
    testBandwidthScott();
    testKernel();
    testKDE();
    testKDECDF1D();
    testKDECDF1D_2();
}