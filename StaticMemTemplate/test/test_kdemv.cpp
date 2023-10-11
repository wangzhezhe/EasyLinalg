#include <include/basic.h>
#include <include/kdemv.h>
#include <assert.h>
#include "./tester/kdepp.h"

using namespace EASYLINALG;
using namespace kdepp;

double TestKDESource()
{
    std::cout << "---TestKDESource" << std::endl;
    // test code comes from
    // https://github.com/duncanmcn/kdepp/tree/master
    // three samples
    Eigen::VectorXf p1(2);
    p1 << 0, 2;
    Eigen::VectorXf p2(2);
    p2 << 0.15, 2.3;
    Eigen::VectorXf p3(2);
    p3 << -0.1, 2.5;
    // std::array<double, 2> p2 = {};
    // std::array<double, 2> p3 = {};

    std::vector<Eigen::VectorXf> data = {p1, p2, p3};

    Kdemv kernel(data);

    // std::array<double, 2> test_point = {0.1, 2.5};
    Eigen::VectorXf test_point(2);
    test_point << 0.1, 2.5;
    double result = kernel.eval(test_point);
    std::cout << "TestKDESource " << result << std::endl;
    return result;
}

void TestKDEMV()
{
    std::cout << "---TestKDEMV" << std::endl;
    Vec<double, 3> attr1;
    attr1[0] = 0;
    attr1[1] = 0.15;
    attr1[2] = -0.1;

    Vec<double, 3> attr2;
    attr2[0] = 2.0;
    attr2[1] = 2.3;
    attr2[2] = 2.5;

    Vec<Vec<double, 3>, 2> inputData;
    inputData[0] = attr1;
    inputData[1] = attr2;

    KDEMV<double, 2, 3> kdemv(inputData);

    Vec<double, 2> sample;
    sample[0] = 0.1;
    sample[1] = 2.5;
    double result = kdemv.evaluatePDV(sample);
    // TODO, continue debug results here
    std::cout << "TestKDEMV " << result << std::endl;

    double srcRes = TestKDESource();
    assert(fabs(result - srcRes) < 0.000001);
}

int main()
{
    TestKDEMV();
}