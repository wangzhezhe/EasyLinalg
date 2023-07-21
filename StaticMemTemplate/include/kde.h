

#ifndef EASY_LIALG_EIGEN
#define EASY_LIALG_EIGEN

#include "./basic.h"
#include "./blas.h"
namespace EASYLINALG
{
    template <typename T, uint Size>
    LIAG_FUNC_MACRO float BandWidthScott(const Vec<T, Size> &input)
    {
        bool sampleStdev = true;
        double stdev = input.GetStdev(sampleStdev);
        int len = input.NUM_COMPONENTS;
        return 3.49 * stdev * pow(len, (-0.333));
    }

    LIAG_FUNC_MACRO float GaussianKernel(float input)
    {
        return (1 / (sqrt(2 * M_PI))) * exp(-0.5 * input * input);
    }

    // input: data vector,dedicated point, bandwidth(default parameter)
    // output kde results at the dedicated point
    template <typename T, uint Size>
    LIAG_FUNC_MACRO double KDE1D(const Vec<T, Size> &inputSample, float dedicatedPoint, float bandWidth = 0.0)
    {
        // choose the bandWidth according to input if it uses default value
        float h = BandWidthScott(inputSample);

        double kde = 0;
        // go through each input value, compuet kde value
        for (int i = 0; i < inputSample.NUM_COMPONENTS; i++)
        {
            double kernelInput = (dedicatedPoint - inputSample[i]) / h;
            kde = kde + GaussianKernel(kernelInput);
        }

        kde = kde / (inputSample.NUM_COMPONENTS * 1.0 * h);

        return kde;
    }

}

#endif