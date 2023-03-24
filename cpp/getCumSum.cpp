#include <Eigen/Dense>

namespace ei = Eigen;

ei::ArrayXXi getCumSum(const ei::ArrayXXi &signal)
{
    int nSamples = signal.rows();
    int nDims = signal.cols();
    ei::MatrixXi out(nSamples + 1, nDims);
    out.row(0).setZero();

    int accumulator;

    for (size_t kDim = 0; kDim < nDims; kDim++)
    {
        accumulator = 0;        
        for (size_t kSample = 0; kSample < nSamples; kSample++)
        {
            accumulator += signal(kSample, kDim);
            out(kSample + 1, kDim) = accumulator;
        }
    }
    return out;
}