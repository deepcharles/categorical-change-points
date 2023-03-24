#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>
#include <Eigen/Dense>

namespace ei = Eigen;

ei::ArrayXXi getCumSum(const ei::ArrayXXi &signal);

ei::ArrayXi getAFPCP(const ei::ArrayXXi &signal, double penalty)
{
    int nSamples = signal.rows();
    int nSymbols = signal.cols();
    ei::ArrayXXi cumSignal = getCumSum(signal);

    // instantiate path vector
    ei::ArrayXi out(nSamples + 1); // path vector
    out(0) = -1;

    // Define the xlogx function
    auto xlogx_func = [](double x)
    { return x == 0 ? 0.0 : x * std::log(x); };
    auto xlogx_ref = std::ref(xlogx_func);

    // temporary variables
    std::vector<size_t> admissibleSetOfCP = {0};         // set of admissible change-point indexes
    std::vector<size_t>::iterator iterAdmissibleSetOfCP; // iterator over the vector of the admissible change-point indexes
    std::vector<double> tempSocVec;                      // current soc (for pruning)
    ei::ArrayXd socVec(nSamples + 1);                    // vector of sums of costs (soc)
    socVec(0) = 0;
    size_t bestCP, kIndex;
    double currentMinSoc, currentSoc, currentCost, segmentLength;
    ei::ArrayXd sumSymbols(nSymbols);

    ei::ArrayXd lastSegment, compareSegment, probaVec;
    size_t compareCP, lastCP;
    double lambda, probaSum, lagrangeVal;
    bool isPruned;
    double infDouble = std::numeric_limits<double>::infinity();

    // helper function to compute sum_i x_i*log(y_i)
    auto getXDotLogY = [](ei::ArrayXd x, ei::ArrayXd y)
    {
        double out = 0.0;
        for (size_t kDim = 0; kDim < x.rows(); kDim++)
        {
            if (x(kDim) > 0)
            {
                out += x(kDim) * log(y(kDim));
            }
        }
        return out;
    };

    // change-point detection for endIndex>1
    for (size_t endIndex = 1; endIndex < nSamples + 1; endIndex++)
    {
        currentMinSoc = infDouble;
        for (size_t lastCP : admissibleSetOfCP)
        {
            // compute the cost on segment lastCP..endIndex
            segmentLength = double(endIndex - lastCP);
            sumSymbols = (cumSignal.row(endIndex) - cumSignal.row(lastCP)).cast<double>();
            currentCost = -segmentLength * (sumSymbols / segmentLength).unaryExpr(xlogx_ref).sum();

            // sum of cost V(0, lastCP) + c(lastCP, endIndex) + beta
            currentSoc = socVec(lastCP) + currentCost + penalty;

            // keep track of SoC for pruning
            tempSocVec.push_back(socVec(lastCP) + currentCost);

            // keep the minimum sum of costs
            if (currentSoc < currentMinSoc)
            {
                currentMinSoc = currentSoc;
                bestCP = lastCP;
            }
        }
        // Record the minimum and argmin
        socVec(endIndex) = currentMinSoc;
        out(endIndex) = bestCP;

        // Pruning
        iterAdmissibleSetOfCP = admissibleSetOfCP.begin();
        compareCP = *iterAdmissibleSetOfCP; // start with the first non-pruned index
        while (iterAdmissibleSetOfCP != admissibleSetOfCP.end())
        {
            isPruned = false;
            lastCP = *iterAdmissibleSetOfCP;

            // compute the cost on segment lastCP..endIndex
            segmentLength = double(endIndex - lastCP);
            sumSymbols = (cumSignal.row(endIndex) - cumSignal.row(lastCP)).cast<double>();
            currentCost = -segmentLength * (sumSymbols / segmentLength).unaryExpr(xlogx_ref).sum();

            if (socVec(lastCP) + currentCost > currentMinSoc) // PELT rule
            {
                iterAdmissibleSetOfCP = admissibleSetOfCP.erase(iterAdmissibleSetOfCP);
                isPruned = true;
            }
            else // not pruned by PELT
            {
                if (compareCP < lastCP)
                {
                    // Find interval for lambda
                    lambda = infDouble;
                    lastSegment = (cumSignal.row(endIndex) - cumSignal.row(lastCP)).cast<double>();
                    compareSegment = (cumSignal.row(lastCP) - cumSignal.row(compareCP)).cast<double>();

                    for (size_t kSymbol = 0; kSymbol < nSymbols; kSymbol++)
                    {
                        if ((lastSegment(kSymbol) == 0) and (compareSegment(kSymbol) > 0))
                        {
                            lambda = infDouble;
                            break;
                        }
                        else if ((lastSegment(kSymbol) > 0) and (compareSegment(kSymbol) > 0))
                        {
                            lambda = std::min(lambda, lastSegment(kSymbol) / compareSegment(kSymbol));
                        }
                    }
                    // Choose a lambda and compute Lagrangian
                    if (lambda != infDouble)
                    {
                        lambda *= 0.75; // go somewhere in the interval [0, lambda_max]
                        probaVec = (lastSegment - lambda * compareSegment);
                        probaSum = probaVec.sum();
                        if (probaSum > 0)
                        {
                            probaVec /= probaSum; // normalize probabilities
                            lagrangeVal = socVec(lastCP) + penalty - getXDotLogY(lastSegment, probaVec);
                            lagrangeVal += lambda * (socVec(lastCP) - socVec(compareCP) + getXDotLogY(compareSegment, probaVec));
                            if (lagrangeVal > currentMinSoc + penalty) // AFP rule
                            {
                                iterAdmissibleSetOfCP = admissibleSetOfCP.erase(iterAdmissibleSetOfCP);
                                isPruned = true;
                            }
                        }
                    }
                }
            }
            if (not isPruned) // not pruned by PELT or AFP
            {
                ++iterAdmissibleSetOfCP;
            }
            compareCP = lastCP;
        }
        admissibleSetOfCP.push_back(endIndex);
    }
    return out;
}