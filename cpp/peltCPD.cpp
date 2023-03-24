#include <vector>
#include <cmath>
#include <functional>
#include <Eigen/Dense>

namespace ei = Eigen;

ei::ArrayXXi getCumSum(const ei::ArrayXXi &signal);

ei::ArrayXi getPeltCP(const ei::ArrayXXi &signal, double penalty)
{
    int nSamples = signal.rows();
    int nSymbols = signal.cols();
    ei::ArrayXXi cumSignal = getCumSum(signal);

    // instantiate path vector
    ei::VectorXi out(nSamples + 1); // path vector
    out(0) = -1;

    // Define the xlogx function
    auto xlogx_func = [](double x)
    { return x == 0 ? 0.0 : x * std::log(x); };
    auto xlogx_ref = std::ref(xlogx_func);
    
    // temporary variables
    std::vector<size_t> admissibleSetOfCP = {0};         // set of admissible change-point indexes
    std::vector<size_t>::iterator iterAdmissibleSetOfCP; // iterator over the vector of the admissible change-point indexes
    std::vector<float> tempSocVec;                       // current soc (for pruning)
    ei::VectorXf socVec(nSamples + 1);                   // vector of sums of costs (soc)
    socVec(0) = 0;
    size_t bestCP, kIndex;
    double currentMinSoc, currentSoc, currentCost, segmentLength;
    ei::ArrayXd sumSymbols(nSymbols);

    // change-point detection for endIndex>1
    for (size_t endIndex = 1; endIndex < nSamples + 1; endIndex++)
    {
        currentMinSoc = std::numeric_limits<float>::infinity();
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
        kIndex = 0;
        while (iterAdmissibleSetOfCP != admissibleSetOfCP.end())
        {
            if (tempSocVec[kIndex] > currentMinSoc) // PELT rule
            {
                iterAdmissibleSetOfCP = admissibleSetOfCP.erase(iterAdmissibleSetOfCP);
            }
            else
            {
                ++iterAdmissibleSetOfCP;
            }
            ++kIndex;
        }
        tempSocVec.clear();
        admissibleSetOfCP.push_back(endIndex);
    }
    return out;
}