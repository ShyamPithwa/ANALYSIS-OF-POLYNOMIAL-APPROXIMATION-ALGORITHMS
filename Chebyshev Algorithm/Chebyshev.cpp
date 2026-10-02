#include "../Headers/factorizations.h"
#include "../Headers/NLA_FUNCTIONS.h"
#include "../Headers/Arnoldi_iteration.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace std;

vector<double> Orthogonal_Polynomial_Approximator(
    const vector<double>& x, const vector<double>& y)
{
    if (x.empty() || x.size() != y.size())
        throw invalid_argument("x and y must have the same non-zero size.");

    const double tolerance = 1e-5;
    const double center = (*min_element(x.begin(), x.end()) +
                           *max_element(x.begin(), x.end())) / 2.0;
    double scale = (*max_element(x.begin(), x.end()) -
                    *min_element(x.begin(), x.end())) / 2.0;
    if (scale == 0.0)
        scale = 1.0;

    Matrix u = vectorToMatrix(y);
    const double dataNorm = tnorm(u);
    const double denominator = dataNorm == 0.0 ? 1.0 : dataNorm;
    double bestResidual = numeric_limits<double>::max();
    double totalCandidateTimeMs = 0.0;
    int iterationCount = 0;
    Polynomial bestPolynomial;
    int bestDegree = 0;

    for (int degree = 0;
         degree < static_cast<int>(x.size()) && degree <= 50;
         ++degree)
    {
        const auto candidateStart = chrono::steady_clock::now();
        const int k = degree + 1;
        Matrix A(static_cast<unsigned>(x.size()), static_cast<unsigned>(k), 0.0);

        // Build the Chebyshev family T_0(t)=1, T_1(t)=t,
        // T_(j+1)(t)=2tT_j(t)-T_(j-1)(t), on the scaled interval [-1,1].
        for (unsigned i = 0; i < x.size(); ++i)
        {
            const double t = (x[i] - center) / scale;
            A(i, 0) = 1.0;
            if (degree >= 1)
                A(i, 1) = t;
            for (int j = 1; j < degree; ++j)
                A(i, j + 1) = 2.0 * t * A(i, j) - A(i, j - 1);
        }

        Householder_QR qr(A);
        Matrix R1 = Sub_Matrix(qr.R, 0, k - 1, 0, k - 1);
        Matrix Qt = qr.Q.transpose();
        Matrix b = Qt * u;
        Matrix rhs = Sub_Matrix(b, 0, k - 1, 0, 0);
        vector<double> coefficients = back_substitution(R1, rhs);

        Matrix residual = MV_multiplication(A, coefficients) - u;
        const double currentResidual = tnorm(residual) / denominator;

        // Convert the Chebyshev series in t back into a single power-basis
        // polynomial in x, combining terms through Polynomial::addPolynomial.
        Polynomial approximatingPolynomial(vector<double>{0.0});
        Polynomial previousBasis(vector<double>{1.0});
        if (degree >= 0)
            approximatingPolynomial.addPolynomial(previousBasis, coefficients[0]);

        if (degree >= 1)
        {
            Polynomial currentBasis(vector<double>{-center / scale, 1.0 / scale});
            approximatingPolynomial.addPolynomial(currentBasis, coefficients[1]);
            for (int j = 1; j < degree; ++j)
            {
                // Multiply currentBasis by t=(x-center)/scale, then apply
                // T_(j+1)(t)=2tT_j(t)-T_(j-1)(t).
                Polynomial nextBasis = currentBasis;
                nextBasis.increasePower();
                nextBasis.scale(1.0 / scale);
                nextBasis.addPolynomial(currentBasis, -center / scale);
                nextBasis.scale(2.0);
                nextBasis.addPolynomial(previousBasis, -1.0);

                approximatingPolynomial.addPolynomial(nextBasis, coefficients[j + 1]);
                previousBasis = currentBasis;
                currentBasis = nextBasis;
            }
        }

        cout << "Degree = " << degree
             << "  Residual Error = " << setprecision(9)
             << currentResidual << endl;
        vector<double> displayCoefficients = approximatingPolynomial.getCoefficients();
        for (double& coefficient : displayCoefficients)
            if (fabs(coefficient) < 5e-7)
                coefficient = 0.0;
        printPolynomial(displayCoefficients);
        const auto candidateEnd = chrono::steady_clock::now();
        const double candidateTimeMs = chrono::duration<double, milli>(candidateEnd - candidateStart).count();
        totalCandidateTimeMs += candidateTimeMs;
        ++iterationCount;
        cout << "Candidate calculation time: "
             << candidateTimeMs
             << " ms" << endl;
        cout << endl;

        if (currentResidual < bestResidual)
        {
            bestResidual = currentResidual;
            bestPolynomial = approximatingPolynomial;
            bestDegree = degree;
        }
        if (currentResidual <= tolerance)
            break;
    }

    const double averageTimePerIterationMs = iterationCount > 0
        ? totalCandidateTimeMs / iterationCount : 0.0;
    cout << "Average time per iteration: " << averageTimePerIterationMs
         << " ms (" << iterationCount << " degree iterations)" << endl;

    cout << endl << "Best Degree    : " << bestDegree << endl
         << "Residual Error : " << setprecision(9) << bestResidual << endl
         << "Best Fit Polynomial" << endl;
    vector<double> powerCoefficients = bestPolynomial.getCoefficients();
    vector<double> displayCoefficients = powerCoefficients;
    for (double& coefficient : displayCoefficients)
        if (fabs(coefficient) < 5e-7)
            coefficient = 0.0;
    printPolynomial(displayCoefficients);
    ofstream resultFile("best_polynomial_results.txt", ios::app);
    if (!resultFile)
        throw runtime_error("Unable to open best_polynomial_results.txt.");
    resultFile << setprecision(17)
               << "Algorithm: Chebyshev-based\n"
               << "Best degree: " << bestDegree << "\n"
               << "Average time per iteration (ms): " << averageTimePerIterationMs << "\n"
               << "Coefficients (ascending powers):";
    for (double coefficient : powerCoefficients)
        resultFile << " " << coefficient;
    resultFile << "\n\n";
    return powerCoefficients;
}

