#include "factorizations.h"
#include "NLA_FUNCTIONS.h"
#include <iostream>
#include <chrono>
#include <fstream>
#include<cfloat>
#include <algorithm>

using namespace std;

Householder_QR::Householder_QR(const Matrix &A)
    : R(A), Q(identity(A.getRows()))
{

    int n = A.getCols();
    int m = A.getRows();

    for (int i = 0; i < std::min(m - 1, n); i++)
    {
        Matrix H = identity(m);
        Matrix A_sub = Sub_Matrix(R, i, m - 1, i, n - 1);

        vector<double> x = access_colunm_vector(A_sub, 0);
        vector<double> e(x.size(), 0.0);
        e[0] = 1.0;

        double y = tnorm(x);
        double z = (x[0] >= 0) ? 1.0 : -1.0;

        vector<double> v = VS_Multiplication(e, y * z);
        v = Vec_Add(v, x);

        double v1 = tnorm(v);
        if (v1 == 0.0)
            continue;

        v = VS_Multiplication(v, 1.0 / v1);
        Matrix H_small = identity(v.size());
        H_small = H_small - (OuterProduct(v, v) * 2.0);

        Matrix temp = VM_multiplication(v, A_sub);
        Matrix update = OuterProduct(v, temp) * 2.0;
        InsertSubMatrix(H, H_small, i, i);
        Q = Q * H;
        A_sub = A_sub - update;
        InsertSubMatrix(R, A_sub, i, i);
    }
}

vector<double> back_substitution(const Matrix &R, const Matrix &b)
{
    int n = R.getCols();
    vector<double> x(n, 0.0);

    for (int i = n - 1; i >= 0; i--)
    {
        double sum = 0.0;

        for (int j = i + 1; j < n; j++)
        {
            sum += R(i, j) * x[j];
        }

        if (fabs(R(i, i)) < 1e-12)
        {
            throw runtime_error("Back substitution failed: zero diagonal element.");
        }

        x[i] = (b(i, 0) - sum) / R(i, i);
    }

    return x;
}

double maximumApproximationError(const vector<double>& x,
                                 const vector<double>& y,
                                 const vector<double>& coefficients)
{
    if (x.size() != y.size())
    {
        throw invalid_argument("x and y must contain the same number of points.");
    }

    double maximumError = 0.0;
    for (size_t i = 0; i < x.size(); ++i)
    {
        double approximation = 0.0;
        for (auto coefficient = coefficients.rbegin();
             coefficient != coefficients.rend(); ++coefficient)
        {
            approximation = approximation * x[i] + *coefficient;
        }

        maximumError = max(maximumError, fabs(approximation - y[i]));
    }

    return maximumError;
}

vector<double> Polynomial_Approximator(const vector<double> &x, const vector<double> &y)
{
    const double tolerance = 1e-5;
    double bestResidual = numeric_limits<double>::max();
    double totalCandidateTimeMs = 0.0;
    int iterationCount = 0;
    int bestDegree = 2;
    vector<double> bestCoeff;
    Matrix u = vectorToMatrix(y);
    vector<double> c;
    const double dataNorm = tnorm(u);
    const double residualDenominator = dataNorm == 0.0 ? 1.0 : dataNorm;
    int degree = 2;
    while (degree < x.size() && degree <= 50)
    {
        const auto candidateStart = chrono::steady_clock::now();
        Matrix A = vandermonde(x, degree);
        Householder_QR qr(A);
        int k = degree + 1;
        Matrix R1 = Sub_Matrix(qr.R, 0, k - 1, 0, k - 1);
        Matrix Qt = qr.Q.transpose();
        Matrix b = Qt * u;
        Matrix rhs = Sub_Matrix(b, 0, k - 1, 0, 0);
        c = back_substitution(R1, rhs);
        double r = 0;
        Matrix z = MV_multiplication(A, c) - u;
        r = tnorm(z) / residualDenominator;
        cout << "Degree = " << degree
             << "  Residual Error = " << setprecision(9) << r << endl;
        printPolynomial(c);
        const auto candidateEnd = chrono::steady_clock::now();
        const double candidateTimeMs = chrono::duration<double, milli>(candidateEnd - candidateStart).count();
        totalCandidateTimeMs += candidateTimeMs;
        ++iterationCount;
        cout << "Candidate calculation time: "
             << candidateTimeMs
             << " ms" << endl;
        cout << endl;
        if (r < bestResidual)
        {
            bestResidual = r;
            bestDegree = degree;
            bestCoeff = c;
        }

        if (r <= tolerance)
        {
            break;
        }

        degree++;
    }
    const double averageTimePerIterationMs = iterationCount > 0
        ? totalCandidateTimeMs / iterationCount : 0.0;
    cout << "Average time per iteration: " << averageTimePerIterationMs
         << " ms (" << iterationCount << " degree iterations)" << endl;
    cout << endl
         << endl;
    cout << "Best Degree    : " << bestDegree << endl;
    cout << "Residual Error : " << setprecision(9) << bestResidual << endl;
    cout << "Maximum Error  : " << setprecision(9)
         << maximumApproximationError(x, y, bestCoeff) << endl;

    ofstream resultFile("best_polynomial_results.txt", ios::app);
    if (!resultFile)
        throw runtime_error("Unable to open best_polynomial_results.txt.");
    resultFile << setprecision(17)
               << "Algorithm: QR-based\n"
               << "Best degree: " << bestDegree << "\n"
               << "Average time per iteration (ms): " << averageTimePerIterationMs << "\n"
               << "Coefficients (ascending powers):";
    for (double coefficient : bestCoeff)
        resultFile << " " << coefficient;
    resultFile << "\n\n";

    return bestCoeff;
}

void printPolynomial(const vector<double> &coeff)
{
    cout << "P(x) = ";

    bool first = true;

    for (int i = 0; i < coeff.size(); i++)
    {
        if (fabs(coeff[i]) < 1e-12)
            continue;

        if (!first)
        {
            if (coeff[i] >= 0)
                cout << " + ";
            else
                cout << " - ";
        }
        else
        {
            if (coeff[i] < 0)
                cout << "-";
            first = false;
        }

        cout << fixed << setprecision(6) << fabs(coeff[i]);

        if (i >= 1)
            cout << "x";

        if (i >= 2)
            cout << "^" << i;
    }

    cout << endl;
}

void readDataPoints(const char *fileName, vector<double> &x, vector<double> &y)
{
    ifstream file(fileName);

    if (!file)
    {
        throw runtime_error("Unable to open file.");
    }

    double xi, yi;

    while (file >> xi >> yi)
    {
        x.push_back(xi);
        y.push_back(yi);
    }

    file.close();

    if (x.empty())
    {
        throw runtime_error("No data found in file.");
    }
}
