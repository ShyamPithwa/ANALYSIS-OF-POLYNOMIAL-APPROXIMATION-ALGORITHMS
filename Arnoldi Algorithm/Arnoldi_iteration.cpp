#include "Arnoldi_iteration.h"
#include "NLA_FUNCTIONS.h"
#include "factorizations.h"
#include <limits>
#include <iostream>
#include <chrono>
#include <fstream>

Arnoldi_iteration::Arnoldi_iteration(Matrix &A, Matrix &q, int k)
{

    int n = A.getRows();

    H = Matrix(k + 1, k, 0.0);
    Q = Matrix(n, k + 1, 0.0);

    // q_0 = q / ||q||
    Matrix q0 = q * (1.0 / tnorm(q));
    set_column(Q, 0, q0);

    for (int i = 0; i < k; ++i)
    {

        Matrix qi = access_colunm(Q, i);

        // v = A q_i
        Matrix v = A * qi;

        // Orthogonalization
        for (int j = 0; j <= i; ++j)
        {

            Matrix qj = access_colunm(Q, j);

            H(j, i) = inner_product(qj, v);

            v = v - qj * H(j, i);
        }

        // h_{i+1,i} = ||v||
        H(i + 1, i) = tnorm(v);

        // Breakdown check
        if (H(i + 1, i) < 1e-12)
        {
            break;
        }

        // q_{i+1} = v / h_{i+1,i}
        Matrix qi_plus_1 = v * (1.0 / H(i + 1, i));

        set_column(Q, i + 1, qi_plus_1);
    }
}

Arnoldi_Krylov_Polynomial::Arnoldi_Krylov_Polynomial(const Matrix &data)
{
    // Determine data orientation and dimensions
    // Expected arrangement:
    // row 0: x1, x2, ..., xn
    // row 1: y1, y2, ..., yn
    int n = 0;
    vector<double> x_vec;
    vector<double> y_vec;

    if (data.getRows() == 2 && data.getCols() >= 1)
    {
        n = data.getCols();
        x_vec.resize(n);
        y_vec.resize(n);
        for (int i = 0; i < n; ++i)
        {
            x_vec[i] = data(0, i);
            y_vec[i] = data(1, i);
        }
    }
    else if (data.getCols() == 2 && data.getRows() >= 1)
    {
        n = data.getRows();
        x_vec.resize(n);
        y_vec.resize(n);
        for (int i = 0; i < n; ++i)
        {
            x_vec[i] = data(i, 0);
            y_vec[i] = data(i, 1);
        }
    }
    else
    {
        n = data.getCols();
        x_vec.resize(n);
        y_vec.resize(n);
        for (int i = 0; i < n; ++i)
        {
            x_vec[i] = data(0, i);
            if (data.getRows() > 1)
            {
                y_vec[i] = data(1, i);
            }
        }
    }

    if (n <= 0)
    {
        effectiveDegree = -1;
        finalPolynomial = Polynomial();
        return;
    }

    // STEP 1 — Construct the Arnoldi problem
    // X = diag(x_1, ..., x_n)
    X = Matrix(n, n, 0.0);
    for (int i = 0; i < n; ++i)
    {
        X(i, i) = x_vec[i];
    }

    // q0 = 1 / sqrt(n) * [1, 1, ..., 1]^T
    Matrix q0(n, 1, 1.0 / std::sqrt(static_cast<double>(n)));

    // Perform Arnoldi iteration ONCE using existing Arnoldi_iteration class
    int max_k = (n > 1) ? std::min(n - 1, 50) : 0;
    Arnoldi_iteration arnoldi(X, q0, max_k);
    this->Q = arnoldi.Q;
    this->H = arnoldi.H;

    // STEP 2 — Generate candidate polynomials
    // p_0(x) = 1
    candidatePolynomials.clear();
    candidatePolynomials.push_back(Polynomial(vector<double>{1.0}));
    vector<double> candidateGenerationTimes(1, 0.0);

    for (int k = 0; k < max_k; ++k)
    {
        double h_next = H(k + 1, k);
        if (h_next < 1e-12)
        {
            // Breakdown occurred; terminate polynomial recurrence
            break;
        }

        // p_{k+1}(x) = (x * p_k(x) - sum_{j=0}^k h_{j,k} * p_j(x)) / h_{k+1,k}
        const auto candidateStart = std::chrono::steady_clock::now();
        Polynomial p_next = candidatePolynomials[k];
        p_next.increasePower();

        for (int j = 0; j <= k; ++j)
        {
            double h_jk = H(j, k);
            p_next.addPolynomial(candidatePolynomials[j], -h_jk);
        }

        p_next.scale(1.0 / h_next);
        candidatePolynomials.push_back(p_next);
        const auto candidateEnd = std::chrono::steady_clock::now();
        candidateGenerationTimes.push_back(
            std::chrono::duration<double, std::milli>(
                candidateEnd - candidateStart).count());
    }

    int numCandidates = static_cast<int>(candidatePolynomials.size());

    // STEP 3 — Calculate approximation coefficients
    // c_j = q_j^T * y
    coefficients.clear();
    coefficients.reserve(numCandidates);
    for (int j = 0; j < numCandidates; ++j)
    {
        vector<double> q_j = access_colunm_vector(Q, j);
        double c_j = Dot_product(q_j, y_vec);
        coefficients.push_back(c_j);
    }

// ==================================================
// STEP 4 + STEP 5
// Evaluate candidates, print every candidate,
// select low degree using residual improvement.
// ==================================================

residuals.clear();

double norm_y = tnorm(y_vec);

vector<double> y_hat(n, 0.0);

// --------------------------------------------------
// sqrt(n)
// --------------------------------------------------

const double sqrt_n =
    std::sqrt(static_cast<double>(n));


// --------------------------------------------------
// Tolerances
// --------------------------------------------------

const double RESIDUAL_THRESHOLD = 1e-5;



// --------------------------------------------------
// Previous residual
// --------------------------------------------------

double previousResidual =
    std::numeric_limits<double>::infinity();


// --------------------------------------------------
// Final selected polynomial
// --------------------------------------------------

int selectedDegree = -1;
double totalCandidateTimeMs = 0.0;
int iterationCount = 0;

Polynomial selectedPolynomial;

double selectedResidual =
    std::numeric_limits<double>::infinity();


// --------------------------------------------------
// Output
// --------------------------------------------------

std::cout << "\n========================================\n";
std::cout << "       ARNOLDI CANDIDATE POLYNOMIALS\n";
std::cout << "========================================\n";


for (int k = 0; k < numCandidates; ++k)
{
    const auto candidateStart = std::chrono::steady_clock::now();
    // ==================================================
    // CURRENT ARNOLDI APPROXIMATION
    //
    // y_hat_k = sum_{j=0}^k c_j q_j
    // ==================================================

    vector<double> q_k =
        access_colunm_vector(Q, k);

    vector<double> c_k_qk =
        VS_Multiplication(q_k, coefficients[k]);

    y_hat =
        Vec_Add(y_hat, c_k_qk);


    // ==================================================
    // CALCULATE RESIDUAL
    // ==================================================

    vector<double> res_vec =
        Vec_Sub(y_vec, y_hat);

    double res_norm =
        tnorm(res_vec);

    double currentResidual =
        (norm_y > 0.0)
        ? (res_norm / norm_y)
        : res_norm;

    residuals.push_back(currentResidual);


    // ==================================================
    // CONSTRUCT CURRENT POLYNOMIAL
    //
    // q_j = p_j / sqrt(n)
    //
    // Therefore:
    //
    // c_j q_j
    // =
    // (c_j / sqrt(n)) p_j
    // ==================================================

    Polynomial currentPolynomial;

    for (int j = 0; j <= k; ++j)
    {
        double polynomialCoefficient =
            coefficients[j] / sqrt_n;

        currentPolynomial.addPolynomial(
            candidatePolynomials[j],
            polynomialCoefficient
        );
    }


    // ==================================================
    // PRINT CURRENT CANDIDATE
    // ==================================================

    std::cout
        << "\n----------------------------------------\n";

    std::cout
        << "Degree: "
        << k
        << "\n";

    std::cout
        << "Polynomial:\n";

    vector<double> currentCoefficients =
        currentPolynomial.getCoefficients();

    printPolynomial(currentCoefficients);


    std::cout
        << std::setprecision(12);

    std::cout
        << "Residual: "
        << currentResidual
        << "\n";
    const auto candidateEnd = std::chrono::steady_clock::now();
    const double candidateTimeMs = candidateGenerationTimes[k]
        + std::chrono::duration<double, std::milli>(
              candidateEnd - candidateStart).count();
    totalCandidateTimeMs += candidateTimeMs;
    ++iterationCount;
    std::cout
        << "Candidate calculation time: "
        << candidateTimeMs
        << " ms\n";


    // ==================================================
    // CALCULATE CHANGE IN RESIDUAL
    // ==================================================

    double residualImprovement =
        std::numeric_limits<double>::infinity();

    if (k > 0)
    {
        residualImprovement =
            std::abs(previousResidual -
                     currentResidual);

        std::cout
            << "Residual improvement: "
            << residualImprovement
            << "\n";
    }


    // ==================================================
    // STOPPING CRITERIA
    // ==================================================

    bool residualThresholdReached =
        (currentResidual <=
         RESIDUAL_THRESHOLD);


    // ==================================================
    // IF STOPPING CRITERION REACHED
    // ==================================================

    if (residualThresholdReached)
    {
        selectedDegree = k;

        selectedPolynomial =
            currentPolynomial;

        selectedResidual =
            currentResidual;


        std::cout
            << "\n>>> TERMINATION CRITERION REACHED <<<\n";


        if (residualThresholdReached)
        {
            std::cout
                << "Reason: Residual below tolerance.\n";
        }
        break;
    }


    // ==================================================
    // UPDATE PREVIOUS RESIDUAL
    // ==================================================

    previousResidual =
        currentResidual;
}


// ==================================================
// FALLBACK
// ==================================================
//
// If no stopping criterion was reached,
// use the final generated candidate.
// ==================================================

if (selectedDegree == -1 &&
    numCandidates > 0)
{
    selectedDegree =
        numCandidates - 1;

    selectedPolynomial =
        candidatePolynomials[selectedDegree];

    selectedResidual =
        residuals[selectedDegree];

    // Reconstruct the polynomial correctly
    // using c_j / sqrt(n)

    selectedPolynomial =
        Polynomial();

    for (int j = 0;
         j <= selectedDegree;
         ++j)
    {
        double polynomialCoefficient =
            coefficients[j] / sqrt_n;

        selectedPolynomial.addPolynomial(
            candidatePolynomials[j],
            polynomialCoefficient
        );
    }
}


// ==================================================
// FINAL SUMMARY
// ==================================================

effectiveDegree =
    selectedDegree;

const double averageTimePerIterationMs = iterationCount > 0
    ? totalCandidateTimeMs / iterationCount : 0.0;

finalPolynomial =
    selectedPolynomial;


std::cout
    << "\n\n========================================\n";

std::cout
    << "          FINAL ARNOLDI SUMMARY\n";

std::cout
    << "========================================\n";

std::cout
    << "Selected degree   : "
    << effectiveDegree
    << "\n";

std::cout
    << std::setprecision(12);

std::cout
    << "Residual           : "
    << selectedResidual
    << "\n";
std::cout
    << "Average time per iteration: "
    << averageTimePerIterationMs
    << " ms (" << iterationCount << " degree iterations)\n";

std::cout
    << "\nSelected polynomial:\n";

vector<double> finalCoefficients =
    finalPolynomial.getCoefficients();

printPolynomial(finalCoefficients);

std::ofstream resultFile("best_polynomial_results.txt", std::ios::app);
if (!resultFile)
    throw std::runtime_error("Unable to open best_polynomial_results.txt.");
resultFile << std::setprecision(17)
           << "Algorithm: Arnoldi-based\n"
           << "Best degree: " << effectiveDegree << "\n"
           << "Average time per iteration (ms): " << averageTimePerIterationMs << "\n"
           << "Coefficients (ascending powers):";
for (double coefficient : finalCoefficients)
    resultFile << " " << coefficient;
resultFile << "\n\n";

std::cout
    << "\n========================================\n";
}
