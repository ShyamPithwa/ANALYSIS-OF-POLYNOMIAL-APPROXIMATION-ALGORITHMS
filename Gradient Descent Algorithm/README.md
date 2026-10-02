# Gradient Descent Polynomial Approximation

This folder fits a polynomial to observed data by minimizing a least-squares objective with the repository's gradient descent optimizer. It solves one optimization problem per degree, starting with a constant and adding powers of x until the fit reaches its residual tolerance or the available degree limit.

## Mathematical formulation

For observations $(x_i,y_i)$, a degree-$d$ polynomial is represented by its coefficient vector

$$
\theta=(\theta_0,\theta_1,\ldots,\theta_d)^T,
\qquad p_d(x)=\sum_{k=0}^{d}\theta_kx^k.
$$

The objective minimized for each degree is the half mean squared error:

$$
J_d(\theta)=\frac{1}{2n}\sum_{i=0}^{n-1}\left(p_d(x_i)-y_i\right)^2.
$$

The factor $1/2$ simplifies the derivative; it does not change the minimizing coefficients. In matrix notation, with a Vandermonde matrix $A_{ik}=x_i^k$, this is $J_d(\theta)=\|A\theta-y\|_2^2/(2n)$. The code evaluates the polynomial by summing powers directly and does not explicitly build or solve this matrix system.

Gradient descent updates the coefficient vector in the negative-gradient direction. Here the gradient is estimated numerically by centered finite differences. At each iteration, an Armijo backtracking line search reduces the step size until the objective decreases sufficiently. The next degree starts from the previous degree's coefficients, with a new coefficient initialized to zero. This warm start preserves the existing polynomial while allowing the optimizer to learn the newly added term.

The objective's residual is reported as

$$
r_d=\frac{\|p_d(x_i)-y_i\|_2}{\|y\|_2},
$$

with the absolute residual norm used when all y values are zero. The fitter retains the candidate with the smallest residual and stops early when a candidate reaches `1e-5`. It considers degrees from 0 through $\min(n-1,50)$, inclusive. If none reaches the tolerance, it returns the best residual candidate encountered. The final output also reports the maximum absolute error on the input points.

## Optimizer behavior

The polynomial fitter sets these optimizer limits:

- Maximum of 10,000 iterations per degree.
- Gradient norm tolerance `1e-8`.
- Objective change tolerance `1e-14`.
- Step norm tolerance `1e-10`.
- Initial step size 1.0, Armijo constant `1e-4`, and backtracking factor 0.5.
- Centering/scaling of the coefficient parameters is disabled (`useScaling = false`).

The optimizer can finish because the gradient is small, the step is small, the objective change is small, the line search fails, or it reaches the iteration limit. The polynomial fitter checks that the returned objective and coefficient vector are finite and throws if they are not. It uses the returned candidate even if the optimizer stopped for a non-convergence status; inspect the residual and fit quality when convergence is difficult.

## Code organization

Implementation files are in this folder; shared interfaces and supporting types are declared in `Headers`.

- `polynomial_main.cpp` loads the root `dataset.txt`, validates that it contains data, calls the polynomial fitter, and converts exceptions into an error message and nonzero exit status.
- `PolynomialObjective.cpp` implements `PolynomialObjective`, which evaluates the half mean squared error and polynomial values, formats a polynomial for output, and contains `Polynomial_Approximator_Gradient_Descent`. The approximator manages the degree loop, warm starts, residual comparisons, timing, summary, and result-file output.
- `GradientDescent.cpp` implements the optimizer loop: it evaluates the objective, computes a numerical gradient, steps opposite that gradient, requests a line-search step, and applies convergence checks.
- `NumericalGradient.cpp` computes each gradient component with centered finite differences.
- `LineSearch.cpp` implements the Armijo backtracking search. It starts at the proposed step size and repeatedly multiplies it by the backtracking factor until the sufficient-decrease condition is met or a limit is reached.
- `Scaling.cpp` provides optional parameter centering and scaling. The polynomial fitter currently leaves this option off.
- `ConvergenceChecker.cpp` contains the gradient, objective-change, and step-size stopping checks.
- `AckleyObjective.cpp` and `gd_tests.cpp` exercise the shared optimizer on other objectives and are not part of the polynomial driver's data-loading path.
- Headers such as `PolynomialObjective.h`, `GradientDescent.h`, `GradientDescentOptions.h`, `NumericalGradient.h`, `LineSearch.h`, `Scaling.h`, `ConvergenceChecker.h`, `ObjectiveFunction.h`, and `OptimizationResult.h` define interfaces, settings, and result/status types.
- `plot_gradient_descent.py` calls the shared root-level plotting script with the gradient-descent algorithm and Gaussian example as defaults.

## Input and output

The configured driver reads `dataset.txt` from the current working directory. The VS Code task runs from the repository root. Each line contains an x value followed by a y value, separated by whitespace; no header is used:

```text
0 1
1 2.5
2 8
3 20
```

The x and y arrays must have the same nonzero length. The fitter supports one data point, in which case it can only attempt a constant model. The candidate degree is capped by the number of observations minus one to avoid fitting an underdetermined polynomial system by construction.

The program prints each candidate's degree, residual, polynomial, and candidate time, followed by the average time per optimizer iteration and the best fit. It appends the selected degree, timing, and coefficients (ascending powers) to `best_polynomial_results.txt` in the current working directory. Existing entries are preserved; repeated runs add more entries.

## Build, run, and plot

Open the repository root in Visual Studio Code and choose **Run Gradient Descent Polynomial Approximation** from the Tasks menu. The task compiles the driver with the optimizer and its helper implementations, then runs it using the root `dataset.txt`.

To plot the result, run `plot_gradient_descent.py` with Python from the repository root. It delegates to `plot_polynomial_fits.py`; the shared script handles plotting options and inputs.

## Numerical considerations

Gradient descent offers an iterative alternative to direct least-squares factorization, but it depends on step selection and convergence settings. Coefficients for high-degree monomial polynomials can have very different scales, and raw powers of large x values can make the objective poorly conditioned. Parameter scaling exists in the shared optimizer but is not enabled by this fitter; scaling or centering x, using a better-conditioned polynomial basis, and checking convergence status can improve challenging fits. Training residual and maximum input-point error measure agreement with the supplied data, not predictive accuracy between or beyond those points.
