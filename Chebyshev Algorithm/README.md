# Chebyshev Polynomial Approximation

This folder fits data with polynomials expressed in the Chebyshev basis. It maps the observed x-range to [-1, 1], builds a Chebyshev design matrix for each candidate degree, and solves each least-squares problem with the shared Householder QR implementation. The final result is converted to ordinary powers of x before it is returned or saved.

## Mathematical method

Given n observations (x_i, y_i), first map each x value to a scaled coordinate:

$$
t_i = \frac{x_i - c}{s}, \qquad
c = \frac{x_{\min}+x_{\max}}{2}, \qquad
s = \frac{x_{\max}-x_{\min}}{2}.
$$

For nonconstant x data, this maps the smallest and largest inputs to -1 and 1. If all x values are equal, the implementation sets s = 1 to avoid division by zero; in that case the scaled inputs are all zero.

The Chebyshev polynomials are defined by

$$
T_0(t)=1, \qquad T_1(t)=t, \qquad
T_{j+1}(t)=2tT_j(t)-T_{j-1}(t).
$$

A degree-d approximation has the form

$$
p_d(x)=\sum_{j=0}^{d} a_j T_j\!\left(\frac{x-c}{s}\right).
$$

For each degree, the code builds a matrix whose entries are

$$
A_{ij}=T_j(t_i), \qquad 0 \le i < n,\quad 0 \le j \le d.
$$

It solves the least-squares problem

$$
\underset{a}{\operatorname{minimize}}\;\|Aa-y\|_2
$$

using Householder QR. With A = QR and Q orthogonal, the coefficient vector is obtained by solving the upper-triangular system

$$
R_{1:d+1,\,1:d+1}a=(Q^Ty)_{1:d+1}
$$

by back substitution. QR avoids explicitly forming A-transpose times A, which would square the condition number.

Chebyshev polynomials are orthogonal on [-1, 1] under the weight (1-t^2)^(-1/2). The discrete data points need not be Chebyshev nodes, so the code still solves a least-squares problem for the actual sample locations. Scaling x and using this basis usually behaves better at higher degrees than directly using raw powers of a large or shifted x range.

## Converting the result to powers of x

The fitting coefficients a_j belong to the scaled Chebyshev basis. The public function returns ordinary power-basis coefficients instead. The code constructs each basis polynomial in x using

$$
\widetilde T_0(x)=1, \qquad
\widetilde T_1(x)=\frac{x-c}{s}, \qquad
\widetilde T_{j+1}(x)=2\frac{x-c}{s}\widetilde T_j(x)-\widetilde T_{j-1}(x),
$$

then combines them as

$$
p_d(x)=\sum_{j=0}^{d}a_j\widetilde T_j(x)
      =\sum_{k=0}^{d}b_kx^k.
$$

The returned vector stores b in ascending-power order: index 0 is the constant term, index 1 multiplies x, etc. The converted power coefficients may be less numerically stable to evaluate at high degrees than the Chebyshev representation; conversion is performed to match the other polynomial fitters' output format.

## Degree selection and residual

The approximator tests degrees from 0 through min(n - 1, 50). For each candidate it computes the relative Euclidean residual

$$
r_d=\frac{\|Aa-y\|_2}{\|y\|_2}.
$$

If y is the zero vector, the denominator is set to 1, leaving the residual as an absolute norm. The fit stops at the first degree with residual at most 1e-5. It retains the smallest-residual candidate seen, so if no candidate reaches the tolerance, the best attempted fit is returned. Candidate output includes the degree, residual, converted polynomial, and calculation time; the summary reports the selected degree, residual, and average time per degree tried.

The selected coefficients and degree are appended to `best_polynomial_results.txt` in the program's current working directory. Values with magnitude below 5e-7 are shown as zero when printing, but the saved and returned coefficients retain their computed values.

## Code organization

- `Chebyshev.cpp` contains `Orthogonal_Polynomial_Approximator`: x scaling, Chebyshev design matrix construction, Householder QR solve, residual and degree selection, conversion to power coefficients, and reporting.
- `orthogonal_main.cpp` is the executable driver. It loads `dataset.txt` through the shared `readDataPoints` helper, invokes the approximator, and reports exceptions.
- `plot_chebyshev.py` delegates to the shared root-level `plot_polynomial_fits.py`, defaulting to the Chebyshev-based Gaussian example.
- Shared QR declarations and data helpers are in `Headers/factorizations.h`; Householder QR and back substitution are implemented in `QR Algorithm/factorizations.cpp`.
- Matrix operations and utility functions are provided by `NLA` and declared in `Headers`.
- `Polynomial` is declared in `Headers/Arnoldi_iteration.h` and implemented in `Arnoldi Algorithm/polynomial.cpp`. The Chebyshev routine uses it to build the converted power-basis polynomial.

## Input, build, and run

The driver reads `dataset.txt` from the current working directory. In the configured VS Code task, the executable runs from the repository root. Each line contains x followed by y, separated by whitespace, with no header:

```text
0 1
1 2.5
2 8
3 20
```

The x and y vectors must have the same nonzero length. Open the repository root in Visual Studio Code and run **Run Orthogonal Polynomial Approximation** from the Tasks menu. Despite the task label, this is the Chebyshev fitter. It builds the driver and shared numerical sources, then runs on the root `dataset.txt`.

To plot this method, run `plot_chebyshev.py` with Python from the repository root. It uses the shared plotting script and its options for choosing inputs and functions.

## Numerical considerations

Scaling the data range and using a Chebyshev basis can improve conditioning, but cannot fix every rank-deficient or badly sampled problem. The QR back substitution code reports an error if a diagonal entry is too close to zero. Very high-degree conversion into powers of x can also create large coefficients through cancellation. The reported residual measures the fit at the supplied observations and does not guarantee accuracy between or outside them.
