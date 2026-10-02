# Householder QR Polynomial Approximation

This folder contains a C++ implementation of polynomial least-squares fitting using Householder QR factorization. It builds a Vandermonde matrix for each candidate degree, factors that matrix into orthogonal and upper-triangular components, and solves for the polynomial coefficients without forming the normal equations.

## What the algorithm computes

Given data points $(x_i,y_i)$, the method seeks a polynomial

$$
p_d(x)=c_0+c_1x+c_2x^2+\cdots+c_dx^d
$$

whose values are close to the observations in the least-squares sense. For degree $d$, the Vandermonde matrix and target vector are

$$
A_{ij}=x_i^j,\qquad
A=\begin{bmatrix}1&x_0&x_0^2&\cdots&x_0^d\\
1&x_1&x_1^2&\cdots&x_1^d\\
\vdots&\vdots&\vdots&&\vdots\\
1&x_{m-1}&x_{m-1}^2&\cdots&x_{m-1}^d\end{bmatrix},\quad
y=\begin{bmatrix}y_0&y_1&\cdots&y_{m-1}\end{bmatrix}^T.
$$

The coefficient vector $c=(c_0,\ldots,c_d)^T$ is chosen to minimize $\|Ac-y\|_2$. Coefficients are stored and returned in ascending-power order: index 0 is the constant term, index 1 multiplies $x$, and so on.

## Householder QR factorization

For a matrix $A$, the factorization constructs an orthogonal matrix $Q$ and an upper-triangular matrix $R$ such that $A=QR$. Orthogonality means $Q^TQ=I$. The implementation starts with `R = A` and `Q = I` and processes one column at a time. At column $k$, it forms a reflector vector from the still-active part of that column and applies

$$
H_k=I-2vv^T
$$

to zero the entries below the diagonal. Each reflector is symmetric and orthogonal. The code applies it to the active submatrix of `R` and accumulates the reflectors in `Q`. Choosing the sign of the column norm to avoid cancellation improves the stability of the reflector construction.

Once $A=QR$, the least-squares problem becomes

$$
Rc=Q^Ty.
$$

Only the first $d+1$ rows of $Q^Ty$ and the leading $(d+1)\times(d+1)$ triangular block of $R$ are needed. `back_substitution` solves this system from the last coefficient to the first. It throws an error if a diagonal entry is smaller than `1e-12`, since that indicates a singular or numerically unusable triangular system.

Using QR avoids the extra conditioning penalty of solving the normal equations $A^TAc=A^Ty$, which square the condition number. The raw-power Vandermonde matrix can still be poorly conditioned, especially at high degrees or for large-magnitude x values.

## Degree selection and reported errors

`Polynomial_Approximator(x, y)` considers degrees starting at 2, increasing by one each time. It stops when either:

- the relative residual is at most `1e-5`;
- the degree reaches the data-supported limit (the loop requires `degree < x.size()`); or
- the degree reaches 50.

For each candidate it computes the relative Euclidean residual

$$
r=\frac{\|Ac-y\|_2}{\|y\|_2}.
$$

When all y values are zero, the implementation uses 1 as the denominator, so the reported residual remains well-defined. It retains the candidate with the smallest relative residual, even if the tolerance is never met. The console output includes each candidate degree and residual, the polynomial, candidate calculation time, average time per degree tried, and the selected fit's maximum absolute pointwise error $\max_i|p(x_i)-y_i|$. The maximum-error calculation evaluates the polynomial with Horner's method.

The result is also appended to `best_polynomial_results.txt` in the process's current working directory. It records the algorithm label, selected degree, average candidate time, and coefficients in ascending-power order. Repeated runs append additional entries.

## Folder contents

- `factorizations.cpp` — Householder QR, triangular back substitution, polynomial fitting, data loading, and reporting routines. It uses the shared matrix and numerical utilities in `NLA` and declarations in `Headers`.
- `plot_qr.py` — small entry point to the repository's shared polynomial-fit plotting script, selecting the QR-based algorithm and the Gaussian example by default.

The executable driver is `NLA/householder_qr_test.cpp`, not a source file in this folder. The dataset and VS Code build tasks are also maintained at the repository root.

## Input data

The driver reads `dataset.txt` from the current working directory. In the repository's configured run task, that is the project root. Each nonempty line must contain two whitespace-separated numbers: x followed by y. There is no header. For example:

```text
0 1
1 2.5
2 8
3 20
```

At least three observations are needed for the degree-2 starting candidate to run. The implementation does not explicitly reject shorter datasets; with fewer than three points, the candidate loop is skipped and no meaningful fit is produced. The number of observations also limits the largest attempted degree so the system is not underdetermined by construction. Duplicate x values or other data that make the Vandermonde columns dependent can still cause a singular-system error.

## Build and run

Open the repository root in Visual Studio Code and run **Run Householder QR Test** from the Tasks menu. That task builds the driver with the shared QR, matrix, and numerical utility sources, then runs it from the repository root so it can find `dataset.txt`.

To plot the fit, run `plot_qr.py` with Python from the repository root. It delegates to `plot_polynomial_fits.py`; use that shared script's options if you want to choose a different function or dataset.

## Numerical considerations

The implementation uses ordinary monomial powers of x. Large x ranges and high degrees can make the Vandermonde matrix ill-conditioned, leading to inaccurate coefficients even when the computed residual is small. Centering and scaling x, or using an orthogonal polynomial basis such as Chebyshev polynomials, can improve behavior for such data. A small residual measures agreement at the supplied points; it does not by itself establish that the polynomial generalizes between or beyond those points.
