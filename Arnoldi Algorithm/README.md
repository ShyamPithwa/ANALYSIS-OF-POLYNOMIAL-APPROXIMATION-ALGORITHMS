# Arnoldi Krylov Polynomial Approximation

This folder implements polynomial least-squares approximation using the Arnoldi process. Instead of constructing and solving a separate Vandermonde least-squares problem for every degree, it builds an orthonormal basis for the sampled polynomial space and uses that basis to generate a sequence of candidate fits.

## Mathematical idea

Suppose the input contains $n$ observations $(x_i,y_i)$. Let

$$
X=\operatorname{diag}(x_0,x_1,\ldots,x_{n-1}),\qquad
q_0=\frac{1}{\sqrt n}(1,1,\ldots,1)^T.
$$

Multiplying a vector of sampled values by $X$ multiplies its value at each data point by that point's x-coordinate. Thus, $Xq_0$ contains the samples of $x/\sqrt n$, and repeated applications of $X$ span the sampled monomials $1,x,x^2,\ldots$. In exact arithmetic, Arnoldi orthogonalization builds an orthonormal basis $q_0,q_1,\ldots$ for these Krylov spaces:

$$
\mathcal K_{k+1}(X,q_0)=\operatorname{span}\{q_0,Xq_0,\ldots,X^kq_0\}.
$$

For each column $q_k$, the iteration computes $v=Xq_k$, subtracts its projections onto the existing basis vectors, and records the projection coefficients in an upper-Hessenberg matrix $H$:

$$
h_{j,k}=q_j^T Xq_k,\qquad
v\leftarrow v-h_{j,k}q_j\quad(j=0,\ldots,k),
$$

$$
h_{k+1,k}=\|v\|_2,\qquad
q_{k+1}=v/h_{k+1,k}.
$$

The resulting relation is $XQ_k=Q_kH_k+h_{k+1,k}q_{k+1}e_k^T$. The columns of $Q$ are orthonormal, so projecting the observation vector onto the first $k+1$ vectors gives the least-squares fit within that sampled polynomial subspace:

$$
\hat y_k=Q_kQ_k^Ty=\sum_{j=0}^{k}(q_j^Ty)q_j.
$$

Because the diagonal matrix $X$ acts on sampled polynomial values, each basis vector can also be represented by a polynomial. If $q_j=p_j(x_i)/\sqrt n$ at every sample, then $p_0(x)=1$, and Arnoldi's relation gives the recurrence

$$
p_{k+1}(x)=\frac{x p_k(x)-\sum_{j=0}^{k}h_{j,k}p_j(x)}{h_{k+1,k}}.
$$

The fitted polynomial corresponding to the projection is therefore

$$
P_k(x)=\sum_{j=0}^{k}\frac{q_j^Ty}{\sqrt n}\,p_j(x).
$$

The division by $\sqrt n$ converts the normalized sample-vector basis $q_j$ back to polynomial values. The coefficients returned by the `Polynomial` class are in ascending powers of x.

## Candidate selection

The implementation evaluates candidates from degree 0 upward. After each added basis vector it computes the relative Euclidean residual

$$
r_k=\frac{\|y-\hat y_k\|_2}{\|y\|_2},
$$

using the absolute residual norm when $\|y\|_2=0$. It stops at the first candidate with residual at most `1e-5`. Arnoldi iteration is capped at $\min(n-1,50)$ steps and can stop earlier if the next-vector norm is less than `1e-12` (breakdown, which means the sampled Krylov space no longer grows). If no candidate meets the residual tolerance, the last generated candidate is selected. This is a residual-tolerance rule, not a minimum-degree penalty or a separate validation procedure.

For each candidate, console output includes its degree, polynomial, residual, and candidate calculation time. The final summary includes the selected degree, residual, average candidate time, and polynomial. The selected result is also appended to `best_polynomial_results.txt` in the program's current working directory.

## Code structure

- `Arnoldi_iteration.cpp` contains the Arnoldi iteration and the `Arnoldi_Krylov_Polynomial` fitting workflow. The constructor reads the data matrix into x and y vectors, creates the diagonal matrix $X$ and initial vector $q_0$, runs Arnoldi once, builds the polynomial basis from the columns of $H$, projects y onto the columns of $Q$, and evaluates candidates.
- `polynomial.cpp` implements the small `Polynomial` value type: increase the power by multiplying by x, add a scaled polynomial, scale all coefficients, evaluate, and retrieve coefficients.
- `arrnoldi_test.cpp` is the runnable driver. It loads `dataset.txt` into a `Matrix` and constructs `Arnoldi_Krylov_Polynomial`; the constructor performs and reports the fit.
- `test.cpp` is another minimal driver for the same workflow.
- `plot_arnoldi.py` delegates plotting to the shared root-level `plot_polynomial_fits.py`, defaulting to the Arnoldi-based Gaussian example.
- `Headers/Arnoldi_iteration.h` declares the iteration, polynomial, and approximator classes. Shared declarations and matrix utilities live in `Headers` and `NLA`.

## Input format

The configured driver reads `dataset.txt` from its current working directory. The repository's **Run Arnoldi Test** VS Code task runs it from the repository root. Use one observation per line, with x then y separated by whitespace and no header:

```text
0 1
1 2.5
2 8
3 20
```

The matrix loader supports the two-column observation layout, and the fitting constructor also handles a 2-by-n matrix with x values in its first row and y values in its second row. Use at least one data point; with a single sample, the only possible candidate is the constant fit.

## Build and run

Open the repository root in Visual Studio Code and run **Run Arnoldi Test** from the Tasks menu. This builds the driver with the Arnoldi implementation, shared polynomial implementation, matrix utilities, and QR support used by shared numerical routines, then runs it against the root `dataset.txt`.

To plot the Arnoldi fit, run `plot_arnoldi.py` with Python from the repository root. It uses the common plotting script; options for selecting other functions or inputs are documented there.

## Numerical notes

Arnoldi creates an orthogonal basis, which helps control the loss of independence that occurs when directly working with high powers of x. The implementation performs one sequential orthogonalization pass and uses a fixed breakdown threshold of `1e-12`; floating-point roundoff can still reduce orthogonality for difficult data. Repeated or tightly clustered x values can cause early breakdown because the sampled polynomial space has limited dimension. Degree is also capped at 50. A low residual describes fit to the supplied samples and does not guarantee good behavior between or outside them.
