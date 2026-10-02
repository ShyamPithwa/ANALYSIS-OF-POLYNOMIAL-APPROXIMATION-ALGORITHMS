# Analysis of Polynomial Approximation Algorithms

A C++17 project comparing four least-squares polynomial fitting approaches on the same two-column dataset:

- **Householder QR:** fits a monomial Vandermonde system using Householder QR and triangular back substitution.
- **Arnoldi:** constructs an orthonormal basis for sampled Krylov subspaces and projects observations onto it.
- **Chebyshev:** maps x to [-1, 1], fits a Chebyshev basis with Householder QR, then returns power-basis coefficients.
- **Gradient descent:** minimizes half mean squared error with numerical gradients and Armijo backtracking.

Each algorithm folder has a README with the mathematical derivation and a guide to its implementation. All methods read the root `dataset.txt` and append their selected degree and coefficients to `best_polynomial_results.txt`.

## Build and run

Install a C++17 compiler such as GCC and open this folder in Visual Studio Code. Select one of the **Run ... Polynomial Approximation** tasks. Each run task builds its executable first and then reads `dataset.txt` from the repository root. You can edit the dataset as whitespace-separated x/y pairs, one observation per line, with no header.

The VS Code tasks require `g++` to be on PATH. They write executables alongside their algorithm source; executables and generated result files are ignored by Git.

## Plotting

Install the Python packages listed in `requirements.txt`, run all four approximation tasks to populate `best_polynomial_results.txt`, then run the shared plotter from the repository root:

```sh
python plot_polynomial_fits.py
```

To plot one method, use `--algorithm` with `QR-based`, `Arnoldi-based`, `Chebyshev-based`, or `Gradient Descent-based`. The script accepts `--data`, `--results`, `--function` (`gaussian`, `exponential`, or `sine`), and optional `--output` arguments. Each algorithm folder also has a small plotting entry point that selects that method by default.

## Shared implementation

`Headers/` contains the matrix, QR, data utility, polynomial, and optimizer interfaces. `NLA/` provides matrix operations and numerical helpers. The four algorithm folders contain their fitters, drivers, and algorithm-specific documentation. `.vscode/tasks.json` contains only build/run tasks for these four polynomial methods.
