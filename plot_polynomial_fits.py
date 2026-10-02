#!/usr/bin/env python3
"""Plot saved polynomial fits against noisy Gaussian source data."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DEFAULT_RESULTS = Path(__file__).with_name("best_polynomial_results.txt")
DEFAULT_DATA = Path(__file__).with_name("dataset.txt")
EXPECTED_ALGORITHMS = (
    "QR-based",
    "Gradient Descent-based",
    "Arnoldi-based",
    "Chebyshev-based",
)


def read_latest_fits(path, required_algorithms=EXPECTED_ALGORITHMS):
    """Return latest fits and ensure the requested algorithms are available."""
    fits = {}
    current_algorithm = None
    current_degree = None

    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = line.strip()
        if line.startswith("Algorithm:"):
            current_algorithm = line.partition(":")[2].strip()
            current_degree = None
        elif line.startswith("Best degree:"):
            try:
                current_degree = int(line.partition(":")[2].strip())
            except ValueError as error:
                raise ValueError(f"Invalid best degree on line {line_number} of {path}.") from error
        elif line.startswith("Coefficients (ascending powers):"):
            if not current_algorithm or current_degree is None:
                raise ValueError(f"Incomplete fit record before line {line_number} of {path}.")
            raw_values = line.partition(":")[2].strip()
            try:
                coefficients = [float(value) for value in raw_values.split()]
            except ValueError as error:
                raise ValueError(f"Invalid coefficient on line {line_number} of {path}.") from error
            if not coefficients:
                raise ValueError(f"No coefficients found on line {line_number} of {path}.")
            if len(coefficients) != current_degree + 1:
                raise ValueError(
                    f"{current_algorithm} declares degree {current_degree} but has "
                    f"{len(coefficients)} coefficients."
                )
            fits[current_algorithm] = (current_degree, coefficients)

    missing = [name for name in required_algorithms if name not in fits]
    if missing:
        raise ValueError(
            "Missing saved fit(s): " + ", ".join(missing)
            + f". Run the corresponding approximation algorithm first; results are read from {path}."
        )
    return fits


def evaluate_polynomial(coefficients, x_values):
    """Evaluate ascending-power coefficients using Horner's method."""
    values = np.zeros_like(x_values, dtype=float)
    for coefficient in reversed(coefficients):
        values = values * x_values + coefficient
    return values


def main(default_algorithm=None, default_function="gaussian"):
    parser = argparse.ArgumentParser(
        description="Plot saved polynomial fits with noisy Gaussian data and a reference curve."
    )
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS,
                        help="fit results text file (default: best_polynomial_results.txt)")
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA,
                        help="two-column x/y data file (default: dataset.txt)")
    parser.add_argument("--output", type=Path,
                        help="optional image path, such as polynomial_fits.png")
    parser.add_argument("--algorithm", choices=EXPECTED_ALGORITHMS,
                        default=default_algorithm,
                        help="plot only this algorithm (default: plot all four)")
    parser.add_argument("--function", choices=("gaussian", "exponential", "sine"),
                        default=default_function,
                        help=f"reference function to display (default: {default_function})")
    args = parser.parse_args()

    if not args.results.is_file():
        parser.error(f"Results file not found: {args.results}")
    if not args.data.is_file():
        parser.error(f"Data file not found: {args.data}")

    try:
        required_algorithms = [args.algorithm] if args.algorithm else EXPECTED_ALGORITHMS
        fits = read_latest_fits(args.results, required_algorithms)
        data = np.loadtxt(args.data, dtype=float)
    except (OSError, ValueError) as error:
        parser.error(str(error))

    if data.ndim != 2 or data.shape[1] < 2 or data.shape[0] == 0:
        parser.error(f"Expected at least two columns (x y) in {args.data}.")
    x_data, y_data = data[:, 0], data[:, 1]
    # Match the displayed interval to the supplied data.
    x_min, x_max = float(np.min(x_data)), float(np.max(x_data))
    x_grid = np.linspace(x_min, x_max, 1600)
    if args.function == "gaussian":
        reference = 0.4 * np.exp(-0.5 * x_grid**2)
        reference_label = "Gaussian reference (height 0.4, σ=1)"
        data_label = "Noisy Gaussian Data"
    elif args.function == "exponential":
        reference = np.exp(x_grid)
        reference_label = "exp(x) reference"
        data_label = "Exponential Data"
    else:
        reference = np.sin(x_grid)
        reference_label = "sin(x) reference"
        data_label = "Sine-Wave Data"

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(x_data, y_data, color="black", marker="o", s=22,
               alpha=0.7, label="Original data", zorder=3)
    ax.plot(x_grid, reference, color="black", linestyle="--",
            linewidth=1.8, label=reference_label)

    colors = ("tab:blue", "tab:orange", "tab:green", "tab:red")
    algorithms_to_plot = (
        [args.algorithm] if args.algorithm else list(EXPECTED_ALGORITHMS)
    )
    color_by_algorithm = dict(zip(EXPECTED_ALGORITHMS, colors))
    for algorithm in algorithms_to_plot:
        color = color_by_algorithm[algorithm]
        degree, coefficients = fits[algorithm]
        y_fit = evaluate_polynomial(coefficients, x_grid)
        ax.plot(x_grid, y_fit, color=color, linewidth=2,
                label=f"{algorithm} (degree {degree})")

    plot_title = (
        f"{args.algorithm} Polynomial Fit to {data_label}"
        if args.algorithm
        else f"Polynomial Fits to {data_label}"
    )
    ax.set_title(plot_title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xlim(x_min, x_max)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()

    if args.output:
        fig.savefig(args.output, dpi=160, bbox_inches="tight")
        print(f"Saved plot to {args.output}")
    plt.show()


if __name__ == "__main__":
    main()
