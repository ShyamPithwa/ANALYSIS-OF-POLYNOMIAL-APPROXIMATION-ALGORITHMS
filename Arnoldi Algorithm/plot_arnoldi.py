import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from plot_polynomial_fits import main


if __name__ == "__main__":
    main(default_algorithm="Arnoldi-based", default_function="gaussian")
