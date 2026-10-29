"""Make the pipeline's bin/ scripts importable as modules in pytest."""

import sys
from pathlib import Path

BIN_DIRECTORY = Path(__file__).parents[2] / "bin"
sys.path.insert(0, str(BIN_DIRECTORY))
