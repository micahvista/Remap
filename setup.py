from setuptools import setup, find_packages
from pathlib import Path

# Simple helper to read the README (optional)
this_dir = Path(__file__).resolve().parent
readme_path = this_dir / "README.md"
if readme_path.exists():
    long_description = readme_path.read_text(encoding="utf-8")
else:
    long_description = "Remap"

setup(
    name="remap",
    version="0.1.0",
    description="Remap / liftover pipeline for long-read alignments.",
    # ...
    packages=find_packages("src"),
    package_dir={"": "src"},
    python_requires=">=3.10",
    install_requires=[
        "biopython>=1.86",
        "cigar>=0.1.3",
        "numba>=0.62.1",
        "numpy>=2.2.0",
        "pysam>=0.23.3",
        "vacmap-index==0.0.2",
    ],
    entry_points={
        "console_scripts": [
            "remap=remap.remap:main",   # NOTE: module is remap/remap.py
        ],
    },
)