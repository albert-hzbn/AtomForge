from setuptools import setup, find_packages

setup(
    name="atomforge-py",
    version="0.1.0",
    description="Python interface for AtomForge: load, edit, and view atomic structures",
    author="AtomForge Contributors",
    python_requires=">=3.8",
    packages=find_packages(),
    install_requires=[],
    extras_require={
        "ase": ["ase>=3.22"],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
)
