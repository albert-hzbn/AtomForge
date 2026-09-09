from setuptools import setup, find_packages

setup(
    name="atomforge-py",
    version="0.1.0",
    description="Atomic structures and native electronic post-processing for AtomForge",
    author="AtomForge Contributors",
    python_requires=">=3.8",
    packages=find_packages(),
    install_requires=[],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
)
