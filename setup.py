from setuptools import setup, find_packages

setup(
    name="first_shell_sdf",
    version="0.1.0",
    packages=find_packages(),
    py_modules=["extract_first_shell"],
    install_requires=[
        "MDAnalysis",
        "numpy",
        "tqdm"
    ],
    entry_points={
        'console_scripts': [
            'first_shell_sdf=extract_first_shell:main',
        ],
    },
    author="Rajib Biswas",
    description="Compute spatial density of first-shell atoms around a selected atom from GROMACS trajectory",
    license="MIT",
    keywords=["MDAnalysis", "GROMACS", "first shell solvation", "spatial density", "dx"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)

