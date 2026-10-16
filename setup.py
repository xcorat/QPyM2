from setuptools import setup, find_packages

setup(
    name='qpym2',
    version='0.5.1',
    description='CUORE M2 analysis',
    author='Sachi Wagaarachchi',
    author_email='sachi@berkeley.edu',
    packages=find_packages(),
    install_requires=[
        'numpy>=2.0',
        'pandas', 
        'h5py',
        'netcdf4',
        'pymc>=5.25',
        'arviz>=0.19',
        'matplotlib',
        'scipy>=1.13',
        'numba>=0.60',
        'PyYAML>=6.0'
    ],
    extras_require={
        'test': ['pytest>=6.0', 'psutil'],
        'dev': ['pytest>=6.0', 'psutil']
    },
    python_requires='>=3.8',
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
)
