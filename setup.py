from distutils.core import setup

setup(name='QPyM2',
      version='0.0.1',
      description='CUORE M2 analysis',
      author='Sachi Wagaarachchi',
      author_email='sachi@berkeley.edu',
      packages=['distutils',
                'numpy',
                'pandas',
                'h5py',
                'sphinx',
                'myst-parser'],
)
