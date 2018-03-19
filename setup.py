from setuptools import setup
from Cython.Build import cythonize

setup(name='advntr',
      version='1.0.0',
      description='A tool for genotyping Variable Number Tandem Repeats (VNTR) from sequence data',
      author='Mehrdad Bakhtiari',
      author_email='mbakhtia@ucsd.edu',
      license='BSD-3-Clause',
      url='https://github.com/mehrdadbakhtiari/adVNTR',
      packages=['advntr', 'pomegranate'],
      package_dir={'advntr': 'advntr', 'advntr.pomegranate': 'pomegranate'},
      provides=["advntr"],
      entry_points={
            'console_scripts': ['advntr=advntr.__main__:main']
      },
      ext_modules=cythonize(["pomegranate/*.pyx"]),
      classifiers=["Environment :: Console",
                   "Intended Audience :: Developers",
                   "Intended Audience :: Science/Research",
                   "Operating System :: OS Independent",
                   "Programming Language :: Python",
                   "Programming Language :: Python :: 2",
                   "Programming Language :: Python :: 3",
                   "Topic :: Scientific/Engineering :: Bio-Informatics"],
      )
