from setuptools import setup

def readme():
    with open('README') as f:
        return f.read()

setup(name='spopdyn',
      version='alpha',
      description='Spatial meta-population dynamics stochastic simulations',
      long_description=readme(),
      url='http://github.com/geeklhem/spopdyn',
      author='Guilhem Doulcier',
      author_email='guilhem.doulcier@ens.fr',
      license='GPLv3',
      packages=['spopdyn'],
      zip_safe=False,
      install_requires=[
          'matplotlib',
          'numpy',
          'scipy',
          'python-libsbml',
          'SALib'
      ],
      test_suite='nose.collector',
      tests_require=['nose'],
      scripts=['bin/spopdyn'],)
