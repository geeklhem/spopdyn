from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='spopdyn',
      version='alpha',
      description='Spatial meta-population dynamics stochastic simulations',
      long_description=readme();
      url='http://github.com/geeklhem/spopdyn',
      author='Guilhem Doulcier',
      author_email='guilhem.doulcier@ens.fr',
      license='GPLv3',
      packages=['spopdyn'],
      zip_safe=False)
