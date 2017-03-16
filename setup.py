"""Setup file."""
from setuptools import setup
from setuptools import find_packages

setup(name='dropSeqPip',
      version='0.21',
      description='A drop-seq pipeline',
      url='http://github.com/hoohm/Drop-seq',
      author='Roelli Patrick',
      author_email='patrick.roelli@gmail.com',
      license='GNU GPL3',
      packages=['dropSeqPip'],
      package_data={'dropSeqPip': ['Rscripts/*.R', 'Snakefiles/*.snake']},
      zip_safe=False,
      install_requires=['snakemake', 'pyyaml'],
      entry_points={
          'console_scripts': ['dropSeqPip = dropSeqPip.__main__:main']}
      )
