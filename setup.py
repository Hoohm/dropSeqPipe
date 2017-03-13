"""Setup file."""
from setuptools import setup
from setuptools import find_packages

setup(name='dropSeqPip',
      version='0.1',
      description='A drop-seq pipeline',
      url='http://github.com/hoohm/Drop-seq',
      author='Roelli Patrick',
      author_email='patrick.roelli@gmail.com',
      license='GNU GPL3',
      packages=['dropSeqPip'],
      package_data={'dropSeqPip':['Rscripts/*.R', 'Snakefiles/*.snake']},
      zip_safe=False,
      #include_package_data=True,
      #data_files=[('Rscripts', ['knee_plot.R', 'species_plot.R']),
                  #('Snakefiles', ['pre_align.snake', 'star_align.snake', 'post_align.snake'])],
      entry_points={
          'console_scripts': ['dropSeqPip = dropSeqPip.__main__:main']}
      )
