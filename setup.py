"""Setup file."""
from setuptools import setup, find_packages
from setuptools.command.install import install


class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        """Command to run after installation is complete."""
        print('before import')
        import rpy2.robjects.packages as rpackages
        from rpy2.robjects.packages import importr
        print('after import')
        utils = rpackages.importr('utils')
        utils.chooseCRANmirror(ind=1)
        packnames = ['ggplot2',
                     'reshape',
                     'gtable',
                     'grid',
                     'stringr',
                     'plyr',
                     'yaml',
                     'gridExtra',
                     'devtools',
                     'ggpubr',
                     'ggseqlogo']
        for package in packnames:
            try:
                importr(package)
            except:
                print('installing {}'.format(package))
                utils.install_packages(package)
        install.run(self)
    

setup(name='dropSeqPipe',
      version='0.23a',
      description='A drop-seq pipeline',
      keywords=['single cell', 'dropseq', 'drop-seq'],
      classifiers = [],
      url='http://github.com/hoohm/dropSeqPipe',
      download_url='http://github.com/hoohm/dropSeqPipe/archive/0.23.tar.gz',
      author='Roelli Patrick',
      author_email='patrick.roelli@gmail.com',
      license='GNU GPL3',
      packages=find_packages(),
      package_data={'': ['*.R',
                        '*.snake',
                        '*.py']},
      zip_safe=False,
      install_requires=['snakemake', 'pyyaml', 'rpy2', 'multiqc'],
      entry_points={
          'console_scripts': ['dropSeqPipe = dropSeqPipe.__main__:main']},
      cmdclass={'install': PostInstallCommand}
      )
