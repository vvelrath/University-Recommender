__author__="Vivekanandh"
__date__ ="$Dec 4, 2014 7:34:08 PM$"

from setuptools import setup,find_packages

setup (
  name = 'DBScanClustering',
  version = '0.1',
  packages = find_packages(),

  # Declare your packages' dependencies here, for eg:
  install_requires=['foo>=3'],

  # Fill in these to make your Egg ready for upload to
  # PyPI
  author = 'Vivekanandh',
  author_email = '',

  summary = 'Just another Python package for the cheese shop',
  url = '',
  license = '',
  long_description= 'Long description of the package',

  # could also include long_description, download_url, classifiers, etc.

  
)