"""
Author(s): Matthew Loper

See LICENCE.txt for licensing and contact information.
"""

from distutils.core import setup

setup(name='chumpy',
      version='0.52',
      #py_modules=['ch', 'ch_ops', 'linalg', 'utils', 'api_compatibility', 'ch_random', 'test_ch', 'test_inner_composition', 'test_linalg'],
      packages = ['chumpy', 'chumpy.test_ch'],
      package_dir = {'chumpy': '.'},
      author='Matthew Loper',
      author_email='matt.loper@gmail.com',
      url='https://github.com/mattloper/chumpy',
      description='chumpy',
      )

