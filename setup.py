"""
Author(s): Matthew Loper

See LICENCE.txt for licensing and contact information.
"""

from distutils.core import setup
import importlib
from pip.req import parse_requirements
from os.path import join, split

req_fname = join(split(__file__)[0], 'requirements.txt')
install_reqs = parse_requirements(req_fname, session=False)
install_requires = [str(ir.req) for ir in install_reqs]

setup(name='chumpy',
    version=importlib.import_module('chumpy').__version__,
    packages = ['chumpy'],
    author='Matthew Loper',
    author_email='matt.loper@gmail.com',
    url='https://github.com/mattloper/chumpy',
    description='chumpy',
    license='MIT',
    install_requires=install_requires,

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Mathematics',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',

        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX :: Linux'
    ],
)

