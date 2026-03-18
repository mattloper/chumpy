"""
Author(s): Matthew Loper

See LICENCE.txt for licensing and contact information.
"""

from setuptools import setup
from runpy import run_path


def get_version():
    namespace = run_path("chumpy/version.py")
    return namespace["version"]


# Read requirements directly from requirements.txt
# Don't use pip's internal APIs - they're not meant for setup.py
with open("requirements.txt") as f:
    install_requires = [
        line.strip() for line in f if line.strip() and not line.startswith("#")
    ]


setup(
    name="chumpy",
    version=get_version(),
    packages=["chumpy"],
    author="Matthew Loper",
    author_email="matt.loper@gmail.com",
    url="https://github.com/mattloper/chumpy",
    description="chumpy",
    license="MIT",
    install_requires=install_requires,
    python_requires=">=2.7",
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 4 - Beta",
        # Indicate who your project is intended for
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Mathematics",
        # Pick your license as you wish (should match "license" above)
        "License :: OSI Approved :: MIT License",
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
    ],
)
