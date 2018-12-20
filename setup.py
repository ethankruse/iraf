import builtins
from setuptools import setup

# Comment and idea taken from numpy.
# This is a bit hackish: we are setting a global variable so that the main
# airaf __init__ can detect if it is being loaded by the setup routine, to
# avoid attempting to load components that aren't built yet.  While ugly, it's
# a lot more robust than what was previously being used.
builtins.__IRAF_SETUP__ = True

import iraf

# XXX: all of this needs to be looked at and reworked and added before release
setup(
    name="iraf",
    version=iraf.__version__,
    author="Ethan Kruse",
    author_email="ethan.a.kruse@gmail.com",
    url="https://github.com/ethankruse/iraf",
    license="MIT",
    description="Astronomical Image Reduction and Analysis Facility",
    long_description=open("README.rst").read(),
    packages=['iraf'],
    install_requires=["numpy", "matplotlib", "scipy", "astropy", "pytest"],
    classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ]
)
