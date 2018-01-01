import builtins
from setuptools import setup

# Hackishly inject a constant into builtins to enable importing of the
# package before all the dependencies are built.
builtins.__IRAF_SETUP__ = True
import iraf

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
    install_requires=["numpy, matplotlib, scipy, astropy, h5py"],
    classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ]
)