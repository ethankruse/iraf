AIRAF
=====
**AIRAF: Astronomical Image Reduction and Analysis Facility**

.. image:: https://travis-ci.com/ethankruse/iraf.svg?token=eppNFbJp1yHxpxgRYHsc&branch=master
    :target: https://travis-ci.com/ethankruse/iraf

AIRAF is a Python implementation of the original `IRAF`_ developed by
the National Optical Astronomy Observatories (NOAO).
The code is open source, well tested, and intended to give identical results
to traditional image reduction and analysis in the original IRAF.

The original goal of this project was to create a native-Python translation 
of the 40-year-old original IRAF package. IRAF was a tremendous resource to 
the astronomy community, and having all observatory data reductions use
the same tool was good for teaching new students, scientific reproducibility,
and productivity.

Now that IRAF has stopped being maintained by NASA and is increasingly impossible to 
install on modern machines, combined with its incredibly antiquated interface,
the data reduction landscape has fractured. Thus nearly all the benefits of
IRAF are gone and the community is stuck piecing together various reduction pipelines
for their own individual case.

I created this package to explore creating a direct Python translation of IRAF, but
with modern graphics and syntax. I downloaded the IRAF source code and taught myself
its custom language to translate the routines into Pythonic reproductions.

Currently, there are 4 working tools (currently nested in the same package structure as the original IRAF):

* ``iraf.plot.implot``, which is an interactive tool identical to the original, but using modern graphics.
* ``iraf.images.imutil.imstat``, which provides stats about the data contents of images.
* ``iraf.noao.imred.ccdred.combine``, which combines multiple images into one (e.g. median flatfield).
* ``iraf.noao.imred.ccdred.ccdproc``, the workhorse of the data reduction pipeline.

Each routine has thorough tests that make sure everything is working and produces 
identical results to using the original IRAF package.

While working through this and writing the tests, I've actually found some bugs 
in IRAF, though none likely to have impacted any scientific results. 


.. _IRAF: https://iraf-community.github.io/
