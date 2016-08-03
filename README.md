SPyFFI
======

[![Build Status](https://travis-ci.org/TESScience/SPyFFI.svg?branch=master)](https://travis-ci.org/TESScience/SPyFFI)

The *Spiffy Python for Full Frame Images* package is a collection of
tools for simulating TESS images.

### Basic usage

Once you have it installed, you should be able to start to play,
either through the command line in `ipython` or through a script, by
creating a `SPyFFI.Observation` object, and then using that
`Observation`'s `create()` method to generate all its exposures.
Observation takes as an input a Python dictionary, which contains the
various input parameters (whether you want a test pattern or stars
drawn from the real sky, whether images should be jittered, whether
stars should be given cartoon light curves, how many exposures of each
cadence to make, etc...).

An introduction to some of the parameters you
may want to change is available in
[`scripts/demonstration.py`](scripts/demonstration.py), or you
could try the following from an `ipython` prompt:

    from SPyFFI.Observation import Observation, default
    o = Observation(default)
    o.create()

### Other Information

For installation instructions, please see [`INSTALL.md`](INSTALL.md)

For additional details, please see the 
[`User Manual`](https://docs.google.com/document/d/1EYwhLq8iRSLVoeTKls7dGEf4LrJ14vA-9UyhfDGhnVA)

### Contributors

SPyFFI was created by Zach Berta-Thompson.

Other contributors include:

  - Jacobi Kosiarek
  - Al Levine
  - Peter Sullivan
  - Deb Woods
  - John Doty
  - Matthew Wampler-Doty