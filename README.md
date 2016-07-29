SPyFFI
======

The Spiffy Python for Full Frame Images package is a collection of tools for simulating TESS images. They were created by Zach Berta-Thompson, with contributions from Al Levine, Peter Sullivan, and Deb Woods.

### basic usage

Once you have it installed, you should be able to start to play, either through the command line in `ipython` or through a script, by creating an "Observation" object, and then using that Observation's create() method to generate all its exposures. Observation takes as an input a Python dictionary, which contains the various input parameters (whether you want a test pattern or stars drawn from the real sky, whether images should be jittered, whether stars should be given cartoon light curves, how many exposures of each cadence to make, etc...). An introduction to some of the parameters you may want to change is available in 'scripts/demonstration.py', or you could try the following from an ipython prompt:

    from SPyFFI.Observation import Observation, default
    o = Observation(default)
    o.create()

### other information

For installation instructions, please see
[`INSTALL.md`](https://github.com/TESScience/SPyFFI/blob/master/INSTALL.md)
