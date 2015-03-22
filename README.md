SPyFFI
======

The Spiffy Python for Full Frame Images package is a collection of tools for simulating TESS images. They were created by Zach Berta-Thompson, with contributions from Al Levine, Peter Sullivan, and Deb Woods.

### basic usage

Once you have it installed, you should be able to start to play, either through the command line in `ipython` or through a script, with something like the following:

    from SPyFFI import Observation
    o = Observation.TestPattern(magnitudes=[10], subarray=10, nexposures=25, random=True)
    o.expose(jitter=True)

More instructions will be posted here as the user interface gets a little bit more finalized.


### other information

For installation instructions, please see
[`INSTALL.md`](https://github.com/zkbt/SPyFFI/blob/master/INSTALL.md)

For a more introduction on how to use this github repository, please see [`git_tutorial.md`](https://github.com/zkbt/SPyFFI/blob/master/git_tutorial.md).
