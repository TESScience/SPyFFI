Installing SPyFFI
=================

These are instructions for installing SPyFFI and its (many) dependencies.

### Python 2.7

If you don't have a relatively modern Python 2.7 distribution, download and install [`anaconda`](https://www.continuum.io/downloads). It's an easy way to get a self-managed Python distribution, where you have easy control over the libraries and can install packages with something like `pip`.

### Installing The Latest Release

To install the latest release, type:

    pip install spyffi

### Installing a Developer Snapshot

To install a developer snapshot, type:

    pip install git+https://github.com/TESScience/SPyFFI.git

### Data Files

When you first import `spyffi`, it will download some files into `.tess/spyffi`.  If you would like these to be imported to somewhere else, set the `SPYFFIDATA` environment variable in your shell.