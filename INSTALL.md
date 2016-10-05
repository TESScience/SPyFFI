Installing SPyFFI
=================

These are instructions for installing SPyFFI and its (many) dependencies.

### Python 2.7

If you don't have a relatively modern Python 2.7 distribution, download and install [`anaconda`](https://www.continuum.io/downloads). It's an easy way to get a self-managed Python distribution, where you have easy control over the libraries and can install packages with something like `pip`.

### Sandboxed Installation

***If you wish to globally install SPyFFI, skip this step.***

It is convenient to make a self-contained sandboxed installation if you are wish to develop SPyFFI or do not otherwise have admin permission on your machine.

To make a sandboxed installation, you will need [`virtualenv`](http://docs.python-guide.org/en/latest/dev/virtualenvs/) installed.  Then type at the command line, assuming you are using the BASH shell:

    virtualenv -p $(which python2.7) spyffi_sandbox
    source ./spyffi_sandbox/bin/activate

Now all `pip` commands (such as the ones below) will install python modules into the `spyffi_sandbox` directory.

### Installing The Latest Release

To install the latest release, type:

    pip install spyffi

### Installing a Developer Snapshot

To install a developer snapshot, type:

    pip install git+https://github.com/TESScience/SPyFFI.git

### Data Files

When you first import `spyffi`, it will download some files into `~/.tess/spyffi`.  If you would like these to be imported to somewhere else, export the `SPYFFIDATA` environment variable in your shell.
