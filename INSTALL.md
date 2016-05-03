Installing SPyFFI
=================

These are instructions for installing SPyFFI and its (many) dependencies. One of these days, I will learn how to package this installation more cleanly!

### Python
If you don't have a relatively modern Python 2.7 distribution, download and install [`anaconda`](). It's an easy way to get a self-managed Python distribution, where you have easy control over the libraries and can install packages with something like `pip`.

### easy Python dependencies
Make sure you have recent versions of the following astronomy-related libraries installed. Using `pip` is an easy way to do this (the `--upgrade` flag will make sure you get the most up-to-date version if you have something already installed):

    pip install matplotlib --upgrade  
    pip install numpy --upgrade  
    pip install astropy --upgrade  
    pip install astroquery --upgrade
    pip install git+https://github.com/ericmandel/pyds9.git#egg=pyds9 --upgrade
    (the last of these allows you to connect to ds9 directly from Python. If you don't have ds9, get it!)

### SpyFFI Python scripts
Change to the directory into which you want to install SPyFFI (e.g. ~/code/). From, here clone this git repository with

`git clone https://username@github.com/zkbt/SPyFFI/`

to create and populate the directory (~/code/SPyFFI), where you replace `username` with your github account name. You must be listed as a collaborator on the SPyFFI repository to have access -- contact Zach (zkbt@mit.edu) to be added.

For an introduction on how to use this github repository more generally, please see [`git_tutorial.md`](https://github.com/zkbt/SPyFFI/blob/master/git_tutorial.md).

### zachopy Python scripts
Install the zachopy toolkit in exactly the same way, with

`git clone https://github.com/zkbt/zachopy/`

You should also run this from ~/code/ (or whatever you choose), so the zachopy tools will be in your PATH.

### set environment variables
Make sure you have the following environment variables set, which I do by modifying my `~/.profile` to read:

    ############################
    # TESS -- the SPyFFI tools #
    ############################

    # where the Python scripts can be found (i.e. the directory containing SPyFFI/ and zachopy/)
    export PYTHONPATH=$HOME/Dropbox/python/:$PYTHONPATH

    # the root directory of SPyFFI (needed for some file-finding by the code)
    export SPYFFIPATH=$HOME/Dropbox/code/SPyFFI/

    # the directory where you want to store TESS data (inputs/intermediates/outputs)
    export SPYFFIDATA=$HOME/Cosmos/Data/TESS/FFIs/

### recompile the cosmic ray code
Go into the cosmical_realistic directory, and run `python setup.py build_ext --inplace` to compile the C extension required to generate cosmic rays. (It's an intenstive process, and Al's cosmic ray code is a super-fast).

### download the SPyFFI inputs
SPyFFI requires a few big files in order to run (primarily the PSF library and a jitterball). Download [SPyFFI_coreinputs.tar.gz](https://www.dropbox.com/s/0e4c2uk34phv4qx/SPyFFI_coreinputs.tar.gz?dl=0), move it to your $SPYFFIDATA directory, run `tar -xvf SPyFFI_inputs.tar.gz` to unpack it.

### asking for help
These instructions worked for Zach, starting with a clean install on the MIT antares cluster, as of 3/21/2015. It was was fairly smooth, but it *probably* will not work on the first try for you. Sorry! Please e-mail me (zkbt@mit.edu) with problems and I'll try to help as quickly as I can!
