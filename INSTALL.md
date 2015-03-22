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


### pyds9 dependency
This allows Python to interface directly with ds9 to display images. (This, of course, assumes you have ds9 install too.) Run this from any directory (and then you can delete the downloaded code if you like):

    wget http://ds9.si.edu/download/pyds9/pyds9-1.7.tar.gz
    tar -xvf pyds9-1.7.tar.gz
    cd pyds9-1.7
    python setup.py install

### SpyFFI Python scripts
Change to the directory into which you want to install SPyFFI (e.g. ~/code/). From, here clone this git repository with

`git clone https://github.com/zkbt/SPyFFI/`

to create and populate the directory (~/code/SPyFFI). You must be listed as a collaborator on the SPyFFI repository to have access -- contact Zach (zkbt@mit.edu) to be added.

For an introduction on how to use this github repository more generally, please see [`git_tutorial.md`](https://github.com/zkbt/SPyFFI/blob/master/git_tutorial.md).

### zachopy Python scripts
Install the zachopy toolkit in exactly the same way, with

`git clone https://github.com/zkbt/zachopy/`

You should also run this from ~/code/ (or whatever you choose), so the zachopy tools will be in your PATH.

### set environment variables
Make sure you have the following environment variables set, which I do by modifying my `/.profile` to read:

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
SPyFFI requires a few big files in order to run (primarily the PSF library and a jitterball). Download and unpack them with the following commands:

    cd $SPYFFIDATA
    wget
    tar -xvf SPyFFI_inputs.tar.gz



### asking for help
These instructions worked for Zach, starting with a clean install on the MIT antares cluster, as of 3/21/2015. It was was fairly smooth, but it *probably* will not work on the first try for you. Please e-mail me (zkbt@mit.edu) with problems and I'll try to help as quickly as I can!
