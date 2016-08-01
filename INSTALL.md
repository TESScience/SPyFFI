Installing SPyFFI
=================

These are instructions for installing SPyFFI and its (many) dependencies. One of these days, I will learn how to package this installation more cleanly!

### Python
If you don't have a relatively modern Python 2.7 distribution, download and install [`anaconda`](). It's an easy way to get a self-managed Python distribution, where you have easy control over the libraries and can install packages with something like `pip`.

### Easy Python dependencies
Make sure you have recent versions of the following astronomy-related libraries installed. Using `pip` is an easy way to do this (the `--upgrade` flag will make sure you get the most up-to-date version if you have something already installed):

    pip install matplotlib --upgrade  
    pip install numpy --upgrade  
    pip install astropy --upgrade  
    pip install astroquery --upgrade
    pip install git+https://github.com/zkbt/zachopy.git#egg=zachopy-0.1.0


### SpyFFI Python scripts
Change to the directory into which you want to install SPyFFI (e.g. ~/code/). From, here clone this git repository with

`git clone https://username@github.com/zkbt/SPyFFI/`

to create and populate the directory (~/code/SPyFFI), where you replace `username` with your github account name. You must be listed as a collaborator on the SPyFFI repository to have access -- contact Zach (zkbt@mit.edu) to be added.

For an introduction on how to use this github repository more generally, please see [`git_tutorial.md`](https://github.com/zkbt/SPyFFI/blob/master/git_tutorial.md).

### Set Environment Variables
Make sure you have the following environment variables set, which I do by modifying my `~/.profile` to read:

    ############################
    # TESS -- the SPyFFI tools #
    ############################

    # the root directory of SPyFFI (needed for some file-finding by the code)
    export SPYFFIPATH=$HOME/Dropbox/code/SPyFFI/

    # the directory where you want to store TESS data (inputs/intermediates/outputs)
    export SPYFFIDATA=$HOME/Cosmos/Data/TESS/FFIs/

### Recompile the cosmic ray code
Go into the cosmical_realistic directory, and run `make` to compile the C extension required to generate cosmic rays. (It's an intenstive process, and Al's cosmic ray code is a super-fast).

### Download the SPyFFI inputs
SPyFFI requires a few big files in order to run (primarily the PSF library and a jitterball). Download [SPyFFI_coreinputs.tar.gz](https://www.dropbox.com/s/0e4c2uk34phv4qx/SPyFFI_coreinputs.tar.gz?dl=0), move it to your $SPYFFIDATA directory, run `tar -xvf SPyFFI_inputs.tar.gz` to unpack it.