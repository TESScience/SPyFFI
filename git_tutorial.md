Using the SPyFFI git repository
===============================

If you don't have it already, [install git](http://git-scm.com/book/en/v2/Getting-Started-Installing-Git) to your local machine.

Things to Do Once
-----------------
To get started, change into the directory in which you'd like to store the code and run

`git clone https://github.com/zkbt/SPyFFI.git`

This will download all the code in the repository, and start a local branch on your computer.

Configure your local git account (using your own name, e-mail, and preferred editor):

`git config --global user.name "Cecilia Payne-Gaposchkin"`

`git config --global user.email cpg@badass.com`

`git config --global core.editor nano`

`git config --global color.ui true`


To Get the Most Recent Code
---------------------------
From within the repository's directory, run
`git pull`
to pull down the most recent version and attempt to merge them into your local repository.

To Upload Your Updates to the Code
----------------------------------
It's a two-step process. First, you need to add and commit your changed files to your local repository. Second, you push your local repository up to be merged into the online one.

To add your modified files to a staging area (getting them ready to be committed), do something like
`git add myupdatedfile.c`
or more generally
`git add *.py`
or even just
`git add .`

Then, to commit these changes, type
`git commit -m "Add snazzy feature to someneatfile.py, and modified documentation."`
where the commit message inside the double quotes summarizes the key changes this commit makes.

Now, you can push your committed changes up to the online github repository for everyone else to use. Type
`git push`
If you look at the [online view of the repository](https://github.com/zkbt/SPyFFI) you should see your changes propagating through.

To Check the Status of Any Changes
----------------------------------
Running `git status` from the command line will tell you what files you've modified, and whether or not they are staged for the next commit.

Runnning `git log` will show you the history of this repository, documented at every commit that's been made.

Any questions whatsoever, please contact Zach at zkbt@mit.edu!
