Notes on how to integrate Al's cosmic ray C code into Python:

Read http://dan.iel.fm/posts/python-c-extensions/. It has everything in it. The summary is below:

1) Modify the cosmic ray generating code, defining a function cosmical() in cosmical.c, which returns a pointer to an array of doubles
2) Make sure the function definition in cosmical.h matches the one in cosmical.c
3) Make sure the variable parsing is correct in _cosmical.c
4) Run `python setup.py build_ext --inplace` from command line.
5) Use `from import cosmical_realistic._cosmical` to import the function for use in Python.
