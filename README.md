# Sedimentology
Collection of functions (largely in Python) for various sedimentological problems. **Very much still a work in progress.**

Much of the code here is either directly adapted from or heavily influenced by the excellent textbook from Allmendinger et al. 2012: Structural geology algorithms: vectors and tensors.
- Find a copy here:  https://doi.org/10.1017/CBO9780511920202
- Find their (probably better written) Matlab routines here: https://www.mathworks.com/academia/books/structural-geology-algorithms-allmendinger.html

Relevant literature cited in files drawing from it - please cite your sources!

## File list
- vectorstats.py: contains functions for calculating vector mean statistics. Created with palaeoflow calculations in mind but probably has other applications out there. Contains the following functions:
  - currayMean
  - plane2pole
  - z2p
  - sph2cart
  - cart2sph
  - fisherMean
  - intersect
- planeintersect.c: a C program taking dip azimuth and dip readings of planes on the command line and calculating their intersection. Probably not hugely useful in and of itself but could be cannibalised into something else (I'd imagine).
