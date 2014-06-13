# User's manual for the SwAMP demo

Using this demo is supposed to be straightforward: one needs only to
open MATLAB, go to the current folder and type 'demo'.

As soon as the demo starts, the compilation process is going to take place:
SwAMP is written in C, and must be compiled using MATLAB's MEX API. If you
have a C compiler in your computer, everything should (hopefully) go smoothly!
We have tested the compilation using gcc in different platforms, but we'd
expect it to work with other compilers as well.

If you have problems, you can try the Python version, which in spite of
being much more slow achieves the same results.

## A few details

- The demo script calls on its turn functions that are on the 'examples'
  folder. By exploring these, one may get a better grasp of how to use
  SwAMP.
  
- SwAMP's source code is located on the 'src' folder; in particular, the
  bulk of the algorithm is contained in the 'src/solvers/amp.c' file. This version
  follows exactly the listings in the paper, and is already optimized to
  work with sparse matrices. Additionally, 3 other versions are present in
  the same folder: 'gamp.c', which implements G-SwAMP; 'amp_dense.c',
  a version that isn't optimized for sparse matrices; and 'amp_alt.c', a slight
  modification of the algorithm that, in spite of reaching the same
  results, sometimes converges faster.
