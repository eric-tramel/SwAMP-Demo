# SwAMP Demo User's Manual
![Artist's Rendition of SwAMP](http://i.imgur.com/29r8l6B.png)

Using this demo is supposed to be straightforward: one needs only to
open Matlab, go to the current folder and run the command `demo`.

When the demo starts, a compilation will take place.
SwAMP is written in C and must be compiled using Matlab's MEX API. If you
have a C compiler on your computer, everything should (hopefully) go smoothly!
We have tested the compilation using `gcc` in different platforms, but we'd
expect it to work with other compilers as well. Make sure to run `mex -setup` if you have no previously used Matlab's MEX feature. 

If you have problems, you can try the Python version which, in spite of
being much slower, achieves the same results.

## Key Reference
A. Manoel, F. Krzakala, E. W. Tramel, L. Zdeborová, 
"Sparse Estimation with the Swept Approximated Message-Passing Algorithm," *arXiv submitted.*

## Contributors to this Repository
* **Andre Manoel,** *original source author* `[andremanoel@gmail.com]`
* **Eric W. Tramel,** *maintainer* `[eric.tramel@gmail.com]`

## A few details

- The demo script calls functions from the the `examples`
  folder. By exploring these, one may get a better grasp of how to use
  SwAMP.

- SwAMP's source code is located on the `src` folder; in particular, the
  bulk of the algorithm is contained in the `src/solvers/amp.c` file. 
  This version
  follows exactly the listings in the paper, and is already optimized to
  work with sparse matrices. Additionally, 3 other versions are present in
  the same folder: 
    * `gamp.c`, which implements G-SwAMP; 
    * `amp_dense.c`, a version that isn't optimized for sparse matrices; 
    * and `amp_alt.c`, a slight modification of the algorithm that, in spite of reaching the same results, sometimes converges faster.

