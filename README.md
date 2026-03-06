# RAGlib: Real Algebraic Geometry Library

`RAGlib` is a package for solving systems of polynomial equations and
inequalities over the real numbers. At the moment, it allows to decide
the existence of real solutions to systems of equations, (strict)
inequalities, and/or inequations. It can also compute at least one
point in each connected component of the sets of real solutions to
such systems of polynomial constraints. 

It is based on computer algebra algorithms (a.k.a. symbolic
computation) using exact computations. 

It is written using the computer algebra system `maple` and is based
on the [msolve](https://msolve.lip6.fr) library for Groebner bases
computations and real root isolation of systems of equations with
finitely many solutions. 

Full documentation and examples can be found
at [https://msolve.lip6.fr](https://msolve.lip6.fr).  
Source code and installation instructions for
[msolve](https://msolve.lip6.fr) can be found at  
[https://github.com/algebraic-solving/msolve](https://github.com/algebraic-solving/msolve).

## Installation instructions

- Install [msolve](https://msolve.lip6.fr)  
- Install the file interface between [msolve](https://msolve.lip6.fr)
  and `maple` which is given here:  
  [https://github.com/algebraic-solving/msolve/blob/master/interfaces/msolve-to-maple-file-interface.mpl](https://github.com/algebraic-solving/msolve/blob/master/interfaces/msolve-to-maple-file-interface.mpl)  
  (note that you may need here to adapt folder names given in lines 25
  to 28 in the above file, depending on how your home directory is
  organized)  
  After the `msolve` package is created and installed in the folder
  `savelibname` (say `/home/<your-login>/libs`), add the following
  line to your `.mapleinit` file
  which should be at the root of your home directory.  
  `libname:=savelibname,libname:`   (e.g.
  `libname="/home/<your-login>/libs`)  
- add the following line to you `.mapleinit` file   
  `kernelopts(includepath=<src_folder>):`  
  where `src_folder` is the string containing the absolute path to the
  folder containing the sources of `RAGlib`.  
- after launching `maple`, just read the file `rag.mm`. 

## Basic usage

At the moment, `RAGlib` provides two main functions:
- `HasRealSolutions(eqs, pos, ineqs)` where eqs, pos, ineqs are lists
    of polynomials encoding equality, positivity and non-vanishing
    constraints respectively.   
    It decides whether the set of real solutions to the input
    constraints is empty or not. In case of emptiness, it returns an
    empty list, otherwise, it returns a list of witness points.  
    All such points are encoded by isolating boxes. 
- `PointsPerComponents(eqs, pos, ineqs)` where eqs, pos, ineqs are lists
    of polynomials encoding equality, positivity and non-vanishing
    constraints respectively.  
    It returns a list of points (encoded with isolating boxes) meeting
    all connected components of the set of real solutions to the input
    polynomial constraints. 

For instance, the call 
``PointsPerComponents([x*y-1], [x^2+y^2-4], []);``

