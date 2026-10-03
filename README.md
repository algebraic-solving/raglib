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

1. Install [msolve](https://msolve.lip6.fr) (make sure it is available in your PATH).
2. Install the file interface package between [msolve](https://msolve.lip6.fr)
  and `maple`:
    - download the source file:
      [https://github.com/algebraic-solving/msolve/blob/master/interfaces/msolve-to-maple-file-interface.mpl](https://github.com/algebraic-solving/msolve/blob/master/interfaces/msolve-to-maple-file-interface.mpl);
    - (optional) choose a custom installation path at lines 25 and 26, in the
      above file (default: `/home/<your-login>/libs`);
    - run `maple msolve-to-maple-file-interface.mpl` from the file's directory.
3. Install `RAGlib`:
    - (optional) open `install.mm` and choose a custom installation path at line 6 (default: `/home/<your-login>/libs`);
    - run `maple -q install.mm` from the repository root.

4. If your maple packages have been installed in the folder at path `savelibname`, make sure you have the following line in your `.mapleinit` file, which should be at the root of your home directory:
`libname:=savelibname,libname:` (e.g. `libname:="/home/<your-login>/libs:",libname:`).

## Basic usage
After lauching maple, import `RAGlib` library as follows:
```maple
> with(RAG);
                    [HasRealSolutions, PointsPerComponents]
```

At the moment, `RAGlib` provides two main functions:
- `RAG[HasRealSolutions](eqs, pos, ineqs)` where eqs, pos, ineqs are lists
    of polynomials encoding equality, positivity and non-vanishing
    constraints respectively.   
    It decides whether the set of real solutions to the input
    constraints is empty or not. In case of emptiness, it returns an
    empty list, otherwise, it returns a list of witness points.  
    All such points are encoded by isolating boxes. 
- `RAG[PointsPerComponents](eqs, pos, ineqs)` where eqs, pos, ineqs are lists
    of polynomials encoding equality, positivity and non-vanishing
    constraints respectively.  
    It returns a list of points (encoded with isolating boxes) meeting
    all connected components of the set of real solutions to the input
    polynomial constraints. 

For instance, the call 

    RAG[PointsPerComponents]([x*y-1], [x^2+y^2-4], []);

returns 

    [[x = [1/3, 1/3], y = [3, 3]], [x = [3, 3], y = [523091811282223396986315785267305534675196287038669542741/1569275433846670190958947355801916604025588861116008628224, 1046183622564446793972631570534611069350392574077339085483/3138550867693340381917894711603833208051177722232017256448]], [x = [214641234886981963472998975199122065933/85070591730234615865843651857942052864, 429282469773963926945997950398244131867/170141183460469231731687303715884105728], y = [539468053742263529310033444669866418093/1361129467683753853853498429727072845824, 1078936107484527058620066889339732836187/2722258935367507707706996859454145691648]], [x = [165065513768260236452534495220057220641/340282366920938463463374607431768211456, 165065513768260236452534495220057220643/340282366920938463463374607431768211456], y = [1402983416631528715794310114511939995141/680564733841876926926749214863536422912, 701491708315764357897155057255969997571/340282366920938463463374607431768211456]]]

whose numerical approximation is 

    [[x = [.3333333333, .3333333333], y = [3., 3.]], [x = [3., 3.], y = [.3333333333, .3333333333]], [x = [2.523095591, 2.523095591], y = [.3963385310, .3963385310]], [x = [.4850839474, .4850839474], y = [2.061498851, 2.061498851]]]

The use of computer algebra algorithms make `RAGlib` capable of
tackling problems which may be quite difficult for numerical methods
because of precision issues.

This can be crucial since some semi-algebraic sets can have very thin
connected components. 

For instance, the call 

    RAG[PointsPerComponents]([], [], [(x-1)^2+y^2-1, (x+1)^2+y^2-1-1/2^32]);

returns 

    [[x = [-4, -4], y = [0, 0]], [x = [-1, -1], y = [-3, -3]], [x = [-1, -1], y = [0, 0]], [x =[-1, -1], y = [3, 3]], [x = [1, 1], y = [-2, -2]], [x = [1, 1], y = [0, 0]], [x = [1, 1], y= [2, 2]], [x = [2, 2], y = [-1, -1]], [x = [2, 2], y = [1, 1]], [x = [3, 3], y = [0, 0]], [x = [1/17179869184, 1/17179869184], y = [0, 0]], [x = [1/17179869184, 1/17179869184], y = [-30798844054400941813321769662360858826861/21778071482940061661655974875633165533184, -7699711013600235453330442415590214706715/5444517870735015415413993718908291383296]], [x = [1/17179869184, 1/17179869184], y = [7699711013600235453330442415590214706715/5444517870735015415413993718908291383296, 30798844054400941813321769662360858826861/21778071482940061661655974875633165533184]], [x = [19807040628566084398385987583/340282366920938463463374607431768211456, 19807040628566084398385987585/340282366920938463463374607431768211456], y = [-680564733848479273802510901698334739553/340282366920938463463374607431768211456, -21267647932764977306328465678072960611/10633823966279326983230456482242756608]], [x = [19807040628566084398385987583/340282366920938463463374607431768211456, 19807040628566084398385987585/340282366920938463463374607431768211456], y = [10633823966485650323098009195830204003/10633823966279326983230456482242756608, 340282366927540810339136294266566528097/340282366920938463463374607431768211456]], [y = [0, 0], x = [-3671508318603584465714016303674699/340282366920938463463374607431768211456, -1835754159301792232857008151837349/170141183460469231731687303715884105728]], [y = [0, 0], x= [1835754159301792232857008151837349/170141183460469231731687303715884105728, 3671508318603584465714016303674699/340282366920938463463374607431768211456]]]

whose numerical approximation is 

    [[x = [-4., -4.], y = [0., 0.]], [x = [-1., -1.], y = [-3., -3.]], [x = [-1., -1.], y = [0.,0.]], [x = [-1., -1.], y = [3., 3.]], [x = [1., 1.], y = [-2., -2.]], [x = [1., 1.], y = [0., 0.]], [x = [1., 1.], y = [2., 2.]], [x = [2., 2.], y = [-1., -1.]], [x = [2., 2.], y = [1., 1.]], [x = [3., 3.], y = [0., 0.]], [x = [.5820766091e-10, .5820766091e-10], y = [0., 0.]], [x = [.5820766091e-10, .5820766091e-10], y = [-1.414213562, -1.414213562]], [x = [.5820766091e-10, .5820766091e-10], y = [1.414213562, 1.414213562]], [x = [.5820766091e-10, .5820766091e-10], y = [-2.000000000, -2.000000000]], [x = [.5820766091e-10, .5820766091e-10], y = [1.000000000, 1.000000000]], [y = [0., 0.], x = [-.1078959322e-4, -.1078959322e-4]], [y = [0., 0.], x = [.1078959322e-4, .1078959322e-4]]]


## Future plans

- Many optimizations remain to be implemented, especially for early
    termination when deciding if a system of polynomial constraints is
    consistent.
- At the moment, some singular situations are not handled by `RAGlib`
    (in that case, it raises an error). These implementations will
    come soon.
- More functionalities will come soon, notably, for computing the
    dimension of real algebraic sets or semi-algebraic sets, one-block
    quantifier elimination and roadmap computations for answering
    connectivity queries. 
- A translation of `RAGlib` in [`AlgebraicSolving.jl`](https://algebraic-solving.github.io/) will be hopefully started soon (see also [https://github.com/algebraic-solving/AlgebraicSolving.jl](https://github.com/algebraic-solving/AlgebraicSolving.jl)). 
