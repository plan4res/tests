# test/PolyhedralFunction

A tester which provides very comprehensive tests for `PolyhedralFunction`
and some tests for any `CDASolver` able to handle `C05Function` in the
objective (such as `BundleSolver`) and any `CDASolver` able to handle
Linear Programs (such as `MILPSolver` and its derived classes
`CPXMILPSolver` , `SCIPMILPSolver` and `GRBMILPSolver`), as well as for 
some of the mechanics of the "core" SMS++ library.

This executable, given the input parameter n, constructs a "random"
`PolyhedralFunction` and put it as the only `Objective` of the
`AbstractBlock` NDOBlock, otherwise "empty" save for the n `ColVariable`
active in the `PolyhedralFunction`, that can have simple bound constraints
imposed on them if the macro `HAVE_CONSTRAINTS` is properly set.

The same `PolyhedralFunction` is represented as a Linear Program in
another `AbstractBlock` (LPBlock) having the same number of `ColVariable`, a
"linear objective" (`FRealObjective` with a `LinearFunction` inside) and
"linear constraints" (`FRowConstraint` with a `LinearFunction` inside).

An appropriate `CDASolver` is attached to NDOBlock, which can be any
`Solver` capable of handling a `C05Function` objective (say,
`BundleSolver`). An appropriate `CDASolver` is attached to LPBlock, which
can be any `Solver` capable of handling Linear Programs (say, some derived
class of `MILPSolver` such as `CPXMILPSolver` , `SCIPMILPSolver` or `GRBMILPSolver`).

After all this is done, the NDOBlock and LPBlock are solved with the
registered `Solver` and the results (termination status and objective
value, if applicable) are compared.

The `PolyhedralFunction` and the LP are then repeatedly randomly modified
"in the same way", and re-solved several times; each time the results
of the two `Solver` are compared.

The usage of the executable is the following:

       ./PolyhedralFunction_test seed [wchg nvar dens #rounds #chng %chng]
       wchg: what to change, coded bit-wise [127]
             0 = add rows, 1 = delete rows 
             2 = modify rows, 3 = modify constants
             4 = change global lower/upper bound
       nvar: number of variables [10]
       dens: rows / variables [4]
       #rounds: how many iterations [40]
       #chng: number changes [10]
       %chng: probability of changing [0.5]

A batch file is provided that runs a not-so-large set of tests with
different sizes and seeds of the random generator; all these passing is a
good sign that no regressions have been done for the tested modules, and
in particular for `PolyhedralFunction`.

A makefile is also provided that builds the executable including the
`BundleSolver` module and all its dependencies, in particular
`MILPSolver` (and, obviously, the core SMS++ library).


## Authors

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa

## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details.
