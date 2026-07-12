(** This file is part of RAGlib (Real Algebraic Geometry Library).
 *
 * RAGlib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * RAGlib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with RAGlib.  If not, see <https://www.gnu.org/licenses/>
 *
 * Authors:
 * Mohab Safey El Din **)

(**
TODO: 
- comment inputs/ouputs of all functions
- secure the use of MSolve package
- management of overlapping intervals should be done inside msolve
  (see ComputeBoundsRegular and ComputeBoundsSingular)
- FindGenericLine can be made significantly more efficient (especially
  in the singular case)
- implement early termination for detecting the emptiness of the
  semi-algebraic set under study
- better choices of constraints in the incremental solving procedure
**)

RAG:=module()
option package, load=ModuleLoad;

export
HasRealSolutions, PointsPerComponents;

local
IsMinimal, MinimalGeneratorsDichotomy, MinimalGenerators, WitnessLift,
MonomialSupport, TableCoeffsSinglePoly, OldTableCoeffs, TableCoeffs,
LiftPolynomials, NewValuesWitness, GeneratorsLift, ElimModSatIntersect,
ElimModSatIntersectLM, ModSatIntersect, ModSatIntersectLM,
ElimSaturateIntersect, SaturateIntersect, ComputeMaximalMinors, IsRegular,
HaveFiniteIntersections, GoodFiberValue_svars, GoodFiberValue,
DegreeTruncate, Increment, NextForms, TestGenericLineDegreeSingular,
TestGenericLineDegreeRegular, HasFiniteCriticalLocus, FindGenericLineRegular,
FindGenericLineSingular, CoeffDeform_eps, CoeffDeform, FindGenericLine,
SmallMidRational, ConstructFibers, HasOverLapCoupleOfIntervals, HasOverLap,
ComputeBoundsRegular, ManageOverLapComputeBoundsSingular,
ElimComputeBoundsSingular, ComputeBoundsSingular, ComputeBounds,
ModularLimits, LimitsDeformedCriticalPoints, CriticalPointsSingular,
SplitSystem_cstr, SplitSystem, CriticalPoints, ExactSolSelection,
GenerateCriticalPointsFamilies, FamCriticalPoints, UnivariateSolveFamily,
ConstantSolveFamily,
AdmissibleSolutions, GenerateDeformedFamilies_eps,
OLDGenerateDeformedFamilies_eps, ConstrainedValues, UnboundedComponents,
DegenerateDeformedSystem, InfiniteBranches, ZeroDimBoundaries, SolveFamily,
SemiAlgebraicSolveIterateOnFamilies, PointsPerComponentsAlgebraic,
SemiAlgebraicSolve, ModuleLoad;

ModuleLoad:=proc()
# The msolve library should installed, compiled and its binary should
# be accessible from your PATH
with(MSolve);
end proc:
$include<multi-modular.mm>: 
$include<rag-subroutines.mm>:
$include<rag-main.mm>: 


end module:
