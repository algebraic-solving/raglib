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

## Multi-modular routines used in RAGlib
## Most of them should be implemented in msolve in a not too far
## future (hopefully)

IsMinimal:=proc(candidates, lm, fc, vars, opts:={})
local newlm:
  newlm:=MSolve:-MSolveGroebnerLM(candidates, fc, vars, opts):
  newlm:=convert(newlm, set):
  return evalb(newlm=lm);
end:

MinimalGeneratorsDichotomy:=proc(gb, sys, a, b, lm, fc, vars, opts:={})
local N, boo;
  if b-a<=1 then 
    boo:=IsMinimal([op(gb[1..b]), op(sys)], lm, fc, vars, opts):
    if boo=true then 
      return b; 
    else 
      return MinimalGeneratorsDichotomy(gb, sys, b, b+1, lm, fc, vars, opts);
    fi;
  fi;
  N:=iquo(a + b, 2);

  boo:=IsMinimal([op(gb[1..N]), op(sys)], lm, fc, vars, opts):
  if boo then 
    return MinimalGeneratorsDichotomy(gb, sys, a, N, lm, fc, vars, opts);
  else 
    return MinimalGeneratorsDichotomy(gb, sys, N, b, lm, fc, vars, opts);
  end if;
end proc:

#gb is a GB for the grevlex order over vars, fc is the characteristic 
#gb is ordered increasingly
#sys is a list of extra polynomials known to belong to the ideal 
#returns the number of elements of gb needed to generate the same ideal 
#as gb starting from the first element of gb 
MinimalGenerators:=proc(gb, sys, vars, fc, opts:={})
local g, lm;
  lm:=map(g->Groebner:-LeadingMonomial(g, tdeg(op(vars))), gb):
  lm:=convert(lm, set):
  return MinimalGeneratorsDichotomy(gb, sys, 1, nops(gb), lm, fc, vars, opts);
end:

##############################################################################
##############################################################################

#newfc is a prime
#witnessmod is a list of witness coefficients (taken modulo modulus)
#newvaluesmod is a list of witness coefficients (taken modulo newfc)
WitnessLift:=proc(newfc, witnessmod, newvaluesmod, modulus)
local newmodulus, i, newwitnessmod, u, witness, lc, monomials:
  newmodulus:=modulus*newfc:
  
  newwitnessmod:=chrem([witnessmod, newvaluesmod], 
                       [modulus, newfc]):

  witness:=[]:
  for i from 1 to nops(witnessmod) do 
    u:=iratrecon(newwitnessmod[i], newmodulus):
    if evalb(u=FAIL) then 
      return [], newwitnessmod, newmodulus;
    else 
      witness:=[op(witness), u]:
    end if;
  end do;
  return witness, newwitnessmod, newmodulus;
end proc:

#sys is a list of polynomials
#lsuport is a list monomial support  
#vars is the list of variables
#
#returns boo, set such that boo=false and set = newsupport of some element in
#sys have a support larger than the corresponding element in lsupport, otherwise
#it returns true, {}
MonomialSupport:=proc(sys, lsupport, vars)
local monomials, lc, m, i, boo, newlsupport:
  monomials:={}:
  newlsupport:=[]:
  boo:=true:
  for i from 1 to nops(sys) do 
    lc:=[coeffs(sys[i], vars, 'm')]:
    m:={m}:
    if m subset convert(lsupport[i], set) then 
      newlsupport:=[op(newlsupport), lsupport[i]]:
    else
      newlsupport:=[op(newlsupport), sort(convert(m, list), 
      (a, b)->Groebner:-TestOrder(a, b, tdeg(op(vars))))]:
      boo:=false;
    end if;
  end do;
  if boo = false then
    return false, newlsupport;
  end if;
  return true, {};
end proc:

#pol is a polynomial
#vars is the list of variables
#support is a list of monomials containing the support of pol
#support is assumed to be ordered increasingly by grevlex(vars)
#returns the list of corresponding coefficients
TableCoeffsSinglePoly:=proc(pol, vars, support)
local i, k, lc, lm, lctable;
  lc:=[coeffs(pol, vars, 'lm')]:
  lm:=sort([lm], (a, b)->Groebner:-TestOrder(a, b, tdeg(op(vars)))):
  lctable:=Array([seq(0, i=1..nops(support))]):
  for i from 1 to nops(lm) do
    member(lm[i], support, 'k');
    lctable[k]:=lc[i];
  end do;
  return convert(lctable, list);
end proc:

OldTableCoeffs:=proc(sys, vars, lsupport, lctables)
local newlctables, i, lc, pol, k;
  newlctables:=Array([seq([], i=1..nops(sys))]);
  #newlctables:=[seq([], i=1..nops(sys))];
  for i from 1 to nops(sys) do 
    pol:=sys[i];
    lc:=TableCoeffsSinglePoly(pol, vars, lsupport[i]);
    newlctables[i]:=[seq([op(lctables[i][k]), lc[k]],k=1..nops(lc))]:
  end do;
  return convert(newlctables, list);
end proc:

TableCoeffs:=proc(sys, vars, lsupport, lctables)
local newlctables, i, lc, pol, k;
  newlctables:=[];
  for i from 1 to nops(sys) do 
    pol:=sys[i];
    lc:=TableCoeffsSinglePoly(pol, vars, lsupport[i]);
    newlctables:=[op(newlctables), [seq([op(lctables[i][k]), lc[k]],k=1..nops(lc))]]:
  end do;
  return convert(newlctables, list);
end proc:

LiftPolynomials:=proc(lctables, support, primetable, modulus, islifted)
local i, j, length, lifted, pol, cc;
  lifted:=[];
  for i from 1 to nops(lctables) do 
    pol:=0:
    if i > islifted then 
      length:=nops(lctables[i]);
      for j from 1 to nops(lctables[i]) do 
        cc:=chrem(lctables[i][j], primetable);
        cc:=iratrecon(cc, modulus);
        if evalb(cc=FAIL) then return lifted; end if;
        pol:=pol+cc*support[i][length-j+1];
      end do;
      lifted:=[op(lifted), pol];
    end if;
  end do;
  return lifted;
end proc:

NewValuesWitness:=proc(newsys, vars, support)
local newvaluesmod, boo, m, i, pol, lc;
  newvaluesmod:=[]:
  for i from 1 to nops(newsys) do 
    pol:=newsys[i]:
    lc:=[coeffs(pol, vars, 'm')]:
    m:=sort([m], (a, b)->Groebner:-TestOrder(a,b,tdeg(op(vars)))):
    boo:=member(support[i][-1], m,'k');
    newvaluesmod:=[op(newvaluesmod),lc[k]]:
  end do;
  return newvaluesmod;
end proc:

GeneratorsLift:=proc(newsys, newfc, vars, modulus, witnessmod, support,
systable, primetable)
local boo, newsupport, i, pol, ls, newvaluesmod, newsystable,
newprimetable, lc, witness, newwitnessmod, newmodulus, k;
  newsystable := [op(systable), newsys]:
  newprimetable := [op(primetable), newfc]:
  boo, newsupport := MonomialSupport(newsys, support, vars):
  if boo = false then 
    return false, false, newsupport, newsystable, newprimetable, 0, 0, [];
  end if;

  newvaluesmod:=NewValuesWitness(newsys, vars, support):

  witness, newwitnessmod, newmodulus := WitnessLift(newfc, witnessmod, newvaluesmod, 
                                        modulus):
  if nops(witness)=0 then 
    #witness coefficients could not be lifted
    return true, false, support, newsystable, newprimetable, newwitnessmod, newmodulus, [];
  end if;
  return true, true, support, newsystable, newprimetable, newwitnessmod, newmodulus, witness;
end proc:

ElimModSatIntersect:=proc(eqs1, pol, eqs2, fc, vars, rag_sep_elem, opts:={})
local sqf, var, gb1, gb2, i, p;
  var:=cat(x,0);
  i:=1:
  while member(var, vars) do 
    var:=cat(u,i);
    i:=i+1:
  end do;
  gb1:=MSolve:-MSolveGroebner([var*pol-1,op(eqs1)], fc, [var, op(vars)],
  {"elim"=1} union opts):
  gb1:=map(pol->if indets(pol) subset indets(vars) then pol fi, gb1):
  gb2:=MSolve:-MSolveGroebner([op(gb1), pol, op(eqs2)], fc, [op(vars), rag_sep_elem], opts union {"elim"=nops(vars), "verb"=0}):
  if gb2=[1] then return [1]; end if;
  gb2:=map(p->if indets(p) = {rag_sep_elem} then p fi, gb2):
  sqf:=Sqrfree(gb2[1]) mod fc;
  gb2:=mul(p[1], p in sqf[2]):
  gb2:=expand(gb2) mod fc;
  return [gb2];
end proc:

ElimModSatIntersectLM:=proc(eqs1, pol, eqs2, fc, vars, rag_sep_elem, opts:={})
local var, gb1, gb2, i, p, sqf;
  var:=cat(x,0);
  i:=1:
  while member(var, vars) do 
    var:=cat(u,i);
    i:=i+1:
  end do;
  gb1:=MSolve:-MSolveGroebner([var*pol-1,op(eqs1)], fc, [var, op(vars)],
  {"elim"=1} union opts):
  gb1:=map(pol->if indets(pol) subset indets(vars) then pol fi, gb1):
  gb2:=MSolve:-MSolveGroebner([op(gb1), pol, op(eqs2)], fc, [op(vars), rag_sep_elem], opts union {"elim"=nops(vars)}):
  if gb2=[1] then return [1]; end if;
  gb2:=map(p->if indets(p) = {rag_sep_elem} then p fi, gb2):
  sqf:=Sqrfree(gb2[1]) mod fc;
  gb2:=mul(p[1], p in sqf[2]):
  gb2:=expand(gb2) mod fc;
  return [Groebner:-LeadingMonomial(gb2, tdeg(rag_sep_elem))];
end proc:


ModSatIntersect:=proc(eqs1, pol, eqs2, fc, vars, opts:={})
local var, gb1, gb2, i;
  var:=cat(x,0);
  i:=1:
  while member(var, vars) do 
    var:=cat(u,i);
    i:=i+1:
  end do;
  gb1:=MSolve:-MSolveGroebner([var*pol-1,op(eqs1)], fc, [var, op(vars)],
  {"elim"=1} union opts):
  gb1:=map(pol->if indets(pol) subset indets(vars) then pol fi, gb1):
  gb2:=MSolve:-MSolveGroebner([op(gb1), pol, op(eqs2)], fc, vars, opts):
  return gb2;
end proc:

ModSatIntersectLM:=proc(eqs1, pol, eqs2, fc, vars, opts:={})
local var, gb1, gb2, i;
  var:=cat(x,0);
  i:=1:
  while member(var, vars) do 
    var:=cat(u,i);
    i:=i+1:
  end do;
  gb1:=MSolve:-MSolveGroebner([var*pol-1,op(eqs1)], fc, [var, op(vars)],
  {"elim"=1} union opts):
  gb1:=map(pol->if indets(pol) subset indets(vars) then pol fi, gb1):
  gb2:=MSolve:-MSolveGroebnerLM([op(gb1), pol, op(eqs2)], fc, vars, opts):
  return gb2;
end proc:


#Does the same of SaturateIntersect but instead of returning a
#truncated grevlex GB of the solutions, it returns an elimination
#ideal
ElimSaturateIntersect:=proc(eqs1, pol, eqs2, vars, rag_sep_elem, opts:={})
local gb, prime, lmbg, N, lm, i, j, lifted, witness, boo2, boo1, newsupport,
lctables, sys, nprimes, witnessmod, boo, systable, primetable, support, modulus,
fc1, nthreads, str, fc, lmgb, newopts, rr, oldlifted, fcinit, oldwitness, islifted, 
newlifted, prevlifted;
  randomize():
  str:=subs(opts, "nthreads");
  if type(str, integer) then 
    nthreads:=str:
  else
    nthreads:=1:
  end if;
  newopts:={"linalg"=42, "nthreads"=nthreads}:
  rr:=rand(2^30..1303905300):
  fcinit:=nextprime(rr()):
  lm:=ElimModSatIntersectLM(eqs1, pol, eqs2, fcinit, vars,
  rag_sep_elem, newopts): 
  fc:=nextprime(2^30):
  while fc = fcinit do 
    fc:=nextprime(fc);
  end do;
  gb:=ElimModSatIntersect(eqs1, pol, eqs2, fc, vars, rag_sep_elem, newopts):
  lmgb:=map(pol->Groebner:-LeadingMonomial(pol, tdeg(op(vars),
  rag_sep_elem)), gb):
  while lmgb <> lm do 
    fc:=nextprime(fc):
    while fc=fcinit do 
      fc:=nextprime(fc);
    end do;
    gb:=ElimModSatIntersect(eqs1, pol, eqs2, fc, vars, rag_sep_elem, newopts):
    lmgb:=map(pol->Groebner:-LeadingMonomial(pol, tdeg(op(vars), rag_sep_elem)), gb):
  end do;
  fc1 := fc:
  sys:=[op(eqs1), op(eqs2)]:
  #N:=MinimalGenerators(gb, sys, vars, fc):
  N:=1:
  printf("[deg=%d]", max(map(degree, gb[1..N])));
  modulus:=fc:
  boo, support:=MonomialSupport(gb[1..N], [seq([],i=1..N)], [op(vars),
  rag_sep_elem]):
  boo:=true:
  primetable:=[fc]:
  systable:=[gb[1..N]]:
  lctables:=[seq([seq([],i=1..nops(support[j]))],j=1..N)]:
  lctables:=TableCoeffs(systable[1], [op(vars), rag_sep_elem], support, lctables):
  if nops(lctables) = 0 then 
    lctables:=[seq([seq([],i=1..nops(support[j]))],j=1..N)]:
  end if;
  witnessmod:=NewValuesWitness(gb[1..N], [op(vars), rag_sep_elem], support):
  nprimes:=1:
  printf("{%d}", nprimes);
  oldlifted:=[]:
  oldwitness:=[]:
  islifted:=0:#largest index of polys which are lifted
  lifted:=[]:
  while boo do 
    fc:=nextprime(fc);
    while fc = fcinit do 
      fc:=nextprime(fc);
    end do:
    gb:=ElimModSatIntersect(eqs1, pol, eqs2, fc, vars, rag_sep_elem, newopts):
    lmgb:=map(pol->Groebner:-LeadingMonomial(pol, tdeg(op(vars),rag_sep_elem)), gb):
    sys:=gb[1..N]:
    if lm=lmgb then 
      boo1, boo2, newsupport, systable, primetable, witnessmod, modulus, witness := 
        GeneratorsLift(sys, fc, [op(vars), rag_sep_elem], modulus, witnessmod, support, systable,
                      primetable):
      nprimes := nprimes + 1;
      printf("{%d}", nprimes);
      if boo1 = false then 
        printf("[!]");
        support := newsupport;
      end if;
      lctables:=TableCoeffs(sys, [op(vars), rag_sep_elem], support, lctables):
      if nops(lctables) = 0 then 
        lctables:=[seq([seq([],i=1..nops(support[j]))],j=1..N)]:
      end if;
      if boo1 = true and boo2 = true then 
        if oldwitness <> witness then 
          boo2:=false:
          oldwitness:=witness:
        else 
          oldwitness:=witness:
        end if;
      end if;
      if boo1 = true and boo2 = true then 
        printf("*");
        newlifted:=LiftPolynomials(lctables, support, primetable, modulus, islifted):
        islifted:=0:
        for i from 1 to min(nops(newlifted), nops(prevlifted)) do
          if newlifted[i]=prevlifted[i] and not(member(newlifted[i],
            lifted)) then 
            lifted:=[op(lifted), newlifted[i]]:
            printf("[%d]", nops(lifted));
            islifted:=islifted + 1;
          else 
            prevlifted:=newlifted;
          end if;
        end do;
        if nops(lifted) = N then 
          return map(numer, lifted);
        end if;
      end if;
    else
      lprint("Bad prime");
    end if;
  end do;
  return lifted;
end proc:

SaturateIntersect:=proc(eqs1, pol, eqs2, vars, opts:={})
local gb, prime, lmbg, N, lm, i, j, lifted, witness, boo2, boo1, newsupport,
lctables, sys, nprimes, witnessmod, boo, systable, primetable, support, modulus,
fc1, nthreads, str, fc, lmgb, newopts, rr, oldlifted, fcinit, oldwitness, islifted, 
newlifted, prevlifted;
  randomize():
  str:=subs(opts, "nthreads");
  if type(str, integer) then 
    nthreads:=str:
  else
    nthreads:=1:
  end if;
  newopts:={"linalg"=42, "nthreads"=nthreads}:
  rr:=rand(2^30..1303905300):
  fcinit:=nextprime(rr()):
  lm:=ModSatIntersectLM(eqs1, pol, eqs2, fcinit, vars, newopts): 
  fc:=nextprime(2^30):
  while fc = fcinit do 
    fc:=nextprime(fc);
  end do;
  gb:=ModSatIntersect(eqs1, pol, eqs2, fc, vars, newopts):
  lmgb:=map(pol->Groebner:-LeadingMonomial(pol, tdeg(op(vars))), gb):
  while lmgb <> lm do 
    fc:=nextprime(fc):
    while fc=fcinit do 
      fc:=nextprime(fc);
    end do;
    gb:=ModSatIntersect(eqs1, pol, eqs2, fc, vars, newopts):
    lmgb:=map(pol->Groebner:-LeadingMonomial(pol, tdeg(op(vars))), gb):
  end do;
  fc1 := fc:
  sys:=[op(eqs1), op(eqs2)]:
  N:=MinimalGenerators(gb, sys, vars, fc):
  printf("[ngens=%d,mdeg=%d]", N, max(map(degree, gb[1..N])));
  modulus:=fc:
  boo, support:=MonomialSupport(gb[1..N], [seq([],i=1..N)], vars):
  lprint("support", support);
  boo:=true:
  primetable:=[fc]:
  systable:=[gb[1..N]]:
  lctables:=[seq([seq([],i=1..nops(support[j]))],j=1..N)]:
  lctables:=TableCoeffs(systable[1], vars, support, lctables):
  if nops(lctables) = 0 then 
    lctables:=[seq([seq([],i=1..nops(support[j]))],j=1..N)]:
  end if;
  witnessmod:=NewValuesWitness(gb[1..N], vars, support):
  nprimes:=1:
  printf("{%d}", nprimes);
  oldlifted:=[]:
  oldwitness:=[]:
  islifted:=0:#largest index of polys which are lifted
  lifted:=[]:
  while boo do 
    fc:=nextprime(fc);
    while fc = fcinit do 
      fc:=nextprime(fc);
    end do:
    gb:=ModSatIntersect(eqs1, pol, eqs2, fc, vars, newopts):
    lmgb:=map(pol->Groebner:-LeadingMonomial(pol, tdeg(op(vars))), gb):
    sys:=gb[1..N]:
    if lm=lmgb then 
      boo1, boo2, newsupport, systable, primetable, witnessmod, modulus, witness := 
        GeneratorsLift(sys, fc, vars, modulus, witnessmod, support, systable,
                      primetable):
      nprimes := nprimes + 1;
      printf("{%d}", nprimes);
      if boo1 = false then 
        printf("[!]");
        support := newsupport;
      end if;
      lctables:=TableCoeffs(sys, vars, support, lctables):
      if nops(lctables) = 0 then 
        lctables:=[seq([seq([],i=1..nops(support[j]))],j=1..N)]:
      end if;
      if boo1 = true and boo2 = true then 
        if oldwitness <> witness then 
          boo2:=false:
          oldwitness:=witness:
        else 
          oldwitness:=witness:
        end if;
      end if;
      if boo1 = true and boo2 = true then 
        printf("*");
        newlifted:=LiftPolynomials(lctables, support, primetable, modulus, islifted):
        islifted:=0:
        for i from 1 to min(nops(newlifted), nops(prevlifted)) do
          if newlifted[i]=prevlifted[i] and not(member(newlifted[i],
            lifted)) then 
            lifted:=[op(lifted), newlifted[i]]:
            printf("[%d]", nops(lifted));
            islifted:=islifted + 1;
          else 
            prevlifted:=newlifted;
          end if;
        end do;
        if nops(lifted) = N then 
          return map(numer, lifted);
        end if;
      end if;
    else
      lprint("Bad prime");
    end if;
  end do;
  return lifted;
end proc:
