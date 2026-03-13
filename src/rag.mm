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
- make it a maple package
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

# The msolve library should installed, compiled and its binary should
# be accessible from your PATH
with(MSolve);

$include<multi-modular.mm>: 

ComputeMaximalMinors:=proc(M)
local i, rows, nr, nc, lc, l, minors, a, b, _l, _pol, _p;
  nr:=LinearAlgebra:-RowDimension(M):
  nc:=LinearAlgebra:-ColumnDimension(M):
  if nr > nc then
    error "Not implemented yet";
  end if;
  lc:=combinat:-choose(nc, nr):
  rows:=[seq(i, i=1..nr)]:
  minors:={seq(LinearAlgebra:-Determinant(LinearAlgebra:-SubMatrix(M,rows,l)),l in lc)}:
  minors:=convert(minors, list):
  minors:=map(_pol->if degree(_pol)>0 then 
                expand(mul(_p, _p in map(_l->_l[1],sqrfree(_pol)[2]))) else _pol
                fi, minors);
  return sort(minors, (a, b)-> degree(a) < degree(b));
end proc:

IsRegular:=proc(F, vars,opts:={})
local rr, hyp, _var, ld, gb, J;
  rr   := rand(1..65520):
  hyp  := add(rr()*_var, _var in vars)+rr():
  J    := convert(linalg:-jacobian(F, vars), Matrix):
  ld   := ComputeMaximalMinors(J):
  ld   := remove(member, ld, [0]);
  gb   := MSolveGroebnerLM([op(F), op(ld), hyp], 0, vars,opts):
  if gb=[1] then
    return true, ld;
  else
    return false, ld;
  end if;
end proc:

HaveFiniteIntersections:=proc(eqs, cstr, vars, opts:={})
local i, gb, rr, hyp, j;
  rr:=rand(1..2^30):
  for i from 1 to nops(cstr) do
    hyp:=add(rr()*vars[j],j=1..nops(vars)) + rr():
    gb:=MSolveGroebnerLM([op(eqs), cstr[i], hyp], 0, vars, opts):
    if gb<>[1] then 
      return false;
    end if;
  end do;
  return true;
end proc:

#To be improved; one should also return values that make positive the
#inequalities
GoodFiberValue_svars:=proc(svars, Inequalities, Inequations)
local i, j, boo, _pol, spos, sineq, spt;
  boo := true;
  j:=0:
  while boo do 
    spt := {seq(svars[i]=j,i=1..nops(svars))}:
    sineq := subs(spt, Inequations);
    spos := map(_pol->if degree(_pol) <= 0 then _pol fi,
                map(_pol->subs(spt, _pol), Inequalities));
    if not(member(0, sineq)) and not(member(0, spos)) then 
      return spt;
    end if;
    j:=j+1;
  end do;
end proc;

GoodFiberValue:=proc(var, hyp, Inequalities, Inequations)
local i, boo, pos, ineq, ls, hypsol;
  i := 0;
  boo := true;
  while boo do 
    hypsol:=solve(hyp=i, var);
    ls := {var = hypsol};
    pos:=subs(ls, Inequalities);
    ineq:=subs(ls, Inequations);
    if not(member(0, pos)) and not(member(0, ineq)) then 
      return ls, hypsol;
    end if;
    i := i+1;
  end do;
end proc;


DegreeTruncate:=proc(gb, fam, lm, vars, fc, mdeg, opts:={})
local deg, ldeg, p, lm2, i, tord, tgb, gb2, mgb, thresh, oldN, boo, tgb1, tgb2,
N;

  ldeg := sort(convert({op(map(degree, gb))}, list)):
#Minimum degree in gb 
  mgb := min(map(degree, gb)):
  tord := tdeg(op(vars)):
  thresh := min(mdeg + (mdeg - max(ldeg)) / 4, mgb + (mdeg - mgb) / 4):
  for i from 1 to nops(ldeg) do 
    deg := ldeg[i]:
    if deg >= thresh then 
      tgb := map(p->if degree(p) <= deg then p fi, gb):
      gb2 := MSolveGroebner([op(tgb), op(fam)], fc, vars,opts):
      oldN := nops(map(p->if degree(p) = deg then p fi, tgb)):
      lm2 := map(p->Groebner:-LeadingMonomial(p, tord), gb2):
      if lm2 = lm then 
        boo:=true:
        tgb1:=map(p->if degree(p) < deg then p fi, tgb):
        tgb2:=map(p->if degree(p) = deg then p fi, tgb):
        N:=nops(tgb2):
        while boo do 
          tgb:=[op(tgb1), op(tgb2[1..iquo(N, 2)])]:
          gb2 := MSolveGroebner([op(tgb), op(fam)], fc, vars,opts):
          lm2 := map(p->Groebner:-LeadingMonomial(p, tord), gb2):
          if lm2 = lm and N > 1 then
            oldN := N:
            N:=min(N, iquo(N, 2)):
          else 
            N:= oldN:
            boo:=false:
          end if;
        end do;
        return N+nops(tgb1);
      end if;
    end if;
  end do;
  return nops(gb);
end:

Increment:=proc(ol,ind)
local i, l;
  l:= ol;
  for i from 1 to nops(ind) do 
    l[ind[i]]:=l[ind[i]]+1;
  end do;
  return l;
end proc:

NextForms:=proc(ll)
local j, n, newll, i, lp;
  n := nops(ll[1]);
  newll := [];
  for i from 1 to n-1 do
    lp := combinat:-choose(n-1, i);
    for j from 1 to nops(lp) do 
      newll := [op(newll), op(map(l->Increment(l, lp[j]), ll))];
    end do;
  end do;
  return newll;
end proc:

# gendeg is empty when it has not been pre-determined (and hence needs to be
# computed)
TestGenericLineDegreeSingular:=proc(F, vars, hyp, singminors, gendeg,opts:={})
local tord, rr, gb2, fc, pol, J, minors, gb, hs, rag_hilb_var, deg, lF, ldeg, i, toadd, p, mdeg, isbounded, j, lhyp;

  if type(subs(opts, "isbounded"), integer) then 
    isbounded:=subs(opts, "isbounded");
  else 
    isbounded:=0:
  end if;
  if nops(F) < nops(vars) then 
    J      := convert(linalg:-jacobian([op(F), hyp], vars), Matrix):
    minors := ComputeMaximalMinors(J):
  else 
    minors := []:
  end if;
  ldeg   := []:
  toadd  := []:
  lF     := [seq(F[i]-F[1],i=2..nops(F))]:
  tord   := tdeg(op(vars));
  for i from 1 to nops(singminors) do 
    rr    := rand(2^30..1303905301):
    fc    := nextprime(rr()):

    pol   := singminors[i]:
    lhyp := [seq(randpoly(vars, degree = 2, dense),j=1..2)]:
    if nops(gendeg) = 0 or (nops(gendeg) > 0 and gendeg[i] > 0) then
      mdeg := max(map(degree, [op(lF), op(minors)])):
      gb    := MSolveGroebner([op(lhyp), 
                rag_sat_var*pol - 1, op(lF), op(minors), op(toadd)], 
                fc, 
                [rag_sat_var, op(vars)],opts union {"elim"=1}):
      if gb<>[1] then 
        deg:=-2: 
      else 
        gb    := MSolveGroebner([rag_sat_var*pol - 1, op(lF), op(minors), op(toadd)], 
               fc, 
               [rag_sat_var, op(vars)],opts union {"elim"=1}):
        gb    := map(p->if not(member(rag_sat_var, indets(p))) then p fi, gb):
        gb2   := MSolveGroebner([op(gb), F[1]], 
                 fc, vars, opts):
        gb2   := map(p->Groebner:-LeadingMonomial(p,tord), gb2);
        hs    := Groebner:-HilbertSeries(gb2, vars, rag_hilb_var);
        if degree(denom(hs)) = 0 then 
          deg   := subs(rag_hilb_var=1, hs):
        else 
          deg := -2;
        end if;
      end if;
    else 
      if gendeg[i] > 0 then 
        deg := gendeg[i]-1;
      else 
        deg := 0;
      end if;
    end if;
    ldeg  := [op(ldeg), deg]:
    toadd := [op(toadd), pol]:
  end do;

  return ldeg, minors;
end proc:

TestGenericLineDegreeRegular:=proc(F, vars, hyp,opts:={})
local J, minors, gb, hs, rag_hilb_var, deg;

  if nops(F) < nops(vars) then 
    J     := convert(linalg:-jacobian([op(F), hyp], vars), Matrix):
    minors:= ComputeMaximalMinors(J):
  else 
    minors := []:
  end if;
  gb    := MSolveGroebnerLM([op(F),op(minors)], 0, vars,opts):
  hs    := Groebner:-HilbertSeries(gb, vars, rag_hilb_var);
  if degree(denom(hs)) = 0 then 
    deg   := subs(rag_hilb_var=1, hs):
  else 
    deg   := -2;
  end if;
  return deg, minors;
end proc:

HasFiniteCriticalLocus:=proc(eqs, F, minors, singminors, vars, opts:={})
local gb, i, v, hyp, pol, rr;
  rr:=rand(1..nextprime(2^30)):
  hyp:=add(rr()*v, v in vars)+rr():
  for i from 1 to nops(singminors) do
    pol:=singminors[i]:
    gb:=MSolveGroebnerLM([rag_sat_var*pol-1, op(eqs), op(F), hyp,
                          op(minors), op(singminors[1..i-1])], 
                          0, [rag_sat_var, op(vars)], opts union
                          {"linalg"=42}):
    if gb<>[1] then 
      return false;
    end if;
  end do;
  return true;
end proc;

FindGenericLineRegular:=proc(eqs, F, lc, singminors, vars,opts:={})
local rr, hyp, J, B, n, v, i, j, minors, deg, gendeg, ll, isbounded,
newll, dF, verb:

  if type(subs(opts, "isbounded"), integer) then 
    isbounded:=subs(opts, "isbounded");
  else 
    isbounded:=0:
  end if;
  if type(subs(opts, "verb"), integer) then 
    verb:=subs(opts, "verb");
  else 
    verb:=0:
  end if;
  dF:=[seq(F[i]+lc[i],i=1..nops(F))]:
################################
#Generic degree
  rr:=rand(1..nextprime(2^30)):
  hyp:=add(rr()*v, v in vars):
  gendeg, minors := TestGenericLineDegreeRegular([op(eqs),op(dF)], vars, hyp, opts); 
  if gendeg = -2 then 
    lprint(eqs, F, lc, singminors, vars);
    error "Generic line should provide finitely many critical points" ; 
  end if;
################################

  for i from 1 to nops(vars) do 
    if verb >= 1 then 
      printf("+");
    end if;
    if not(member(vars[i], F)) then 
      hyp        := vars[i]:
      deg, minors := TestGenericLineDegreeRegular([op(eqs),op(dF)], vars, hyp, opts):
      if deg=gendeg or (deg >= 0 and isbounded > 0) then
        if HasFiniteCriticalLocus(eqs, F, minors, singminors, vars,
          opts) then 
          return hyp, minors, gendeg;
        end if;
      fi;
    end if;
  od;

  B:=1:
  n:=nops(vars):
  ll:=[[seq(1, i=1..n)]]:
  while true do 
    for i from 1 to nops(ll) do 
      if verb>= 1 then 
        printf("+");
      end if;
      hyp         := add(ll[i][j]*vars[j], j=1..n):
      deg, minors := TestGenericLineDegreeRegular([op(eqs),op(dF)], vars, hyp, opts):
      if deg=gendeg or (deg >= 0 and isbounded > 0)  then
        if HasFiniteCriticalLocus(eqs, F, minors, singminors, vars,
          opts) then 
          return hyp, minors, gendeg;
        end if;
      quit;
      end if;
    end do:
    newll := NextForms(ll);
    ll:=remove(member, newll, ll):
  od:

end proc:

FindGenericLineSingular:=proc(F, vars, singminors, opts:={})
local rr, hyp, J, B, n, _var, i, j, minors, deg, gendeg, ll, isbounded, newll:

  if type(subs(opts, "isbounded"), integer) then 
    isbounded:=subs(opts, "isbounded");
  else 
    isbounded:=0:
  end if;
################################
#Generic degree
  rr:=rand(1..65520):
  hyp:=add(rr()*_var, _var in vars):
  if isbounded = 0 then 
  gendeg, minors := TestGenericLineDegreeSingular(F, vars, hyp, singminors, [], opts); 
  if member(-2, gendeg) then
    error"Generic line should give finitely many critical points";
  end if;
  else 
    gendeg:=[seq(1, i=1..nops(singminors))]:
  end if;
################################
#This is necessary because when F contains algebraically dependent polynomials
#the next loop will return the first variable (see example 
#  F = [x1,x2,x1-x2],vars=[x1,x2,x3])

  for i from 1 to nops(vars) do 
    printf("+");
    if not(member(vars[i], F)) then 
      hyp         := vars[i]:
      deg, minors := TestGenericLineDegreeSingular(F, vars, hyp, singminors,
                     gendeg, opts):
      if deg=gendeg or (nops(convert(map(sign, deg),set))=1 and isbounded > 0) then
        return hyp, minors, gendeg;
      fi;
    end if;
  od;
  B:=1:
  n:=nops(vars):
  ll:=[[seq(1, i=1..n)]]:
  while true do 
    for i from 1 to nops(ll) do 
      printf("+");
      hyp         := add(ll[i][j]*vars[j], j=1..n):
      deg, minors := TestGenericLineDegreeSingular(F, vars, hyp, singminors,
      gendeg, opts):
      if deg=gendeg or (nops(convert(map(sign, deg),set))=1 and isbounded > 0) then
        return hyp, minors, gendeg;
      fi;
    od:
    newll := NextForms(ll);
    ll:=remove(member, newll, ll):
  end do:

end proc:


CoeffDeform_eps:=proc(eqs, F, vars, eps, cstr, opts:={})
local rr, gb, hs, deg, rag_hilb_var, i, gb0, dim, newlc, lc, k, allvars:
  rr:=rand(1..2^30):
  allvars:=[op(vars), eps]:
  gb0:=MSolveGroebnerLM([op(eqs), seq(F[i]+rr()*eps, i=1..nops(F))], 0, allvars, 
            opts union {"linalg"=42}):
  hs    := Groebner:-HilbertSeries(gb0, allvars, rag_hilb_var);
  deg   := subs(rag_hilb_var=1, numer(hs)):
  dim   := degree(denom(hs));
  lc:=[[seq(1, i=1..nops(F))]]:
  while true do 
    for k from 1 to nops(lc) do 
      gb:=MSolveGroebnerLM([op(eqs), 
              seq(F[i]+lc[k][i]*eps, i=1..nops(F))], 0, allvars, 
              opts union {"linalg"=42}):
      hs    := Groebner:-HilbertSeries(gb, allvars, rag_hilb_var);
      if gb=gb0 then 
        if HaveFiniteIntersections([op(eqs), 
                                    seq(F[i]+lc[k][i]*eps, i=1..nops(F))], 
                                    cstr, [op(vars), eps], opts) then 
          return lc[k];
        end if;
      end if;
    end do;
    newlc:=NextForms(lc):
    lc:=remove(member, newlc, lc):
  end do;
end proc:


CoeffDeform:=proc(eqs, F, singminors, vars, opts:={})
local rr, gb, hs, deg, rag_hilb_var, i, lhyp, hyp, gbsing, newlc, lc, k:
  rr:=rand(1..2^30):
  lhyp:=[seq(add(rr()*vars[i],i=1..nops(vars)),
        nops(vars)-nops(F)-nops(eqs))]:
  gb:=MSolveGroebnerLM([op(eqs), seq(F[i]+rr(), i=1..nops(F)), op(lhyp)], 0, vars, 
            opts union {"linalg"=42}):
  hs    := Groebner:-HilbertSeries(gb, vars, rag_hilb_var);
  if degree(denom(hs)) = 0 then 
    deg   := subs(rag_hilb_var=1, hs):
  else 
    lprint(args);
    error "Bug in CoeffDeform";
  end if;
  lc:=[[seq(1, i=1..nops(F))]]:
  hyp:=add(rr()*vars[i],i=1..nops(vars))+rr():
  while true do 
    for k from 1 to nops(lc) do 
      gb:=MSolveGroebnerLM([op(eqs), 
              seq(F[i]+lc[k][i], i=1..nops(F)), op(lhyp)], 0, vars, 
              opts union {"linalg"=42}):
      gbsing:=MSolveGroebnerLM([op(eqs), op(singminors), 
              seq(F[i]+lc[k][i], i=1..nops(F)), hyp], 0, vars, 
              opts union {"linalg"=42}):
      hs    := Groebner:-HilbertSeries(gb, vars, rag_hilb_var);
      if degree(denom(hs))=0 and subs(rag_hilb_var=1, hs) = deg and
        gbsing = [1] then 
        return lc[k];
      end if;
    end do;
    newlc:=NextForms(lc):
    lc:=remove(member, newlc, lc):
  end do;
end proc:

FindGenericLine:=proc(eqs, F, vars, opts:={})
local boo, singminors, verb, gendeg, hyp, minors, i, lc;
  
  if type(subs(opts, "verb"), integer) then 
    verb:=subs(opts, "verb");
  else 
    verb:=0:
  end if;
  boo, singminors := IsRegular([op(eqs), op(F)], vars, opts):
  
  lc:=CoeffDeform(eqs, F, singminors, vars, opts):
  hyp, minors, gendeg := FindGenericLineRegular(eqs, F, lc,
                         singminors, vars, opts);
  if boo = true then
    if verb >= 1 then printf("[R]"); end if;
    return hyp, minors, gendeg, true, lc;
  else 
    if verb >= 1 then printf("[S]"); end if;
    return hyp, [minors, singminors], gendeg, false, lc;
  end if;
end proc:

SmallMidRational:=proc(rr1, rr2)
local c, mid, i;
  if rr1 >= rr2 then 
    lprint(args);
    error "Bug detected";
  end if;
  mid:=(rr1+rr2)/2:
  numtheory:-cfrac(mid,20,'c');
  for i from 1 to nops(c)-1 do 
    if rr1 < c[i] and c[i] < rr2 then 
      return c[i];
    end if;
  end do;
  return mid;
end proc;

ConstructFibers:=proc(ll, hyp, cstr)
local i, mid, res, vvar, ls, val;
#decides if some fiber cancels constraints it should not
#in that case, another fiber should be chosen
#the way it is done is not that optimal but the overhead should be
#negligeible (and that situation is rather rare)
#smallest fiber 
  vvar:=indets(hyp)[1]:
  val:=floor(ll[1][1])-1:
  ls:={vvar=solve(hyp-val)}:
  if member(0, subs(ls, cstr)) then
    return ConstructFibers([[val, val], op(ll)], hyp, cstr);
  end if;
#largest fiber 
  val:=ceil(ll[nops(ll)][2])+1:
  ls:={vvar=solve(hyp-val)}:
  if member(0, subs(ls, cstr)) then
    return ConstructFibers([op(ll), [val, val]], hyp, cstr);
  end if;
#remaining ones 
  for i from 1 to nops(ll)-1 do 
    val:=(ll[i][2]+ll[i+1][1])/2:
    ls:={vvar=solve(hyp-val)}:
    if member(0, subs(ls, cstr)) then
      return ConstructFibers(sort([[val, val], op(ll)], (a, b)->a[2] <
      b[1]), hyp, cstr);
    end if;
  end do;

#now, we can construct the fibers
  res:=[floor(ll[1][1]-1)]:
  for i from 1 to nops(ll)-1 do 
    mid:=SmallMidRational(ll[i][2], ll[i+1][1]):
    res:=[op(res), mid]:
  end do;
  res:=[op(res), ceil(ll[nops(ll)][2]+1)]:
end proc;

HasOverLapCoupleOfIntervals:=proc(l1, l2)
  return (evalb(l1[2] >= l2[1]))
end proc;

#Assumes that the list of intervals _list has been sorted
HasOverLap:=proc(_list)
local i, boo;
  if nops(_list)<=1 then
     return false;
  end if;
  for i from 1 to nops(_list)-1 do
    boo:=HasOverLapCoupleOfIntervals(_list[i], _list[i+1]);
    if boo=true then
       return boo;
    end if;
  end do;
  return false; 
end proc;

ComputeBoundsRegular:=proc(Equations, Fam, Positive, NotNull, vars, 
                    hyp, minors, gendeg, opts:={})
local vvar, ls, rd, sols, rr, lF, lhyp, gb, i, j, verb;

    if type(subs(opts, "verb"), integer) then 
      verb:=subs(opts, "verb");
    else 
      verb:=0:
    end if;

    rd:=rand(1..65521):
    if gendeg = 0 then 
      sols := [0, []]:
    else 
      if nops(Equations) + nops(Fam) < nops(vars) then 
        sols := MSolveRealRoots([op(Equations), op(Fam), op(minors),
                hyp-rag_sep_elem], [op(vars), rag_sep_elem],
                [op(Positive), op(NotNull)], opts);
      else 
        sols := MSolveRealRoots([op(Equations), op(Fam), op(minors),
                hyp-rag_sep_elem], [op(vars), rag_sep_elem],
                [], opts);
      end if;
      if sols[1] >0 then 
        lprint(args);
        error "Degenerate case in ComputeBounds: to be implemented (bounds)";
      fi;
      sols := [0, AdmissibleSolutions(sols, nops(Positive))];
    end if;
    if nops(sols[2]) > 0 then 
      rr := convert(map(s->subs(s, rag_sep_elem), sols[2]), set);
      rr := sort(convert(rr, list), (a, b)->a[2] < b[1]);
      if HasOverLap(rr) = false then 
        rr := ConstructFibers(rr, hyp, [op(Positive), op(NotNull)]);
      else 
        if verb>= 1 then printf("[Overlap Regular]"); end if;
        gb := MSolveGroebner([op(Equations), op(Fam), op(minors),
              hyp-rag_sep_elem], 0, [op(vars), rag_sep_elem],
              opts union
              {"elim"=nops(vars)});
        sols := MSolveRealRoots(gb, [rag_sep_elem], []):
        rr := convert(map(s->subs(s, rag_sep_elem), sols[2]), set);
        rr := sort(convert(rr, list), (a, b)->a[2] < b[1]);
        rr := ConstructFibers(rr, hyp, [op(Positive), op(NotNull)]);
      end if;
      return hyp, rr;
    else
      i:=0:
      vvar:=indets(hyp)[1]:
      ls:={vvar=solve(hyp=i,vvar)}:
      lF:=[op(Equations),op(map(_p->_p+rd(),Fam)),hyp-i]:
      lhyp:=[seq(randpoly(vars,degree=1,dense),j=1..nops(vars)-nops(lF)+1)]:
      gb:=MSolve:-MSolveGroebnerLM([op(lhyp),op(lF)],0,vars,opts):
      while gb<>[1] or member(0, subs(ls, [op(Positive),
        op(NotNull), op(Fam)])) do
        i:=i+1:
        lF:=[op(Equations),op(map(_p->_p+rd(),Fam)),hyp-i]:
        lhyp:=[seq(randpoly(vars,degree=1,dense),j=1..nops(vars)-nops(lF)+1)]:
        gb:=MSolve:-MSolveGroebnerLM([op(lhyp),op(lF)],0,vars,opts):
        ls:={vvar=solve(hyp=i,vvar)}:
      od:
      return hyp, [i];
    end if;
end proc:

ManageOverLapComputeBoundsSingular:=proc(Equations, Fam, singminors, minors, 
              gendeg, lgb, hyp, Positive, NotNull, vars, opts:={})
local i, toadd, upol, _l, squpol, rr, sols, pol, gb1, gb2;
  toadd:=[]:
  upol:=1:
  for i from 1 to nops(singminors) do 
      pol := singminors[i]:
      if gendeg[i] = 0 then 
      else
        gb1:=MSolveGroebner([op(lgb[i]), op(Equations), op(minors),
                                op(Fam), op(toadd),
                                hyp-rag_sep_elem], 0, [op(vars),
                                rag_sep_elem], opts union
                                {"elim"=nops(vars)}):
        upol:=upol*gb1[1]:
        gb2:=MSolveGroebner([rag_sat_var*pol-1, op(Equations), op(minors),
                                op(Fam), op(toadd), hyp-rag_sep_elem], 0,
                                [rag_sat_var, op(vars), rag_sep_elem], 
                                opts union {"elim"=nops(vars)+1}):
        upol:=upol*gb2[1]:
      end if;
      toadd := [op(toadd), pol]:
  end do;
  squpol:=map(_l->_l[1], sqrfree(upol)[2]):
  if nops(squpol) = 0 then 
    lprint(args);
    error "Bug detected (ManageOverLapComputeBoundsSingular)";
  end if;
  upol:=mul(_l, _l in squpol):
  sols:=MSolveRealRoots([upol], [rag_sep_elem], []):
  sols:=sols[2]:
  rr := convert(map(s->subs(s, rag_sep_elem), sols), set);
  rr := sort(convert(rr, list), (a, b)->a[2] < b[1]);
  rr := ConstructFibers(rr, hyp, [op(Positive), op(NotNull)]);
  return rr;
end proc;

ElimComputeBoundsSingular:=proc(Equations, Fam, Positive, NotNull, vars, 
                    hyp, vminors, gendeg, opts:={})
local rd, minors, singminors, i, pol, gb, nsols, OldDigits, toadd, lgb, sols, 
lF, nsols2, rr, j, lhyp, k, vvar, ls, verb;

    if type(subs(opts, "verb"), integer) then 
      verb:=subs(opts, "verb");
    else 
      verb:=0:
    end if;

    rd:=rand(1..65521):
    singminors := vminors[2];

    minors     := vminors[1]:
    sols       := []:
    toadd      := []:
    lF         := [seq(Fam[i]-Fam[1], i=2..nops(Fam))]:
    lgb := []:
    for i from 1 to nops(singminors) do
      pol := singminors[i]:
      if gendeg[i] = 0 then 
        nsols := [0, []]:
        lgb:=[op(lgb), [1]]:
      else
        OldDigits:=Digits:
        Digits:=max(10,
        max(map(degree,map(expand,[op(Fam),op(Equations)])))+iquo(max(map(ilog2,[seq(coeffs(i), i in
        map(expand, [op(Fam), op(Equations)]))])),2));

        gb:=ElimSaturateIntersect([op(Equations), op(minors), op(lF), 
                               op(toadd)], 
                                pol, 
                                [pol, Fam[1], op(singminors), hyp-rag_sep_elem], 
                                vars, rag_sep_elem, opts):
        lgb:=[op(lgb), gb]:
        nsols := [0, []];
        nsols2:=MSolveRealRoots([rag_sat_var*pol-1, op(Equations), op(minors),
                                op(Fam), op(toadd), hyp-rag_sep_elem], 
                                [rag_sat_var, op(vars), rag_sep_elem], 
                                [], opts):
        if nops(nsols2) < 2 then 
          lprint(args);
          error "nsols2 should have cardinality 2 (2)";
        else 
          nsols2 := [0, AdmissibleSolutions(nsols2, nops(Positive))];
          nsols2:=[0, map(s->map(c->if indets(c) subset indets([op(vars),
                        rag_sep_elem]) then c
                        fi,s), nsols2[2])]:
        end if;
        Digits:=OldDigits:
        nsols:=[0, [op(nsols[2]), op(nsols2[2])]]:
      end if;
      if nops(nsols) < 2 then
        lprint(args);
        error "nsols should have cardinality 2";
      end if;
      toadd := [op(toadd), pol]:
      sols:=[op(sols), op(nsols[2])];
    end do;
    if nops(sols) > 0 then
      rr := convert(map(s->subs(s, rag_sep_elem), sols), set);
      rr := sort(convert(rr, list), (a, b)->a[2] < b[1]);
      if HasOverLap(rr) = false then 
        rr := ConstructFibers(rr, hyp, [op(Positive), op(NotNull)]);
      else 
        if verb>=1 then printf("Overlap Singular"); end if;
        rr := ManageOverLapComputeBoundsSingular(Equations, Fam, singminors, minors, 
              gendeg, lgb, hyp, Positive, NotNull, vars, opts): 
      end if;
      return hyp, rr;
    else
      j:=0:
      vvar:=indets(hyp)[1]:
      if max(gendeg) = 0 then 
        ls:={vvar=solve(hyp=j,vvar)}:
        while member(0,subs(ls, [op(Fam),op(Positive),op(NotNull)]))
          do 
            j:=j+1:
            ls:={vvar=solve(hyp=j,vvar)}:
          end do;
        return hyp, [j];
      end if;
      lF:=[op(Equations),op(map(_p->_p+rd(),Fam)),hyp-j]:
      lhyp:=[seq(randpoly(vars,degree=1,dense),k=1..nops(vars)-nops(lF)+1)]:
      gb:=MSolve:-MSolveGroebnerLM([op(lhyp),op(lF)],0,vars,opts):
      ls:={vvar=solve(hyp=j,vvar)}:
      while gb<>[1] or member(0, subs(ls, [op(Positive),
        op(NotNull), op(Fam)])) do
        j:=j+1:
        lF:=[op(Equations),op(map(_p->_p+rd(),Fam)),hyp-j]:
        lhyp:=[seq(randpoly(vars,degree=1,dense),k=1..nops(vars)-nops(lF)+1)]:
        gb:=MSolve:-MSolveGroebnerLM([op(lhyp),op(lF)],0,vars,opts):
        ls:={vvar=solve(hyp=j,vvar)}:
      od:
      return hyp, [j];
  end if;
end proc;

ComputeBoundsSingular:=proc(Equations, Fam, Positive, NotNull, vars, 
                    hyp, vminors, gendeg, lc, opts:={})
local rd, minors, singminors, i, pol, gb, nsols, OldDigits, toadd, lgb, sols, 
lF, nsols2, rr, j, lhyp, k, vvar, ls;
    rd:=rand(1..65521):
    minors     := vminors[1]:
    singminors := remove(member, vminors[2], [0]);
    singminors := remove(member, singminors, minors):

    sols       := []:
    toadd      := []:
    lF         := [seq(lc[1]*Fam[i]-lc[i]*Fam[1], i=2..nops(Fam))]:
    lgb := []:
    for i from 1 to nops(singminors) do
      pol := singminors[i]:
      if gendeg[i] = 0 then 
        nsols := [0, []]:
        lgb:=[op(lgb), [1]]:
      else
        OldDigits:=Digits:
        Digits:=max(10,
        max(map(degree,map(expand,[op(Fam),op(Equations)])))+iquo(max(map(ilog2,[seq(coeffs(i), i in
        map(expand, [op(Fam), op(Equations)]))])),2));

        gb:=SaturateIntersect([op(Equations), op(minors), op(lF), op(toadd)], 
                                pol, 
                                [pol, Fam[1]], vars, opts):
        lgb:=[op(lgb), gb]:
        if nops(Equations) + nops(Fam) < nops(vars) then 
          nsols:=MSolveRealRoots([op(gb), op(Equations), op(minors),
                                  op(Fam), op(toadd),
                                  hyp-rag_sep_elem], [op(vars),
                                  rag_sep_elem], [op(Positive),
                                  op(NotNull)], opts):
        else 
          nsols:=MSolveRealRoots([op(gb), op(Equations), op(minors),
                                  op(Fam), op(toadd),
                                  hyp-rag_sep_elem], [op(vars),
                                  rag_sep_elem], [], opts):
        end if;
        if nops(nsols) < 2 then 
          lprint(args);
          error "nsols should have cardinality 2 (1)";
        end if;
        nsols := [0, AdmissibleSolutions(nsols, nops(Positive))];
        if nops(Equations) + nops(Fam) < nops(vars) then 
          nsols2:=MSolveRealRoots([rag_sat_var*pol-1, op(Equations), op(minors),
                                op(Fam), op(toadd), hyp-rag_sep_elem], 
                                [rag_sat_var, op(vars), rag_sep_elem], 
                                [op(Positive), op(NotNull)], opts):
        else 
          nsols2:=MSolveRealRoots([rag_sat_var*pol-1, op(Equations), op(minors),
                                op(Fam), op(toadd), hyp-rag_sep_elem], 
                                [rag_sat_var, op(vars), rag_sep_elem], 
                                [], opts):
        end if;
        if nops(nsols2) < 2 then 
          #lprint(args);
          error "nsols2 should have cardinality 2 (2)";
        else 
          nsols2 := [0, AdmissibleSolutions(nsols2, nops(Positive))];
          nsols2:=[0, map(s->map(c->if indets(c) subset indets([op(vars),
                        rag_sep_elem]) then c
                        fi,s), nsols2[2])]:
        end if;
        Digits:=OldDigits:
        nsols:=[0, [op(nsols[2]), op(nsols2[2])]]:
      end if;
      if nops(nsols) < 2 then
        lprint(args);
        error "nsols should have cardinality 2";
      end if;
      toadd := [op(toadd), pol]:
      sols:=[op(sols), op(nsols[2])];
    end do;
    if nops(sols) > 0 then
      rr := convert(map(s->subs(s, rag_sep_elem), sols), set);
      rr := sort(convert(rr, list), (a, b)->a[2] < b[1]);
      if HasOverLap(rr) = false then 
        rr := ConstructFibers(rr, hyp, [op(Positive), op(NotNull)]);
      else 
        if verb >= 1 then printf("Overlap Singular"); end if;
        rr := ManageOverLapComputeBoundsSingular(Equations, Fam, singminors, minors, 
              gendeg, lgb, hyp, Positive, NotNull, vars, opts): 
      end if;
      return hyp, rr;
    else
      j:=0:
      vvar:=indets(hyp)[1]:
      if max(gendeg) = 0 then 
        ls:={vvar=solve(hyp=j,vvar)}:
        while member(0,subs(ls, [op(Fam),op(Positive),op(NotNull)]))
          do 
            j:=j+1:
            ls:={vvar=solve(hyp=j,vvar)}:
          end do;
        return hyp, [j];
      end if;
      lF:=[op(Equations),op(map(_p->_p+rd(),Fam)),hyp-j]:
      lhyp:=[seq(randpoly(vars,degree=1,dense),k=1..nops(vars)-nops(lF)+1)]:
      gb:=MSolve:-MSolveGroebnerLM([op(lhyp),op(lF)],0,vars,opts):
      ls:={vvar=solve(hyp=j,vvar)}:
      while gb<>[1] or member(0, subs(ls, [op(Positive),
        op(NotNull), op(Fam)])) do
        j:=j+1:
        lF:=[op(Equations),op(map(_p->_p+rd(),Fam)),hyp-j]:
        lhyp:=[seq(randpoly(vars,degree=1,dense),k=1..nops(vars)-nops(lF)+1)]:
        gb:=MSolve:-MSolveGroebnerLM([op(lhyp),op(lF)],0,vars,opts):
        ls:={vvar=solve(hyp=j,vvar)}:
      od:
      return hyp, [j];
  end if;
end proc;

ComputeBounds:=proc(Equations, Fam, Positive, NotNull, vars, opts:={})
local boo, i, pol, nsols, gb, hyp, rr, s, sols,j, k, gendeg, singminors, 
sysminors, nsols2, toadd, minors, OldDigits, lF, lhyp, lc;

  hyp, minors, gendeg, boo, lc := FindGenericLine(Equations, Fam, vars, opts):
#Regular case 
  if boo then 
    return ComputeBoundsRegular(Equations, Fam, Positive, NotNull, vars, 
                    hyp, minors, gendeg, opts);
  else
    return ComputeBoundsSingular(Equations, Fam, Positive, NotNull, vars, 
                    hyp, minors, gendeg, lc, opts);
    end if;
end proc:

ModularLimits:=proc(Equations, pol, vars, eps, opts:={})
local p, rrfc, gb, fc, mdeg, gb2, tord, trunc_deg;
  rrfc := rand(2^30..1303905301):
  fc   := nextprime(rrfc()):
  mdeg := max(map(degree, [pol, op(Equations)])):
  gb := MSolveGroebner([rag_sat_var * pol - 1, op(Equations)], fc, 
        [rag_sat_var, op(vars)], opts union {"elim"=1});
  gb := map(p->if not(member(rag_sat_var, indets(p))) then p fi, gb):
  gb2 := MSolveGroebner([eps, op(gb)], fc, vars, opts):
  if gb2=[1] then return false, 0;
  else
    tord := tdeg(op(vars)):
    gb2 := map(p->Groebner:-LeadingMonomial(p, tord), gb2):
    trunc_deg := DegreeTruncate(gb, [eps], gb2, vars, fc, mdeg, opts):
    return true, trunc_deg;
  end if;
end proc:

LimitsDeformedCriticalPoints:=proc(Equations, Fam, Inequalities,
    Inequations, pol1, pol2, vars, opts:={})
local v, dFam, lf, pol, boo, N, sols, nsols, hyp, toadd, rr, i, 
J, JS, minors, sminors, rag_sat_var, gb, a, b, hyp1, hyp2, J1, J2, lminors, j;
  rr   := rand(1..65520):
  lf   := [seq(add(rr()*v, v in vars), i=1..nops(Fam))];
  dFam := [seq(Fam[i]+rag_eps_var*lf[i],i=1..nops(Fam))];
  hyp1  := add(rr()*v, v in vars):
  hyp2  := add(rr()*v, v in vars):

  J1   := convert(linalg:-jacobian([op(Equations), op(dFam), hyp1], vars), Matrix):
  J2   := convert(linalg:-jacobian([op(Equations), op(dFam), hyp2], vars), Matrix):
  JS    := convert(linalg:-jacobian([op(Equations), op(dFam)], vars), Matrix):

  sminors := ComputeMaximalMinors(JS);
  if nops(Equations)+nops(Fam) < nops(vars) then 
  minors  := [op({op(ComputeMaximalMinors(J1)), op(ComputeMaximalMinors(J2))})]:
  else 
  minors:= [];
  end if;
  sminors := sort(remove(member, sminors, minors),(a,b)->degree(a)<=degree(b)):

  toadd := []:
  sols  := []:
  lminors:=SplitSystem(minors):
  for i from 1 to nops(sminors) do
    pol   := sminors[i];
    for j from 1 to nops(lminors) do 
      boo, N := ModularLimits([op(Equations), op(dFam), op(toadd), op(lminors[j])],
                pol*rag_eps_var, [rag_eps_var, op(vars)], rag_eps_var, opts);
      if boo = true then
        gb    := MSolveGroebner([op(Equations), op(dFam), 
                             rag_sat_var * pol*rag_eps_var - 1, 
                             op(lminors[j]), op(toadd)], 0, [rag_sat_var, rag_eps_var, op(vars)], opts union {"trunc"=N, "elim"=1});
        nsols := MSolveRealRoots( [rag_sat_var * Fam[1] - 1, op(gb), pol, op(toadd)], 
                          [rag_sat_var, op(vars), rag_eps_var],
                          [pol1, op(Inequalities), pol2, op(Inequations)], opts);
        if nops(nsols) < 2 then 
          lprint(nsols);
          error "nsols should have cardinality 2";
        end if;
        nsols:=AdmissibleSolutions(nsols, nops(Inequalities)+1);
        lprint(evalf(nsols));
        nsols := map(sol->map(c-> if member(lhs(c), vars) then c fi, sol),nsols);
      else 
        nsols := []:
      end if;
      toadd :=[op(toadd), pol]:
      sols  :=[op(sols), op(nsols)];

    end do;
  end do;
  return sols;
end proc:

CriticalPointsSingular:=proc(Equations, Fam, Inequalities, Inequations, pol,
                             pol2, vars, minors, opts:={})
local i, toremove, sols, positive, nnull, cstr, j, gb, hs, np, rag_sat_var1, rag_sat_var2;

#One first tries to remove from the set of singular points vanishing loci of
#constraints to see if we get finitely many points (turned out to be useful in
#many applications).
  cstr :=[op(Inequalities), op(Inequations)]:
  toremove:=[]:
  for i from 1 to nops(cstr) do 
    gb := MSolve:-MSolveGroebnerLM([op(Equations),
                  seq(Fam[j]-Fam[1],j=1..nops(Fam)), op(minors),
                  rag_sat_var1*Fam[1]-1, rag_sat_var2*cstr[1]-1], 0, 
                  [rag_sat_var2, rag_sat_var1, op(vars)], {"elim"=2} union opts):
    hs:=Groebner:-HilbertSeries(gb, vars, rag_hilb_var):
    if degree(denom(hs))=0 then 
      toremove:=[cstr[i]]:
      break;
    end if;
  end do;
  if nops(toremove) > 0 then 
  sols:=MSolveRealRoots([op(Equations), seq(Fam[i]-Fam[1],i=2..nops(Fam)), op(minors), 
                         rag_sat_var1*Fam[1]-1, rag_sat_var2*toremove[1]-1], 
                        [rag_sat_var2, rag_sat_var1, op(vars)], 
                        [pol, op(Inequalities), pol2, op(Inequations)], opts):
  sols:=AdmissibleSolutions(sols, nops(Inequalities) + 1);
  sols := map(_p->map(_c->if member(lhs(_c), vars) then _c else fi, _p), sols);
  return sols;
  end if;
  return LimitsDeformedCriticalPoints(Equations, Fam, Inequalities,
                                      Inequations, pol, pol2, vars, opts);

end proc:

SplitSystem_cstr:=proc(F, cstr)
local i, pol, _l, lsys, sys, newlsys, lf, _p, vars, degone;
  if nops(F)=0 then return F; end if;
  if nops(F)=1 then 
    pol:=F[1]:
    if degree(pol)<=0 then return [[pol]]:
    else 
      lf:=map(_l->_l[1], factors(pol)[2]):
      lf:=map(_l->if not(member(true, map(_p->divide(_p, _l), cstr))) then [_l] fi,lf):
      return lf;
    end if;
  end if;
  pol:=F[1]:
  if degree(pol)>0 then 
    lf:=map(_l->primpart(_l[1]),factors(pol)[2]);
    lf:=map(_l-> if not(member(true, map(_p->divide(_p, _l), cstr))) then
    _l fi, lf):
  else 
    lf:=[pol];
  end if;
  lsys:=convert(convert(map(_l->convert(_l, set), SplitSystem(F[2..-1])), set),list);
  lsys:=map(_l->op(map(_p->[op(_l), _p], lf)), lsys);
  vars:=[op(indets([op(F), op(cstr)]))]:
  newlsys:=[]:
  for i from 1 to nops(lsys) do 
    sys:=lsys[i]:
    degone:=map(_p->if degree(_p)=1 then _p fi, sys):
    if not(member(0, map(_p->Groebner:-NormalForm(_p, degone, tdeg(op(vars))),
      cstr))) then 
      newlsys:=[op(newlsys), sys]:
    end if;
  end do:
  return newlsys;
end:

SplitSystem:=proc(F)
local i, pol, _l, lsys, lf;
  if nops(F)=0 then return F; end if;
  if nops(F)=1 then 
    pol:=F[1]:
    if degree(pol)<=0 then return [[pol]]:
    else 
      lf:=map(_l->[_l[1]],factors(pol)[2]):
      return lf;
    end if;
  end if;
  pol:=F[1]:
  if degree(pol)>0 then 
    lf:=map(_l->primpart(_l[1]),factors(pol)[2]);
  else 
    lf:=[pol];
  end if;
  lsys:=convert(convert(map(_l->convert(_l, set), SplitSystem(F[2..-1])), set),list);
  return map(_l->op(map(_p-> if not(member(_p, _l)) then [op(_l), _p] fi, lf)), lsys);
end:

CriticalPoints:=proc(Equations, Fam, Inequalities, Inequations, vars, opts:={})
local i, newsols, j, pol, sols, positive, J, minors, sysminors, a, b, pos_pol,
nz_pol, pol2, np, cstr;
  J      := convert(linalg:-jacobian([op(Equations), op(Fam)], vars), Matrix):
  minors := ComputeMaximalMinors(J):

  positive := remove(member, Inequalities, Fam);
  pos_pol := sort(select(member, Inequalities, Fam), 
             (a, b)->degree(a) <= degree(b)):
  if nops(pos_pol)>0 then 
    pol:=pos_pol[1]:
  else 
    pol:=1:
  end if;
  nz_pol := select(member, Inequations, [op(Fam),op(map(a->-a,Fam))]):
  if nops(nz_pol)>0 then 
    pol2:=nz_pol[1]:
  else 
    pol2:=1:
  end if;
  np := nops(positive);
  cstr:=[op(remove(member, Inequalities, Fam)), op(remove(member,
  Inequations, [op(Fam), op(map(a->-a, Fam))]))]:
  sysminors:=SplitSystem_cstr(minors, cstr):
  sysminors:=map(sys->if not(member(0,map(degree, sys))) then sys fi,
  sysminors):
  sysminors:=map(sys-> remove(member, sys, [0]), sysminors);
  sysminors:=map(sys-> if nops(convert(sys, set) intersect convert(Inequations, set))=0 then sys fi, sysminors);
  sysminors:=map(sys-> if nops(convert(sys, set) intersect convert(Inequalities, set))=0 then sys fi, sysminors);
  sols:=[]:
  for j from 1 to nops(sysminors) do 
    newsols:=MSolveRealRoots([op(Equations), seq(Fam[i]-Fam[1],i=2..nops(Fam)), 
                    rag_sat_var*Fam[1]-1, op(sysminors[j])], [rag_sat_var, op(vars)], 
                    [pol, op(positive), 
                    pol2, op(remove(member, Inequations, [op(Fam),op(map(a->-a, Fam))]))], opts):

    if newsols = [1] then
      newsols:= CriticalPointsSingular(Equations, Fam, positive,
      nz_pol, pol, pol2, vars, sysminors[j], opts);
    else 
      newsols:=AdmissibleSolutions(newsols, np + 1);
    end if;
    newsols:=
          map(pt->map(_c->if member(lhs(_c),vars) then _c fi, pt), newsols):
    sols:=[op(sols), op(newsols)]:
  end do;
  return sols;
end proc:

#Warning: ce n'est correct que sur les solutions exactes
ExactSolSelection:=proc(sols, Inequalities, Inequations, vars)
local osol, i, sol, badsols;
  badsols:=[];
  for i from 1 to nops(sols) do 
    osol := sols[i];
    sol := map(c -> lhs(c) = (rhs(c)[1] + rhs(c)[2]) / 2, osol);
    if nops(map(p->if subs(sol,p)<=0 then 1 fi, Inequalities)) > 0 or
       nops(map(p->if subs(sol,p)=0 then 1 fi, Inequations)) > 0 then
       badsols:=[op(badsols), osol];
    end if;
  end do;
  return remove(member,sols,badsols);
end proc:

GenerateCriticalPointsFamilies:=proc(FamPositive, FamNotNull)
local lFam, j, _l, lFam1, lFam2, pol;
  lFam:=[FamPositive]:
  for j from 1 to nops(FamNotNull) do 
    pol:=FamNotNull[j];
    lFam1:=map(_l->[op(_l), pol], lFam);
    lFam2:=map(_l->if nops(_l)>0 then [op(_l), -pol] fi, lFam);
    lFam:=[op(lFam1),op(lFam2)]:
  end do;
  return lFam;
end proc:

FamCriticalPoints:=proc(Equations, FamPositive, FamNotNull, 
                     Inequalities, Inequations, vars, opts:={})
local i, Families, Fam, sols, a, b;
  sols:=[]:
  Families:=GenerateCriticalPointsFamilies(FamPositive, FamNotNull);
  for i from 1 to nops(Families) do
    Fam:=sort(Families[i],(a, b)->degree(a)<=degree(b)):
    sols:=[op(sols), 
           op(CriticalPoints(Equations, Fam, Inequalities, Inequations,
           vars, opts))]:
  end do;
  return sols;
end proc:

UnivariateSolveFamily:=proc(Equations, Fam, Inequalities, Inequations, vars)
local g, upol, f, uroots, i, sq, sols, q, newpol, p, mid;
  lprint(args);
  g:=0:
  for i from 1 to nops(Equations) do 
    g:=gcd(g, Equations[i]):
  end do;
  if degree(g)=0 then return []; end if;

  if nops(Equations) > 0 then 
  for i  from 1 to nops(Inequalities) do
    q:=gcd(g, Inequalities[i]):
    if degree(q)>0 then
      g:=numer(normal(g/q)):
    end if;
  end do;
  if degree(g)=0 then return []; end if;

  for i  from 1 to nops(Inequations) do
    q:=gcd(g, Inequations[i]):
    if degree(q)>0 then
      g:=numer(normal(g/q)):
    end if;
    g:=gcd(g, Inequations[i]):
  end do;
  end if;
  if degree(g)=0 then return []; end if;

  #Now, we search for real roots of g which satisfy the constraints.
  if degree(g) > 0 then 
    sols:=MSolveRealRoots([g], [op(indets(g))], [op(Inequalities),
                          op(Inequations)]):
    sols:=AdmissibleSolutions(sols, nops(Inequalities)):
    return sols;
  end if;

  #Case where g = 0
  upol:=mul(f, f in Fam):
  for i from 1 to nops(Equations) do 
    upol:=gcd(upol, Equations[i]):
  od:
  sq:=sqrfree(upol)[2]:
  if nops(sq)>0 then 
    upol:=mul(p[1], p in sq);
  else 
    upol:=1:
  end if;
  newpol:=1:
  for i from 1 to nops(Equations) do 
    newpol:=expand(newpol*Equations[i]):
    if degree(newpol)>1 then 
      sq:=sqrfree(newpol)[2]:
      newpol:=mul(p[1], p in sq);
    end if;
  end do;
  for i from 1 to nops(Fam) do 
    newpol:=expand(newpol*Fam[i]):
    if degree(newpol)>1 then 
      sq:=sqrfree(newpol)[2]:
      newpol:=mul(p[1], p in sq);
    end if;
  end do;
  for i from 1 to nops(Inequalities) do 
    newpol:=expand(newpol*Inequalities[i]):
    if degree(newpol)>1 then 
      sq:=sqrfree(newpol)[2]:
      newpol:=mul(p[1], p in sq);
    end if;
  end do;
  for i from 1 to nops(Inequations) do
    newpol:=expand(newpol*Inequations[i]):
    if degree(newpol)>1 then 
      sq:=sqrfree(newpol)[2]:
      newpol:=mul(p[1], p in sq);
    end if;
  end do;
  newpol:=expand(newpol);
  uroots:=MSolveRealRoots([newpol], vars, [])[2]:
  uroots:=map(p->rhs(p[1]), uroots);
  uroots:=sort(uroots, (a, b)->a[2]<b[1]);
  if HasOverLap(uroots) then 
    lprint(newpol);
    error "Bug in msolve";
  end if;
  sols:=[]:
  for i from 1 to nops(uroots)-1 do 
    mid:=SmallMidRational(uroots[i][2], uroots[i+1][1]);
    lprint(evalf(mid));
    if not(member(-1, map(sign, subs(vars[1]=mid, Inequalities))))
      then 
      sols:=[op(sols), [vars[1]=[mid,mid]]];
    end if;
  end do;
  if nops(uroots) > 0 then 
    mid:=floor(uroots[1][1]-1);
    if not(member(-1, map(sign, subs(vars[1]=mid, Inequalities))))
       then 
      sols:=[[vars[1]=[mid,mid]],op(sols)];
    end if;
    mid:=ceil(uroots[-1][2]+1);
    if not(member(-1, map(sign, subs(vars[1]=mid, Inequalities))))
       then 
      sols:=[op(sols), [vars[1]=[mid, mid]]];
    end if;
  else 
    if not(member(-1, map(sign, subs(vars[1]=0, Inequalities)))) then
      sols:=[[vars[1]=[0,0]]]:
    end if;
  end if;
  return sols;
end proc;


AdmissibleSolutions:=proc(tsols, np)
local j, sgn, sols;
  if nops(tsols)=2 then 
    return tsols[2];
  end if;
  sols:=[];
  if nops(tsols[2]) > 0 and nops(tsols) > 2 then 
    for j from 1 to nops(tsols[2]) do 
      sgn := map(op, tsols[3][j][1..np]);
      if not(member(-1, map(sign, sgn))) and not(member(0, sgn)) and
         not(member(0,  map(op, tsols[3][j][np+1..-1]))) then 
        sols:=[op(sols), tsols[2][j]];
      end if;
    end do;
  end if;
  return sols;
end proc: 

#Generates families to solve
GenerateDeformedFamilies_eps:=proc(eqs, FamPositive, FamNotNull, 
                                   vars, eps, cstr, opts:={})
local i, lsys, deform, j, deform1, deform2, pol, lc, sys;
  deform:=[FamPositive]:
  for i from 1 to nops(FamNotNull) do 
    pol:=FamNotNull[i]:
    deform1:=map(_l->[op(_l), pol], deform):
    deform2:=map(_l->if nops(_l)>0 then [op(_l), -pol] fi, deform):
    deform:=[op(deform1), op(deform2)]:
  end do;
  lsys:=[]:
  for i from 1 to nops(deform) do 
    lc:=CoeffDeform_eps(eqs, deform[i], vars, eps, cstr, opts);
    sys:=[op(eqs), seq(deform[i][j]-lc[j]*eps,j=1..nops(deform[i]))]:
    lsys:=[op(lsys), sys]:
  end do;
  return lsys;
end proc:

#Generates families to solve
OLDGenerateDeformedFamilies_eps:=proc(FamPositive, FamNotNull, lcpos,
                                   lcineq, eps)
local i, deform, deform1, deform2, pol;
  deform:=[[seq(FamPositive[i]-lcpos[i]*eps, i=1..nops(lcpos))]]:
  #deform:=[map(_p->_p-eps, FamPositive)]:
  for i from 1 to nops(FamNotNull) do 
    pol:=FamNotNull[i]:
    deform1:=map(_l->[op(_l), pol-lcineq[i]*eps], deform):
    deform2:=map(_l->if nops(_l)>0 then [op(_l), pol+lcineq[i]*eps] fi, deform):
    deform:=[op(deform1), op(deform2)]:
  end do;
  return deform;
end proc:

ConstrainedValues := proc(Equations, FamPositive, FamNotNull, vars, Inequalities, Inequations, opts:={})
local tsols, sgn, npos, i, eps, sols, _l, deform, mvars, rr, j, k,
newtsols, lc; 
  deform := GenerateDeformedFamilies_eps(Equations,
                                         FamPositive, FamNotNull, 
                                         vars, eps);
  sols:=[]:
  npos := nops(Inequalities);
  rr:=rand(1..65520):
  for i from 1 to nops(deform) do 
    if indets(deform[i])<>indets(vars) then
      mvars:=indets(vars) minus indets(deform[i]);
    else
      mvars:=[]:
    end if;
    tsols := MSolveRealRoots([op(mvars), op(Equations), op(deform[i])], 
             [eps, op(vars)], [op(Inequalities), op(Inequations)], opts);
    if tsols[1] > 0 or tsols[1]=-1 then 
      tsols := MSolveRealRoots(
               [_l*eps-1,op(mvars), op(Equations), 
                op(deform[i])], 
               [_l, eps, op(vars)], [op(Inequalities), op(Inequations)], opts);
      if tsols[1] > 0 then 
        lprint(args);
        error "Deformed systems are not 0-dim";
      end if;
    end if;
    if nops(Inequalities)>0 then 
      tsols := AdmissibleSolutions(tsols, npos);
    else 
      if nops(Inequations) > 0 then 
        newtsols:=[]:
        for j from 1 to nops(tsols[2]) do 
          for k from 1 to nops(tsols[3][j]) do 
            if tsols[3][j][k][1] = 0 and tsols[3][j][k][2] = 0 then 
              break;
            end if;
            if tsols[3][j][k][1] = 0 and tsols[3][j][k][2] <> 0 then 
              lprint([[_l*eps-1,op(mvars), op(Equations), op(map(_p->coeff(_p, eps, 0)-rr()*eps, deform[i]))], 
               [_l, eps, op(vars)], [op(Inequalities), op(Inequations)], opts]);
              error "Bug in MSolveRealRoots to be investigated";
            end if;
            if tsols[3][j][k][1] <> 0 and tsols[3][j][k][2] = 0 then 
              lprint([[_l*eps-1,op(mvars), op(Equations), op(map(_p->coeff(_p, eps, 0)-rr()*eps, deform[i]))], 
               [_l, eps, op(vars)], [op(Inequalities), op(Inequations)], opts]);
              error "Bug in MSolveRealRoots to be investigated";
            end if;
          end do;
        end do;
        tsols:=newtsols:
      else 
        tsols:=tsols[2]:
      end if;
    end if;
    sols := [op(sols), op(tsols)];
  end do;
  if nops(FamPositive) > 0 then 
    sols := map(_p-> if subs(_p, eps)[1]>0 then _p else fi, sols);
  end if;
  sols := map(_p->map(_c->if member(lhs(_c), vars) then _c else fi, _p), sols);
  return sols;
end proc;


UnboundedComponents:=proc(Equations, FamPositive, FamNotNull, Inequalities, Inequations, vars, opts:={})
local cstr, i, verb, isempty, hyp, bounds, v, hypsol, ls,
NewEquations, NewFamPositive, NewFamNotNull, NewInequations,
NewInequalities, newvars, tsols, Fam, sols, Positive, NotNull;

  if type(subs(opts, "verb"), integer) then 
    verb:=subs(opts, "verb");
  else 
    verb:=0:
  end if;
  if type(subs(opts, "isempty"), integer) then 
    isempty:=subs(opts, "isempty"):
  else 
    isempty:=0:
  end if;
  
  sols:=[]:
  Fam:=[op(FamPositive), op(FamNotNull)]:
  Fam:=sort(Fam, (a, b)->degree(a)<=degree(b)):
  Positive:=remove(member, Inequalities, FamPositive): 
  NotNull:=remove(member, Inequations, FamNotNull):
  if verb>=1 then printf("[b:") end if;
  hyp, bounds := ComputeBounds(Equations, Fam, Positive, NotNull, vars, opts);
  if verb>=1 then printf(" -> %a, %a", hyp, bounds); end if;
  if verb>=1 then printf("]") end if;
  for i from 1 to nops(bounds) do
      v               := indets(hyp)[1];
      hypsol          := solve(hyp-bounds[i], v);
      ls              := {v = hypsol}:
      NewEquations    := map(expand, subs(ls, Equations));
      NewFamPositive  := map(expand, subs(ls, FamPositive));
      NewFamNotNull   := map(expand, subs(ls, FamNotNull));
      NewInequations  := map(expand, subs(ls, Inequations));
      NewInequalities := map(expand, subs(ls, Inequalities));
      newvars:=remove(member, vars, [v]);
      tsols:=SolveFamily(NewEquations, NewFamPositive, NewFamNotNull, 
                         NewInequalities, NewInequations,
             newvars, opts);
      if degree(hypsol) <= 0 then 
        sols := [op(sols), op(map(s->[v = [hypsol, hypsol], op(s)], tsols))];
      else
        sols := [op(sols), op(map(s->[v = [subs(map(c->lhs(c)=rhs(c)[1], s), hypsol), 
                                           subs(map(c->lhs(c)=rhs(c)[2], s), hypsol)], 
                                      op(s)], tsols))];
      end if;
      if nops(sols) > 0 and isempty > 0 then 
        return sols;
      end if;
  end do;
  return sols;
end proc;

DegenerateDeformedSystem:=proc(sys, ld, Inequalities, Inequations, vars, eps, opts)
local gb, sols;
  lprint(args);
  #debug(SaturateIntersect, ModSatIntersectLM,MinimalGeneratorsDichotomy);
  gb:=SaturateIntersect(sys, ld[1], ld, [op(vars), eps], opts):
  lprint(gb);
  sols:=MSolveRealRoots([op(gb), op(sys), ld[1]], [op(vars), eps],
        [op(Inequalities), eps, op(Inequations)], opts):
  lprint(sols);
  quit;
  sols:=AdmissibleSolutions(sols, nops(Inequalities)):
  sols:=map(_p->map(_c->if member(lhs(_c), vars) then _c fi, _p), sols):
  return sols;
end proc:

ZeroDimBoundaries:=proc(Equations, FamPositive, FamNotNull,
                        Inequalities, Inequations, vars, opts:={})
local verb, eps, lsys, emin, delta, J, i, j, sols, lsols;

  if type(subs(opts, "verb"), integer) then 
    verb:=subs(opts, "verb");
  else 
    verb:=0:
  end if;
  
  lsys:=GenerateDeformedFamilies_eps(Equations, FamPositive,
  FamNotNull, vars, eps, {op(Inequalities), op(Inequations)}):
  lprint(lsys);
  J:=convert(linalg:-jacobian([op(Equations), op(FamPositive), op(FamNotNull)],
              vars), Matrix);
  delta:=ComputeMaximalMinors(J):
  lsols:=[]:
  for i from 1 to nops(lsys) do
    sols:=MSolveRealRoots([rag_sat_var*eps-1,op(lsys[i]),op(delta)], 
          [rag_sat_var, eps, op(vars)],
          [op(Inequalities), eps, op(Inequations)], opts):
    if sols[1]>0 then 
      debug(DegenerateDeformedSystem);
      sols:=DegenerateDeformedSystem(lsys[i], delta, Inequalities,
            Inequations, vars, eps, opts):
    end if;
    lprint(evalf(sols));
    if nops(FamPositive)>0 then 
      sols:=AdmissibleSolutions(sols, nops(Inequalities)+1);
    else 
      sols:=AdmissibleSolutions(sols, nops(Inequalities));
    end if;
    lprint(evalf(sols));
    sols:=map(_p->map(_c->if member(lhs(_c), vars) then _c fi, _p), sols):
    lsols:=[op(lsols), op(sols)]:
  
    emin:=2:
    for j from 1 to nops(Inequalities) do 
      sols:=MSolveRealRoots([rag_sat_var*eps-1,op(lsys[i]),Inequalities[j]], 
            [rag_sat_var, eps, op(vars)],
            [], opts):
      if sols[1]>0 then 
        lprint(args);
        error "Bug in ZeroDimBoundaries";
      end if;
      if nops(sols[2])>0 then 
        emin:=min(map(abs, map(op, map(_p->subs(_p, eps), sols[2]))));
      end if;
    end do;
    for j from 1 to nops(Inequations) do 
      sols:=MSolveRealRoots([rag_sat_var*eps-1,op(lsys[i]),Inequations[j]], 
            [rag_sat_var, eps, op(vars)],
            [], opts):
      if sols[1]>0 then 
        lprint(args);
        error "Bug in ZeroDimBoundaries";
      end if;
      lprint("sols", evalf(sols));
      if nops(sols[2])>0 then 
        emin:=min(map(abs, map(op, map(_p->subs(_p, eps), sols[2]))));
      end if;
    end do;

    if emin=2 then 
      sols:=MSolveRealRoots(subs(eps=emin,lsys[i]), vars,
                            [op(Inequalities), op(Inequations)], opts):
      while sols[1]>0 do 
        if verb >= 1 then printf("*"); end if;
        emin:=emin/2:
        sols:=MSolveRealRoots(subs(eps=emin,lsys[i]), vars,
                            [op(Inequalities), op(Inequations)], opts):
      end do:
    else 
      emin:=emin/2:
      sols:=MSolveRealRoots(subs(eps=emin,lsys[i]), vars,
                            [op(Inequalities), op(Inequations)], opts):
      while sols[1]>0 do 
        if verb >= 1 then printf("*"); end if;
        emin:=emin/2:
        sols:=MSolveRealRoots(subs(eps=emin,lsys[i]), vars,
                            [op(Inequalities), op(Inequations)], opts):
      end do:
    end if;
    sols:=AdmissibleSolutions(sols, nops(Inequalities));
    lprint(evalf(sols));

    lsols:=[op(lsols), op(sols)]:
  end do;

  return lsols;
end proc;

SolveFamily:=proc(Equations, FamPositive, FamNotNull, Inequalities, Inequations, vars, opts:={})
local i, ls, cp, hyp, bounds, sols, tsols, a, b, verb, Fam, cstr, NewFamPositive,
NewFamNotNull, isempty, newsols, isbounded;

  if type(subs(opts, "verb"), integer) then 
    verb:=subs(opts, "verb");
  else 
    verb:=0:
  end if;
  if type(subs(opts, "isempty"), integer) then 
    isempty:=subs(opts, "isempty"):
  else 
    isempty:=0:
  end if;
  if type(subs(opts, "isbounded"), integer) then 
    isbounded:=subs(opts, "isbounded"):
  else 
    isbounded:=0:
  end if;


  if member(0, [op(FamPositive), op(Inequalities)]) or 
     member(0, [op(FamNotNull), op(Inequations)]) then 
    return [];
  end if;
  if member(-1,
            map(_pol->if degree(_pol)<=0 then sign(_pol) fi, 
            [op(FamPositive), op(Inequalities)])) then 
    return [];
  end if;

  Fam:=[op(FamPositive), op(FamNotNull)]:
  Fam:=sort(Fam, (a, b)->degree(a)<=degree(b)):
  sols:=[];
  if nops(vars) = 1 then 
    return UnivariateSolveFamily(Equations, FamPositive, Inequalities, Inequations, vars); 
  end if;
  if nops(Equations) = nops(vars) then 
    sols:=MSolveRealRoots(Equations, vars, [op(FamPositive),
                          op(Inequalities), op(FamNotNull),
                          op(Inequations)], opts);
    if sols[1]>0 then 
      error "Bug detected";
    end if;
    sols:=AdmissibleSolutions(sols,
    nops(FamPositive)+nops(Inequalities));
    return sols;
  end if;

  if nops(Equations) + nops(Fam) < nops(vars) then 
    sols := UnboundedComponents(Equations, FamPositive, FamNotNull, 
              Inequalities, Inequations, vars, opts):
    if nops(sols) > 0 and isempty >0 then 
      return sols;
    fi;
  end if;
  if nops(Equations) + nops(Fam) = nops(vars) then 
    sols:=ZeroDimBoundaries(Equations, FamPositive, FamNotNull,
    Inequalities, Inequations, vars, opts):
    if nops(sols) > 0 and isempty > 0 then 
      return sols;
    fi;
  end if;
  if 0=1 and nops(Equations) + nops(Fam) > nops(vars) then 
    if verb>=1 then printf("[cn"); end if;
    try
      sols := ConstrainedValues(Equations, FamPositive, FamNotNull, vars, remove(member,
               Inequalities, FamPositive), remove(member, Inequations, FamNotNull), opts);
    catch:
      lprint(FamPositive, FamNotNull, vars, opts);
      error "Problem when calling ConstrainedValues";
    end:
    if verb>=1 then printf("]"); end if;
    if nops(sols) > 0 and isempty >0 then 
      return sols;
    fi;
  end if;
  return sols;
end proc:

SemiAlgebraicSolveIterateOnFamilies:=proc(Equations, Families, Inequalities, Inequations, vars, opts)
local i, Fam, sols, newsols, verb, st, pos, nonzero, isempty;

  if type(subs(opts, "verb"), integer) then 
    verb:=subs(opts, "verb");
  else 
    verb:=0:
  end if;
  if type(subs(opts, "isempty"), integer) then 
    isempty:=subs(opts, "isempty"):
  else 
    isempty:=0:
  end if;

  sols := []:
  st:=time[real]():
  for i from 1 to nops(Families) do 
    Fam  := Families[i]:
    pos := map(idx->Inequalities[idx], Fam[1]);
    nonzero := map(idx->Inequations[idx], Fam[2]):
    if nops(indets([op(Equations), op(Inequalities), op(Inequations)])) = nops(vars) then
      if verb>=1 then 
        printf("<");
      end if;
      newsols:=SolveFamily(Equations, pos, nonzero, 
                           Inequalities, Inequations, vars, opts):
      if verb>=1 then 
        printf(":->[ns=%d]:",nops(newsols));
      end if:
      if verb >=1 and (time[real]()-st>=10 or i = nops(Families)) then
        printf("{%.2f%%}",100.0*i/nops(Families));
        st:=time[real]():
      end if; 

      if verb>=1 then 
        printf(">");
      end if;
      sols := [op(sols), 
               op(newsols)
              ];
      if isempty > 0 and nops(sols) > 0 then 
        return sols; 
      end if;
    end if;
  end do;
  if verb >=1 then 
    printf("++{%d}++",nops(sols));
  end if;
  return sols;
end proc:

PointsPerComponentsAlgebraic:=proc(Equations, Inequalities, Inequations, opts:={})
local vars, boo, singminors, svars, spt, i, sineq, sols, _sol, spos,
newsols, ls, var, hypsol, isempty, hyp, minors, gendeg, _pol; 
  if type(subs(opts, "isempty"), integer) then 
    isempty:=subs(opts, "isempty"):
  else 
    isempty:=0:
  end if;
  if member(0, Inequations) or member(0, Inequalities) or 
      member(-1, map(sign, map(_pol->if degree(_pol)=0 then _pol fi, Inequalities))) then 
    return [];
  end if; 
  vars:=[op(indets(Equations))]:
# Check that Equations is regular enough
  boo, singminors := IsRegular(Equations, vars):
  if boo = false then 
    error "Singular case not implemented yet";
  end if;
  svars:=remove(member, [op(indets(Inequalities) union indets(Inequations))], vars):
  spt:=GoodFiberValue_svars(svars, Inequalities, Inequations);
  spos:=subs(spt, Inequalities):
  sineq:=subs(spt, Inequations):

  if nops(vars) = nops(Equations) then 
    if member(-1, map(sign, map(_pol-> if degree(_pol)<=0 then _pol fi, spos))) then 
      return [];
    end if;
    sols := MSolveRealRoots(Equations, vars, [op(spos), op(sineq)],
      opts):
    sols := AdmissibleSolutions(sols, nops(spos)):
    sols:=map(_sol->[op(_sol), op(spt)], sols):
    return sols;
  end if;

# Find a good choice of projection
  hyp, minors, gendeg := FindGenericLineRegular(Equations, [], [], singminors,  
                                                vars, opts):

# Compute the associated critical points and selects those satisfying
# the constraints
  if gendeg = 0 then 
    return [0, []];
  end if;
  sols := MSolveRealRoots([op(Equations), op(minors)], vars,
          [op(spos), op(sineq)], opts);
  sols:=AdmissibleSolutions(sols, nops(spos));
  sols:=map(_sol->[op(_sol), op(spt)], sols):

  if nops(sols) > 0 and isempty > 0 then 
     return sols;
  end if;

# Recursive call on an arbitrary fiber
  var := indets(hyp)[1];
  ls, hypsol := GoodFiberValue(var, hyp, Inequalities, Inequations);
  newsols := PointsPerComponentsAlgebraic(subs(ls, Equations),
    subs(ls, Inequalities), subs(ls, Inequations), opts);
  if degree(hypsol) <= 0 then 
     newsols := map(s->[var = [hypsol, hypsol], op(s)], newsols);
  else
     newsols := map(s->[var = [subs(map(c->lhs(c)=rhs(c)[1], s), hypsol), 
                                        subs(map(c->lhs(c)=rhs(c)[2], s), hypsol)], 
                                   op(s)], newsols);
  end if;
  newsols:=map(_sol->[op(_sol), op(spt)], newsols):
  if nops(newsols) > 0 and isempty > 0 then 
     return newsols;
  end if;

  return [op(sols), op(newsols)];
end proc;

SemiAlgebraicSolve:=proc(Equations, Inequalities, Inequations, opts:={})
local newsols, _toStudy, l, i, pol, boo, vars, _Studied, singminors, pt, _c,
midsols, verb, sols, lsigns, _l, nc, isempty;

  if type(subs(opts, "verb"), integer) then 
    verb:=subs(opts, "verb");
  else 
    verb:=0:
  end if;

  if type(subs(opts, "isempty"), integer) then 
    isempty:=subs(opts, "isempty");
  else 
    isempty:=0:
  end if;

  vars:=[op(indets([op(Equations), op(Inequalities), op(Inequations)]))];
  sols := [];
  lsigns:={}:
  if nops(Equations) > 0 then
    sols := PointsPerComponentsAlgebraic(Equations, Inequalities, Inequations, opts);
    if isempty>=1 and nops(sols) > 0 then
      return sols;
    end if;
  end if;

  nc := nops(Inequalities) + nops(Inequations);
  _Studied := [[[],[]]]:
  newsols := []:

#TODO: improve the choice of Inequalities and introduce criteria to
#reduce the combinatorial complexity  
  for i from 1 to nops(Inequalities) do 
    pol := Inequalities[i];
    _toStudy := map(l->if nops(l[1])+nops(l[2]) <= nops(vars) then [[op(l[1]), i],[]] fi, _Studied):

    if verb>=1 then 
      printf("\nDealing with (positivity) constraint %a / %a of degree %d\n", 
             i, nc, degree(pol));
      printf("Number of systems to study is %d\n", nops(_toStudy));
    end if;

    newsols:=SemiAlgebraicSolveIterateOnFamilies(Equations, _toStudy, Inequalities,
          Inequations, vars, opts);

    midsols:=map(pt->map(_c->lhs(_c)=(rhs(_c)[1]+rhs(_c)[2])/2, pt), newsols):
    lsigns:=lsigns union convert(map(l->map(sign, l), map(pt->subs(pt,Inequations), midsols)), set):
    if verb >= 1 and nops(Inequations)>0 then 
     printf("Signs of inequations at computed points: %a\n", lsigns);
    end if;
    sols:=[op(sols), op(newsols)]:
    _Studied := [op(_Studied), op(_toStudy)]:
    if isempty>=1 and nops(sols) > 0 then
      return sols;
    end if;
  end do;

#TODO: improve the choice of Inequations and introduce criteria to
#reduce the combinatorial complexity  
  for i from 1 to nops(Inequations) do 
    pol := Inequations[i]:
    _toStudy := map(l->if nops(l[1])+nops(l[2]) <= nops(vars) then [l[1],[op(l[2]), i]] fi, _Studied):
    if verb>=1 then 
      printf("\nDealing with (non-zero) constraint %a / %a of degree %d\n", 
              nops(Inequalities) + i, nc, degree(pol));
      printf("Number of systems to study is %d\n", nops(_toStudy));
    end if;

    newsols:=SemiAlgebraicSolveIterateOnFamilies(Equations, _toStudy, Inequalities,
          Inequations, vars, opts);
    midsols:=map(pt->map(_c->lhs(_c)=(rhs(_c)[1]+rhs(_c)[2])/2, pt), newsols):
    lsigns:=lsigns union convert(map(l->map(sign, l), map(pt->subs(pt,Inequations), midsols)), set):

    if verb >= 1 and nops(Inequations)>0 then 
     printf("Signs of inequations at computed points: %a\n", lsigns);
    end if;
    
    if nops(map(_l->_l[i], lsigns)) = 1 then 
#Inequations[i] has no sign change, it can be removed
      lprint("Removes non-zero constraint", Inequations[i]);
      _toStudy := map(_l->if not(member(Inequations[i], _l)) then _l fi,
      _toStudy);
    end if;
    sols:=convert(convert([op(sols), op(newsols)], set), list):
    if isempty>=1 and nops(sols) > 0 then
      return sols;
    end if;

    _Studied := [op(_Studied), op(_toStudy)]:
  end do;
  return sols;
end proc:


HasRealSolutions:=proc(eqs, pos, ineqs, opts:={})
local isempty, _l, newopts;

  if type(subs(opts, "isempty"), integer) then 
    isempty:=subs(opts, "isempty");
  else 
    isempty:=0:
  end if;

  if isempty = 0 then 
    return SemiAlgebraicSolve(eqs, pos, ineqs, opts union {"isempty"=1});
  else 
    newopts:=map(_l->if lhs(_l)="isempty" then lhs(_l)=1 else _l fi, opts);
    return SemiAlgebraicSolve(eqs, pos, ineqs, newopts);
  end if;
end proc;

PointsPerComponents:=proc(eqs, pos, ineqs, opts:={})
local isempty, _l, newopts;

  if type(subs(opts, "isempty"), integer) then 
    isempty:=subs(opts, "isempty");
  else 
    isempty:=0:
  end if;

  if isempty = 0 then 
    return SemiAlgebraicSolve(eqs, pos, ineqs, opts);
  else 
    newopts:=map(_l->if lhs(_l)="isempty" then lhs(_l)=0 else _l fi, opts);
    return SemiAlgebraicSolve(eqs, pos, ineqs, newopts);
  end if;
end proc;
