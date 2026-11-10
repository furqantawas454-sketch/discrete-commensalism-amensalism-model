(* ::Package:: *)
(* ============================================================================
   01_fixed_points.wl   --  Sections 3 and 4 of the paper
   (a) symbolic check of the fixed point equations used in Theorem 4.1
   (b) numerical fixed points of Examples 8.1-8.3, compared with a direct NSolve
   Run:  wolframscript -file 01_fixed_points.wl   (from this folder)
   ============================================================================ *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "00_model.wl"}]];

Print["\n=== (a) Symbolic checks for Theorem 4.1 ==="];
Clear[x, y, z, a1, b1, a2, b2, a3, b3, be, al, c, m, n, xs, L];
ps = {a1, b1, a2, b2, a3, b3, be, al, c, m, n};

(* Proposition 3.1: G(X) = X  <=>  h F(X) = 0 *)
Print["G(X) - X == h F(X): ",
  Simplify[mapG[{x, y, z}, ps, h] - {x, y, z} - h fieldF[{x, y, z}, ps]] === {0, 0, 0}];

(* Case z = 0, x = a1/b1: quadratic y^2 + A2 y + A3 = 0 (Eq. 4.4) *)
X1 = a1/b1;
lhsE2 = Expand[-((a2 - b2 y)(m X1 + n y + 1) + c X1)/(b2 n)];
rhsE2 = y^2 + (a1 m/(b1 n) + 1/n - a2/b2) y - (a1 a2 m/(b1 b2 n) + a2/(b2 n) + a1 c/(b1 b2 n));
Print["E2 quadratic matches Eq. (4.4): ", Simplify[lhsE2 - rhsE2] === 0];

(* Case z = a3/(b3 x): first equation times -x/b1 gives the cubic phi(x) (Eq. 4.2) *)
cubicCheck = Simplify[(x (a1 - b1 x) + be a3/(b3 x)) (-x/b1) - (x^3 - a1/b1 x^2 - a3 be/(b1 b3))];
Print["Cubic phi(x) matches Eq. (4.2): ", cubicCheck === 0];

(* Discriminant of phi is -4 p^3 q - 27 q^2 < 0  (Remark 4) *)
Print["Discriminant of x^3 - p x^2 - q: ", Factor[Discriminant[x^3 - pp x^2 - qq, x]]];

(* Coexistence: psi(y) = b2 n y^2 + B y - C (Eq. 4.7) and C = (m x* + 1) Delta *)
lhsE4 = Expand[-((L - b2 y)(m xs + n y + 1) + c xs)];
rhsE4 = b2 n y^2 + (b2 (m xs + 1) - n L) y - (L (m xs + 1) + c xs);
Print["E4 quadratic matches Eq. (4.7): ", Simplify[lhsE4 - rhsE4] === 0];
Print["C == (m x* + 1) Delta: ", Simplify[(L (m xs + 1) + c xs) - (m xs + 1) (L + c xs/(m xs + 1))] === 0];

Print["\n=== (b) Fixed points of the three examples ==="];
Do[
  Module[{p = examples[name], fp, direct, eqs, xx, yy, zz},
   fp = fixedPoints[p];
   eqs = Numerator[Together[#]] == 0 & /@ fieldF[{xx, yy, zz}, p];
   direct = Sort[Select[{xx, yy, zz} /. NSolve[eqs, {xx, yy, zz}, Reals],
       #[[1]] > 0 && #[[2]] >= 0 && #[[3]] >= 0 &]];
   Print["\n", name, "   Delta = ", NumberForm[fp["Delta"], 6]];
   Print["  E1 = ", NumberForm[N[fp["E1"]], 6], "   E2 = ", NumberForm[N[fp["E2"]], 6]];
   Print["  E3 = ", NumberForm[N[fp["E3"]], 6], "   E4 = ", NumberForm[N[fp["E4"]], 6]];
   Print["  formulas agree with direct NSolve: ",
     Max[Abs[Flatten[Sort[N[{fp["E1"], fp["E2"], fp["E3"], fp["E4"]}]] - direct]]] < 10^-8]],
  {name, Keys[examples]}];
