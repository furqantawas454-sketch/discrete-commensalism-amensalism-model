(* ::Package:: *)
(* ============================================================================
   02_stability.wl   --  Section 5 of the paper
   (a) symbolic checks: structure of the Jacobian, factorization (5.2),
       trace/determinant (5.6) of the (x,z)-block, y-eigenvalue at E4,
       eigenvalues at E1 and E2
   (b) stability thresholds of Examples 8.1-8.3 and a numerical test of
       Theorem 5.4 just below and just above h_E4
   Run:  wolframscript -file 02_stability.wl
   ============================================================================ *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "00_model.wl"}]];

Print["\n=== (a) Symbolic checks for Section 5 ==="];
Clear[x, y, z, a1, b1, a2, b2, a3, b3, be, al, c, m, n, h, lam];
ps = {a1, b1, a2, b2, a3, b3, be, al, c, m, n};
J = jacG[{x, y, z}, ps, h];
Print["J12 == 0 and J32 == 0: ", Simplify[{J[[1, 2]], J[[3, 2]]}] === {0, 0}];
Print["Characteristic polynomial factors as in (5.2): ",
  Simplify[Det[lam IdentityMatrix[3] - J] -
     (lam - J[[2, 2]]) ((lam - J[[1, 1]]) (lam - J[[3, 3]]) - J[[1, 3]] J[[3, 1]])] === 0];

(* (x,z)-block at E3/E4: use z = a3/(b3 x) and a1 = b1 x - beta z/x *)
A = {{a1 - 2 b1 x, be}, {-a3 z/x^2, a3/x - 2 b3 z}} /. z -> a3/(b3 x) /. a1 -> b1 x - be a3/(b3 x^2);
Print["trace(A) == -(b1 x* + beta z*/x* + a3/x*): ",
  Simplify[Tr[A] + (b1 x + be (a3/(b3 x))/x + a3/x)] === 0];
Print["det(A) == a3 b1 + 2 a3 beta z*/x*^2: ",
  Simplify[Det[A] - (a3 b1 + 2 a3 be (a3/(b3 x))/x^2)] === 0];

(* y-eigenvalue at E4: eliminate a2 with the fixed point relation *)
J22E4 = Simplify[J[[2, 2]] /. a2 -> b2 y - c x/(m x + n y + 1) + al z];
Print["J22(E4) == 1 - h rho4: ",
  Simplify[J22E4 - (1 - h y (b2 + c n x/(m x + n y + 1)^2))] === 0];

(* eigenvalues at E1 *)
eigE1 = Eigenvalues[J /. {x -> a1/b1, y -> 0, z -> 0}];
Print["Eigenvalues at E1: ", Simplify[eigE1]];

(* Lemma 5.2: roots of F(-1) = D h^2 + 2 T h + 4 as a quadratic in h *)
Clear[TT, DD];
Print["Roots of F(-1) in h: ", Solve[DD hh^2 + 2 TT hh + 4 == 0, hh]];

Print["\n=== (b) Stability thresholds of the examples ==="];
Do[
  Module[{p = examples[name], sd, E4, hE4, mods},
   sd = stabilityData[p]; E4 = fixedPoints[p]["E4"]; hE4 = sd["hE4"];
   Print["\n", name];
   Print["  T = ", NumberForm[sd["T"], 6], ",  D = ", NumberForm[sd["D"], 6],
     ",  T^2-4D = ", NumberForm[sd["T^2-4D"], 6], ",  T^2/D = ", NumberForm[sd["T^2/D"], 6]];
   Print["  h* = ", NumberForm[sd["hstar"], 6], ",  rho4 = ", NumberForm[sd["rho4"], 6],
     ",  2/rho4 = ", NumberForm[sd["hF"], 6], ",  h_E4 = ", NumberForm[hE4, 6]];
   Print["  E4 loses stability through: ", sd["route"]];
   mods[hh_] := Sort[Abs[Eigenvalues[jacG[E4, p, hh]]]];
   Print["  |eigenvalues| at 0.98 h_E4: ", NumberForm[mods[0.98 hE4], 5], "  (all < 1: sink)"];
   Print["  |eigenvalues| at 1.02 h_E4: ", NumberForm[mods[1.02 hE4], 5], "  (one > 1: unstable)"]],
  {name, Keys[examples]}];
