(* ::Package:: *)
(* ============================================================================
   04_chaos_control.wl   --  Section 7 and the control parts of Section 8
   (a) symbolic check of Lemma 7.1: B G(X) + (1-B) X = X + (B h) F(X)
   (b) for each example: B* = h_E4/h, eigenvalue moduli of J_B(E4) at
       B* -/+ 0.02, and direct iteration of the controlled map
   Run:  wolframscript -file 04_chaos_control.wl
   ============================================================================ *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "00_model.wl"}]];

Print["\n=== (a) Lemma 7.1 ==="];
Clear[x, y, z, a1, b1, a2, b2, a3, b3, be, al, c, m, n, h, B];
ps = {a1, b1, a2, b2, a3, b3, be, al, c, m, n};
Print["B G(X) + (1-B) X == X + (B h) F(X): ",
  Simplify[B mapG[{x, y, z}, ps, h] + (1 - B) {x, y, z} - mapG[{x, y, z}, ps, B h]] === {0, 0, 0}];

controlledMap[X_, p_, h_, Bc_] := Bc mapG[X, p, h] + (1 - Bc) X;

Print["\n=== (b) Control of the three examples ==="];
controlCases = {{"Example 8.1", ex1, 1.40}, {"Example 8.2", ex2, 0.74}, {"Example 8.3", ex3, 0.91}};
Do[
  Module[{name = cs[[1]], p = cs[[2]], h0 = cs[[3]], E4, Bs, mods, tail},
   E4 = fixedPoints[p]["E4"];
   Bs = stabilityData[p]["hE4"]/h0;
   mods[Bc_] := Sort[Abs[Eigenvalues[(1 - Bc) IdentityMatrix[3] + Bc jacG[E4, p, h0]]]];
   tail[Bc_] := Module[{orb = NestList[controlledMap[#, p, h0, Bc] &, E4 + 0.01, 20000]},
      Max[Norm /@ ((# - E4) & /@ Take[orb, -100])]];
   Print["\n", name, "  (h = ", h0, ")"];
   Print["  B* = h_E4/h = ", NumberForm[Bs, 6]];
   Print["  B = B* - 0.02: |eig J_B(E4)| = ", NumberForm[mods[Bs - 0.02], 5],
     ",  distance of orbit from E4 after 20000 steps = ", ScientificForm[tail[Bs - 0.02], 3]];
   Print["  B = B* + 0.02: |eig J_B(E4)| = ", NumberForm[mods[Bs + 0.02], 5],
     ",  distance of orbit from E4 after 20000 steps = ", NumberForm[tail[Bs + 0.02], 4]]],
  {cs, controlCases}];
