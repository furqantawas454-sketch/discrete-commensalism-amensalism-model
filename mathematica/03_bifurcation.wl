(* ::Package:: *)
(* ============================================================================
   03_bifurcation.wl   --  Section 6 of the paper
   (a) symbolic check of the flip coefficients a, b in (6.7) and of b > 0
   (b) validation of the Neimark-Sacker formula (6.10) on the delayed
       logistic map, whose bifurcation at r = 2 is known to be supercritical
   (c) all conditions of Theorems 6.1-6.3 for Examples 8.1-8.3, and the
       subcritical case of Remark 9
   Run:  wolframscript -file 03_bifurcation.wl
   ============================================================================ *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "00_model.wl"}]];

(* ---------- normal form coefficients ---------- *)

(* flip in the y-direction, Theorem 6.1 / Eq. (6.7) *)
flipYcoef[p_] := Module[{a1, b1, a2, b2, a3, b3, be, al, c, m, n, xs, ys, zs, K, S, rho4, hF, a, b},
   {a1, b1, a2, b2, a3, b3, be, al, c, m, n} = p; {xs, ys, zs} = fixedPoints[p]["E4"];
   K = c xs; S = m xs + n ys + 1; rho4 = ys (b2 + K n/S^2); hF = 2/rho4;
   a = hF (-b2 - K n/S^2 + ys K n^2/S^3);
   b = hF (K n^2/S^3 - ys K n^3/S^4);
   <|"hF" -> hF, "a" -> a, "b" -> b, "lF" -> a^2 + b,
     "eig M at hF" -> Eigenvalues[planarForms[p, hF][[1]]], "dg'/dh" -> -rho4|>];

(* flip in the (x,z)-block, Theorem 6.2 / Eq. (6.8) *)
flipBcoef[p_] := Module[{sd, h1, M, Bf, Cf, q, pp, rho4, T, Dd},
   sd = stabilityData[p]; T = sd["T"]; Dd = sd["D"];
   h1 = (-T - Sqrt[T^2 - 4 Dd])/Dd; {M, Bf, Cf} = planarForms[p, h1];
   q = NullSpace[M + IdentityMatrix[2]][[1]];
   pp = NullSpace[Transpose[M] + IdentityMatrix[2]][[1]]; pp = pp/(pp . q);
   <|"h1" -> h1, "eig M at h1" -> Eigenvalues[M], "lambda_y at h1" -> 1 - h1 sd["rho4"],
     "mu1 (transversality)" -> (T - Sqrt[T^2 - 4 Dd])/2,
     "lB" -> 1/6 pp . Cf[q, q, q] - 1/2 pp . Bf[q, Inverse[M - IdentityMatrix[2]] . Bf[q, q]]|>];

(* first Lyapunov coefficient of a planar map at a Neimark-Sacker point, Eq. (6.10) *)
lyapunovNS[M_, Bf_, Cf_] := Module[{ev, es, ev2, es2, lam, q, p0, th, ip},
   {ev, es} = Eigensystem[M]; lam = First[Select[ev, Im[#] > 0 &]];
   q = es[[Position[ev, lam][[1, 1]]]];
   {ev2, es2} = Eigensystem[Transpose[M]];
   p0 = es2[[Position[ev2, Conjugate[lam]][[1, 1]]]];
   ip[u_, v_] := Conjugate[u] . v; p0 = p0/Conjugate[ip[p0, q]];
   th = Arg[lam];
   1/2 Re[Exp[-I th] (ip[p0, Cf[q, q, Conjugate[q]]] +
        2 ip[p0, Bf[q, Inverse[IdentityMatrix[2] - M] . Bf[q, Conjugate[q]]]] +
        ip[p0, Bf[Conjugate[q], Inverse[Exp[2 I th] IdentityMatrix[2] - M] . Bf[q, q]]])]];

nsCoef[p_] := Module[{sd, T, Dd, h, M, Bf, Cf, lam},
   sd = stabilityData[p]; T = sd["T"]; Dd = sd["D"]; h = -T/Dd;
   {M, Bf, Cf} = planarForms[p, h]; lam = Eigenvalues[M];
   <|"hNS" -> h, "eig M at hNS" -> lam, "|lambda|" -> Abs[lam], "theta0" -> Abs[Arg[lam[[1]]]],
     "T^2/D (must differ from 2 and 3)" -> T^2/Dd, "d|lambda|/dh = -T/2" -> -T/2,
     "lambda_y at hNS" -> 1 - h sd["rho4"], "lNS" -> lyapunovNS[M, Bf, Cf]|>];

(* ---------- (a) symbolic check of Eq. (6.7) ---------- *)
Print["\n=== (a) Symbolic check of the y-flip coefficients (6.7) ==="];
Clear[y, h, L0, K, Q, n, b2, ys];
phi[y_] := L0 - b2 y + K/(Q + n y);              (* phi in Eq. (6.1), Q = m x* + 1, K = c x* *)
g[y_] := y + h y phi[y];
relation = L0 -> b2 ys - K/(Q + n ys);            (* phi(y*) = 0 *)
S0 = Q + n ys;
Print["g'(y*) == 1 - h rho4: ",
  Simplify[(D[g[y], y] /. y -> ys /. relation) - (1 - h ys (b2 + K n/S0^2))] === 0];
Print["a == g''(y*)/2 matches (6.7): ",
  Simplify[(D[g[y], {y, 2}]/2 /. y -> ys) - h (-b2 - K n/S0^2 + ys K n^2/S0^3)] === 0];
Print["b == g'''(y*)/6 matches (6.7): ",
  Simplify[(D[g[y], {y, 3}]/6 /. y -> ys) - h (K n^2/S0^3 - ys K n^3/S0^4)] === 0];
Print["b == h K n^2 Q / S^4  (hence b > 0): ",
  Simplify[h (K n^2/S0^3 - ys K n^3/S0^4) - h K n^2 Q/S0^4] === 0];

(* ---------- (b) validation of the Neimark-Sacker formula ---------- *)
Print["\n=== (b) Validation of (6.10) on the delayed logistic map ==="];
Module[{r = 2, A0, B0, C0},
  A0 = N[D[{yy, r yy (1 - xx)}, {{xx, yy}}] /. {xx -> 1/2, yy -> 1/2}];
  B0[u_, v_] := {0, -r (u[[1]] v[[2]] + u[[2]] v[[1]])};
  C0[u_, v_, w_] := {0, 0};
  Print["(x,y) -> (y, 2 y (1 - x)) at (1/2,1/2): l_NS = ", lyapunovNS[A0, B0, C0],
    "   (negative = supercritical, as known)"]];

(* ---------- (c) the three examples ---------- *)
Print["\n=== (c) Conditions of Theorems 6.1-6.3 ==="];
Print["\nExample 8.1 (Theorem 6.1, flip in y-direction):"];
Print["  h* = ", stabilityData[ex1]["hstar"], "  (must exceed hF)"];
KeyValueMap[Print["  ", #1, " = ", #2] &, flipYcoef[ex1]];

Print["\nExample 8.2 (Theorem 6.2, flip in (x,z)-block):"];
Print["  hF = ", stabilityData[ex2]["hF"], "  (must exceed h1)"];
KeyValueMap[Print["  ", #1, " = ", #2] &, flipBcoef[ex2]];

Print["\nExample 8.3 (Theorem 6.3, Neimark-Sacker):"];
Print["  hF = ", stabilityData[ex3]["hF"], "  (must exceed hNS)"];
KeyValueMap[Print["  ", #1, " = ", #2] &, nsCoef[ex3]];

Print["\nRemark 9 (subcritical Neimark-Sacker):"];
exSub = {0.2, 3., 2., 3., 3., 3., 1., 0.2, 1.5, 1., 3.};
Print["  hF = ", stabilityData[exSub]["hF"]];
KeyValueMap[Print["  ", #1, " = ", #2] &, nsCoef[exSub]];
