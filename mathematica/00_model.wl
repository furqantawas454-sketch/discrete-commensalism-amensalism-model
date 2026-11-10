(* ::Package:: *)
(* ============================================================================
   00_model.wl
   Shared definitions for the Euler-discretized three-species
   commensalism-amensalism model with Beddington-DeAngelis response.

   Discrete model (Section 3 of the paper):
     x_{k+1} = x + h [ x (a1 - b1 x) + beta z ]
     y_{k+1} = y + h y [ a2 - b2 y + c x/(m x + n y + 1) - alpha z ]
     z_{k+1} = z + h z [ a3/x - b3 z ]

   Parameter vector order everywhere: {a1, b1, a2, b2, a3, b3, beta, alpha, c, m, n}
   This file is loaded by every other script with Get[...].
   ============================================================================ *)

(* ---- the three parameter sets used in Section 8 ---- *)
ex1 = {1.2, 0.6, 0.5, 0.4, 0.9, 0.7, 0.3, 0.6, 2.0, 0.5, 0.4};   (* Example 8.1: flip in y-direction *)
ex2 = {3.0, 1.5, 1.0, 0.5, 1.5, 0.2, 0.2, 0.2, 1.0, 0.2, 0.2};   (* Example 8.2: flip in (x,z)-block *)
ex3 = {2.0, 3.0, 0.5, 0.2, 1.5, 3.0, 0.5, 0.5, 3.0, 1.0, 2.0};   (* Example 8.3: Neimark-Sacker   *)
examples = <|"Example 8.1" -> ex1, "Example 8.2" -> ex2, "Example 8.3" -> ex3|>;

(* ---- vector field of the continuous model (Eq. 3.1) and the Euler map (Eq. 3.2) ---- *)
fieldF[{x_, y_, z_}, {a1_, b1_, a2_, b2_, a3_, b3_, be_, al_, c_, m_, n_}] :=
  {x (a1 - b1 x) + be z,
   y (a2 - b2 y + c x/(m x + n y + 1) - al z),
   z (a3/x - b3 z)};

mapG[X_, p_, h_] := X + h fieldF[X, p];

(* ---- Jacobian of the map (Eq. 5.1) ---- *)
jacG[X_, p_, h_] := Module[{xx, yy, zz},
   D[mapG[{xx, yy, zz}, p, h], {{xx, yy, zz}}] /. Thread[{xx, yy, zz} -> X]];

(* ---- fixed points (Theorem 4.1) ---- *)
xStar[p_] := Module[{a1, b1, a3, b3, be, t},
   {a1, b1, a3, b3, be} = p[[{1, 2, 5, 6, 7}]];
   Max[t /. NSolve[t^3 - a1/b1 t^2 - a3 be/(b1 b3) == 0, t, Reals]]];

fixedPoints[p_] := Module[{a1, b1, a2, b2, a3, b3, be, al, c, m, n, x1, A2, A3, y2, xs, zs, L, Bq, Cq, Delta, ys},
   {a1, b1, a2, b2, a3, b3, be, al, c, m, n} = p;
   x1 = a1/b1;
   A2 = a1 m/(b1 n) + 1/n - a2/b2;
   A3 = -(a1 a2 m/(b1 b2 n) + a2/(b2 n) + a1 c/(b1 b2 n));
   y2 = (-A2 + Sqrt[A2^2 - 4 A3])/2;                               (* Eq. 4.4 *)
   xs = xStar[p]; zs = a3/(b3 xs);                                  (* Eq. 4.2, 4.5 *)
   L = a2 - al zs;
   Delta = L + c xs/(m xs + 1);                                     (* Eq. 4.3 *)
   Bq = b2 (m xs + 1) - n L; Cq = L (m xs + 1) + c xs;
   ys = If[Delta > 0, (-Bq + Sqrt[Bq^2 + 4 b2 n Cq])/(2 b2 n), None];  (* Eq. 4.6 *)
   <|"E1" -> {x1, 0, 0}, "E2" -> {x1, y2, 0}, "E3" -> {xs, 0, zs},
     "E4" -> If[ys === None, None, {xs, ys, zs}], "Delta" -> Delta|>];

(* ---- stability quantities (Section 5) ---- *)
stabilityData[p_] := Module[{a1, b1, a2, b2, a3, b3, be, al, c, m, n, fp, xs, ys, zs, T, Dd, S, rho4, hstar, h2},
   {a1, b1, a2, b2, a3, b3, be, al, c, m, n} = p;
   fp = fixedPoints[p]; {xs, ys, zs} = fp["E4"];
   T = -(b1 xs + be zs/xs + a3/xs);                                 (* Eq. 5.6 *)
   Dd = a3 b1 + 2 a3 be zs/xs^2;
   S = m xs + n ys + 1;
   rho4 = ys (b2 + c n xs/S^2);                                     (* Eq. 5.9 *)
   hstar = If[T^2 >= 4 Dd, (-T - Sqrt[T^2 - 4 Dd])/Dd, -T/Dd];      (* Eq. 5.8 *)
   h2 = If[T^2 >= 4 Dd, (-T + Sqrt[T^2 - 4 Dd])/Dd, None];
   <|"T" -> T, "D" -> Dd, "T^2-4D" -> T^2 - 4 Dd, "T^2/D" -> T^2/Dd,
     "hstar" -> hstar, "h2" -> h2, "rho4" -> rho4, "hF" -> 2/rho4,
     "hE4" -> Min[hstar, 2/rho4],
     "route" -> Which[2/rho4 < hstar, "flip in y-direction",
                      T^2 > 4 Dd, "flip in (x,z)-block", True, "Neimark-Sacker"]|>];

(* ---- multilinear forms of the planar (x,z)-map at (x*, z*)  (Eqs. 6.4-6.5) ---- *)
planarForms[p_, h_] := Module[{a1, b1, a3, b3, be, xs, zs, M, Bf, Cf},
   {a1, b1, a3, b3, be} = p[[{1, 2, 5, 6, 7}]];
   xs = xStar[p]; zs = a3/(b3 xs);
   M = IdentityMatrix[2] + h {{-(b1 xs + be zs/xs), be}, {-a3 zs/xs^2, -a3/xs}};
   Bf[u_, v_] := {-2 h b1 u[[1]] v[[1]],
      2 h a3 zs/xs^3 u[[1]] v[[1]] - h a3/xs^2 (u[[1]] v[[2]] + u[[2]] v[[1]]) - 2 h b3 u[[2]] v[[2]]};
   Cf[u_, v_, w_] := {0,
      -6 h a3 zs/xs^4 u[[1]] v[[1]] w[[1]] +
       2 h a3/xs^3 (u[[1]] v[[1]] w[[2]] + u[[1]] v[[2]] w[[1]] + u[[2]] v[[1]] w[[1]])};
   {M, Bf, Cf}];

Print["00_model.wl loaded: examples ", Keys[examples]];
