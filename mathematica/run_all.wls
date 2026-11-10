#!/usr/bin/env wolframscript
(* run_all.wls -- runs the four Mathematica scripts in order.
   Usage:  wolframscript -file run_all.wls   (from this folder) *)
dir = DirectoryName[$InputFileName];
Do[
  Print["\n", StringRepeat["=", 70], "\nRunning ", f, "\n", StringRepeat["=", 70]];
  Get[FileNameJoin[{dir, f}]],
  {f, {"01_fixed_points.wl", "02_stability.wl", "03_bifurcation.wl", "04_chaos_control.wl"}}];
