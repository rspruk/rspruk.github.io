splitsc 0.9.6 for Stata 15 or later. Time-split honest inference for synthetic controls with imperfect fit.

Files
  splitsc.ado           the command (ado program and Mata functions)
  splitsc.sthlp         help file
  splitsc_cert.do       certification against the reference values of the paper
  splitsc_examples.do   California, West Germany (with Theorem 3 and without 1990), Basque Country
  splitsc_graphs.do     all figures for the three applications, exported to figures/
  prop99.csv germany.csv basque.csv   the three application panels (unit, year, y)
  splitsc.pkg stata.toc  for net install

Install
  net install splitsc, from("C:/path/to/splitsc") replace
  net get splitsc, from("C:/path/to/splitsc")      // copies the do-files and data to the current folder
or copy splitsc.ado and splitsc.sthlp to your PERSONAL ado folder. Then: help splitsc

First run, from the folder with the do-files and the csv files
  do splitsc_cert.do       // every line should print "ok"
  do splitsc_examples.do
  do splitsc_graphs.do     // figures, written to the subfolder figures

Not implemented: the curvature-optimal anchor (needs a second-order cone programme) and the competing methods.
