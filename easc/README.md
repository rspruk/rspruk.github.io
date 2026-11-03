# easc: exposure-augmented synthetic control for Stata

Implements the method of Spruk (2026), *Synthetic Control under Out-of-Span Common Shocks: Diagnosis, Exposure Balance and Correction*.

## Install

From Stata (requires version 16 or later):

    net install easc, from("https://rspruk.github.io/easc/") replace
    net get easc, from("https://rspruk.github.io/easc/")     // test_easc.do and the Iran data

or copy `easc.ado` and `easc.sthlp` into your PERSONAL ado directory (type `sysdir` to find it).

## Use

    help easc
    do test_easc.do       // reproduces the Iran application of the paper

## Citation

Spruk, R. (2026). Synthetic Control under Out-of-Span Common Shocks: Diagnosis, Exposure Balance and Correction. Working paper, University of Ljubljana.
