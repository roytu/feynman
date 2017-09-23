# feynman

![Electron self-energy correction](http://jefferywinkler.com/beystandmod106.gif)

This software uses the electron self-energy correction as a demo, and the output of the calculation can be seen in [reports/report.pdf](reports/report.pdf).


# Overview

Program that takes the integrals from QFT scattering amplitudes and simplifies them with some funny tricks.  Currently supports only QED.  WIP

Steps:
1. Write out the amplitude, using regularized vertices and propagators.
2. Use Feynman's trick to write all denominators as integrals.
3. Perform integrals over four-momenta.
4. Take highest order (in the cutoff momenta).
5. Evaluate spin matrices.

