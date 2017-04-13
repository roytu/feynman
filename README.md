# feynman

Program that takes the integrals from QFT scattering amplitudes and simplifies them with some funny tricks.  Currently supports only QED.  WIP

Steps:
1. Write out the amplitude, using regularized vertices and propagators.
2. Use Feynman's trick to write all denominators as integrals.
3. Perform integrals over four-momenta
4. Take highest order (in the cutoff momenta)
5. 
