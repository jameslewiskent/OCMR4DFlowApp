The code is provided to reproduce the Atrial 4D Flow Analyis published by the authors.

The code was developed by Aaron T. Hess, Jack J. Miller and James L. Kent with input from Marco Spartera and Antonio Stracquadanio.


The code has been tested in Matlab R2021a.


Third party software downloaded from the internet was used extenstively in this project, these include:

unwrap
==============
https://github.com/mloecher/4dflow-lapunwrap

It was code that accompanied the paper 
https://onlinelibrary.wiley.com/doi/full/10.1002/jmri.25045

All credit goes to the authors of the paper 

Michael Loecher PhD  Eric Schrauben PhD  Kevin M. Johnson PhD  Oliver Wieben PhD



DivFreeWavelet
=====================
https://people.eecs.berkeley.edu/~mlustig/Software.html

Frank Ong, Martin Uecker, Umar Tariq, Albert Hsiao, Marcus T Alley, Shreyas S. Vasanawala and Michael Lustig , Robust 4D flow denoising using divergence-free wavelet transform, Magnetic Resonance in Medicine, 2014 Published on-line DOI: 10.1002/mrm.25176

see DivFreeWavelet/README.m

this may need compiling (or can be disabled in processing)



Rodgers_dcm_read
========================
this code was provided by Prof Christopher Rodgers
https://github.com/OXSAtoolbox/OXSA
https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0185356

dcm4che
====================
version dcm4che-2.0.24-bin
https://www.dcm4che.org/


XOM
================
by Elliotte Rusty Harold
see cmr42read\XOM\README.txt