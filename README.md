# GI-LEARNING LIBRARY

`GI-LEARNING LIBRARY` is a C++ framework for [grammatical inference](https://en.wikipedia.org/wiki/Grammar_induction).

A full description is available in the following paper: [GI-learning: an optimized framework for grammatical
inference.](http://www.dicgim.unipa.it/networks/ndslab/pdf/0130.pdf), CompSysTech16.

Grammar induction algorithms:
- RPNI
- EDSM
- Blue*
- L*

Recently, it was added the W-method algorithm into the library, a tenicque to empirically evaluate similarity between formal languages. 

Easy configuration instructions:
```
1. Go inside the project folder "GI-learning"

2. make distclean
   autoreconf -i
   ./configure
   make
   cd src/
   ./giLearning ../examples/examples_big.txt ../examples/lstar.txt
```

If you obtain some error from the "configure" command, you must resolve it installing the necessary libraries or tools.
