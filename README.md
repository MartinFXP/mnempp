# mnempp

Files:
------

mnem.cpp contains the code for the main function to be build into an executeable binary file.

functions.h contains the mnem algorithm and helper functions.

Install:
--------

Use your favourite compiler, e.g.,

```
g++ mnem.cpp -o mnem
```

Usage:
------

The data file needs to be tab separated and the first row/line corresponds
to the S-genes named from 0 to n-1 with n as the number of S-genes, e.g, for n=4

```
0	0	0	1	1	2	2	2	2	3	3
0.523	1.432	-5.3	9.4	1.5	-0.3	-7.6	9.1	3.5	-6.3	1.1
-4.5    1.2     -3.1    1.7     2.8     -8.9    -1.3    -5.3    3.6     0.53    -0.14
...
```

The columns do not have to be ordered by S-gene names.

If you do not have suitable data, you can use the compiled binary to simulate
data. Use '-s 1' flag (default: 0) to simulate data and '-sout filename' to
define the text file containing the ground truth networks and mixture weights.

```
mnem -h # for help

mnem -n 4 -m 100 -l 1000 -K 2 -s 1 -out data.txt -sout gtn.txt
```

Inference is done similar, but remember to place the data.txt in front of
the '-in' flag instead of the '-out' flag.

```
mnem -n 4 -m 100 -l 1000 -K 2 -in data.txt -out res.txt
```

Several runs are suggested, depending on noise and complexity of the ground truth. If
the number of components K is unknown, several different Ks are suggested. However, the
log likelihood cannot be compared directly. Some form of Information Criterion like BIC
or AIC (not implemented) is necessary to select the best K.

DISCLAIMER:
-----------

This program is fully functional and computes a mixture of directed graphs based on
the M&NEM algorithm (Pirkl & Beerenwinkel. 2018).

However, this program is written in basic C++ as a programming exercise and should
only be used for testing purposes (feedback is appreciated). For real data analysis
I recommend the more thoroughly tested R/Bioconductor package mnem (https://github.com/cbg-ethz/mnem).

## References:

Pirkl, M., Beerenwinkel, N.; Single cell network analysis with a mixture
of Nested Effects Models, Bioinformatics, Volume 34, Issue 17, 1 September
2018,
Pages i964-i971, https://doi.org/10.1093/bioinformatics/bty602.
