# SVIME
A C++ library that implements Stochastic Variational Inference for Motif Elicitation (SVIME). It discovers an unbounded number of motifs over DNA sequences from FASTA files and produces their logos.

## Installation
Simply copy the source and header files to the src folder of your project.

## Dependencies
C++
* Boost
* Eigen
* OpenMP

Python3
* Matplotlib
* Numpy
* Pandas

## Tutorial
We'll find Oct4 binding motifs in DNA sequences overlapping Oct4 ChIP-seq peaks from [1]. Download the 'oct4\_sorted.fa' file from https://github.com/tahmidmehdi/svime/tree/master/data. This file contains the binding sites. 
1. Create a project with a source file and include the following files:
```python
#include "svime.h"
#include "distribution.h"
#include "util.h"
#include "asa103.hpp"
#include "processFasta.h"
#include <omp.h>
#include <boost/foreach.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <Eigen/StdVector>
```
2. In your main function, create a mapping of chromosomes to sizes (number of 15-mers in the chromosome based on the FASTA file) with the `faToMatrix` function. This stores sequences and genomic coordinates to .txt files in a specified output directory.
```python
std::map<std::string, int> chrSizes = faToMatrix("/path/to/oct4_sorted.fa", 15, "/path/to/output");
```
3. Create an array of parameters for the step-size function of SVI. The function is described in [2]. The first element of the array is the tau parameter and the second is the kappa parameter. Then, create a pointer to tau.
```python
float step[2] = {0, 0.5};
float* stepPtr = step;
```
4. Create a `svime` object. The arguments are described in the next section.
```python
svime model = svime(15, 1, 1, 20, stepPtr, 1000, 10, 4, 42);
```
5. Fit the model \& find motifs. Check /path/to/output/results for logos.
```python
svime::variationalDist q = model.fit_predict("/path/to/output", chrSizes, NULL);
```

## Classes
### svime
Implements SVI for Dirichlet Process Mixture of Product-Multinomials [3].

Argument | Data type | Description
:---: | :---: | :---
window | int | required. Length of motifs.
alpha | float | required. The alpha parameter for the Dirichlet Process. Determines how precisely the model should look for motifs. Higher values will create more motifs.
epochs | int | required. The maximum number of epochs.
max\_clusters | int | required. The maximum number of motifs the model can create.
step\_pars | float* | required. Array of parameters for step-size function.
batch\_size | int | optional (default: 1000). Number of window-mers in each batch.
tol | float | optional (default: 0.001). The algorithm stops when the difference between evidence lower bounds (ELBOs) in 2 consecutive iterations is less than tol.
n\_jobs | int | optional (default: 1). The number of threads to use.
random\_state | int | optional (default: 42). Determines the initial clusters and ensures reproducibility.

`fit_predict(outDir, chrSizes, hyperparameters = NULL)`

Argument | Data type | Description
:---: | :---: | :---
outDir | string | required. Output directory.
chrSizes | map<string, int> | required. A mapping of chromosomes to their sizes.
hyperparameters | psm* | optional (default: all concentrations are set to 1). A position score matrix (psm struct) of prior concentration parameters for each base and position.

## References
[1] Kopp, W. and Schulte-Sasse, R. (2017). Unsupervised learning of dna sequence features using a convolutional restricted boltzmann machine. bioRxiv.

[2] Hoffman, M. D. et al. (2013). Stochastic variational inference. The Journal of Machine Learning Research, 14(1), 1303-1347.

[3] Dunson, D. B. and Xing, C. (2009). Nonparametric bayes modeling of multivariate categorical data. Journal of the American Statistical Association, 104, 1042-1051.

