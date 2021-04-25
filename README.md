# NGWP.jl

<!-- [![Build Status](https://travis-ci.com/haotian127/NGWP.jl.svg?branch=master)](https://travis-ci.com/haotian127/NGWP.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/haotian127/NGWP.jl?svg=true)](https://ci.appveyor.com/project/haotian127/NGWP-jl)
[![Coverage](https://codecov.io/gh/haotian127/NGWP.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/haotian127/NGWP.jl)
[![Coverage](https://coveralls.io/repos/github/haotian127/NGWP.jl/badge.svg?branch=master)](https://coveralls.io/github/haotian127/NGWP.jl?branch=master) -->

**WARNING**: `NGWP.jl` is deprecated and will be part of [`MultiscaleGraphSignalTransforms.jl`](https://github.com/UCD4IDS/MultiscaleGraphSignalTransforms.jl) from April 24, 2021.

NGWP.jl is the repository for the [Natural Graph Wavelet Packet Dictionaries](https://link.springer.com/article/10.1007/s00041-021-09832-3).

<!-- ## SETUP

To install the NGWP.jl (Natural Graph Wavelet Packets), run
```julia
]
(@v1.6) pkg> add https://gitlab.com/UCD4IDS/MTSG.jl.git
(@v1.6) pkg> add https://github.com/UCD4IDS/NGWP.jl.git
(@v1.6) pkg> test NGWP
using NGWP
``` -->

## EXAMPLE

* The scripts for reproducing the figures in [2] can be found at `MultiscaleGraphSignalTransforms/test/paperscripts/NGWP_JFAA2021/scripts`.
* The scripts for reproducing the figures in [3] can be found at `MultiscaleGraphSignalTransforms/test/dissertations/htli/scripts`.

## REFERENCES

1. H. Li and N. Saito, [Metrics of graph Laplacian eigenvectors](https://www.math.ucdavis.edu/~saito/publications/metgraphlap.html), in Wavelets and Sparsity XVIII (D. Van De Ville, M. Papadakis, and Y. M. Lu, eds.), Proc. SPIE 11138, Paper #111381K, 2019.
2. C. Alexander, H. Li and N. Saito, [Natural graph wavelet packet dictionaries](https://link.springer.com/article/10.1007/s00041-021-09832-3), J. Fourier Anal. Appl., vol. 27, Article \#41, 2021.
3. H. Li, Natural Graph Wavelet Dictionaries: Methods and Applications, Ph.D. dissertation, University of California, Davis, May. 2021.
