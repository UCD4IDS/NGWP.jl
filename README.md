# NGWP.jl

<!-- [![Build Status](https://travis-ci.com/haotian127/NGWP.jl.svg?branch=master)](https://travis-ci.com/haotian127/NGWP.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/haotian127/NGWP.jl?svg=true)](https://ci.appveyor.com/project/haotian127/NGWP-jl)
[![Coverage](https://codecov.io/gh/haotian127/NGWP.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/haotian127/NGWP.jl)
[![Coverage](https://coveralls.io/repos/github/haotian127/NGWP.jl/badge.svg?branch=master)](https://coveralls.io/github/haotian127/NGWP.jl?branch=master) -->

This is the repository for the [Natural Graph Wavelet Packet Dictionaries](https://arxiv.org/abs/2009.09020).

## SETUP

To install the NGWP.jl (Natural Graph Wavelet Packets), run
```julia
]
(@v1.6) pkg> add https://gitlab.com/BoundaryValueProblems/MTSG.jl.git
(@v1.6) pkg> add https://github.com/haotian127/NGWP.jl.git
(@v1.6) pkg> test NGWP
using NGWP
```

## EXAMPLE

See [`test/paperscript`](https://github.com/haotian127/NGWP.jl/tree/master/test/paperscript) and [`dissertation/scripts`](https://github.com/haotian127/NGWP.jl/tree/master/dissertation/scripts).

## REFERENCES

1. H. Li and N. Saito, [Metrics of graph Laplacian eigenvectors](https://www.math.ucdavis.edu/~saito/publications/metgraphlap.html), in Wavelets and Sparsity XVIII (D. Van De Ville, M. Papadakis, and Y. M. Lu, eds.), Proc. SPIE 11138, Paper #111381K, 2019.
2. C. Alexander, H. Li and N. Saito, [Natural Graph Wavelet Packet Dictionaries](https://arxiv.org/abs/2009.09020), 2021.
