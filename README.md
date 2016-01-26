# ParSpMatVec.jl
Shared-memory implementation of parallel sparse matrix vector product in Julia

## Installation
To install on a unix machine, follow these steps
```
Pkg.clone("https://github.com/lruthotto/ParSpMatVec.jl")
Pkg.build("ParSpMatVec")
Pkg.test("ParSpMatVec")

```

## Usage
Currently, we do not overload the matrix vector product in `Base` (this might be added in the future). Let `A` be a sparse matrix, `alpha` and `beta` floating point numbers and `x` and `y` be real- or complex values vectors of appropriate size. Then, the following commands are equivalent 
```
nproc = 4; # choose number of OMP threads
yt = copy(y)
y = beta*y + alpha * A*x
ParSpMatVec.A_mul_B!( alpha, A, x, beta, yt, nproc)
```
Similarly, for the transpose matrix-vector product:
```
yt= copy(y)
y = beta*y + alpha * A'*x
ParSpMatVec.Ac_mul_B!( alpha, A, x, beta, yt, nproc)
```

The last input, `nproc`, determines how many OpenMP threads are used. Note that, due to the compressed column storage, products with the adjoint of `A` are expected to scale better. 

## ToDo

A few things to do:
- [ ] automatic build on Windows

