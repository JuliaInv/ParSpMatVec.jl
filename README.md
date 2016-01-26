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
yt = copy(y)
y = beta*y + alpha * A*x
ParSpMatVec.A_mul_B!( alpha, A, x, beta, yt, 4)
```
Similarly, for the transpose matrix-vector product:
```
yt= copy(y)
y = beta*y + alpha * A'*x
ParSpMatVec.Ac_mul_B!( alpha, A, x, beta, yt, 4 )
```

## ToDo

A few things to do:
-[ ] automatic build on Windows

