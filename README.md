# LU factorization

Test xGETRF (LU factorization) with vrious BLAS libraries.
All programs solve an NxN matrix and compute the numerical errors.

## Outputs 
```
1st col. : N
2nd col. : Rel. error defined as || Ax - b ||_2 / || b ||_2
3rd col. : HPL error defined as || Ax - b ||_oo / ( eps * ( || A ||_oo * || x ||_oo + || b ||_oo ) * n )
4th col. : Execution time
```
