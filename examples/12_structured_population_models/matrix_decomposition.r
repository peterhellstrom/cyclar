A.a <- 2
A.b <- 1
A.c <- 1
A.d <- 2 

A <- matrix(c(A.a,A.b,A.c,A.d),ncol=2,nrow=2,byrow=TRUE)

det(A)
trace(A)

eigen(A)

# Trace
trc <- sum(diag(A))
trc
A.a + A.d

# Determinant
det(A)
A.a*A.d - A.b*A.c

# Characteristic polynomial for a 2x2 matrix
# lambda^2 - lambda*alpha + beta
# alpha = trace
# beta = determinant

# Eigenvalues, 2x2 matrix
# Solve characteristic polynomial
(trc + sqrt(trc^2 - 4*det(A))) / 2
(trc - sqrt(trc^2 - 4*det(A))) / 2

# Characteristic polynomial, solve with polyroot
polyroot(c(det(A),-trc,1))
Re(polyroot(c(det(A),-trc,1)))
sort(Re(polyroot(c(det(A),-trc,1))), decreasing=TRUE)

z <- polyroot(c(det(A),-trc,1))
is.real(z)
is.complex(z)
Re(z)
Im(z)

# Find eigenvectors
# solve(a,b)
# a = coefficients of the system, matrix
# b = right hand side of the system (RHS)
solve(A,matrix(c(3,3),nrow=2))
