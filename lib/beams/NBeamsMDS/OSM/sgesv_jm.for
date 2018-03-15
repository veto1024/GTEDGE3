      subroutine SGESV_NR (n,nrhs,a,lda,ipiv,b,ldb,info)

      implicit none
      integer info, lda, ldb, n, nrhs
      integer ipiv(*)
      real a(lda,*), b(ldb,*), d

c/    This routine solves the system of linear equations with 
c/    multiple right hand sides of the form A*X = B where A is
c/    an n by n coefficient matrix, B is an n by nrhs matrix of
c/    right hand sides which, on output, holds the solution 
c/    matrix X. The routine is written in such a way as to be 
c/    interchangeable with the SGESV routine from the LAPACK 
c/    package, but it uses routines from Numerical Recipes.
c/    Written by John Mandrekas, GIT, 3-10-95

c/    Parameters:
c/    n     : Number of linear equations
c/    a     : On entry, the nxn coefficient matrix. On exit, a = LU
c/    lda   : Leading dimension of a
c/    ipiv  : permutation array (row i was interchanged with ipiv(i))
c/    b     : On entry, the n by nrhs matrix B. On exit, if info = 0,
c/            the solution matrix X
c/    info  : If 0 --> OK!

c/    Local variables
      integer j, ntty

      data ntty /6/

c/    First, call LUDCMP to calculate the LU decomposition of A:

      call ludcmp(a, n, lda, ipiv, d, info)

c/    Exit with message if LU failed:

      if (info.NE.0) then
         write (ntty,'(1x, A28)') 'The LU decomposition failed!'
         return
      endif

c/    Call the routine LUBKSB for each right hand side to perform
c/    the forward and back substitutions:

      do j = 1, nrhs
        call lubksb(a, n, lda, ipiv, b(1,j))
      enddo

      return
      end

      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL a(np,np),b(n)

c/    Solves the set of n linear equations A*X = B. Here A is input
c/    not as the matrix A but its LU decomposition determined by the
c/    routine ludcmp. This routine just does the forward and
c/    back substitutions. This routine does not modify a, n, np, and indx.

c/    Adapted from Numerical Recipes by John Mandrekas, GIT

c/    Parameters:
c/    a(1:n,1:n)   : The LU decomposition of the matrix A
c/    n            : The number of equations
c/    np           : The leading (physical) dimension of A
c/    indx(1:n)    : The permutation vector returned by LUDCMP
c/    b(1:n)       : On input, the RHS vector B. On output, the
c/                   solution vector X

      INTEGER i,ii,j,ll
      REAL sum

      ii=0
      do i = 1, n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do j = ii, i-1
            sum=sum-a(i,j)*b(j)
          enddo
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
      enddo

c/    Do the back substitution:

      do i = n, 1, -1
        sum=b(i)
        do j = i+1, n
          sum=sum-a(i,j)*b(j)
        enddo
        b(i)=sum/a(i,i)
      enddo

      return
      END


      SUBROUTINE ludcmp(a,n,np,indx,d,ifail)
      implicit none
      INTEGER n,np,indx(n),ifail,NMAX
      REAL d,a(np,np),TINY
      
c/    Largest expected n and a small number:
      PARAMETER (NMAX=500,TINY=1.0e-20)

c/    Adapted from "Numerical Recipes" by John Mandrekas, GIT
c/    Given a matrix A(1:n,1:n) with physical dimension np by np,
c/    this routine replaces it by the LU decomposition of a rowwise
c/    permutation of itself. 

c/    a     : On input, the n by n matrix. On output, the LU decomposition
c/            of a arranged in the usual way (see Eq. 2.3.14 of Numerical
c/            Recipes)
c/    n     : The order of the matrix a
c/    np    : The physical dimension of a
c/    indx  : The array indx(1:n) is an output vector that records the 
c/            row permutations effected by the partial pivoting.
c/    d     : +1 or -1 depending whether the number of permutations was
c/            even or odd.
c/    ifail : If ifail < 0, then we have a singular matrix
      
      INTEGER i,imax,j,k
      REAL aamax,dum,sum,vv(NMAX)

      ifail = 0
      d=1.

c/    Loop over rows to get the implicit scaling information:

      do i=1,n
        aamax=0.
        do j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
        enddo
        if (aamax.eq.0.) then  ! singular matrix in ludcmp
          ifail = -1
          return
        endif
        vv(i)=1./aamax
      enddo

c/    This is the loop over columns of Crout's method:
      do j=1,n
        do i=1,j-1
          sum=a(i,j)
          do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
        enddo

c/    Initialize for search of largest pivot element:

        aamax=0.
        do i=j,n
          sum=a(i,j)
          do k=1,j-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
        enddo

c/    Decide whether to interchange rows:

        if (j.ne.imax)then
          do k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          enddo
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do i=j+1,n
            a(i,j)=a(i,j)*dum
          enddo
        endif
      enddo
      return
      END
