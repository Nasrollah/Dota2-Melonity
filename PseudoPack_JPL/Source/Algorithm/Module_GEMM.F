#if defined (LOCAL_GEMM)

#include "pseudopack.h"


! ------------------------------------------------------------------------
!  Purpose
!  =======
!
!  GEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!              TRANSA = 'T' or 't',  op( A ) = A'.
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!              TRANSB = 'T' or 't',  op( B ) = B'.
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  MAX( 1, m ), otherwise  LDA must be at
!           least  MAX( 1, k ).
!           Unchanged on exit.
!
!  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  MAX( 1, k ), otherwise  LDB must be at
!           least  MAX( 1, n ).
!           Unchanged on exit.
!
!  BETA   - REAL            .
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - REAL             array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           MAX( 1, m ).
!           Unchanged on exit.
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!  -- Modified on 7-7-96
!     Wai-Sun Don, Div. of Applied Math., Brown University
!       for the PseudoPack
!
! ------------------------------------------------------------------------


MODULE General_Matrix_Multiply

implicit NONE

PRIVATE

PUBLIC  :: GEMM
PUBLIC  :: L_Same
PUBLIC  :: XERBLA 

CONTAINS

  Subroutine GEMM (Trans_A, Trans_B, M, N, K, alpha, A,LDA, B,LDB, beta, C,LDC)

  integer  :: i, j, l, M, N, K, LDA, LDB, LDC
  integer  :: N_Col_A, N_Row_A, N_Row_B
  logical  :: Not_A, Not_B
  REALTYPE :: alpha, beta, Temp

  character(LEN=1) :: Trans_A, Trans_B

  REALTYPE, dimension(LDA,*) :: A
  REALTYPE, dimension(LDB,*) :: B
  REALTYPE, dimension(LDC,*) :: C

#ifdef DEBUG
  integer :: Info = 0
#endif

  Not_A = L_Same (Trans_A,'N')
  if (Not_A) then
    N_Row_A = M ; N_Col_A = K
  else
    N_Row_A = K ; N_Col_A = M
  endif

  Not_B = L_Same (Trans_B,'N')
  if (Not_B) then
    N_Row_B = K
  else
    N_Row_B = N
  endif

#ifdef DEBUG
  if ((.NOT. Not_A) .AND. (.NOT. L_Same (Trans_A,'C'))  &
                    .AND. (.NOT. L_Same (Trans_A,'T'))) Info = 1
  if (Info == 1) then ; write (6 ,101) ; write (99,101) ; endif
 101 format (1x/1x,'FATAL ERROR : GEMM' / &
                   '              Trans_A (N, C, T) = ',a/)

  if ((.NOT. Not_B) .AND. (.NOT. L_Same (Trans_B,'C'))  &
                    .AND. (.NOT. L_Same (Trans_B,'T'))) Info = 2
  if (Info == 2) then ; write (6 ,102) ; write (99,102) ; endif
 102 format (1x/1x,'FATAL ERROR : GEMM' / &
                   '              Trans_B (N, C, T) = ',a/)

  if (M <  0)                Info = 3 
  if (Info == 3) then ; write (6 ,103) M ; write (99,103) M ; endif
 103 format (1x/1x,'FATAL ERROR : GEMM' / &
                1x,'              M = ',i5,' < 0'/)

  if (N <  0)                Info = 4
  if (Info == 4) then ; write (6 ,104) N ; write (99,104) N ; endif
 104 format (1x/1x,'FATAL ERROR : GEMM' / &
                1x,'              N = ',i5,' < 0'/)

  if (K <  0)                Info = 5
  if (Info == 5) then ; write (6 ,105) K ; write (99,105) K ; endif
 105 format (1x/1x,'FATAL ERROR : GEMM' / &
                1x,'              K = ',i5,' < 0'/)

  if (LDA <  MAX(1,N_Row_A)) Info = 6
  if (Info == 6) then ; write (6 ,106) LDA, N_Row_A 
                        write (99,106) LDA, N_Row_A ; endif
 106 format (1x/1x,'FATAL ERROR : GEMM' / &
                1x,'              LDA, N_Row_A = ',2i5/ &
                1x,'              LDA < Max(1, N_Row_A)'/)

  if (LDB <  MAX(1,N_Row_B)) Info = 7
  if (Info == 7) then ; write (6 ,107) LDB, N_Row_B 
                        write (99,107) LDB, N_Row_B ; endif
 107 format (1x/1x,'FATAL ERROR : GEMM' / &
                1x,'              LDB, N_Row_B = ',2i5/ &
                1x,'              LDB < Max(1, N_Row_B)'/)

  if (LDC <  MAX(1,M      )) Info = 8
  if (Info == 8) then ; write (6 ,108) LDC, M
                        write (99,108) LDC, M ; endif
 108 format (1x/1x,'FATAL ERROR : GEMM' / &
                1x,'              LDC, M       = ',2i5/ &
                1x,'              LDC < Max(1, M)      '/)

  if (Info /= 0) STOP
#endif

  if ((M == 0) .OR. (N == 0) .OR.  &
      (((alpha == ZERO) .OR. (K == 0)) .AND. (beta == ONE)) ) RETURN
 
  if (alpha == ZERO) then  ! alpha = 0
    if (beta == ZERO) then
      do j = 1, N
        do i = 1, M
          C(i,j) = ZERO
        enddo
      enddo
    else
      do j = 1, N
        do i = 1, M
          C(i,j) = beta*C(i,j)
        enddo
      enddo
   endif

   RETURN
  endif
 
  if (Not_B) then
    if (Not_A) then
 
      do j = 1, N  !       Form  C := alpha*A*B + beta*C.

        if (beta == ZERO) then
          do i = 1, M
            C(i,j) = ZERO
          enddo
        else if (beta /= ONE) then
          do i = 1, M
            C(i,j) = beta*C(i,j)
          enddo
        endif

        do l = 1, K
          if (B(l,j) /= ZERO) then
            Temp = alpha*B(l,j)
            do i = 1, M
              C(i,j) = C(i,j) + Temp*A(i,l)
            enddo
          endif
        enddo

      enddo

    else

      do j = 1, N  !       Form  C := alpha*A'*B + beta*C
        do i = 1, M

          Temp = ZERO
          do l = 1, K
            Temp = Temp + A(l,i)*B(l,j)
          enddo

          if (beta == ZERO) then
            C(i,j) = alpha*Temp
          else
            C(i,j) = alpha*Temp + beta*C(i,j)
          endif

        enddo
      enddo

    endif

  else

    if (Not_A) then

      do j = 1, N  !       Form  C := alpha*A*B' + beta*C

        if (beta == ZERO) then
          do i = 1, M
            C(i,j) = ZERO
          enddo
        else if (beta /= ONE) then
          do i = 1, M
            C(i,j) = beta*C(i,j)
          enddo
        endif

        do l = 1, K
          if (B(j,l) /= ZERO) then
            Temp = alpha*B(j,l)
            do i = 1, M
              C(i,j) = C(i,j) + Temp*A(i,l)
            enddo
          endif
        enddo

      enddo

    else

      do j = 1, N  !       Form  C := alpha*A'*B' + beta*C
        do i = 1, M

          Temp = ZERO
          do l = 1, K
            Temp = Temp + A(l,i)*B(j,l)
          enddo

          if (beta == ZERO) then
            C(i,j) = alpha*Temp
          else
            C(i,j) = alpha*Temp + beta*C(i,j)
          endif

        enddo
      enddo

    endif

  endif
 
  END Subroutine GEMM 
!
! ==================================================================
!
  Logical FUNCTION L_Same (CA, CB)

  integer  :: Int_A, Int_B, Z_Code
  character :: CA, CB

  L_Same = CA == CB
  if (L_Same) RETURN

  Z_Code = ICHAR('Z')
  Int_A  = ICHAR(CA)
  Int_B  = ICHAR(CB)

  SELECT CASE (Z_Code)
    CASE (90 , 122)
      if (Int_A >= 97 .AND. Int_A <= 122) Int_A = Int_A - 32
      if (Int_B >= 97 .AND. Int_B <= 122) Int_B = Int_B - 32

    CASE (169, 233)
      if (Int_A >= 129 .AND. Int_A <= 137 .OR. &
          Int_A >= 145 .AND. Int_A <= 153 .OR. &
          Int_A >= 162 .AND. Int_A <= 169) Int_A = Int_A + 64

      if (Int_B >= 129 .AND. Int_B <= 137 .OR. &
          Int_B >= 145 .AND. Int_B <= 153 .OR. &
          Int_B >= 162 .AND. Int_B <= 169) Int_B = Int_B + 64

    CASE (218, 250)
      if (Int_A >= 225 .AND. Int_A <= 250) Int_A = Int_A - 32
      if (Int_B >= 225 .AND. Int_B <= 250) Int_B = Int_B - 32

  END SELECT

  L_Same = Int_A == Int_B

  END FUNCTION L_Same 
!
! ==================================================================
!
  Subroutine XERBLA (Subroutine_Name, Info)

  integer ::  Info
  character(LEN=*) :: Subroutine_Name

  write (6 ,100) TRIM(Subroutine_Name), Info
  write (99,100) TRIM(Subroutine_Name), Info

  STOP

 100 format (1x/1x,'FATAL ERROR : Illegal argument for Subroutine ',(a)/ &
                1x,'              parameter number : ',i2/)

  END Subroutine XERBLA 

END MODULE General_Matrix_Multiply

#else

MODULE General_Matrix_Multiply

implicit NONE

END MODULE General_Matrix_Multiply

#endif
