#ifndef MACRO_DEFINITION_H
#define MACRO_DEFINITION_H

!===========================================================================
!                         Version Rio 2000
!
! AUTHOR  : Wai-Sun Don and Bruno Costa 
!         
!           Some useful macro definitions.
!
!===========================================================================
!                             Precision

!#if defined (SINGLE_PRECISION) || defined (CRAY)
#if defined (SINGLE_PRECISION) 
# define REALTYPE	REAL

# define CMPXTYPE	COMPLEX
# define CPLXTYPE	COMPLEX
# define CMLXTYPE	COMPLEX
# define CMPLTYPE	COMPLEX

# define ZERO		0.0e0
# define ONE		1.0e0
# define TWO		2.0e0
# define THREE		3.0e0
# define FOUR		4.0e0
# define FIVE		5.0e0
# define SIX		6.0e0
# define SEVEN		7.0e0
# define EIGHT		8.0e0
# define NINE 		9.0e0
# define TEN		10.0e0
# define ELEVEN		11.0e0
# define TWELVE		12.0e0
# define THIRDTEEN	13.0e0
# define SIXTEEN	16.0e0
# define NINTY  	90.0e0
# define HUNDRED	100.0e0
# define THOUSAND	1000.0e0
# define HALF		0.50e0
# define QUARTER	0.250e0
# define SIXTH		0.166666666666666666666e0
# define NINE_TEN	0.9e0
# define HUNDREDTH	0.01e0
# define SQRT_HALF	0.7071067811865476e0
# define PI		3.14159265358979323844e0
# define CIRCLE		360.0e0
# define HALF_CIRCLE	180.0e0

#else
# define REALTYPE	DOUBLE PRECISION

!# define CMPXTYPE	DOUBLE COMPLEX
!# define CPLXTYPE	DOUBLE COMPLEX
!# define CMLXTYPE	DOUBLE COMPLEX
!# define CMPLTYPE	DOUBLE COMPLEX

# define CMPXTYPE	COMPLEX(8) 
# define CPLXTYPE	COMPLEX(8) 
# define CMLXTYPE	COMPLEX(8)
# define CMPLTYPE	COMPLEX(8)

# define ZERO		0.0d0
# define ONE		1.0d0
# define TWO		2.0d0
# define THREE		3.0d0
# define FOUR		4.0d0
# define FIVE		5.0d0
# define SIX		6.0d0
# define SEVEN		7.0d0
# define EIGHT		8.0d0
# define NINE 		9.0d0
# define TEN		10.0d0
# define ELEVEN		11.0d0
# define TWELVE		12.0d0
# define THIRDTEEN	13.0d0
# define SIXTEEN	16.0d0
# define NINTY  	90.0d0
# define HUNDRED	100.0d0
# define THOUSAND	1000.0d0
# define HALF		0.50d0
# define QUARTER	0.250d0
# define SIXTH		0.166666666666666666666d0
# define NINE_TEN	0.9d0
# define HUNDREDTH	0.01d0
# define SQRT_HALF	0.7071067811865476d0
# define PI		3.14159265358979323844d0
# define CIRCLE		360.0d0
# define HALF_CIRCLE	180.0d0

#endif

#define R_2_D		HALF_CIRCLE/PI
#define D_2_R		PI/HALF_CIRCLE

#define ONE_MILLISECOND ONE_SECOND/THOUSAND
#define ONE_SECOND      1
#define ONE_MINUTE      60*ONE_SECOND
#define ONE_HOUR        60*ONE_MINUTE
#define ONE_DAY         24*ONE_HOUR
#define ONE_WEEK         7*ONE_DAY
#define ONE_YEAR        52*ONE_WEEK

#define CHAR_N(n)	ACHAR(48+(n))
#define CHAR_A(n)	ACHAR(64+(n))

#define RING		ACHAR(7)

!
!===========================================================================
!
!                           Parallel Routines (MPI, PESSL)

#if defined (PARALLEL_MPI)

# if defined (SINGLE_PRECISION)
#  define MPI_REALTYPE         MPI_REAL
# else
#  define MPI_REALTYPE         MPI_DOUBLE_PRECISION
# endif

#endif

#if defined (PARALLEL_MPI)
#  define TIMER                MPI_WTIME()
#else
#  define TIMER                CPU_TIME()
#endif

#endif 

#include "grid.h"
