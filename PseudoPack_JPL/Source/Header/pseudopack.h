#ifndef PSEUDOPACK_H
#define PSEUDOPACK_H

!===========================================================================
!                           Version Rio 2000
!
! AUTHOR  : Wai-Sun Don and Bruno Costa 
!
!       This file contains macros (and some structure definitions)
!       for the spectral derivative library project Version Rio.
!
!===========================================================================
!  n is the Total number of points 
!===========================================================================
!                 General Definition

#define HEADER_SIZE 		100
#define MAX_DERIVATIVE		4

#define METRIC_SIZE(N_Derv)	(((N_Derv)*((N_Derv)+1))/2)

#define _MD			MAX_DERIVATIVE
#define _MS			METRIC_SIZE(_MD)
 
#define _G(n)			_H(n+1)
#define _H(n)			(((n)+1)/2)
#define _Q(n)			((_G(n)+1)/2)
#define _R(n)			((_H(n)+1)/2)

!===========================================================================
!
#define _I_LIBRARY			1
#define _I_PRECISION			2
#define _I_CASE     			3
#define _I_RANK             		4 
#define _I_INDEX    			5
#define _I_METHOD			6
#define _I_POINT_TYPE			7
#define _I_MAX_ORDER_OF_D		8 
#define _I_ALGORITHM			9
#define _I_SYMMETRY			10
#define _I_N_PNTS			11
#define _I_N_VECTORS			12
#define _I_BLOCKSIZE			13
#define _I_X_LEFT 			14
#define _I_X_RIGHT			15
#define _I_MAP_F			16
#define _I_MAP				17
#define _I_MANUAL			18
#define _I_ALPHA			19
#define _I_BETA 			20
#define _I_SMOOTH_1			21
#define _I_SMOOTH_2			22
#define _I_EOD_CODE 			23
#define _I_EPSILON  			24
#define _I_STORAGE   			25
#define _I_MATRIX_DIMENSION		26
#define _I_MATRIX_START     		27
#define _I_ODD_MATRIX_START    		28
#define _I_SIGMA_START     		29
#define _I_METRIC_START     		30
#define _I_CFT_WSAVE_START  		31
#define _I_CFT_QWAVE_START  		32
#define _I_NUMBER_THREADS     		33
#define _I_MPI_DISTRIBUTED    		34
#define _I_N_REMAIN           		35
#define _I_N_0                		36
#define _I_N_STAR             		37
#define _I_N_BLOCK            		38
#define _I_NN                 		39
#define _I_END				HEADER_SIZE
 
!===========================================================================
!       Data Structure for Operator D, S, I

#define OPERATOR_LIBRARY(D)		int(D(_I_LIBRARY))
#define OPERATOR_PRECISION(D)		int(D(_I_PRECISION))
#define OPERATOR_CASE(D)		int(D(_I_CASE))
#define OPERATOR_RANK(D)	      	int(D(_I_RANK))
#define OPERATOR_INDEX(D)		int(D(_I_INDEX))
#define OPERATOR_METHOD(D)		int(D(_I_METHOD))
#define OPERATOR_POINT_TYPE(D)		int(D(_I_POINT_TYPE))
#define OPERATOR_MAX_ORDER_OF_D(D)	int(D(_I_MAX_ORDER_OF_D))
#define OPERATOR_ALGORITHM(D)		int(D(_I_ALGORITHM))
#define OPERATOR_SYMMETRY(D)		int(D(_I_SYMMETRY))
#define OPERATOR_N_PNTS(D)		int(D(_I_N_PNTS))
#define OPERATOR_N_VECTORS(D)		int(D(_I_N_VECTORS))
#define OPERATOR_BLOCKSIZE(D)		int(D(_I_BLOCKSIZE))
#define OPERATOR_X_LEFT(D)		   (D(_I_X_LEFT))
#define OPERATOR_X_RIGHT(D)		   (D(_I_X_RIGHT))
#define OPERATOR_MAP_F(D)		int(D(_I_MAP_F))
#define OPERATOR_MAP(D)			int(D(_I_MAP))
#define OPERATOR_MANUAL(D)		int(D(_I_MANUAL))
#define OPERATOR_ALPHA(D)		   (D(_I_ALPHA))
#define OPERATOR_BETA(D)		   (D(_I_BETA)) 
#define OPERATOR_SMOOTH_1(D)		int(D(_I_SMOOTH_1))
#define OPERATOR_SMOOTH_2(D)		int(D(_I_SMOOTH_2))
#define OPERATOR_EOD_CODE(D)		int(D(_I_EOD_CODE))
#define OPERATOR_EPSILON(D)		   (D(_I_EPSILON))
#define OPERATOR_STORAGE(D)		int(D(_I_STORAGE))
#define OPERATOR_MATRIX_DIMENSION(D)	int(D(_I_MATRIX_DIMENSION))
#define OPERATOR_MATRIX_START(D)	int(D(_I_MATRIX_START))
#define OPERATOR_ODD_MATRIX_START(D)	int(D(_I_ODD_MATRIX_START))
#define OPERATOR_SIGMA_START(D)		int(D(_I_SIGMA_START))
#define OPERATOR_METRIC_START(D)	int(D(_I_METRIC_START))
#define OPERATOR_CFT_WSAVE_START(D)	int(D(_I_CFT_WSAVE_START))
#define OPERATOR_CFT_QWAVE_START(D)	int(D(_I_CFT_QWAVE_START))
#define OPERATOR_MPI_DISTRIBUTED(D)	int(D(_I_MPI_DISTRIBUTED))
#define OPERATOR_N_REMAIN(D)		int(D(_I_N_REMAIN))    
#define OPERATOR_N_0(D)			int(D(_I_N_0))    
#define OPERATOR_N_STAR(D)		int(D(_I_N_STAR))    
#define OPERATOR_N_BLOCK(D)		int(D(_I_N_BLOCK))
#define OPERATOR_NN(D)			int(D(_I_NN))    
#define OPERATOR_NUMBER_THREADS(D)	int(D(_I_NUMBER_THREADS))
 
!===========================================================================
!                  Matrix-Matrix Multiply Algorithm

! For the Fourier/Chebyshev/Legendre Method, n = N_Pnts

#define MXM_FULL_SIZE(n)		(n)*(n)
#define MXM_AUX_FULL_SIZE(n)		HEADER_SIZE+_MD*MXM_FULL_SIZE(n)+_MS*(n)

#define MXM_MATRIX_START(n)		HEADER_SIZE                     +1
#define MXM_METRIC_START(n)		HEADER_SIZE+_MD*MXM_FULL_SIZE(n)+1

!===========================================================================
!                  Even-Odd Decomposition Algorithm

! For the Fourier            Method, n = N_Pnts+1
! For the Chebyshev/Legendre Method, n = N_Pnts

#define EOD_FULL_SIZE(n)		_H(n)*_H(n)*2
#define EOD_HALF_SIZE(n)		_H(n)*_H(n)

#define EOD_AUX_FULL_SIZE(n)		HEADER_SIZE+_MD*EOD_FULL_SIZE(n)+_MS*(n)
#define EOD_AUX_HALF_SIZE(n)		HEADER_SIZE+_MD*EOD_FULL_SIZE(n)+_MS*_H(n)

#define EOD_MATRIX_START_EVEN(n)	HEADER_SIZE+1
#define EOD_MATRIX_START_ODD(n)		HEADER_SIZE+_MD*EOD_HALF_SIZE(n)+1
#define EOD_METRIC_START(n)		HEADER_SIZE+_MD*EOD_FULL_SIZE(n)+1
 
!===========================================================================
!                Transform-Recursion-Transform Algorithm.

#if defined (ESSL)
!                   --- IBM ESSL Library Call ---

!                   --- Fast Cosine Transform Macro ---
# if defined (SINGLE_PRECISION)
#  define SINF				ssinf
#  define COSF				scosf
#  define RCFT				srcft
#  define CRFT				scrft
# else
#  define SINF				dsinf
#  define COSF				dcosf
#  define RCFT				drcft
#  define CRFT				dcrft
# endif

!                   --- Trig_Table Macros ---
# if defined (SINGLE_PRECISION)
#  define N_CFT_WSAVE(n)		40000
#  define N_CFT_WSAVE_FULL(n)		40000
#  define N_CFT_WSAVE_HALF(n)		40000
# else
#  define N_CFT_WSAVE(n)		35000
#  define N_CFT_WSAVE_FULL(n)		35000
#  define N_CFT_WSAVE_HALF(n)		35000
# endif

# if defined (SINGLE_PRECISION)
#  define N_CFT_QWAVE_FULL(n)		25000
#  define N_CFT_QWAVE_HALF(n)		25000
# else
#  define N_CFT_QWAVE_FULL(n)		22000
#  define N_CFT_QWAVE_HALF(n)		22000
# endif

#  define N_WSAVE_CFT_HALF(n)           N_CFT_WSAVE_HALF(n)
#  define N_QWAVE_CFT_HALF(n)           N_CFT_QWAVE_HALF(n)

#  define N_WSAVE_CFT(n)		N_CFT_WSAVE_HALF(n)
#  define N_QWAVE_CFT(n)		N_CFT_QWAVE_HALF(n)
#  define N_WSAVE_HALF(n)		N_WSAVE_CFT(n)+N_QWAVE_CFT(n)
#  define N_QWAVE_HALF(n)		_H(n)
#  define C_TABLE_HALF(n)		N_WSAVE_HALF(n)+_H(n)

#  define N_WSAVE_CFT_FULL(n)           N_CFT_WSAVE_FULL(n)
#  define N_WSAVE_FULL(n)		N_CFT_WSAVE_FULL(n)
#  define C_TABLE_FULL(n)		N_WSAVE_FULL(n)
#endif
 
!---------------------------------------------------------------------------
!                 --- NetLib VFFTLIB Call ---

#if defined (VFFT) 
!                   --- Trig_Table Macros ---
# define N_WSAVE_CFT(n)			_H(n)
# define N_WSAVE_HALF(n)		3*_H(n)+15
# define N_QWAVE_HALF(n)		3*_H(n)+15
# define C_TABLE_HALF(n)		2*(3*_H(n)+15)

# define N_WSAVE_FULL(n)		3*(n)+15
# define C_TABLE_FULL(n)		N_WSAVE_FULL(n)
#endif
 
!---------------------------------------------------------------------------
!                 --- Cray Scientific Library Call ---

#if defined (CRAY)
!                   --- Trig_Table Macros ---
# define N_WSAVE_CFT(n)			_R(n)
# define N_WSAVE_HALF(n)		_R(n)+2*_H(n)+19
# define N_QWAVE_HALF(n)		_H(n)
# define C_TABLE_HALF(n)		_R(n)+3*_H(n)+19

# define N_WSAVE_FULL(n)		_H(n)+2*(n)+19
# define C_TABLE_FULL(n)		N_WSAVE_FULL(n)
#endif
 
!---------------------------------------------------------------------------
!   --- Starting index for Trig Table, Sigma, and Metric Vectors ---
 
# define CFT_WSAVE_START_FULL(n)	HEADER_SIZE+1
# define CFT_SIGMA_START_FULL(n)	HEADER_SIZE+C_TABLE_FULL(n)+1
# define CFT_METRIC_START_FULL(n)	CFT_SIGMA_START_FULL(n)+(n)

# define CFT_WSAVE_START_HALF(n)	HEADER_SIZE+1
# define CFT_QSAVE_START_HALF(n)	HEADER_SIZE+N_WSAVE_CFT(n)+1
# define CFT_QWAVE_START_HALF(n)	HEADER_SIZE+N_WSAVE_HALF(n)+1
# define CFT_SIGMA_START_HALF(n)	HEADER_SIZE+C_TABLE_HALF(n)+1
# define CFT_METRIC_START_HALF(n)	CFT_SIGMA_START_HALF(n)+(n)

# define CFT_WSAVE_START(n)		CFT_WSAVE_START_FULL(n)
# define CFT_QWAVE_START(n)		CFT_QWAVE_START_HALF(n)
# define CFT_SIGMA_START(n)		CFT_SIGMA_START_FULL(n)
# define CFT_METRIC_START(n)		CFT_METRIC_START_FULL(n)
 
!---------------------------------------------------------------------------
!               --- Aux_Diff and Aux_Smth Macros ---

# define C_CFT_AUX_XXXX_FULL(n) 	HEADER_SIZE+C_TABLE_FULL(n)+(n)
# define C_CFT_AUX_XXXX_HALF(n) 	HEADER_SIZE+C_TABLE_HALF(n)+(n)

# define C_CFT_AUX_DIFF_FULL(n)		C_CFT_AUX_XXXX_FULL(n)+_MS*(n)
# define C_CFT_AUX_SMTH_FULL(n)		C_CFT_AUX_XXXX_FULL(n)

# define C_CFT_AUX_DIFF_HALF(n)		C_CFT_AUX_XXXX_HALF(n)+_MS*_H(n)
# define C_CFT_AUX_SMTH_HALF(n)		C_CFT_AUX_XXXX_HALF(n)

!===========================================================================
!                        FOURIER METHOD
!                Transform-ik-Transform Algorithm.

#if defined (ESSL)
!               --- IBM ESSL Library Call ---

!               --- Trig_Table Macros ---
# if defined (SINGLE_PRECISION)
#  define N_FFT_FWAVE(n)		25000
# else
#  define N_FFT_FWAVE(n)		22000
# endif

# if defined (SINGLE_PRECISION)
#  define N_FFT_CWAVE(n)		40000
#  define N_FFT_SWAVE(n)		60000
# else
#  define N_FFT_CWAVE(n)		35000
#  define N_FFT_SWAVE(n)		50000
# endif

#  define N_CWAVE(n)			N_FFT_CWAVE(n)
#  define N_SWAVE(n)			N_FFT_SWAVE(n)
#  define N_FWAVE_HALF(n)		N_CWAVE(n)+N_SWAVE(n)
#  define F_TABLE_HALF(n)		N_FWAVE_HALF(n)

#  define N_FWAVE(n)			N_FFT_FWAVE(n)
#  define N_BWAVE(n)			N_FFT_FWAVE(n)
#  define N_FWAVE_FULL(n)		N_FWAVE(n)+N_BWAVE(n)
#  define F_TABLE_FULL(n)		N_FWAVE_FULL(n)
#endif
 
!---------------------------------------------------------------------------
!                 --- NetLib VFFTLIB Call ---

#if defined (VFFT) 

!              --- Trig_Table Macros ---
# define N_CWAVE(n)			3*_G(n)+15
# define N_SWAVE(n)			2*_G(n)+_Q(n)+15
# define F_FWAVE_HALF(n)		N_CWAVE(n)+N_SWAVE(n)
# define F_TABLE_HALF(n)		5*_G(n)+_Q(n)+30

# define N_FWAVE(n)			(n)+15
# define N_FWAVE_FULL(n)		N_FWAVE(n)
# define F_TABLE_FULL(n)		N_FWAVE_FULL(n)
#endif
 
!---------------------------------------------------------------------------
!                 --- Cray Scientific Library Call ---

#if defined (CRAY)

!                   --- Trig_Table Macros ---
# define N_CWAVE(n)			_Q(n)+2*_G(n)+19
# define N_SWAVE(n)			_Q(n)+2*_G(n)+19
# define F_TABLE_HALF(n)		2*(_Q(n)+2*_G(n)+19)

# define N_FWAVE(n)			2*(n)+19
# define N_FWAVE_FULL(n)		N_FWAVE(n)
# define F_TABLE_FULL(n)		N_FWAVE_FULL(n)
#endif
 
!---------------------------------------------------------------------------
!   --- Starting index for Trig Table, Sigma, and Metric Vectors ---

# define FFT_FWAVE_START_FULL(n)	HEADER_SIZE+1
# define FFT_BWAVE_START_FULL(n)	HEADER_SIZE+N_FWAVE(n)+1
# define FFT_SIGMA_START_FULL(n)	HEADER_SIZE+F_TABLE_FULL(n)+1
# define FFT_METRIC_START_FULL(n)	FFT_SIGMA_START_FULL(n)+_G(n)

# define FFT_CWAVE_START_HALF(n)	HEADER_SIZE+1
# define FFT_SWAVE_START_HALF(n)	HEADER_SIZE+N_CWAVE(n)+1
# define FFT_SIGMA_START_HALF(n)	HEADER_SIZE+F_TABLE_HALF(n)+1
# define FFT_METRIC_START_HALF(n)	FFT_SIGMA_START_HALF(n)+_G(n)

!---------------------------------------------------------------------------
!             --- Aux_Diff and Aux_Smth Macros ---

# define F_FFT_AUX_XXXX_FULL(n) 	HEADER_SIZE+F_TABLE_FULL(n)+_G(n)
# define F_FFT_AUX_XXXX_HALF(n) 	HEADER_SIZE+F_TABLE_HALF(n)+_G(n)

# define F_FFT_AUX_DIFF_FULL(n)		F_FFT_AUX_XXXX_FULL(n)+_MS*(n)
# define F_FFT_AUX_SMTH_FULL(n)		F_FFT_AUX_XXXX_FULL(n)

# define F_FFT_AUX_DIFF_HALF(n)		F_FFT_AUX_XXXX_HALF(n)+_MS*_G(n)
# define F_FFT_AUX_SMTH_HALF(n)		F_FFT_AUX_XXXX_HALF(n)
 
!===========================================================================
!                       Macros for debugging

#if defined (DEBUG)
# ifndef DONT_USE_STOP
#  define ASSERT(cond, string) if (.not.(cond)) stop string
# else
#  define ASSERT(cond, string) if (.not.(cond)) write (99,*) __LINE__, string
# endif
#else
#  define ASSERT(cond, string)
#endif

!===========================================================================
!                             Precision

!#if defined (SINGLE_PRECISION) || defined (CRAY)
#if defined (SINGLE_PRECISION) 
# define REALTYPE	REAL
# define CMPXTYPE	COMPLEX 
# define CLPXTYPE	COMPLEX 
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

# define GEMV		SGEMV
# define GEMM		SGEMM
# define GEMUL		SGEMUL
# define GAMMA		Gamma
#else
# define REALTYPE	DOUBLE PRECISION
# define CMPXTYPE	COMPLEX(8) 
# define CLPXTYPE	COMPLEX(8)
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

# define GEMV		DGEMV
# define GEMM		DGEMM
# define GEMUL		DGEMUL
# define GAMMA		DGamma
#endif

#define R_2_D		HALF_CIRCLE/PI
#define D_2_R		PI/HALF_CIRCLE
#define TWOPI		TWO*PI

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

#if defined (PESSL) || defined (SCALAPACK)
# if defined (SINGLE_PRECISION)
# else
#  define PGEMM                PDGEMM
#  define PTRAN                PDTRAN
# endif
#endif

#if defined (PARALLEL_MPI)
#  define TIMER                MPI_WTIME()
#else
#  define TIMER                CPU_TIME()
#endif

#if defined (ESSL)
# if defined (SINGLE_PRECISION)
#  define RCFT2		       SRCFT2
#  define CRFT2		       SCRFT2
# else
#  define RCFT2		       DRCFT2
#  define CRFT2		       DCRFT2
# endif
#endif

#endif
