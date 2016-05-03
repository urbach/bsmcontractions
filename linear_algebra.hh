// ********************



// linear_algebra.hh

// Author: Marc Wagner
// Date: September 2007



// ********************



#ifndef __LINEAR_ALGEBRA_HH__

#define __LINEAR_ALGEBRA_HH__



// ********************


#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


using namespace std;

// ********************



// A complex number.

typedef struct
{
  double re;
  double im;
} complex;

// c = 0.

inline void co_eq_zero(complex *c)
{
  c->re = 0.0;
  c->im = 0.0;
}

// c1 = c1 + c2.

inline void co_pl_eq_co(complex *c1, const complex *c2)
{
  c1->re += c2->re;
  c1->im += c2->im;
}

// c1 = c1 - c2.

inline void co_mi_eq_co(complex *c1, const complex *c2)
{
  c1->re -= c2->re;
  c1->im -= c2->im;
}

// c = c / d.

inline void co_di_eq_re(complex *c, double d)
{
  c->re /= d;
  c->im /= d;
}

// c1 = c2 * c3.

inline void co_eq_co_ti_co(complex *c1, const complex *c2, const complex *c3)
{
  c1->re = c2->re * c3->re - c2->im * c3->im;
  c1->im = c2->re * c3->im + c2->im * c3->re;
}

// c1 = c2^\dagger * c3.

inline void co_eq_co_dag_ti_co(complex *c1, const complex *c2, const complex *c3)
{
  c1->re = c2->re * c3->re + c2->im * c3->im;
  c1->im = c2->re * c3->im - c2->im * c3->re;
}

// c1 = c1 * c2.

inline void co_ti_eq_co(complex *c1, const complex *c2)
{
  double c1_re = c1->re * c2->re - c1->im * c2->im;
  c1->im = c1->re * c2->im + c1->im * c2->re;
  c1->re = c1_re;
}


// ********************



// The gamma matrices.

// Standard choice (as in gwc):

// gamma_0:

// |  0  0 -1  0 |
// |  0  0  0 -1 |
// | -1  0  0  0 |
// |  0 -1  0  0 |

// gamma_1:

// |  0  0  0 -i |
// |  0  0 -i  0 |
// |  0 +i  0  0 |
// | +i  0  0  0 |

// gamma_2:

// |  0  0  0 -1 |
// |  0  0 +1  0 |
// |  0 +1  0  0 |
// | -1  0  0  0 |

// gamma_3:

// |  0  0 -i  0 |
// |  0  0  0 +i |
// | +i  0  0  0 |
// |  0 -i  0  0 |

// Permutation of the eight elements of a spinor (re0, im0, re1, im1, ...).

const int gamma_permutation[16][24] = {

  {12, 13, 14, 15, 16, 17,   18, 19, 20, 21, 22, 23,    0,  1,  2,  3,  4,  5,    6,  7,  8,  9, 10, 11},  // gamma_0 (index 0)
  {19, 18, 21, 20, 23, 22,   13, 12, 15, 14, 17, 16,    7,  6,  9,  8, 11, 10,    1,  0,  3,  2,  5,  4},  // gamma_1 (index 1)
  {18, 19, 20, 21, 22, 23,   12, 13, 14, 15, 16, 17,    6,  7,  8,  9, 10, 11,    0,  1,  2,  3,  4,  5},  // gamma_2 (index 2)
  {13, 12, 15, 14, 17, 16,   19, 18, 21, 20, 23, 22,    1,  0,  3,  2,  5,  4,    7,  6,  9,  8, 11, 10},  // gamma_3 (index 3)

   {0,  1,  2,  3,  4,  5,    6,  7,  8,  9, 10, 11,   12, 13, 14, 15, 16, 17,   18, 19, 20, 21, 22, 23},  // id (index 4)
   {0,  1,  2,  3,  4,  5,    6,  7,  8,  9, 10, 11,   12, 13, 14, 15, 16, 17,   18, 19, 20, 21, 22, 23},  // gamma_5 (index 5)

  {12, 13, 14, 15, 16, 17,   18, 19, 20, 21, 22, 23,    0,  1,  2,  3,  4,  5,    6,  7,  8,  9, 10, 11},  // gamma_0 gamma_5 (index 6)
  {19, 18, 21, 20, 23, 22,   13, 12, 15, 14, 17, 16,    7,  6,  9,  8, 11, 10,    1,  0,  3,  2,  5,  4},  // gamma_1 gamma_5(index 7)
  {18, 19, 20, 21, 22, 23,   12, 13, 14, 15, 16, 17,    6,  7,  8,  9, 10, 11,    0,  1,  2,  3,  4,  5},  // gamma_2 gamma_5 (index 8)
  {13, 12, 15, 14, 17, 16,   19, 18, 21, 20, 23, 22,    1,  0,  3,  2,  5,  4,    7,  6,  9,  8, 11, 10},  // gamma_3 gamma_5 (index 9)

  {7, 6, 9, 8, 11, 10,   1, 0, 3, 2, 5, 4,   19, 18, 21, 20, 23, 22,   13, 12, 15, 14, 17, 16},  // sigma_01 (index 10)
  {6, 7, 8, 9, 10, 11,   0, 1, 2, 3, 4, 5,   18, 19, 20, 21, 22, 23,   12, 13, 14, 15, 16, 17},  // sigma_02 (index 11)
  {1, 0, 3, 2, 5, 4,   7, 6, 9, 8, 11, 10,   13, 12, 15, 14, 17, 16,   19, 18, 21, 20, 23, 22},  // sigma_03 (index 12)
  {1, 0, 3, 2, 5, 4,   7, 6, 9, 8, 11, 10,   13, 12, 15, 14, 17, 16,   19, 18, 21, 20, 23, 22},  // sigma_12 (index 13)
  {6, 7, 8, 9, 10, 11,   0, 1, 2, 3, 4, 5,   18, 19, 20, 21, 22, 23,   12, 13, 14, 15, 16, 17},  // sigma_13 (index 14)
  {7, 6, 9, 8, 11, 10,   1, 0, 3, 2, 5, 4,   19, 18, 21, 20, 23, 22,   13, 12, 15, 14, 17, 16}   // sigma_23 (index 15)

};

const int gamma_sign[16][24] = {

  {-1, -1, -1, -1, -1, -1,   -1, -1, -1, -1, -1, -1,   -1, -1, -1, -1, -1, -1,   -1, -1, -1, -1, -1, -1},  // gamma_0 (index 0)
  {+1, -1, +1, -1, +1, -1,   +1, -1, +1, -1, +1, -1,   -1, +1, -1, +1, -1, +1,   -1, +1, -1, +1, -1, +1},  // gamma_1 (index 1)
  {-1, -1, -1, -1, -1, -1,   +1, +1, +1, +1, +1, +1,   +1, +1, +1, +1, +1, +1,   -1, -1, -1, -1, -1, -1},  // gamma_2 (index 2)
  {+1, -1, +1, -1, +1, -1,   -1, +1, -1, +1, -1, +1,   -1, +1, -1, +1, -1, +1,   +1, -1, +1, -1, +1, -1},  // gamma_3 (index 3)

  {+1, +1, +1, +1, +1, +1,   +1, +1, +1, +1, +1, +1,   +1, +1, +1, +1, +1, +1,   +1, +1, +1, +1, +1, +1},  // id (index 4)
  {+1, +1, +1, +1, +1, +1,   +1, +1, +1, +1, +1, +1,   -1, -1, -1, -1, -1, -1,   -1, -1, -1, -1, -1, -1},  // gamma_5 (index 5)
  {+1, +1, +1, +1, +1, +1,   +1, +1, +1, +1, +1, +1,   -1, -1, -1, -1, -1, -1,   -1, -1, -1, -1, -1, -1},  // gamma_0 gamma_5 (index 6)
  {-1, +1, -1, +1, -1, +1,   -1, +1, -1, +1, -1, +1,   -1, +1, -1, +1, -1, +1,   -1, +1, -1, +1, -1, +1},  // gamma_1 gamma_5 (index 7)
  {+1, +1, +1, +1, +1, +1,   -1, -1, -1, -1, -1, -1,   +1, +1, +1, +1, +1, +1,   -1, -1, -1, -1, -1, -1},  // gamma_2 gamma_5 (index 8)
  {-1, +1, -1, +1, -1, +1,   +1, -1, +1, -1, +1, -1,   -1, +1, -1, +1, -1, +1,   +1, -1, +1, -1, +1, -1},  // gamma_3 gamma_5 (index 9)

  {+1, -1, +1, -1, +1, -1,   +1, -1, +1, -1, +1, -1,   -1, +1, -1, +1, -1, +1,   -1, +1, -1, +1, -1, +1},  // sigma_01 (index 10)
  {-1, -1, -1, -1, -1, -1,   +1, +1, +1, +1, +1, +1,   +1, +1, +1, +1, +1, +1,   -1, -1, -1, -1, -1, -1},  // sigma_02 (index 11)
  {+1, -1, +1, -1, +1, -1,   -1, +1, -1, +1, -1, +1,   -1, +1, -1, +1, -1, +1,   +1, -1, +1, -1, +1, -1},  // sigma_03 (index 12)
  {-1, +1, -1, +1, -1, +1,   +1, -1, +1, -1, +1, -1,   -1, +1, -1, +1, -1, +1,   +1, -1, +1, -1, +1, -1},  // sigma_12 (index 13)
  {-1, -1, -1, -1, -1, -1,   +1, +1, +1, +1, +1, +1,   -1, -1, -1, -1, -1, -1,   +1, +1, +1, +1, +1, +1},  // sigma_13 (index 14)
  {-1, +1, -1, +1, -1, +1,   -1, +1, -1, +1, -1, +1,   -1, +1, -1, +1, -1, +1,   -1, +1, -1, +1, -1, +1}   // sigma_23 (index 15)

};


//UKQCD gammas
// gamma_0 diagonal

// | +1  0  0  0 |
// |  0 +1  0  0 |
// |  0  0 -1  0 |
// |  0  0  0 -1 |


// gamma_5 =  g1 g2 g3 g0:

// |  0  0 +1  0 |
// |  0  0  0 +1 |
// | +1  0  0  0 |
// |  0 +1  0  0 |

// gamma_1:

// |  0  0  0 +i |
// |  0  0 +i  0 |
// |  0 -i  0  0 |
// | -i  0  0  0 |

// gamma_2:

// |  0  0  0 +1 |
// |  0  0 -1  0 |
// |  0 -1  0  0 |
// | +1  0  0  0 |

// gamma_3:

// |  0  0 +i  0 |
// |  0  0  0 -i |
// | -i  0  0  0 |
// |  0 +i  0  0 |

const int ukqcd_permutation[6][24] = {

  { 0,  1,  2,  3,  4,  5,    6,  7,  8,  9, 10, 11,   12, 13, 14, 15, 16, 17,   18, 19, 20, 21, 22, 23},  // gamma_0 (index 0)

  {19, 18, 21, 20, 23, 22,   13, 12, 15, 14, 17, 16,    7,  6,  9,  8, 11, 10,    1,  0,  3,  2,  5,  4},  // gamma_1 (index 1)
  {18, 19, 20, 21, 22, 23,   12, 13, 14, 15, 16, 17,    6,  7,  8,  9, 10, 11,    0,  1,  2,  3,  4,  5},  // gamma_2 (index 2)
  {13, 12, 15, 14, 17, 16,   19, 18, 21, 20, 23, 22,    1,  0,  3,  2,  5,  4,    7,  6,  9,  8, 11, 10},  // gamma_3 (index 3)

  { 0,  1,  2,  3,  4,  5,    6,  7,  8,  9, 10, 11,   12, 13, 14, 15, 16, 17,   18, 19, 20, 21, 22, 23},  // id (index 4)

  {12, 13, 14, 15, 16, 17,   18, 19, 20, 21, 22, 23,    0,  1,  2,  3,  4,  5,    6,  7,  8,  9, 10, 11}  // gamma_5 (index 5)

};

const int ukqcd_sign[6][24] = {
  {+1, +1, +1, +1, +1, +1,   +1, +1, +1, +1, +1, +1,   -1, -1, -1, -1, -1, -1,   -1, -1, -1, -1, -1, -1},  // gamma_0 (index 5)

  {-1, +1, -1, +1, -1, +1,   -1, +1, -1, +1, -1, +1,   +1, -1, +1, -1, +1, -1,   +1, -1, +1, -1, +1, -1},  // gamma_1 (index 1)
  {+1, +1, +1, +1, +1, +1,   -1, -1, -1, -1, -1, -1,   -1, -1, -1, -1, -1, -1,   +1, +1, +1, +1, +1, +1},  // gamma_2 (index 2)
  {-1, +1, -1, +1, -1, +1,   +1, -1, +1, -1, +1, -1,   +1, -1, +1, -1, +1, -1,   -1, +1, -1, +1, -1, +1},  // gamma_3 (index 3)

  {+1, +1, +1, +1, +1, +1,   +1, +1, +1, +1, +1, +1,   +1, +1, +1, +1, +1, +1,   +1, +1, +1, +1, +1, +1},  // id (index 4)
  {+1, +1, +1, +1, +1, +1,   +1, +1, +1, +1, +1, +1,   +1, +1, +1, +1, +1, +1,   +1, +1, +1, +1, +1, +1}  // gamma_5 (index 0)

};



// ********************



// A = 0.
inline void cm_eq_zero(double *A);

// A = 1.
inline void cm_eq_id(double *A);

// A = B.
inline void cm_eq_cm(double *A, const double *B);

// A = B^\dagger.
inline void cm_eq_cm_dag(double *A, const double *B);

// c = det(A).
inline void co_eq_det_cm(complex *c, const double *A);

// c = tr(A).
inline void co_eq_tr_cm(complex *c, const double *A);

// A = Proj_SU3(A).
inline void cm_proj(double *A);

// Prints A.
inline void cm_fprintf(const double *A, FILE *file);

// A = A + B.
inline void cm_pl_eq_cm(double *A, const double *B);

// A = B * d.
inline void cm_eq_cm_ti_re(double *A, const double *B, double d);

// A = A * d.
inline void cm_ti_eq_re(double *A, double d);

// A = B * c.
inline void cm_eq_cm_ti_co(double *A, const double *B, const complex *c);

// A = B * C.
inline void cm_eq_cm_ti_cm(double *A, const double *B, const double *C);

// A = B^\dagger * C.
inline void cm_eq_cm_dag_ti_cm(double *A, const double *B, const double *C);

// A = B * C^\dagger.
inline void cm_eq_cm_ti_cm_dag(double *A, const double *B, const double *C);

// A = B^\dagger * C^\dagger.
inline void cm_eq_cm_dag_ti_cm_dag(double *A, const double *B, const double *C);

// A += B * C.
inline void cm_pl_eq_cm_ti_cm(double *A, const double *B, const double *C);

// A += B^\dagger * C.
inline void cm_pl_eq_cm_dag_ti_cm(double *A, const double *B, const double *C);

// A += B * C^\dagger.
inline void cm_pl_eq_cm_ti_cm_dag(double *A, const double *B, const double *C);

// A += B^\dagger * C^\dagger.
inline void cm_pl_eq_cm_dag_ti_cm_dag(double *A, const double *B, const double *C);

// A^{a b} = t_A^a s_A^b.
inline void cm_eq_fv_ti_fv(double *A, const double *s, const double *t);

// c = \epsilon^{a b c} \epsilon^{d e f} A^{a d} B^{b c} C^{e f}.
inline void co_eq_eps_eps_cm_cm_cm(complex *c, const double *A, const double *B, const double *C);



// ********************



// Prints a gamma matrix.
inline void gamma_fprintf(int gamma_index, FILE *file);
// same, but for permutation array per and sign array sig
inline void gamma_fprintf(FILE *file, int const per[], int const sig[]);


// ********************



// s = 0.
inline void fv_eq_zero(double *s);

// s = -s.
inline void fv_mi(double *s);

// s = t.
inline void fv_eq_fv(double *s, const double *t);

// s = t^\dagger.
inline void fv_eq_fv_dag(double *s, const double *t);

// Prints s.
inline void fv_fprintf(const double *s, FILE *file);

// c = s^dagger * t.
inline void co_eq_fv_dag_ti_fv(complex *c, const double *s, const double *t);

// c += s^dagger * t.
inline void co_pl_eq_fv_dag_ti_fv(complex *c, const double *s, const double *t);

// s = t + u.
inline void fv_eq_fv_pl_fv(double *s, const double *t, const double *u);
inline void fv_eq_fv_pl_fv(float *s, const float *t, const float *u);

// s = t - u.
inline void fv_eq_fv_mi_fv(double *s, const double *t, const double *u);

// s = s + t.
inline void fv_pl_eq_fv(double *s, const double *t);

// s = s - t.
inline void fv_mi_eq_fv(double *s, const double *t);

// s = t * d.
inline void fv_eq_fv_ti_re(double *s, const double *t, double d);

// s = s * d.
inline void fv_ti_eq_re(double *s, double d);
inline void fv_ti_eq_re(float *s, double d);

// s = t * i * d.
inline void fv_eq_fv_ti_im(double *s, const double *t, double d);

// s = s * i * d.
inline void fv_ti_eq_im(double *s, double d);

// s = t * c.
inline void fv_eq_fv_ti_co(double *s, const double *t, const complex *c);

// s = A * t.
inline void fv_eq_cm_ti_fv(double *s, const double *A, const double *t);

// s = A^\dagger * t.
inline void fv_eq_cm_dag_ti_fv(double *s, const double *A, const double *t);

// s = s + A * t.
inline void fv_pl_eq_cm_ti_fv(double *s, const double *A, const double *t);

// s = s - A * t.
inline void fv_mi_eq_cm_ti_fv(double *s, const double *A, const double *t);

// s = s + A^\dagger * t.
inline void fv_pl_eq_cm_dag_ti_fv(double *s, const double *A, const double *t);

// s = s - A^\dagger * t.
inline void fv_mi_eq_cm_dag_ti_fv(double *s, const double *A, const double *t);

// s = gamma * t.
inline void fv_eq_gamma_ti_fv(double *s, int gamma_index, const double *t);
inline void fv_eq_gamma_ti_fv(float *s, int gamma_index, const float *t);
// the same, but with permutation and sign array as parameters.
inline void fv_eq_gamma_ti_fv(double *s, const double *t, int const per[], int const sig[]);


// multiplies two gamma matrices
inline void gamma_eq_gamma_ti_gamma(int per[], int sig[], 
				    int const per1[], int const per2[], int const sig1[], int const sig2[]);

inline void i_times_gamma(int per[], int sig[]);
inline void minus_i_times_gamma(int per[], int sig[]);
inline void minus_gamma(int per[], int sig[]);
inline int rotate_etmc_ukqcd(double * const field, const int T, const int L);
inline int rotate_etmc_ukqcd(float * const field, const int T, const int L);

// ********************



// A = 0.

inline void cm_eq_zero(double *A)
{
  A[ 0] = 0.0;
  A[ 1] = 0.0;
  A[ 2] = 0.0;
  A[ 3] = 0.0;
  A[ 4] = 0.0;
  A[ 5] = 0.0;

  A[ 6] = 0.0;
  A[ 7] = 0.0;
  A[ 8] = 0.0;
  A[ 9] = 0.0;
  A[10] = 0.0;
  A[11] = 0.0;

  A[12] = 0.0;
  A[13] = 0.0;
  A[14] = 0.0;
  A[15] = 0.0;
  A[16] = 0.0;
  A[17] = 0.0;
}



// ********************



// A = 1.

inline void cm_eq_id(double *A)
{
  A[ 0] = 1.0;
  A[ 1] = 0.0;
  A[ 2] = 0.0;
  A[ 3] = 0.0;
  A[ 4] = 0.0;
  A[ 5] = 0.0;

  A[ 6] = 0.0;
  A[ 7] = 0.0;
  A[ 8] = 1.0;
  A[ 9] = 0.0;
  A[10] = 0.0;
  A[11] = 0.0;

  A[12] = 0.0;
  A[13] = 0.0;
  A[14] = 0.0;
  A[15] = 0.0;
  A[16] = 1.0;
  A[17] = 0.0;
}



// ********************



// A = B.

inline void cm_eq_cm(double *A, const double *B)
{
  A[ 0] = B[ 0];
  A[ 1] = B[ 1];
  A[ 2] = B[ 2];
  A[ 3] = B[ 3];
  A[ 4] = B[ 4];
  A[ 5] = B[ 5];

  A[ 6] = B[ 6];
  A[ 7] = B[ 7];
  A[ 8] = B[ 8];
  A[ 9] = B[ 9];
  A[10] = B[10];
  A[11] = B[11];

  A[12] = B[12];
  A[13] = B[13];
  A[14] = B[14];
  A[15] = B[15];
  A[16] = B[16];
  A[17] = B[17];
}



// ********************



// A = B^\dagger.

inline void cm_eq_cm_dag(double *A, const double *B)
{
  A[ 0] =  B[ 0];
  A[ 1] = -B[ 1];
  A[ 2] =  B[ 6];
  A[ 3] = -B[ 7];
  A[ 4] =  B[12];
  A[ 5] = -B[13];

  A[ 6] =  B[ 2];
  A[ 7] = -B[ 3];
  A[ 8] =  B[ 8];
  A[ 9] = -B[ 9];
  A[10] =  B[14];
  A[11] = -B[15];

  A[12] =  B[ 4];
  A[13] = -B[ 5];
  A[14] =  B[10];
  A[15] = -B[11];
  A[16] =  B[16];
  A[17] = -B[17];
}



// ********************



// c = det(A).

inline void co_eq_det_cm(complex *c, const double *A)
{
  c->re =
    (A[ 0]*A[ 8]*A[16] - A[ 0]*A[ 9]*A[17] - A[ 1]*A[ 8]*A[17] - A[ 1]*A[ 9]*A[16]) +
    (A[ 2]*A[10]*A[12] - A[ 2]*A[11]*A[13] - A[ 3]*A[10]*A[13] - A[ 3]*A[11]*A[12]) +
    (A[ 4]*A[ 6]*A[14] - A[ 4]*A[ 7]*A[15] - A[ 5]*A[ 6]*A[15] - A[ 5]*A[ 7]*A[14]) -
    (A[12]*A[ 8]*A[ 4] - A[12]*A[ 9]*A[ 5] - A[13]*A[ 8]*A[ 5] - A[13]*A[ 9]*A[ 4]) -
    (A[14]*A[10]*A[ 0] - A[14]*A[11]*A[ 1] - A[15]*A[10]*A[ 1] - A[15]*A[11]*A[ 0]) -
    (A[16]*A[ 6]*A[ 2] - A[16]*A[ 7]*A[ 3] - A[17]*A[ 6]*A[ 3] - A[17]*A[ 7]*A[ 2]);

  c->im =
    (A[ 1]*A[ 8]*A[16] + A[ 0]*A[ 9]*A[16] + A[ 0]*A[ 8]*A[17] - A[ 1]*A[ 9]*A[17]) +
    (A[ 3]*A[10]*A[12] + A[ 2]*A[11]*A[12] + A[ 2]*A[10]*A[13] - A[ 3]*A[11]*A[13]) +
    (A[ 5]*A[ 6]*A[14] + A[ 4]*A[ 7]*A[14] + A[ 4]*A[ 6]*A[15] - A[ 5]*A[ 7]*A[15]) -
    (A[13]*A[ 8]*A[ 4] + A[12]*A[ 9]*A[ 4] + A[12]*A[ 8]*A[ 5] - A[13]*A[ 9]*A[ 5]) -
    (A[15]*A[10]*A[ 0] + A[14]*A[11]*A[ 0] + A[14]*A[10]*A[ 1] - A[15]*A[11]*A[ 1]) -
    (A[17]*A[ 6]*A[ 2] + A[16]*A[ 7]*A[ 2] + A[16]*A[ 6]*A[ 3] - A[17]*A[ 7]*A[ 3]);
}



// ********************



// c = tr(A).

inline void co_eq_tr_cm(complex *c, const double *A)
{
  c->re = A[ 0] + A[ 8] + A[16];
  c->im = A[ 1] + A[ 9] + A[17];
}



// ********************



#ifdef F_
#define _F(s) s##_
#else
#define _F(s) s
#endif

extern "C"
{
  int _F(ilaenv)(int *ispec, char name[], char opts[], int *n1, int *n2, int *n3, int *n4);

  void _F(zheev)(char *jobz, char *uplo, int *n, double a[], int *lda, double w[], double work[], int *lwork, double *rwork, int *info);
}



// Computes the eigenvalues and the eigenvectors of a hermitian 3x3 matrix.

// n = size of the matrix
// M = hermitean 3x3 matrix (only the upper triangle is needed);

//   M[ 0] + i M[ 1]   M[ 6] + i M[ 7]   M[12] + i M[13]
//   xxx               M[ 8] + i M[ 9]   M[14] + i M[15]
//   xxx               xxx               M[16] + i M[17]

//   at the end of the computation the eigenvectors are stored here

//   v[0] = (lambda[ 0] + i lambda[ 1]   lambda[2] + i lambda[3]   ...
//   v[1] = (lambda[ 6] + i lambda[ 7]   ...
//   v[2] = (lambda[12] + i lambda[13]   ...

// lambda = eigenvalues sorted in ascending order

inline void EV_Hermitian_3x3_Matrix(double *M, double *lambda)
{
  int n = 3;

  //   int one = 1;
  //   int m_one = -1;

  // Call ilaenv_(... to determine the optimal size for work.

  // fprintf(stderr, "    ilaenv...\n");

  // !!!!!!!!!!
  // !!!!!!!!!!
  // !!!!!!!!!!
  // !!!!!!!!!!
  // !!!!!!!!!!
  int lwork = 102;
  // (2 + _F(ilaenv)(&one, "zhetrd", "VU", &n, &m_one, &m_one, &m_one)) * 3;
  // !!!!!!!!!!
  // !!!!!!!!!!
  // !!!!!!!!!!
  // !!!!!!!!!!
  // !!!!!!!!!!

  // fprintf(stderr, "    ilaenv...end (lwork = %d)\n", lwork);

  double *work = (double *)malloc(lwork * 2 * sizeof(double));

  double rwork[9];

  int info;

  _F(zheev)("V", "U", &n, M, &n, lambda, work, &lwork, rwork, &info);

  free(work);
}



// A = Proj_SU3(A).

// Projects a color matrix on SU(3).

// A' = A / sqrt(A^\dagger A)

// P_{SU(3)}(A) = A' / det(A')^{1/3}

inline void cm_proj(double *A)
{
  double d1;
  double M1[18], M2[18], M3[18];


  // Compute A^\dagger A.

  cm_eq_cm_dag_ti_cm(M1, A, A);


  // Compute the eigenvalues and the eigenvectors of A^\dagger A.

  // Transpose A^\dagger A (this is needed to call a Fortan/Lapack function to
  // compute the eigenvectors and eigenvalues).

  M1[ 1] = -M1[ 1];
  M1[ 3] = -M1[ 3];
  M1[ 5] = -M1[ 5];

  M1[ 7] = -M1[ 7];
  M1[ 9] = -M1[ 9];
  M1[11] = -M1[11];

  M1[13] = -M1[13];
  M1[15] = -M1[15];
  M1[17] = -M1[17];

  double lambda[3];

  EV_Hermitian_3x3_Matrix(M1, lambda);

  // fprintf(stderr, "lambda = (%+6.3f  , %+6.3f , %+6.3f).\n",
  //   lambda[0], lambda[1], lambda[2]);

  if(lambda[0] <= 0.000000000001 || lambda[1] <= 0.000000000001 ||
     lambda[2] <= 0.000000000001)
    {
      fprintf(stderr, "lambda = (%+6.3f  , %+6.3f , %+6.3f).\n",
	      lambda[0], lambda[1], lambda[2]);

      fprintf(stderr, "Error: inline void SU3_proj(...\n");
      exit(EXIT_FAILURE);
    }


  // Compute "T^\dagger".

  M1[ 1] = -M1[ 1];
  M1[ 3] = -M1[ 3];
  M1[ 5] = -M1[ 5];

  M1[ 7] = -M1[ 7];
  M1[ 9] = -M1[ 9];
  M1[11] = -M1[11];

  M1[13] = -M1[13];
  M1[15] = -M1[15];
  M1[17] = -M1[17];

  // SU3_fprintf(M1, stderr);


  // Compute "T D^{-1/2}".

  cm_eq_cm_dag(M2, M1);

  d1 = 1.0 / sqrt(lambda[0]);
  M2[ 0] *= d1;
  M2[ 1] *= d1;
  M2[ 6] *= d1;
  M2[ 7] *= d1;
  M2[12] *= d1;
  M2[13] *= d1;

  d1 = 1.0 / sqrt(lambda[1]);
  M2[ 2] *= d1;
  M2[ 3] *= d1;
  M2[ 8] *= d1;
  M2[ 9] *= d1;
  M2[14] *= d1;
  M2[15] *= d1;

  d1 = 1.0 / sqrt(lambda[2]);
  M2[ 4] *= d1;
  M2[ 5] *= d1;
  M2[10] *= d1;
  M2[11] *= d1;
  M2[16] *= d1;
  M2[17] *= d1;


  // Compute "T D^{-1/2} T^\dagger".

  cm_eq_cm_ti_cm(M3, M2, M1);

  // SU3_fprintf(M3, stderr);


  // Compute A'.

  cm_eq_cm_ti_cm(M1, A, M3);

  // SU3_fprintf(M1, stderr);


  // Divide by det(A')^{1/3}.

  complex det;
  co_eq_det_cm(&det, M1);

  double phi = atan2(det.im, det.re) / 3.0;

  complex de1_3_cc;
  de1_3_cc.re = +cos(phi);
  de1_3_cc.im = -sin(phi);

  cm_eq_cm_ti_co(A, M1, &de1_3_cc);

  // SU3_fprintf(A, stderr);
}



// ********************



// Prints A.

inline void cm_fprintf(const double *A, FILE *file)
{
  fprintf(file, "%+9.6f %+9.6f I   %+9.6f %+9.6f I   %+9.6f %+9.6f I\n",
	  A[ 0], A[ 1], A[ 2], A[ 3], A[ 4], A[ 5]);

  fprintf(file, "%+9.6f %+9.6f I   %+9.6f %+9.6f I   %+9.6f %+9.6f I\n",
	  A[ 6], A[ 7], A[ 8], A[ 9], A[10], A[11]);

  fprintf(file, "%+9.6f %+9.6f I   %+9.6f %+9.6f I   %+9.6f %+9.6f I\n",
	  A[12], A[13], A[14], A[15], A[16], A[17]);
}



// ********************



// A = A + B.

inline void cm_pl_eq_cm(double *A, const double *B)
{
  A[ 0] += B[ 0];
  A[ 1] += B[ 1];
  A[ 2] += B[ 2];
  A[ 3] += B[ 3];
  A[ 4] += B[ 4];
  A[ 5] += B[ 5];

  A[ 6] += B[ 6];
  A[ 7] += B[ 7];
  A[ 8] += B[ 8];
  A[ 9] += B[ 9];
  A[10] += B[10];
  A[11] += B[11];

  A[12] += B[12];
  A[13] += B[13];
  A[14] += B[14];
  A[15] += B[15];
  A[16] += B[16];
  A[17] += B[17];
}



// ********************



// A = B * d.

inline void cm_eq_cm_ti_re(double *A, const double *B, double d)
{
  A[ 0] = B[ 0] * d;
  A[ 1] = B[ 1] * d;
  A[ 2] = B[ 2] * d;
  A[ 3] = B[ 3] * d;
  A[ 4] = B[ 4] * d;
  A[ 5] = B[ 5] * d;

  A[ 6] = B[ 6] * d;
  A[ 7] = B[ 7] * d;
  A[ 8] = B[ 8] * d;
  A[ 9] = B[ 9] * d;
  A[10] = B[10] * d;
  A[11] = B[11] * d;

  A[12] = B[12] * d;
  A[13] = B[13] * d;
  A[14] = B[14] * d;
  A[15] = B[15] * d;
  A[16] = B[16] * d;
  A[17] = B[17] * d;
}



// ********************



// A = A * d.

inline void cm_ti_eq_re(double *A, double d)
{
  A[ 0] *= d;
  A[ 1] *= d;
  A[ 2] *= d;
  A[ 3] *= d;
  A[ 4] *= d;
  A[ 5] *= d;

  A[ 6] *= d;
  A[ 7] *= d;
  A[ 8] *= d;
  A[ 9] *= d;
  A[10] *= d;
  A[11] *= d;

  A[12] *= d;
  A[13] *= d;
  A[14] *= d;
  A[15] *= d;
  A[16] *= d;
  A[17] *= d;
}



// ********************



// A = B * c.

inline void cm_eq_cm_ti_co(double *A, const double *B, const complex *c)
{
  double re = c->re;
  double im = c->im;

  A[ 0] = B[ 0] * re - B[ 1] * im;
  A[ 1] = B[ 1] * re + B[ 0] * im;
  A[ 2] = B[ 2] * re - B[ 3] * im;
  A[ 3] = B[ 3] * re + B[ 2] * im;
  A[ 4] = B[ 4] * re - B[ 5] * im;
  A[ 5] = B[ 5] * re + B[ 4] * im;

  A[ 6] = B[ 6] * re - B[ 7] * im;
  A[ 7] = B[ 7] * re + B[ 6] * im;
  A[ 8] = B[ 8] * re - B[ 9] * im;
  A[ 9] = B[ 9] * re + B[ 8] * im;
  A[10] = B[10] * re - B[11] * im;
  A[11] = B[11] * re + B[10] * im;

  A[12] = B[12] * re - B[13] * im;
  A[13] = B[13] * re + B[12] * im;
  A[14] = B[14] * re - B[15] * im;
  A[15] = B[15] * re + B[14] * im;
  A[16] = B[16] * re - B[17] * im;
  A[17] = B[17] * re + B[16] * im;
}



// ********************



// A = B * C.

inline void cm_eq_cm_ti_cm(double *A, const double *B, const double *C)
{
  // A00 = B00*C00 + B01*C10 + B02*C20
  A[ 0] =
    B[ 0]*C[ 0] - B[ 1]*C[ 1] +
    B[ 2]*C[ 6] - B[ 3]*C[ 7] +
    B[ 4]*C[12] - B[ 5]*C[13];
  A[ 1] =
    B[ 0]*C[ 1] + B[ 1]*C[ 0] +
    B[ 2]*C[ 7] + B[ 3]*C[ 6] +
    B[ 4]*C[13] + B[ 5]*C[12];

  // A01 = B00*C01 + B01*C11 + B02*C21
  A[ 2] =
    B[ 0]*C[ 2] - B[ 1]*C[ 3] +
    B[ 2]*C[ 8] - B[ 3]*C[ 9] +
    B[ 4]*C[14] - B[ 5]*C[15];
  A[ 3] =
    B[ 0]*C[ 3] + B[ 1]*C[ 2] +
    B[ 2]*C[ 9] + B[ 3]*C[ 8] +
    B[ 4]*C[15] + B[ 5]*C[14];

  // A02 = B00*C02 + B01*C12 + B02*C22
  A[ 4] =
    B[ 0]*C[ 4] - B[ 1]*C[ 5] +
    B[ 2]*C[10] - B[ 3]*C[11] +
    B[ 4]*C[16] - B[ 5]*C[17];
  A[ 5] =
    B[ 0]*C[ 5] + B[ 1]*C[ 4] +
    B[ 2]*C[11] + B[ 3]*C[10] +
    B[ 4]*C[17] + B[ 5]*C[16];


  // A10 = B10*C00 + B11*C10 + B12*C20
  A[ 6] =
    B[ 6]*C[ 0] - B[ 7]*C[ 1] +
    B[ 8]*C[ 6] - B[ 9]*C[ 7] +
    B[10]*C[12] - B[11]*C[13];
  A[ 7] =
    B[ 6]*C[ 1] + B[ 7]*C[ 0] +
    B[ 8]*C[ 7] + B[ 9]*C[ 6] +
    B[10]*C[13] + B[11]*C[12];

  // A11 = B10*C01 + B11*C11 + B12*C21
  A[ 8] =
    B[ 6]*C[ 2] - B[ 7]*C[ 3] +
    B[ 8]*C[ 8] - B[ 9]*C[ 9] +
    B[10]*C[14] - B[11]*C[15];
  A[ 9] =
    B[ 6]*C[ 3] + B[ 7]*C[ 2] +
    B[ 8]*C[ 9] + B[ 9]*C[ 8] +
    B[10]*C[15] + B[11]*C[14];

  // A12 = B10*C02 + B11*C12 + B12*C22
  A[10] =
    B[ 6]*C[ 4] - B[ 7]*C[ 5] +
    B[ 8]*C[10] - B[ 9]*C[11] +
    B[10]*C[16] - B[11]*C[17];
  A[11] =
    B[ 6]*C[ 5] + B[ 7]*C[ 4] +
    B[ 8]*C[11] + B[ 9]*C[10] +
    B[10]*C[17] + B[11]*C[16];


  // A20 = B20*C00 + B21*C10 + B22*C20
  A[12] =
    B[12]*C[ 0] - B[13]*C[ 1] +
    B[14]*C[ 6] - B[15]*C[ 7] +
    B[16]*C[12] - B[17]*C[13];
  A[13] =
    B[12]*C[ 1] + B[13]*C[ 0] +
    B[14]*C[ 7] + B[15]*C[ 6] +
    B[16]*C[13] + B[17]*C[12];

  // A21 = B20*C01 + B21*C11 + B22*C21
  A[14] =
    B[12]*C[ 2] - B[13]*C[ 3] +
    B[14]*C[ 8] - B[15]*C[ 9] +
    B[16]*C[14] - B[17]*C[15];
  A[15] =
    B[12]*C[ 3] + B[13]*C[ 2] +
    B[14]*C[ 9] + B[15]*C[ 8] +
    B[16]*C[15] + B[17]*C[14];

  // A22 = B20*C02 + B21*C12 + B22*C22
  A[16] =
    B[12]*C[ 4] - B[13]*C[ 5] +
    B[14]*C[10] - B[15]*C[11] +
    B[16]*C[16] - B[17]*C[17];
  A[17] =
    B[12]*C[ 5] + B[13]*C[ 4] +
    B[14]*C[11] + B[15]*C[10] +
    B[16]*C[17] + B[17]*C[16];
}



// ********************

// A += B * C.

inline void cm_pl_eq_cm_ti_cm(double *A, const double *B, const double *C)
{
  // A00 += B00*C00 + B01*C10 + B02*C20
  A[ 0] +=
    B[ 0]*C[ 0] - B[ 1]*C[ 1] +
    B[ 2]*C[ 6] - B[ 3]*C[ 7] +
    B[ 4]*C[12] - B[ 5]*C[13];
  A[ 1] +=
    B[ 0]*C[ 1] + B[ 1]*C[ 0] +
    B[ 2]*C[ 7] + B[ 3]*C[ 6] +
    B[ 4]*C[13] + B[ 5]*C[12];

  // A01 += B00*C01 + B01*C11 + B02*C21
  A[ 2] +=
    B[ 0]*C[ 2] - B[ 1]*C[ 3] +
    B[ 2]*C[ 8] - B[ 3]*C[ 9] +
    B[ 4]*C[14] - B[ 5]*C[15];
  A[ 3] +=
    B[ 0]*C[ 3] + B[ 1]*C[ 2] +
    B[ 2]*C[ 9] + B[ 3]*C[ 8] +
    B[ 4]*C[15] + B[ 5]*C[14];

  // A02 += B00*C02 + B01*C12 + B02*C22
  A[ 4] +=
    B[ 0]*C[ 4] - B[ 1]*C[ 5] +
    B[ 2]*C[10] - B[ 3]*C[11] +
    B[ 4]*C[16] - B[ 5]*C[17];
  A[ 5] +=
    B[ 0]*C[ 5] + B[ 1]*C[ 4] +
    B[ 2]*C[11] + B[ 3]*C[10] +
    B[ 4]*C[17] + B[ 5]*C[16];


  // A10 += B10*C00 + B11*C10 + B12*C20
  A[ 6] +=
    B[ 6]*C[ 0] - B[ 7]*C[ 1] +
    B[ 8]*C[ 6] - B[ 9]*C[ 7] +
    B[10]*C[12] - B[11]*C[13];
  A[ 7] +=
    B[ 6]*C[ 1] + B[ 7]*C[ 0] +
    B[ 8]*C[ 7] + B[ 9]*C[ 6] +
    B[10]*C[13] + B[11]*C[12];

  // A11 += B10*C01 + B11*C11 + B12*C21
  A[ 8] +=
    B[ 6]*C[ 2] - B[ 7]*C[ 3] +
    B[ 8]*C[ 8] - B[ 9]*C[ 9] +
    B[10]*C[14] - B[11]*C[15];
  A[ 9] +=
    B[ 6]*C[ 3] + B[ 7]*C[ 2] +
    B[ 8]*C[ 9] + B[ 9]*C[ 8] +
    B[10]*C[15] + B[11]*C[14];

  // A12 += B10*C02 + B11*C12 + B12*C22
  A[10] +=
    B[ 6]*C[ 4] - B[ 7]*C[ 5] +
    B[ 8]*C[10] - B[ 9]*C[11] +
    B[10]*C[16] - B[11]*C[17];
  A[11] +=
    B[ 6]*C[ 5] + B[ 7]*C[ 4] +
    B[ 8]*C[11] + B[ 9]*C[10] +
    B[10]*C[17] + B[11]*C[16];


  // A20 += B20*C00 + B21*C10 + B22*C20
  A[12] +=
    B[12]*C[ 0] - B[13]*C[ 1] +
    B[14]*C[ 6] - B[15]*C[ 7] +
    B[16]*C[12] - B[17]*C[13];
  A[13] +=
    B[12]*C[ 1] + B[13]*C[ 0] +
    B[14]*C[ 7] + B[15]*C[ 6] +
    B[16]*C[13] + B[17]*C[12];

  // A21 += B20*C01 + B21*C11 + B22*C21
  A[14] +=
    B[12]*C[ 2] - B[13]*C[ 3] +
    B[14]*C[ 8] - B[15]*C[ 9] +
    B[16]*C[14] - B[17]*C[15];
  A[15] +=
    B[12]*C[ 3] + B[13]*C[ 2] +
    B[14]*C[ 9] + B[15]*C[ 8] +
    B[16]*C[15] + B[17]*C[14];

  // A22 += B20*C02 + B21*C12 + B22*C22
  A[16] +=
    B[12]*C[ 4] - B[13]*C[ 5] +
    B[14]*C[10] - B[15]*C[11] +
    B[16]*C[16] - B[17]*C[17];
  A[17] +=
    B[12]*C[ 5] + B[13]*C[ 4] +
    B[14]*C[11] + B[15]*C[10] +
    B[16]*C[17] + B[17]*C[16];
}



// ********************



// A = B^\dagger * C.

inline void cm_eq_cm_dag_ti_cm(double *A, const double *B, const double *C)
{
  // A00 = B00*C00 + B01*C10 + B02*C20
  A[ 0] =
    B[ 0]*C[ 0] + B[ 1]*C[ 1] +
    B[ 6]*C[ 6] + B[ 7]*C[ 7] +
    B[12]*C[12] + B[13]*C[13];
  A[ 1] =
    B[ 0]*C[ 1] - B[ 1]*C[ 0] +
    B[ 6]*C[ 7] - B[ 7]*C[ 6] +
    B[12]*C[13] - B[13]*C[12];

  // A01 = B00*C01 + B01*C11 + B02*C21
  A[ 2] =
    B[ 0]*C[ 2] + B[ 1]*C[ 3] +
    B[ 6]*C[ 8] + B[ 7]*C[ 9] +
    B[12]*C[14] + B[13]*C[15];
  A[ 3] =
    B[ 0]*C[ 3] - B[ 1]*C[ 2] +
    B[ 6]*C[ 9] - B[ 7]*C[ 8] +
    B[12]*C[15] - B[13]*C[14];

  // A02 = B00*C02 + B01*C12 + B02*C22
  A[ 4] =
    B[ 0]*C[ 4] + B[ 1]*C[ 5] +
    B[ 6]*C[10] + B[ 7]*C[11] +
    B[12]*C[16] + B[13]*C[17];
  A[ 5] =
    B[ 0]*C[ 5] - B[ 1]*C[ 4] +
    B[ 6]*C[11] - B[ 7]*C[10] +
    B[12]*C[17] - B[13]*C[16];


  // A10 = B10*C00 + B11*C10 + B12*C20
  A[ 6] =
    B[ 2]*C[ 0] + B[ 3]*C[ 1] +
    B[ 8]*C[ 6] + B[ 9]*C[ 7] +
    B[14]*C[12] + B[15]*C[13];
  A[ 7] =
    B[ 2]*C[ 1] - B[ 3]*C[ 0] +
    B[ 8]*C[ 7] - B[ 9]*C[ 6] +
    B[14]*C[13] - B[15]*C[12];

  // A11 = B10*C01 + B11*C11 + B12*C21
  A[ 8] =
    B[ 2]*C[ 2] + B[ 3]*C[ 3] +
    B[ 8]*C[ 8] + B[ 9]*C[ 9] +
    B[14]*C[14] + B[15]*C[15];
  A[ 9] =
    B[ 2]*C[ 3] - B[ 3]*C[ 2] +
    B[ 8]*C[ 9] - B[ 9]*C[ 8] +
    B[14]*C[15] - B[15]*C[14];

  // A12 = B10*C02 + B11*C12 + B12*C22
  A[10] =
    B[ 2]*C[ 4] + B[ 3]*C[ 5] +
    B[ 8]*C[10] + B[ 9]*C[11] +
    B[14]*C[16] + B[15]*C[17];
  A[11] =
    B[ 2]*C[ 5] - B[ 3]*C[ 4] +
    B[ 8]*C[11] - B[ 9]*C[10] +
    B[14]*C[17] - B[15]*C[16];


  // A20 = B20*C00 + B21*C10 + B22*C20
  A[12] =
    B[ 4]*C[ 0] + B[ 5]*C[ 1] +
    B[10]*C[ 6] + B[11]*C[ 7] +
    B[16]*C[12] + B[17]*C[13];
  A[13] =
    B[ 4]*C[ 1] - B[ 5]*C[ 0] +
    B[10]*C[ 7] - B[11]*C[ 6] +
    B[16]*C[13] - B[17]*C[12];

  // A21 = B20*C01 + B21*C11 + B22*C21
  A[14] =
    B[ 4]*C[ 2] + B[ 5]*C[ 3] +
    B[10]*C[ 8] + B[11]*C[ 9] +
    B[16]*C[14] + B[17]*C[15];
  A[15] =
    B[ 4]*C[ 3] - B[ 5]*C[ 2] +
    B[10]*C[ 9] - B[11]*C[ 8] +
    B[16]*C[15] - B[17]*C[14];

  // A22 = B20*C02 + B21*C12 + B22*C22
  A[16] =
    B[ 4]*C[ 4] + B[ 5]*C[ 5] +
    B[10]*C[10] + B[11]*C[11] +
    B[16]*C[16] + B[17]*C[17];
  A[17] =
    B[ 4]*C[ 5] - B[ 5]*C[ 4] +
    B[10]*C[11] - B[11]*C[10] +
    B[16]*C[17] - B[17]*C[16];
}


// ********************



// A += B^\dagger * C.

inline void cm_pl_eq_cm_dag_ti_cm(double *A, const double *B, const double *C)
{
  // A00 += B00*C00 + B01*C10 + B02*C20
  A[ 0] +=
    B[ 0]*C[ 0] + B[ 1]*C[ 1] +
    B[ 6]*C[ 6] + B[ 7]*C[ 7] +
    B[12]*C[12] + B[13]*C[13];
  A[ 1] +=
    B[ 0]*C[ 1] - B[ 1]*C[ 0] +
    B[ 6]*C[ 7] - B[ 7]*C[ 6] +
    B[12]*C[13] - B[13]*C[12];

  // A01 += B00*C01 + B01*C11 + B02*C21
  A[ 2] +=
    B[ 0]*C[ 2] + B[ 1]*C[ 3] +
    B[ 6]*C[ 8] + B[ 7]*C[ 9] +
    B[12]*C[14] + B[13]*C[15];
  A[ 3] +=
    B[ 0]*C[ 3] - B[ 1]*C[ 2] +
    B[ 6]*C[ 9] - B[ 7]*C[ 8] +
    B[12]*C[15] - B[13]*C[14];

  // A02 += B00*C02 + B01*C12 + B02*C22
  A[ 4] +=
    B[ 0]*C[ 4] + B[ 1]*C[ 5] +
    B[ 6]*C[10] + B[ 7]*C[11] +
    B[12]*C[16] + B[13]*C[17];
  A[ 5] +=
    B[ 0]*C[ 5] - B[ 1]*C[ 4] +
    B[ 6]*C[11] - B[ 7]*C[10] +
    B[12]*C[17] - B[13]*C[16];


  // A10 += B10*C00 + B11*C10 + B12*C20
  A[ 6] +=
    B[ 2]*C[ 0] + B[ 3]*C[ 1] +
    B[ 8]*C[ 6] + B[ 9]*C[ 7] +
    B[14]*C[12] + B[15]*C[13];
  A[ 7] +=
    B[ 2]*C[ 1] - B[ 3]*C[ 0] +
    B[ 8]*C[ 7] - B[ 9]*C[ 6] +
    B[14]*C[13] - B[15]*C[12];

  // A11 += B10*C01 + B11*C11 + B12*C21
  A[ 8] +=
    B[ 2]*C[ 2] + B[ 3]*C[ 3] +
    B[ 8]*C[ 8] + B[ 9]*C[ 9] +
    B[14]*C[14] + B[15]*C[15];
  A[ 9] +=
    B[ 2]*C[ 3] - B[ 3]*C[ 2] +
    B[ 8]*C[ 9] - B[ 9]*C[ 8] +
    B[14]*C[15] - B[15]*C[14];

  // A12 += B10*C02 + B11*C12 + B12*C22
  A[10] +=
    B[ 2]*C[ 4] + B[ 3]*C[ 5] +
    B[ 8]*C[10] + B[ 9]*C[11] +
    B[14]*C[16] + B[15]*C[17];
  A[11] +=
    B[ 2]*C[ 5] - B[ 3]*C[ 4] +
    B[ 8]*C[11] - B[ 9]*C[10] +
    B[14]*C[17] - B[15]*C[16];


  // A20 += B20*C00 + B21*C10 + B22*C20
  A[12] +=
    B[ 4]*C[ 0] + B[ 5]*C[ 1] +
    B[10]*C[ 6] + B[11]*C[ 7] +
    B[16]*C[12] + B[17]*C[13];
  A[13] +=
    B[ 4]*C[ 1] - B[ 5]*C[ 0] +
    B[10]*C[ 7] - B[11]*C[ 6] +
    B[16]*C[13] - B[17]*C[12];

  // A21 += B20*C01 + B21*C11 + B22*C21
  A[14] +=
    B[ 4]*C[ 2] + B[ 5]*C[ 3] +
    B[10]*C[ 8] + B[11]*C[ 9] +
    B[16]*C[14] + B[17]*C[15];
  A[15] +=
    B[ 4]*C[ 3] - B[ 5]*C[ 2] +
    B[10]*C[ 9] - B[11]*C[ 8] +
    B[16]*C[15] - B[17]*C[14];

  // A22 += B20*C02 + B21*C12 + B22*C22
  A[16] +=
    B[ 4]*C[ 4] + B[ 5]*C[ 5] +
    B[10]*C[10] + B[11]*C[11] +
    B[16]*C[16] + B[17]*C[17];
  A[17] +=
    B[ 4]*C[ 5] - B[ 5]*C[ 4] +
    B[10]*C[11] - B[11]*C[10] +
    B[16]*C[17] - B[17]*C[16];
}



// ********************



// A = B * C^\dagger.

inline void cm_eq_cm_ti_cm_dag(double *A, const double *B, const double *C)
{
  // A00 = B00*C00 + B01*C10 + B02*C20
  A[ 0] =
    B[ 0]*C[ 0] + B[ 1]*C[ 1] +
    B[ 2]*C[ 2] + B[ 3]*C[ 3] +
    B[ 4]*C[ 4] + B[ 5]*C[ 5];
  A[ 1] =
    -B[ 0]*C[ 1] + B[ 1]*C[ 0] +
    -B[ 2]*C[ 3] + B[ 3]*C[ 2] +
    -B[ 4]*C[ 5] + B[ 5]*C[ 4];

  // A01 = B00*C01 + B01*C11 + B02*C21
  A[ 2] =
    B[ 0]*C[ 6] + B[ 1]*C[ 7] +
    B[ 2]*C[ 8] + B[ 3]*C[ 9] +
    B[ 4]*C[10] + B[ 5]*C[11];
  A[ 3] =
    -B[ 0]*C[ 7] + B[ 1]*C[ 6] +
    -B[ 2]*C[ 9] + B[ 3]*C[ 8] +
    -B[ 4]*C[11] + B[ 5]*C[10];

  // A02 = B00*C02 + B01*C12 + B02*C22
  A[ 4] =
    B[ 0]*C[12] + B[ 1]*C[13] +
    B[ 2]*C[14] + B[ 3]*C[15] +
    B[ 4]*C[16] + B[ 5]*C[17];
  A[ 5] =
    -B[ 0]*C[13] + B[ 1]*C[12] +
    -B[ 2]*C[15] + B[ 3]*C[14] +
    -B[ 4]*C[17] + B[ 5]*C[16];


  // A10 = B10*C00 + B11*C10 + B12*C20
  A[ 6] =
    B[ 6]*C[ 0] + B[ 7]*C[ 1] +
    B[ 8]*C[ 2] + B[ 9]*C[ 3] +
    B[10]*C[ 4] + B[11]*C[ 5];
  A[ 7] =
    -B[ 6]*C[ 1] + B[ 7]*C[ 0] +
    -B[ 8]*C[ 3] + B[ 9]*C[ 2] +
    -B[10]*C[ 5] + B[11]*C[ 4];

  // A11 = B10*C01 + B11*C11 + B12*C21
  A[ 8] =
    B[ 6]*C[ 6] + B[ 7]*C[ 7] +
    B[ 8]*C[ 8] + B[ 9]*C[ 9] +
    B[10]*C[10] + B[11]*C[11];
  A[ 9] =
    -B[ 6]*C[ 7] + B[ 7]*C[ 6] +
    -B[ 8]*C[ 9] + B[ 9]*C[ 8] +
    -B[10]*C[11] + B[11]*C[10];

  // A12 = B10*C02 + B11*C12 + B12*C22
  A[10] =
    B[ 6]*C[12] + B[ 7]*C[13] +
    B[ 8]*C[14] + B[ 9]*C[15] +
    B[10]*C[16] + B[11]*C[17];
  A[11] =
    -B[ 6]*C[13] + B[ 7]*C[12] +
    -B[ 8]*C[15] + B[ 9]*C[14] +
    -B[10]*C[17] + B[11]*C[16];


  // A20 = B20*C00 + B21*C10 + B22*C20
  A[12] =
    B[12]*C[ 0] + B[13]*C[ 1] +
    B[14]*C[ 2] + B[15]*C[ 3] +
    B[16]*C[ 4] + B[17]*C[ 5];
  A[13] =
    -B[12]*C[ 1] + B[13]*C[ 0] +
    -B[14]*C[ 3] + B[15]*C[ 2] +
    -B[16]*C[ 5] + B[17]*C[ 4];

  // A21 = B20*C01 + B21*C11 + B22*C21
  A[14] =
    B[12]*C[ 6] + B[13]*C[ 7] +
    B[14]*C[ 8] + B[15]*C[ 9] +
    B[16]*C[10] + B[17]*C[11];
  A[15] =
    -B[12]*C[ 7] + B[13]*C[ 6] +
    -B[14]*C[ 9] + B[15]*C[ 8] +
    -B[16]*C[11] + B[17]*C[10];

  // A22 = B20*C02 + B21*C12 + B22*C22
  A[16] =
    B[12]*C[12] + B[13]*C[13] +
    B[14]*C[14] + B[15]*C[15] +
    B[16]*C[16] + B[17]*C[17];
  A[17] =
    -B[12]*C[13] + B[13]*C[12] +
    -B[14]*C[15] + B[15]*C[14] +
    -B[16]*C[17] + B[17]*C[16];
}

// ********************



// A += B * C^\dagger.

inline void cm_pl_eq_cm_ti_cm_dag(double *A, const double *B, const double *C)
{
  // A00 += B00*C00 + B01*C10 + B02*C20
  A[ 0] +=
    B[ 0]*C[ 0] + B[ 1]*C[ 1] +
    B[ 2]*C[ 2] + B[ 3]*C[ 3] +
    B[ 4]*C[ 4] + B[ 5]*C[ 5];
  A[ 1] +=
    -B[ 0]*C[ 1] + B[ 1]*C[ 0] +
    -B[ 2]*C[ 3] + B[ 3]*C[ 2] +
    -B[ 4]*C[ 5] + B[ 5]*C[ 4];

  // A01 += B00*C01 + B01*C11 + B02*C21
  A[ 2] +=
    B[ 0]*C[ 6] + B[ 1]*C[ 7] +
    B[ 2]*C[ 8] + B[ 3]*C[ 9] +
    B[ 4]*C[10] + B[ 5]*C[11];
  A[ 3] +=
    -B[ 0]*C[ 7] + B[ 1]*C[ 6] +
    -B[ 2]*C[ 9] + B[ 3]*C[ 8] +
    -B[ 4]*C[11] + B[ 5]*C[10];

  // A02 += B00*C02 + B01*C12 + B02*C22
  A[ 4] +=
    B[ 0]*C[12] + B[ 1]*C[13] +
    B[ 2]*C[14] + B[ 3]*C[15] +
    B[ 4]*C[16] + B[ 5]*C[17];
  A[ 5] +=
    -B[ 0]*C[13] + B[ 1]*C[12] +
    -B[ 2]*C[15] + B[ 3]*C[14] +
    -B[ 4]*C[17] + B[ 5]*C[16];


  // A10 += B10*C00 + B11*C10 + B12*C20
  A[ 6] +=
    B[ 6]*C[ 0] + B[ 7]*C[ 1] +
    B[ 8]*C[ 2] + B[ 9]*C[ 3] +
    B[10]*C[ 4] + B[11]*C[ 5];
  A[ 7] +=
    -B[ 6]*C[ 1] + B[ 7]*C[ 0] +
    -B[ 8]*C[ 3] + B[ 9]*C[ 2] +
    -B[10]*C[ 5] + B[11]*C[ 4];

  // A11 += B10*C01 + B11*C11 + B12*C21
  A[ 8] +=
    B[ 6]*C[ 6] + B[ 7]*C[ 7] +
    B[ 8]*C[ 8] + B[ 9]*C[ 9] +
    B[10]*C[10] + B[11]*C[11];
  A[ 9] +=
    -B[ 6]*C[ 7] + B[ 7]*C[ 6] +
    -B[ 8]*C[ 9] + B[ 9]*C[ 8] +
    -B[10]*C[11] + B[11]*C[10];

  // A12 += B10*C02 + B11*C12 + B12*C22
  A[10] +=
    B[ 6]*C[12] + B[ 7]*C[13] +
    B[ 8]*C[14] + B[ 9]*C[15] +
    B[10]*C[16] + B[11]*C[17];
  A[11] +=
    -B[ 6]*C[13] + B[ 7]*C[12] +
    -B[ 8]*C[15] + B[ 9]*C[14] +
    -B[10]*C[17] + B[11]*C[16];


  // A20 += B20*C00 + B21*C10 + B22*C20
  A[12] +=
    B[12]*C[ 0] + B[13]*C[ 1] +
    B[14]*C[ 2] + B[15]*C[ 3] +
    B[16]*C[ 4] + B[17]*C[ 5];
  A[13] +=
    -B[12]*C[ 1] + B[13]*C[ 0] +
    -B[14]*C[ 3] + B[15]*C[ 2] +
    -B[16]*C[ 5] + B[17]*C[ 4];

  // A21 += B20*C01 + B21*C11 + B22*C21
  A[14] +=
    B[12]*C[ 6] + B[13]*C[ 7] +
    B[14]*C[ 8] + B[15]*C[ 9] +
    B[16]*C[10] + B[17]*C[11];
  A[15] +=
    -B[12]*C[ 7] + B[13]*C[ 6] +
    -B[14]*C[ 9] + B[15]*C[ 8] +
    -B[16]*C[11] + B[17]*C[10];

  // A22 += B20*C02 + B21*C12 + B22*C22
  A[16] +=
    B[12]*C[12] + B[13]*C[13] +
    B[14]*C[14] + B[15]*C[15] +
    B[16]*C[16] + B[17]*C[17];
  A[17] +=
    -B[12]*C[13] + B[13]*C[12] +
    -B[14]*C[15] + B[15]*C[14] +
    -B[16]*C[17] + B[17]*C[16];
}



// ********************



// A = B^\dagger * C^\dagger.

inline void cm_eq_cm_dag_ti_cm_dag(double *A, const double *B, const double *C)
{
  // A00 = B00*C00 + B01*C10 + B02*C20
  A[ 0] =
    B[ 0]*C[ 0] - B[ 1]*C[ 1] +
    B[ 6]*C[ 2] - B[ 7]*C[ 3] +
    B[12]*C[ 4] - B[13]*C[ 5];
  A[ 1] =
    -B[ 0]*C[ 1] - B[ 1]*C[ 0] +
    -B[ 6]*C[ 3] - B[ 7]*C[ 2] +
    -B[12]*C[ 5] - B[13]*C[ 4];

  // A01 = B00*C01 + B01*C11 + B02*C21
  A[ 2] =
    B[ 0]*C[ 6] - B[ 1]*C[ 7] +
    B[ 6]*C[ 8] - B[ 7]*C[ 9] +
    B[12]*C[10] - B[13]*C[11];
  A[ 3] =
    -B[ 0]*C[ 7] - B[ 1]*C[ 6] +
    -B[ 6]*C[ 9] - B[ 7]*C[ 8] +
    -B[12]*C[11] - B[13]*C[10];

  // A02 = B00*C02 + B01*C12 + B02*C22
  A[ 4] =
    B[ 0]*C[12] - B[ 1]*C[13] +
    B[ 6]*C[14] - B[ 7]*C[15] +
    B[12]*C[16] - B[13]*C[17];
  A[ 5] =
    -B[ 0]*C[13] - B[ 1]*C[12] +
    -B[ 6]*C[15] - B[ 7]*C[14] +
    -B[12]*C[17] - B[13]*C[16];


  // A10 = B10*C00 + B11*C10 + B12*C20
  A[ 6] =
    B[ 2]*C[ 0] - B[ 3]*C[ 1] +
    B[ 8]*C[ 2] - B[ 9]*C[ 3] +
    B[14]*C[ 4] - B[15]*C[ 5];
  A[ 7] =
    -B[ 2]*C[ 1] - B[ 3]*C[ 0] +
    -B[ 8]*C[ 3] - B[ 9]*C[ 2] +
    -B[14]*C[ 5] - B[15]*C[ 4];

  // A11 = B10*C01 + B11*C11 + B12*C21
  A[ 8] =
    B[ 2]*C[ 6] - B[ 3]*C[ 7] +
    B[ 8]*C[ 8] - B[ 9]*C[ 9] +
    B[14]*C[10] - B[15]*C[11];
  A[ 9] =
    -B[ 2]*C[ 7] - B[ 3]*C[ 6] +
    -B[ 8]*C[ 9] - B[ 9]*C[ 8] +
    -B[14]*C[11] - B[15]*C[10];

  // A12 = B10*C02 + B11*C12 + B12*C22
  A[10] =
    B[ 2]*C[12] - B[ 3]*C[13] +
    B[ 8]*C[14] - B[ 9]*C[15] +
    B[14]*C[16] - B[15]*C[17];
  A[11] =
    -B[ 2]*C[13] - B[ 3]*C[12] +
    -B[ 8]*C[15] - B[ 9]*C[14] +
    -B[14]*C[17] - B[15]*C[16];


  // A20 = B20*C00 + B21*C10 + B22*C20
  A[12] =
    B[ 4]*C[ 0] - B[ 5]*C[ 1] +
    B[10]*C[ 2] - B[11]*C[ 3] +
    B[16]*C[ 4] - B[17]*C[ 5];
  A[13] =
    -B[ 4]*C[ 1] - B[ 5]*C[ 0] +
    -B[10]*C[ 3] - B[11]*C[ 2] +
    -B[16]*C[ 5] - B[17]*C[ 4];

  // A21 = B20*C01 + B21*C11 + B22*C21
  A[14] =
    B[ 4]*C[ 6] - B[ 5]*C[ 7] +
    B[10]*C[ 8] - B[11]*C[ 9] +
    B[16]*C[10] - B[17]*C[11];
  A[15] =
    -B[ 4]*C[ 7] - B[ 5]*C[ 6] +
    -B[10]*C[ 9] - B[11]*C[ 8] +
    -B[16]*C[11] - B[17]*C[10];

  // A22 = B20*C02 + B21*C12 + B22*C22
  A[16] =
    B[ 4]*C[12] - B[ 5]*C[13] +
    B[10]*C[14] - B[11]*C[15] +
    B[16]*C[16] - B[17]*C[17];
  A[17] =
    -B[ 4]*C[13] - B[ 5]*C[12] +
    -B[10]*C[15] - B[11]*C[14] +
    -B[16]*C[17] - B[17]*C[16];
}


// ********************



// A += B^\dagger * C^\dagger.

inline void cm_pl_eq_cm_dag_ti_cm_dag(double *A, const double *B, const double *C)
{
  // A00 += B00*C00 + B01*C10 + B02*C20
  A[ 0] +=
    B[ 0]*C[ 0] - B[ 1]*C[ 1] +
    B[ 6]*C[ 2] - B[ 7]*C[ 3] +
    B[12]*C[ 4] - B[13]*C[ 5];
  A[ 1] +=
    -B[ 0]*C[ 1] - B[ 1]*C[ 0] +
    -B[ 6]*C[ 3] - B[ 7]*C[ 2] +
    -B[12]*C[ 5] - B[13]*C[ 4];

  // A01 += B00*C01 + B01*C11 + B02*C21
  A[ 2] +=
    B[ 0]*C[ 6] - B[ 1]*C[ 7] +
    B[ 6]*C[ 8] - B[ 7]*C[ 9] +
    B[12]*C[10] - B[13]*C[11];
  A[ 3] +=
    -B[ 0]*C[ 7] - B[ 1]*C[ 6] +
    -B[ 6]*C[ 9] - B[ 7]*C[ 8] +
    -B[12]*C[11] - B[13]*C[10];

  // A02 += B00*C02 + B01*C12 + B02*C22
  A[ 4] +=
    B[ 0]*C[12] - B[ 1]*C[13] +
    B[ 6]*C[14] - B[ 7]*C[15] +
    B[12]*C[16] - B[13]*C[17];
  A[ 5] +=
    -B[ 0]*C[13] - B[ 1]*C[12] +
    -B[ 6]*C[15] - B[ 7]*C[14] +
    -B[12]*C[17] - B[13]*C[16];


  // A10 += B10*C00 + B11*C10 + B12*C20
  A[ 6] +=
    B[ 2]*C[ 0] - B[ 3]*C[ 1] +
    B[ 8]*C[ 2] - B[ 9]*C[ 3] +
    B[14]*C[ 4] - B[15]*C[ 5];
  A[ 7] +=
    -B[ 2]*C[ 1] - B[ 3]*C[ 0] +
    -B[ 8]*C[ 3] - B[ 9]*C[ 2] +
    -B[14]*C[ 5] - B[15]*C[ 4];

  // A11 += B10*C01 + B11*C11 + B12*C21
  A[ 8] +=
    B[ 2]*C[ 6] - B[ 3]*C[ 7] +
    B[ 8]*C[ 8] - B[ 9]*C[ 9] +
    B[14]*C[10] - B[15]*C[11];
  A[ 9] +=
    -B[ 2]*C[ 7] - B[ 3]*C[ 6] +
    -B[ 8]*C[ 9] - B[ 9]*C[ 8] +
    -B[14]*C[11] - B[15]*C[10];

  // A12 += B10*C02 + B11*C12 + B12*C22
  A[10] +=
    B[ 2]*C[12] - B[ 3]*C[13] +
    B[ 8]*C[14] - B[ 9]*C[15] +
    B[14]*C[16] - B[15]*C[17];
  A[11] +=
    -B[ 2]*C[13] - B[ 3]*C[12] +
    -B[ 8]*C[15] - B[ 9]*C[14] +
    -B[14]*C[17] - B[15]*C[16];


  // A20 += B20*C00 + B21*C10 + B22*C20
  A[12] +=
    B[ 4]*C[ 0] - B[ 5]*C[ 1] +
    B[10]*C[ 2] - B[11]*C[ 3] +
    B[16]*C[ 4] - B[17]*C[ 5];
  A[13] +=
    -B[ 4]*C[ 1] - B[ 5]*C[ 0] +
    -B[10]*C[ 3] - B[11]*C[ 2] +
    -B[16]*C[ 5] - B[17]*C[ 4];

  // A21 += B20*C01 + B21*C11 + B22*C21
  A[14] +=
    B[ 4]*C[ 6] - B[ 5]*C[ 7] +
    B[10]*C[ 8] - B[11]*C[ 9] +
    B[16]*C[10] - B[17]*C[11];
  A[15] +=
    -B[ 4]*C[ 7] - B[ 5]*C[ 6] +
    -B[10]*C[ 9] - B[11]*C[ 8] +
    -B[16]*C[11] - B[17]*C[10];

  // A22 += B20*C02 + B21*C12 + B22*C22
  A[16] +=
    B[ 4]*C[12] - B[ 5]*C[13] +
    B[10]*C[14] - B[11]*C[15] +
    B[16]*C[16] - B[17]*C[17];
  A[17] +=
    -B[ 4]*C[13] - B[ 5]*C[12] +
    -B[10]*C[15] - B[11]*C[14] +
    -B[16]*C[17] - B[17]*C[16];
}



// ********************



// A^{a b} = t_A^a s_A^b.

inline void cm_eq_fv_ti_fv(double *A, const double *s, const double *t)
{
  /*
  // a simple but slow code to cross check

  int i1, i2, i3;

  for(i1 = 0; i1 < 3; i1++)
    {
      for(i2 = 0; i2 < 3; i2++)
	{
	  A[(i1*3+i2)*2+0] = 0.0;
	  A[(i1*3+i2)*2+1] = 0.0;

	  for(i3 = 0; i3 < 4; i3++)
	    {
	      A[(i1*3+i2)*2+0] += s[(i3*3+i1)*2+0] * t[(i3*3+i2)*2+0] - s[(i3*3+i1)*2+1] * t[(i3*3+i2)*2+1];
	      A[(i1*3+i2)*2+1] += s[(i3*3+i1)*2+0] * t[(i3*3+i2)*2+1] + s[(i3*3+i1)*2+1] * t[(i3*3+i2)*2+0];
	    }
	}
    }

  */

  // A^{1 1}

  A[ 0] =
    s[ 0]*t[ 0] - s[ 1]*t[ 1] +
    s[ 6]*t[ 6] - s[ 7]*t[ 7] +
    s[12]*t[12] - s[13]*t[13] +
    s[18]*t[18] - s[19]*t[19];

  A[ 1] =
    s[ 0]*t[ 1] + s[ 1]*t[ 0] +
    s[ 6]*t[ 7] + s[ 7]*t[ 6] +
    s[12]*t[13] + s[13]*t[12] +
    s[18]*t[19] + s[19]*t[18];

  // A^{1 2}

  A[ 2] =
    s[ 0]*t[ 2] - s[ 1]*t[ 3] +
    s[ 6]*t[ 8] - s[ 7]*t[ 9] +
    s[12]*t[14] - s[13]*t[15] +
    s[18]*t[20] - s[19]*t[21];

  A[ 3] =
    s[ 0]*t[ 3] + s[ 1]*t[ 2] +
    s[ 6]*t[ 9] + s[ 7]*t[ 8] +
    s[12]*t[15] + s[13]*t[14] +
    s[18]*t[21] + s[19]*t[20];

  // A^{1 3}

  A[ 4] =
    s[ 0]*t[ 4] - s[ 1]*t[ 5] +
    s[ 6]*t[10] - s[ 7]*t[11] +
    s[12]*t[16] - s[13]*t[17] +
    s[18]*t[22] - s[19]*t[23];

  A[ 5] =
    s[ 0]*t[ 5] + s[ 1]*t[ 4] +
    s[ 6]*t[11] + s[ 7]*t[10] +
    s[12]*t[17] + s[13]*t[16] +
    s[18]*t[23] + s[19]*t[22];

  // ***

  // A^{2 1}

  A[ 6] =
    s[ 2]*t[ 0] - s[ 3]*t[ 1] +
    s[ 8]*t[ 6] - s[ 9]*t[ 7] +
    s[14]*t[12] - s[15]*t[13] +
    s[20]*t[18] - s[21]*t[19];

  A[ 7] =
    s[ 2]*t[ 1] + s[ 3]*t[ 0] +
    s[ 8]*t[ 7] + s[ 9]*t[ 6] +
    s[14]*t[13] + s[15]*t[12] +
    s[20]*t[19] + s[21]*t[18];

  // A^{2 2}

  A[ 8] =
    s[ 2]*t[ 2] - s[ 3]*t[ 3] +
    s[ 8]*t[ 8] - s[ 9]*t[ 9] +
    s[14]*t[14] - s[15]*t[15] +
    s[20]*t[20] - s[21]*t[21];

  A[ 9] =
    s[ 2]*t[ 3] + s[ 3]*t[ 2] +
    s[ 8]*t[ 9] + s[ 9]*t[ 8] +
    s[14]*t[15] + s[15]*t[14] +
    s[20]*t[21] + s[21]*t[20];

  // A^{2 3}

  A[10] =
    s[ 2]*t[ 4] - s[ 3]*t[ 5] +
    s[ 8]*t[10] - s[ 9]*t[11] +
    s[14]*t[16] - s[15]*t[17] +
    s[20]*t[22] - s[21]*t[23];

  A[11] =
    s[ 2]*t[ 5] + s[ 3]*t[ 4] +
    s[ 8]*t[11] + s[ 9]*t[10] +
    s[14]*t[17] + s[15]*t[16] +
    s[20]*t[23] + s[21]*t[22];

  // ***

  // A^{3 1}

  A[12] =
    s[ 4]*t[ 0] - s[ 5]*t[ 1] +
    s[10]*t[ 6] - s[11]*t[ 7] +
    s[16]*t[12] - s[17]*t[13] +
    s[22]*t[18] - s[23]*t[19];

  A[13] =
    s[ 4]*t[ 1] + s[ 5]*t[ 0] +
    s[10]*t[ 7] + s[11]*t[ 6] +
    s[16]*t[13] + s[17]*t[12] +
    s[22]*t[19] + s[23]*t[18];

  // A^{3 2}

  A[14] =
    s[ 4]*t[ 2] - s[ 5]*t[ 3] +
    s[10]*t[ 8] - s[11]*t[ 9] +
    s[16]*t[14] - s[17]*t[15] +
    s[22]*t[20] - s[23]*t[21];

  A[15] =
    s[ 4]*t[ 3] + s[ 5]*t[ 2] +
    s[10]*t[ 9] + s[11]*t[ 8] +
    s[16]*t[15] + s[17]*t[14] +
    s[22]*t[21] + s[23]*t[20];

  // A^{3 3}

  A[16] =
    s[ 4]*t[ 4] - s[ 5]*t[ 5] +
    s[10]*t[10] - s[11]*t[11] +
    s[16]*t[16] - s[17]*t[17] +
    s[22]*t[22] - s[23]*t[23];

  A[17] =
    s[ 4]*t[ 5] + s[ 5]*t[ 4] +
    s[10]*t[11] + s[11]*t[10] +
    s[16]*t[17] + s[17]*t[16] +
    s[22]*t[23] + s[23]*t[22];
}



// ********************



// c = \epsilon^{a b c} \epsilon^{d e f} A^{a d} B^{b c} C^{e f}.

inline void co_eq_eps_eps_cm_cm_cm(complex *c, const double *A, const double *B, const double *C)
{
  /*
  // a simple but slow code to cross check

  c->re = 0.0;
  c->im = 0.0;

  int i1, i2, i3, i4, i5, i6;

  for(i1 = 0; i1 < 3; i1++)
    {
      for(i2 = 0; i2 < 3; i2++)
	{
	  for(i3 = 0; i3 < 3; i3++)
	    {
	      for(i4 = 0; i4 < 3; i4++)
		{
		  for(i5 = 0; i5 < 3; i5++)
		    {
		      for(i6 = 0; i6 < 3; i6++)
			{
			  double sign1 = 0.0;

			  if((i1 == 0 && i2 == 1 && i3 == 2) ||
			     (i1 == 1 && i2 == 2 && i3 == 0) ||
			     (i1 == 2 && i2 == 0 && i3 == 1))
			    sign1 = +1.0;

			  if((i1 == 2 && i2 == 1 && i3 == 0) ||
			     (i1 == 1 && i2 == 0 && i3 == 2) ||
			     (i1 == 0 && i2 == 2 && i3 == 1))
			    sign1 = -1.0;

			  double sign2 = 0.0;

			  if((i4 == 0 && i5 == 1 && i6 == 2) ||
			     (i4 == 1 && i5 == 2 && i6 == 0) ||
			     (i4 == 2 && i5 == 0 && i6 == 1))
			    sign2 = +1.0;

			  if((i4 == 2 && i5 == 1 && i6 == 0) ||
			     (i4 == 1 && i5 == 0 && i6 == 2) ||
			     (i4 == 0 && i5 == 2 && i6 == 1))
			    sign2 = -1.0;

			  c->re += sign1 * sign2 * A[(i1*3+i4)*2+0] * B[(i2*3+i3)*2+0] * C[(i5*3+i6)*2+0];
			  c->re -= sign1 * sign2 * A[(i1*3+i4)*2+0] * B[(i2*3+i3)*2+1] * C[(i5*3+i6)*2+1];
			  c->re -= sign1 * sign2 * A[(i1*3+i4)*2+1] * B[(i2*3+i3)*2+0] * C[(i5*3+i6)*2+1];
			  c->re -= sign1 * sign2 * A[(i1*3+i4)*2+1] * B[(i2*3+i3)*2+1] * C[(i5*3+i6)*2+0];

			  c->im -= sign1 * sign2 * A[(i1*3+i4)*2+1] * B[(i2*3+i3)*2+1] * C[(i5*3+i6)*2+1];
			  c->im += sign1 * sign2 * A[(i1*3+i4)*2+1] * B[(i2*3+i3)*2+0] * C[(i5*3+i6)*2+0];
			  c->im += sign1 * sign2 * A[(i1*3+i4)*2+0] * B[(i2*3+i3)*2+1] * C[(i5*3+i6)*2+0];
			  c->im += sign1 * sign2 * A[(i1*3+i4)*2+0] * B[(i2*3+i3)*2+0] * C[(i5*3+i6)*2+1];
			}
		    }
		}
	    }
	}
    }
  */

  c->re =


    // 1{23} 1{23}
    +A[ 0] * (B[10]-B[14]) * (C[10]-C[14])
    -A[ 0] * (B[11]-B[15]) * (C[11]-C[15])
    -A[ 1] * (B[10]-B[14]) * (C[11]-C[15])
    -A[ 1] * (B[11]-B[15]) * (C[10]-C[14])

    // 1{23} 2{31}
    +A[ 2] * (B[10]-B[14]) * (C[12]-C[ 4])
    -A[ 2] * (B[11]-B[15]) * (C[13]-C[ 5])
    -A[ 3] * (B[10]-B[14]) * (C[13]-C[ 5])
    -A[ 3] * (B[11]-B[15]) * (C[12]-C[ 4])

    // 1{23} 3{12}
    +A[ 4] * (B[10]-B[14]) * (C[ 2]-C[ 6])
    -A[ 4] * (B[11]-B[15]) * (C[ 3]-C[ 7])
    -A[ 5] * (B[10]-B[14]) * (C[ 3]-C[ 7])
    -A[ 5] * (B[11]-B[15]) * (C[ 2]-C[ 6])


    // 2{31} 1{23}
    +A[ 6] * (B[12]-B[ 4]) * (C[10]-C[14])
    -A[ 6] * (B[13]-B[ 5]) * (C[11]-C[15])
    -A[ 7] * (B[12]-B[ 4]) * (C[11]-C[15])
    -A[ 7] * (B[13]-B[ 5]) * (C[10]-C[14])

    // 2{31} 2{31}
    +A[ 8] * (B[12]-B[ 4]) * (C[12]-C[ 4])
    -A[ 8] * (B[13]-B[ 5]) * (C[13]-C[ 5])
    -A[ 9] * (B[12]-B[ 4]) * (C[13]-C[ 5])
    -A[ 9] * (B[13]-B[ 5]) * (C[12]-C[ 4])

    // 2{31} 3{12}
    +A[10] * (B[12]-B[ 4]) * (C[ 2]-C[ 6])
    -A[10] * (B[13]-B[ 5]) * (C[ 3]-C[ 7])
    -A[11] * (B[12]-B[ 4]) * (C[ 3]-C[ 7])
    -A[11] * (B[13]-B[ 5]) * (C[ 2]-C[ 6])


    // 3{12} 1{23}
    +A[12] * (B[ 2]-B[ 6]) * (C[10]-C[14])
    -A[12] * (B[ 3]-B[ 7]) * (C[11]-C[15])
    -A[13] * (B[ 2]-B[ 6]) * (C[11]-C[15])
    -A[13] * (B[ 3]-B[ 7]) * (C[10]-C[14])

    // 3{12} 2{31}
    +A[14] * (B[ 2]-B[ 6]) * (C[12]-C[ 4])
    -A[14] * (B[ 3]-B[ 7]) * (C[13]-C[ 5])
    -A[15] * (B[ 2]-B[ 6]) * (C[13]-C[ 5])
    -A[15] * (B[ 3]-B[ 7]) * (C[12]-C[ 4])

    // 3{12} 3{12}
    +A[16] * (B[ 2]-B[ 6]) * (C[ 2]-C[ 6])
    -A[16] * (B[ 3]-B[ 7]) * (C[ 3]-C[ 7])
    -A[17] * (B[ 2]-B[ 6]) * (C[ 3]-C[ 7])
    -A[17] * (B[ 3]-B[ 7]) * (C[ 2]-C[ 6]);

  // *****

  c->im =


    // 1{23} 1{23}
    -A[ 1] * (B[11]-B[15]) * (C[11]-C[15])
    +A[ 1] * (B[10]-B[14]) * (C[10]-C[14])
    +A[ 0] * (B[11]-B[15]) * (C[10]-C[14])
    +A[ 0] * (B[10]-B[14]) * (C[11]-C[15])

    // 1{23} 2{31}
    -A[ 3] * (B[11]-B[15]) * (C[13]-C[ 5])
    +A[ 3] * (B[10]-B[14]) * (C[12]-C[ 4])
    +A[ 2] * (B[11]-B[15]) * (C[12]-C[ 4])
    +A[ 2] * (B[10]-B[14]) * (C[13]-C[ 5])

    // 1{23} 3{12}
    -A[ 5] * (B[11]-B[15]) * (C[ 3]-C[ 7])
    +A[ 5] * (B[10]-B[14]) * (C[ 2]-C[ 6])
    +A[ 4] * (B[11]-B[15]) * (C[ 2]-C[ 6])
    +A[ 4] * (B[10]-B[14]) * (C[ 3]-C[ 7])


    // 2{31} 1{23}
    -A[ 7] * (B[13]-B[ 5]) * (C[11]-C[15])
    +A[ 7] * (B[12]-B[ 4]) * (C[10]-C[14])
    +A[ 6] * (B[13]-B[ 5]) * (C[10]-C[14])
    +A[ 6] * (B[12]-B[ 4]) * (C[11]-C[15])

    // 2{31} 2{31}
    -A[ 9] * (B[13]-B[ 5]) * (C[13]-C[ 5])
    +A[ 9] * (B[12]-B[ 4]) * (C[12]-C[ 4])
    +A[ 8] * (B[13]-B[ 5]) * (C[12]-C[ 4])
    +A[ 8] * (B[12]-B[ 4]) * (C[13]-C[ 5])

    // 2{31} 3{12}
    -A[11] * (B[13]-B[ 5]) * (C[ 3]-C[ 7])
    +A[11] * (B[12]-B[ 4]) * (C[ 2]-C[ 6])
    +A[10] * (B[13]-B[ 5]) * (C[ 2]-C[ 6])
    +A[10] * (B[12]-B[ 4]) * (C[ 3]-C[ 7])


    // 3{12} 1{23}
    -A[13] * (B[ 3]-B[ 7]) * (C[11]-C[15])
    +A[13] * (B[ 2]-B[ 6]) * (C[10]-C[14])
    +A[12] * (B[ 3]-B[ 7]) * (C[10]-C[14])
    +A[12] * (B[ 2]-B[ 6]) * (C[11]-C[15])

    // 3{12} 2{31}
    -A[15] * (B[ 3]-B[ 7]) * (C[13]-C[ 5])
    +A[15] * (B[ 2]-B[ 6]) * (C[12]-C[ 4])
    +A[14] * (B[ 3]-B[ 7]) * (C[12]-C[ 4])
    +A[14] * (B[ 2]-B[ 6]) * (C[13]-C[ 5])

    // 3{12} 3{12}
    -A[17] * (B[ 3]-B[ 7]) * (C[ 3]-C[ 7])
    +A[17] * (B[ 2]-B[ 6]) * (C[ 2]-C[ 6])
    +A[16] * (B[ 3]-B[ 7]) * (C[ 2]-C[ 6])
    +A[16] * (B[ 2]-B[ 6]) * (C[ 3]-C[ 7]);
}



// ********************
// ********************
// ********************
// ********************
// ********************



// Prints a gamma matrix.

inline void gamma_fprintf(int gamma_index, FILE *file)
{
  complex c1;
  int i1, i2;


  double s[24], t[24];

  for(i1 = 0; i1 < 4; i1++)
    {
      for(i2 = 0; i2 < 4; i2++)
	{
	  fv_eq_zero(t);
	  t[i2*6] = 1.0;

	  fv_eq_gamma_ti_fv(s, gamma_index, t);

	  fv_eq_zero(t);
	  t[i1*6] = 1.0;

	  co_eq_fv_dag_ti_fv(&c1, t, s);

	  fprintf(file, "%+.1f %+.1f   ", c1.re, c1.im);
	}

      fprintf(stderr, "\n");
    }
}

inline void gamma_fprintf(FILE *file, int const per[], int const sig[])
{
  complex c1;
  int i1, i2;


  double s[24], t[24];

  for(i1 = 0; i1 < 4; i1++)
    {
      for(i2 = 0; i2 < 4; i2++)
	{
	  fv_eq_zero(t);
	  t[i2*6] = 1.0;

	  fv_eq_gamma_ti_fv(s, t, per, sig);

	  fv_eq_zero(t);
	  t[i1*6] = 1.0;

	  co_eq_fv_dag_ti_fv(&c1, t, s);

	  fprintf(file, "%+.1f %+.1f   ", c1.re, c1.im);
	}

      fprintf(stderr, "\n");
    }
}



// ********************
// ********************
// ********************
// ********************
// ********************



// s = 0.

inline void fv_eq_zero(double *s)
{

  s[ 0] = 0.0;
  s[ 1] = 0.0;
  s[ 2] = 0.0;
  s[ 3] = 0.0;
  s[ 4] = 0.0;
  s[ 5] = 0.0;

  s[ 6] = 0.0;
  s[ 7] = 0.0;
  s[ 8] = 0.0;
  s[ 9] = 0.0;
  s[10] = 0.0;
  s[11] = 0.0;

  s[12] = 0.0;
  s[13] = 0.0;
  s[14] = 0.0;
  s[15] = 0.0;
  s[16] = 0.0;
  s[17] = 0.0;

  s[18] = 0.0;
  s[19] = 0.0;
  s[20] = 0.0;
  s[21] = 0.0;
  s[22] = 0.0;
  s[23] = 0.0;
}



// ********************



// s = -s.

inline void fv_mi(double *s)
{
  s[ 0] = -s[ 0];
  s[ 1] = -s[ 1];
  s[ 2] = -s[ 2];
  s[ 3] = -s[ 3];
  s[ 4] = -s[ 4];
  s[ 5] = -s[ 5];

  s[ 6] = -s[ 6];
  s[ 7] = -s[ 7];
  s[ 8] = -s[ 8];
  s[ 9] = -s[ 9];
  s[10] = -s[10];
  s[11] = -s[11];

  s[12] = -s[12];
  s[13] = -s[13];
  s[14] = -s[14];
  s[15] = -s[15];
  s[16] = -s[16];
  s[17] = -s[17];

  s[18] = -s[18];
  s[19] = -s[19];
  s[20] = -s[20];
  s[21] = -s[21];
  s[22] = -s[22];
  s[23] = -s[23];
}



// ********************



// s = t.

inline void fv_eq_fv(double *s, const double *t)
{
  s[ 0] = t[ 0];
  s[ 1] = t[ 1];
  s[ 2] = t[ 2];
  s[ 3] = t[ 3];
  s[ 4] = t[ 4];
  s[ 5] = t[ 5];

  s[ 6] = t[ 6];
  s[ 7] = t[ 7];
  s[ 8] = t[ 8];
  s[ 9] = t[ 9];
  s[10] = t[10];
  s[11] = t[11];

  s[12] = t[12];
  s[13] = t[13];
  s[14] = t[14];
  s[15] = t[15];
  s[16] = t[16];
  s[17] = t[17];

  s[18] = t[18];
  s[19] = t[19];
  s[20] = t[20];
  s[21] = t[21];
  s[22] = t[22];
  s[23] = t[23];
}




// ********************



// s = t^\dagger.

inline void fv_eq_fv_dag(double *s, const double *t)
{
  s[ 0] = t[ 0];
  s[ 1] = -t[ 1];
  s[ 2] = t[ 2];
  s[ 3] = -t[ 3];
  s[ 4] = t[ 4];
  s[ 5] = -t[ 5];

  s[ 6] = t[ 6];
  s[ 7] = -t[ 7];
  s[ 8] = t[ 8];
  s[ 9] = -t[ 9];
  s[10] = t[10];
  s[11] = -t[11];

  s[12] = t[12];
  s[13] = -t[13];
  s[14] = t[14];
  s[15] = -t[15];
  s[16] = t[16];
  s[17] = -t[17];

  s[18] = t[18];
  s[19] = -t[19];
  s[20] = t[20];
  s[21] = -t[21];
  s[22] = t[22];
  s[23] = -t[23];
}



// ********************



// Prints s.

inline void fv_fprintf(const double *s, FILE *file)
{
  fprintf(file, "%+9.6f %+9.6f I   %+9.6f %+9.6f I   %+9.6f %+9.6f I\n",
	  s[ 0], s[ 1], s[ 2], s[ 3], s[ 4], s[ 5]);

  fprintf(file, "%+9.6f %+9.6f I   %+9.6f %+9.6f I   %+9.6f %+9.6f I\n",
	  s[ 6], s[ 7], s[ 8], s[ 9], s[10], s[11]);

  fprintf(file, "%+9.6f %+9.6f I   %+9.6f %+9.6f I   %+9.6f %+9.6f I\n",
	  s[12], s[13], s[14], s[15], s[16], s[17]);

  fprintf(file, "%+9.6f %+9.6f I   %+9.6f %+9.6f I   %+9.6f %+9.6f I\n",
	  s[18], s[19], s[20], s[21], s[22], s[23]);
}



// ********************



// c = s^dagger * t.

inline void co_eq_fv_dag_ti_fv(complex *c, const double *s, const double *t)
{
  c->re = 

    s[ 0]*t[ 0] + s[ 1]*t[ 1] +
    s[ 2]*t[ 2] + s[ 3]*t[ 3] + 
    s[ 4]*t[ 4] + s[ 5]*t[ 5] + 

    s[ 6]*t[ 6] + s[ 7]*t[ 7] + 
    s[ 8]*t[ 8] + s[ 9]*t[ 9] + 
    s[10]*t[10] + s[11]*t[11] + 

    s[12]*t[12] + s[13]*t[13] + 
    s[14]*t[14] + s[15]*t[15] + 
    s[16]*t[16] + s[17]*t[17] + 

    s[18]*t[18] + s[19]*t[19] + 
    s[20]*t[20] + s[21]*t[21] + 
    s[22]*t[22] + s[23]*t[23];

  c->im =

    s[ 0]*t[ 1] - s[ 1]*t[ 0] +
    s[ 2]*t[ 3] - s[ 3]*t[ 2] + 
    s[ 4]*t[ 5] - s[ 5]*t[ 4] + 

    s[ 6]*t[ 7] - s[ 7]*t[ 6] + 
    s[ 8]*t[ 9] - s[ 9]*t[ 8] + 
    s[10]*t[11] - s[11]*t[10] + 

    s[12]*t[13] - s[13]*t[12] + 
    s[14]*t[15] - s[15]*t[14] + 
    s[16]*t[17] - s[17]*t[16] + 

    s[18]*t[19] - s[19]*t[18] + 
    s[20]*t[21] - s[21]*t[20] + 
    s[22]*t[23] - s[23]*t[22];
}

// c += s^dagger * t.

inline void co_pl_eq_fv_dag_ti_fv(complex *c, const double *s, const double *t)
{
  c->re +=

    s[ 0]*t[ 0] + s[ 1]*t[ 1] +
    s[ 2]*t[ 2] + s[ 3]*t[ 3] + 
    s[ 4]*t[ 4] + s[ 5]*t[ 5] + 

    s[ 6]*t[ 6] + s[ 7]*t[ 7] + 
    s[ 8]*t[ 8] + s[ 9]*t[ 9] + 
    s[10]*t[10] + s[11]*t[11] + 

    s[12]*t[12] + s[13]*t[13] + 
    s[14]*t[14] + s[15]*t[15] + 
    s[16]*t[16] + s[17]*t[17] + 

    s[18]*t[18] + s[19]*t[19] + 
    s[20]*t[20] + s[21]*t[21] + 
    s[22]*t[22] + s[23]*t[23];

  c->im +=

    s[ 0]*t[ 1] - s[ 1]*t[ 0] +
    s[ 2]*t[ 3] - s[ 3]*t[ 2] + 
    s[ 4]*t[ 5] - s[ 5]*t[ 4] + 

    s[ 6]*t[ 7] - s[ 7]*t[ 6] + 
    s[ 8]*t[ 9] - s[ 9]*t[ 8] + 
    s[10]*t[11] - s[11]*t[10] + 

    s[12]*t[13] - s[13]*t[12] + 
    s[14]*t[15] - s[15]*t[14] + 
    s[16]*t[17] - s[17]*t[16] + 

    s[18]*t[19] - s[19]*t[18] + 
    s[20]*t[21] - s[21]*t[20] + 
    s[22]*t[23] - s[23]*t[22];
}



// ********************



// s = t + u.

inline void fv_eq_fv_pl_fv(double *s, const double *t, const double *u)
{
  s[ 0] = t[ 0] + u[ 0];
  s[ 1] = t[ 1] + u[ 1];
  s[ 2] = t[ 2] + u[ 2];
  s[ 3] = t[ 3] + u[ 3];
  s[ 4] = t[ 4] + u[ 4];
  s[ 5] = t[ 5] + u[ 5];

  s[ 6] = t[ 6] + u[ 6];
  s[ 7] = t[ 7] + u[ 7];
  s[ 8] = t[ 8] + u[ 8];
  s[ 9] = t[ 9] + u[ 9];
  s[10] = t[10] + u[10];
  s[11] = t[11] + u[11];

  s[12] = t[12] + u[12];
  s[13] = t[13] + u[13];
  s[14] = t[14] + u[14];
  s[15] = t[15] + u[15];
  s[16] = t[16] + u[16];
  s[17] = t[17] + u[17];

  s[18] = t[18] + u[18];
  s[19] = t[19] + u[19];
  s[20] = t[20] + u[20];
  s[21] = t[21] + u[21];
  s[22] = t[22] + u[22];
  s[23] = t[23] + u[23];
}

inline void fv_eq_fv_pl_fv(float *s, const float *t, const float *u)
{
  s[ 0] = t[ 0] + u[ 0];
  s[ 1] = t[ 1] + u[ 1];
  s[ 2] = t[ 2] + u[ 2];
  s[ 3] = t[ 3] + u[ 3];
  s[ 4] = t[ 4] + u[ 4];
  s[ 5] = t[ 5] + u[ 5];

  s[ 6] = t[ 6] + u[ 6];
  s[ 7] = t[ 7] + u[ 7];
  s[ 8] = t[ 8] + u[ 8];
  s[ 9] = t[ 9] + u[ 9];
  s[10] = t[10] + u[10];
  s[11] = t[11] + u[11];

  s[12] = t[12] + u[12];
  s[13] = t[13] + u[13];
  s[14] = t[14] + u[14];
  s[15] = t[15] + u[15];
  s[16] = t[16] + u[16];
  s[17] = t[17] + u[17];

  s[18] = t[18] + u[18];
  s[19] = t[19] + u[19];
  s[20] = t[20] + u[20];
  s[21] = t[21] + u[21];
  s[22] = t[22] + u[22];
  s[23] = t[23] + u[23];
}



// ********************



// s = t - u.

inline void fv_eq_fv_mi_fv(double *s, const double *t, const double *u)
{
  s[ 0] = t[ 0] - u[ 0];
  s[ 1] = t[ 1] - u[ 1];
  s[ 2] = t[ 2] - u[ 2];
  s[ 3] = t[ 3] - u[ 3];
  s[ 4] = t[ 4] - u[ 4];
  s[ 5] = t[ 5] - u[ 5];

  s[ 6] = t[ 6] - u[ 6];
  s[ 7] = t[ 7] - u[ 7];
  s[ 8] = t[ 8] - u[ 8];
  s[ 9] = t[ 9] - u[ 9];
  s[10] = t[10] - u[10];
  s[11] = t[11] - u[11];

  s[12] = t[12] - u[12];
  s[13] = t[13] - u[13];
  s[14] = t[14] - u[14];
  s[15] = t[15] - u[15];
  s[16] = t[16] - u[16];
  s[17] = t[17] - u[17];

  s[18] = t[18] - u[18];
  s[19] = t[19] - u[19];
  s[20] = t[20] - u[20];
  s[21] = t[21] - u[21];
  s[22] = t[22] - u[22];
  s[23] = t[23] - u[23];
}



// ********************



// s = s + t.

inline void fv_pl_eq_fv(double *s, const double *t)
{
  s[ 0] += t[ 0];
  s[ 1] += t[ 1];
  s[ 2] += t[ 2];
  s[ 3] += t[ 3];
  s[ 4] += t[ 4];
  s[ 5] += t[ 5];

  s[ 6] += t[ 6];
  s[ 7] += t[ 7];
  s[ 8] += t[ 8];
  s[ 9] += t[ 9];
  s[10] += t[10];
  s[11] += t[11];

  s[12] += t[12];
  s[13] += t[13];
  s[14] += t[14];
  s[15] += t[15];
  s[16] += t[16];
  s[17] += t[17];

  s[18] += t[18];
  s[19] += t[19];
  s[20] += t[20];
  s[21] += t[21];
  s[22] += t[22];
  s[23] += t[23];
}



// ********************



// s = s - t.

inline void fv_mi_eq_fv(double *s, const double *t)
{
  s[ 0] -= t[ 0];
  s[ 1] -= t[ 1];
  s[ 2] -= t[ 2];
  s[ 3] -= t[ 3];
  s[ 4] -= t[ 4];
  s[ 5] -= t[ 5];

  s[ 6] -= t[ 6];
  s[ 7] -= t[ 7];
  s[ 8] -= t[ 8];
  s[ 9] -= t[ 9];
  s[10] -= t[10];
  s[11] -= t[11];

  s[12] -= t[12];
  s[13] -= t[13];
  s[14] -= t[14];
  s[15] -= t[15];
  s[16] -= t[16];
  s[17] -= t[17];

  s[18] -= t[18];
  s[19] -= t[19];
  s[20] -= t[20];
  s[21] -= t[21];
  s[22] -= t[22];
  s[23] -= t[23];
}



// ********************



// s = t * d.

inline void fv_eq_fv_ti_re(double *s, const double *t, double d)
{
  s[ 0] = t[ 0] * d;
  s[ 1] = t[ 1] * d;
  s[ 2] = t[ 2] * d;
  s[ 3] = t[ 3] * d;
  s[ 4] = t[ 4] * d;
  s[ 5] = t[ 5] * d;

  s[ 6] = t[ 6] * d;
  s[ 7] = t[ 7] * d;
  s[ 8] = t[ 8] * d;
  s[ 9] = t[ 9] * d;
  s[10] = t[10] * d;
  s[11] = t[11] * d;

  s[12] = t[12] * d;
  s[13] = t[13] * d;
  s[14] = t[14] * d;
  s[15] = t[15] * d;
  s[16] = t[16] * d;
  s[17] = t[17] * d;

  s[18] = t[18] * d;
  s[19] = t[19] * d;
  s[20] = t[20] * d;
  s[21] = t[21] * d;
  s[22] = t[22] * d;
  s[23] = t[23] * d;
}



// ********************



// s = s * d.

inline void fv_ti_eq_re(double *s, double d)
{
  s[ 0] *= d;
  s[ 1] *= d;
  s[ 2] *= d;
  s[ 3] *= d;
  s[ 4] *= d;
  s[ 5] *= d;

  s[ 6] *= d;
  s[ 7] *= d;
  s[ 8] *= d;
  s[ 9] *= d;
  s[10] *= d;
  s[11] *= d;

  s[12] *= d;
  s[13] *= d;
  s[14] *= d;
  s[15] *= d;
  s[16] *= d;
  s[17] *= d;

  s[18] *= d;
  s[19] *= d;
  s[20] *= d;
  s[21] *= d;
  s[22] *= d;
  s[23] *= d;
}

inline void fv_ti_eq_re(float *s, double d)
{
  s[ 0] *= d;
  s[ 1] *= d;
  s[ 2] *= d;
  s[ 3] *= d;
  s[ 4] *= d;
  s[ 5] *= d;

  s[ 6] *= d;
  s[ 7] *= d;
  s[ 8] *= d;
  s[ 9] *= d;
  s[10] *= d;
  s[11] *= d;

  s[12] *= d;
  s[13] *= d;
  s[14] *= d;
  s[15] *= d;
  s[16] *= d;
  s[17] *= d;

  s[18] *= d;
  s[19] *= d;
  s[20] *= d;
  s[21] *= d;
  s[22] *= d;
  s[23] *= d;
}



// ********************



// s = t * i * d.

inline void fv_eq_fv_ti_im(double *s, const double *t, double d)
{
  s[ 0] = -t[ 1] * d;
  s[ 1] =  t[ 0] * d;
  s[ 2] = -t[ 3] * d;
  s[ 3] =  t[ 2] * d;
  s[ 4] = -t[ 5] * d;
  s[ 5] =  t[ 4] * d;

  s[ 6] = -t[ 7] * d;
  s[ 7] =  t[ 6] * d;
  s[ 8] = -t[ 9] * d;
  s[ 9] =  t[ 8] * d;
  s[10] = -t[11] * d;
  s[11] =  t[10] * d;

  s[12] = -t[13] * d;
  s[13] =  t[12] * d;
  s[14] = -t[15] * d;
  s[15] =  t[14] * d;
  s[16] = -t[17] * d;
  s[17] =  t[16] * d;

  s[18] = -t[19] * d;
  s[19] =  t[18] * d;
  s[20] = -t[21] * d;
  s[21] =  t[20] * d;
  s[22] = -t[23] * d;
  s[23] =  t[22] * d;
}



// ********************



// s = s * i * d.

inline void fv_ti_eq_im(double *s, double d)
{
  double s_tmp;


  s_tmp = s[ 0];
  s[ 0] = -s[ 1]*d;
  s[ 1] = s_tmp*d;

  s_tmp = s[ 2];
  s[ 2] = -s[ 3]*d;
  s[ 3] = s_tmp*d;

  s_tmp = s[ 4];
  s[ 4] = -s[ 5]*d;
  s[ 5] = s_tmp*d;


  s_tmp = s[ 6];
  s[ 6] = -s[ 7]*d;
  s[ 7] = s_tmp*d;

  s_tmp = s[ 8];
  s[ 8] = -s[ 9]*d;
  s[ 9] = s_tmp*d;

  s_tmp = s[10];
  s[10] = -s[11]*d;
  s[11] = s_tmp*d;


  s_tmp = s[12];
  s[12] = -s[13]*d;
  s[13] = s_tmp*d;

  s_tmp = s[14];
  s[14] = -s[15]*d;
  s[15] = s_tmp*d;

  s_tmp = s[16];
  s[16] = -s[17]*d;
  s[17] = s_tmp*d;


  s_tmp = s[18];
  s[18] = -s[19]*d;
  s[19] = s_tmp*d;

  s_tmp = s[20];
  s[20] = -s[21]*d;
  s[21] = s_tmp*d;

  s_tmp = s[22];
  s[22] = -s[23]*d;
  s[23] = s_tmp*d;
}



// ********************



// s = t * c.

inline void fv_eq_fv_ti_co(double *s, const double *t, const complex *c)
{
  double re = c->re;
  double im = c->im;

  s[ 0] = t[ 0] * re - t[ 1] * im;
  s[ 1] = t[ 1] * re + t[ 0] * im;
  s[ 2] = t[ 2] * re - t[ 3] * im;
  s[ 3] = t[ 3] * re + t[ 2] * im;
  s[ 4] = t[ 4] * re - t[ 5] * im;
  s[ 5] = t[ 5] * re + t[ 4] * im;

  s[ 6] = t[ 6] * re - t[ 7] * im;
  s[ 7] = t[ 7] * re + t[ 6] * im;
  s[ 8] = t[ 8] * re - t[ 9] * im;
  s[ 9] = t[ 9] * re + t[ 8] * im;
  s[10] = t[10] * re - t[11] * im;
  s[11] = t[11] * re + t[10] * im;

  s[12] = t[12] * re - t[13] * im;
  s[13] = t[13] * re + t[12] * im;
  s[14] = t[14] * re - t[15] * im;
  s[15] = t[15] * re + t[14] * im;
  s[16] = t[16] * re - t[17] * im;
  s[17] = t[17] * re + t[16] * im;

  s[18] = t[18] * re - t[19] * im;
  s[19] = t[19] * re + t[18] * im;
  s[20] = t[20] * re - t[21] * im;
  s[21] = t[21] * re + t[20] * im;
  s[22] = t[22] * re - t[23] * im;
  s[23] = t[23] * re + t[22] * im;
}



// ********************



// s = A * t.

inline void fv_eq_cm_ti_fv(double *s, const double *A, const double *t)
{
  // s0 = U00*t0 + U01*t1 + U02*t2
  s[ 0] =
    A[ 0]*t[ 0] - A[ 1]*t[ 1] +
    A[ 2]*t[ 2] - A[ 3]*t[ 3] +
    A[ 4]*t[ 4] - A[ 5]*t[ 5];
  s[ 1] =
    A[ 0]*t[ 1] + A[ 1]*t[ 0] +
    A[ 2]*t[ 3] + A[ 3]*t[ 2] +
    A[ 4]*t[ 5] + A[ 5]*t[ 4];

  // s1 = U10*t0 + U11*t1 + U12*t2
  s[ 2] =
    A[ 6]*t[ 0] - A[ 7]*t[ 1] +
    A[ 8]*t[ 2] - A[ 9]*t[ 3] +
    A[10]*t[ 4] - A[11]*t[ 5];
  s[ 3] =
    A[ 6]*t[ 1] + A[ 7]*t[ 0] +
    A[ 8]*t[ 3] + A[ 9]*t[ 2] +
    A[10]*t[ 5] + A[11]*t[ 4];

  // s2 = U20*t0 + U21*t1 + U22*t2
  s[ 4] =
    A[12]*t[ 0] - A[13]*t[ 1] +
    A[14]*t[ 2] - A[15]*t[ 3] +
    A[16]*t[ 4] - A[17]*t[ 5];
  s[ 5] =
    A[12]*t[ 1] + A[13]*t[ 0] +
    A[14]*t[ 3] + A[15]*t[ 2] +
    A[16]*t[ 5] + A[17]*t[ 4];


  // s3 = U00*t3 + U01*t4 + U02*t5
  s[ 6] =
    A[ 0]*t[ 6] - A[ 1]*t[ 7] +
    A[ 2]*t[ 8] - A[ 3]*t[ 9] +
    A[ 4]*t[10] - A[ 5]*t[11];
  s[ 7] =
    A[ 0]*t[ 7] + A[ 1]*t[ 6] +
    A[ 2]*t[ 9] + A[ 3]*t[ 8] +
    A[ 4]*t[11] + A[ 5]*t[10];

  // s4 = U10*t3 + U11*t4 + U12*t5
  s[ 8] =
    A[ 6]*t[ 6] - A[ 7]*t[ 7] +
    A[ 8]*t[ 8] - A[ 9]*t[ 9] +
    A[10]*t[10] - A[11]*t[11];
  s[ 9] =
    A[ 6]*t[ 7] + A[ 7]*t[ 6] +
    A[ 8]*t[ 9] + A[ 9]*t[ 8] +
    A[10]*t[11] + A[11]*t[10];

  // s5 = U20*t3 + U21*t4 + U22*t5
  s[10] =
    A[12]*t[ 6] - A[13]*t[ 7] +
    A[14]*t[ 8] - A[15]*t[ 9] +
    A[16]*t[10] - A[17]*t[11];
  s[11] =
    A[12]*t[ 7] + A[13]*t[ 6] +
    A[14]*t[ 9] + A[15]*t[ 8] +
    A[16]*t[11] + A[17]*t[10];


  // s6 = U00*t6 + U01*t7 + U02*t8
  s[12] =
    A[ 0]*t[12] - A[ 1]*t[13] +
    A[ 2]*t[14] - A[ 3]*t[15] +
    A[ 4]*t[16] - A[ 5]*t[17];
  s[13] =
    A[ 0]*t[13] + A[ 1]*t[12] +
    A[ 2]*t[15] + A[ 3]*t[14] +
    A[ 4]*t[17] + A[ 5]*t[16];

  // s7 = U10*t6 + U11*t7 + U12*t8
  s[14] =
    A[ 6]*t[12] - A[ 7]*t[13] +
    A[ 8]*t[14] - A[ 9]*t[15] +
    A[10]*t[16] - A[11]*t[17];
  s[15] =
    A[ 6]*t[13] + A[ 7]*t[12] +
    A[ 8]*t[15] + A[ 9]*t[14] +
    A[10]*t[17] + A[11]*t[16];

  // s8 = U20*t6 + U21*t7 + U22*t8
  s[16] =
    A[12]*t[12] - A[13]*t[13] +
    A[14]*t[14] - A[15]*t[15] +
    A[16]*t[16] - A[17]*t[17];
  s[17] =
    A[12]*t[13] + A[13]*t[12] +
    A[14]*t[15] + A[15]*t[14] +
    A[16]*t[17] + A[17]*t[16];


  // s9 = U00*t9 + U01*t10 + U02*t11
  s[18] =
    A[ 0]*t[18] - A[ 1]*t[19] +
    A[ 2]*t[20] - A[ 3]*t[21] +
    A[ 4]*t[22] - A[ 5]*t[23];
  s[19] =
    A[ 0]*t[19] + A[ 1]*t[18] +
    A[ 2]*t[21] + A[ 3]*t[20] +
    A[ 4]*t[23] + A[ 5]*t[22];

  // s10 = U10*t9 + U11*t10 + U12*t11
  s[20] =
    A[ 6]*t[18] - A[ 7]*t[19] +
    A[ 8]*t[20] - A[ 9]*t[21] +
    A[10]*t[22] - A[11]*t[23];
  s[21] =
    A[ 6]*t[19] + A[ 7]*t[18] +
    A[ 8]*t[21] + A[ 9]*t[20] +
    A[10]*t[23] + A[11]*t[22];

  // s11 = U20*t9 + U21*t10 + U22*t11
  s[22] =
    A[12]*t[18] - A[13]*t[19] +
    A[14]*t[20] - A[15]*t[21] +
    A[16]*t[22] - A[17]*t[23];
  s[23] =
    A[12]*t[19] + A[13]*t[18] +
    A[14]*t[21] + A[15]*t[20] +
    A[16]*t[23] + A[17]*t[22];
}



// ********************



// s = A^\dagger * t.

inline void fv_eq_cm_dag_ti_fv(double *s, const double *A, const double *t)
{
  // s0 = U00*t0 + U01*t1 + U02*t2
  s[ 0] =
    A[ 0]*t[ 0] + A[ 1]*t[ 1] +
    A[ 6]*t[ 2] + A[ 7]*t[ 3] +
    A[12]*t[ 4] + A[13]*t[ 5];
  s[ 1] =
    A[ 0]*t[ 1] - A[ 1]*t[ 0] +
    A[ 6]*t[ 3] - A[ 7]*t[ 2] +
    A[12]*t[ 5] - A[13]*t[ 4];

  // s1 = U10*t0 + U11*t1 + U12*t2
  s[ 2] =
    A[ 2]*t[ 0] + A[ 3]*t[ 1] +
    A[ 8]*t[ 2] + A[ 9]*t[ 3] +
    A[14]*t[ 4] + A[15]*t[ 5];
  s[ 3] =
    A[ 2]*t[ 1] - A[ 3]*t[ 0] +
    A[ 8]*t[ 3] - A[ 9]*t[ 2] +
    A[14]*t[ 5] - A[15]*t[ 4];

  // s2 = U20*t0 + U21*t1 + U22*t2
  s[ 4] =
    A[ 4]*t[ 0] + A[ 5]*t[ 1] +
    A[10]*t[ 2] + A[11]*t[ 3] +
    A[16]*t[ 4] + A[17]*t[ 5];
  s[ 5] =
    A[ 4]*t[ 1] - A[ 5]*t[ 0] +
    A[10]*t[ 3] - A[11]*t[ 2] +
    A[16]*t[ 5] - A[17]*t[ 4];


  // s3 = U00*t3 + U01*t4 + U02*t5
  s[ 6] =
    A[ 0]*t[ 6] + A[ 1]*t[ 7] +
    A[ 6]*t[ 8] + A[ 7]*t[ 9] +
    A[12]*t[10] + A[13]*t[11];
  s[ 7] =
    A[ 0]*t[ 7] - A[ 1]*t[ 6] +
    A[ 6]*t[ 9] - A[ 7]*t[ 8] +
    A[12]*t[11] - A[13]*t[10];

  // s4 = U10*t3 + U11*t4 + U12*t5
  s[ 8] =
    A[ 2]*t[ 6] + A[ 3]*t[ 7] +
    A[ 8]*t[ 8] + A[ 9]*t[ 9] +
    A[14]*t[10] + A[15]*t[11];
  s[ 9] =
    A[ 2]*t[ 7] - A[ 3]*t[ 6] +
    A[ 8]*t[ 9] - A[ 9]*t[ 8] +
    A[14]*t[11] - A[15]*t[10];

  // s5 = U20*t3 + U21*t4 + U22*t5
  s[10] =
    A[ 4]*t[ 6] + A[ 5]*t[ 7] +
    A[10]*t[ 8] + A[11]*t[ 9] +
    A[16]*t[10] + A[17]*t[11];
  s[11] =
    A[ 4]*t[ 7] - A[ 5]*t[ 6] +
    A[10]*t[ 9] - A[11]*t[ 8] +
    A[16]*t[11] - A[17]*t[10];


  // s6 = U00*t6 + U01*t7 + U02*t8
  s[12] =
    A[ 0]*t[12] + A[ 1]*t[13] +
    A[ 6]*t[14] + A[ 7]*t[15] +
    A[12]*t[16] + A[13]*t[17];
  s[13] =
    A[ 0]*t[13] - A[ 1]*t[12] +
    A[ 6]*t[15] - A[ 7]*t[14] +
    A[12]*t[17] - A[13]*t[16];

  // s7 = U10*t6 + U11*t7 + U12*t8
  s[14] =
    A[ 2]*t[12] + A[ 3]*t[13] +
    A[ 8]*t[14] + A[ 9]*t[15] +
    A[14]*t[16] + A[15]*t[17];
  s[15] =
    A[ 2]*t[13] - A[ 3]*t[12] +
    A[ 8]*t[15] - A[ 9]*t[14] +
    A[14]*t[17] - A[15]*t[16];

  // s8 = U20*t6 + U21*t7 + U22*t8
  s[16] =
    A[ 4]*t[12] + A[ 5]*t[13] +
    A[10]*t[14] + A[11]*t[15] +
    A[16]*t[16] + A[17]*t[17];
  s[17] =
    A[ 4]*t[13] - A[ 5]*t[12] +
    A[10]*t[15] - A[11]*t[14] +
    A[16]*t[17] - A[17]*t[16];


  // s9 = U00*t9 + U01*t10 + U02*t11
  s[18] =
    A[ 0]*t[18] + A[ 1]*t[19] +
    A[ 6]*t[20] + A[ 7]*t[21] +
    A[12]*t[22] + A[13]*t[23];
  s[19] =
    A[ 0]*t[19] - A[ 1]*t[18] +
    A[ 6]*t[21] - A[ 7]*t[20] +
    A[12]*t[23] - A[13]*t[22];

  // s10 = U10*t9 + U11*t10 + U12*t11
  s[20] =
    A[ 2]*t[18] + A[ 3]*t[19] +
    A[ 8]*t[20] + A[ 9]*t[21] +
    A[14]*t[22] + A[15]*t[23];
  s[21] =
    A[ 2]*t[19] - A[ 3]*t[18] +
    A[ 8]*t[21] - A[ 9]*t[20] +
    A[14]*t[23] - A[15]*t[22];

  // s11 = U20*t9 + U21*t10 + U22*t11
  s[22] =
    A[ 4]*t[18] + A[ 5]*t[19] +
    A[10]*t[20] + A[11]*t[21] +
    A[16]*t[22] + A[17]*t[23];
  s[23] =
    A[ 4]*t[19] - A[ 5]*t[18] +
    A[10]*t[21] - A[11]*t[20] +
    A[16]*t[23] - A[17]*t[22];
}



// ********************



// s = s + A * t.

inline void fv_pl_eq_cm_ti_fv(double *s, const double *A, const double *t)
{
  // s0 = s0 + U00*t0 + U01*t1 + U02*t2
  s[ 0] +=
    A[ 0]*t[ 0] - A[ 1]*t[ 1] +
    A[ 2]*t[ 2] - A[ 3]*t[ 3] +
    A[ 4]*t[ 4] - A[ 5]*t[ 5];
  s[ 1] +=
    A[ 0]*t[ 1] + A[ 1]*t[ 0] +
    A[ 2]*t[ 3] + A[ 3]*t[ 2] +
    A[ 4]*t[ 5] + A[ 5]*t[ 4];

  // s1 = s1 + U10*t0 + U11*t1 + U12*t2
  s[ 2] +=
    A[ 6]*t[ 0] - A[ 7]*t[ 1] +
    A[ 8]*t[ 2] - A[ 9]*t[ 3] +
    A[10]*t[ 4] - A[11]*t[ 5];
  s[ 3] +=
    A[ 6]*t[ 1] + A[ 7]*t[ 0] +
    A[ 8]*t[ 3] + A[ 9]*t[ 2] +
    A[10]*t[ 5] + A[11]*t[ 4];

  // s2 = s2 + U20*t0 + U21*t1 + U22*t2
  s[ 4] +=
    A[12]*t[ 0] - A[13]*t[ 1] +
    A[14]*t[ 2] - A[15]*t[ 3] +
    A[16]*t[ 4] - A[17]*t[ 5];
  s[ 5] +=
    A[12]*t[ 1] + A[13]*t[ 0] +
    A[14]*t[ 3] + A[15]*t[ 2] +
    A[16]*t[ 5] + A[17]*t[ 4];


  // s3 = s3 + U00*t3 + U01*t4 + U02*t5
  s[ 6] +=
    A[ 0]*t[ 6] - A[ 1]*t[ 7] +
    A[ 2]*t[ 8] - A[ 3]*t[ 9] +
    A[ 4]*t[10] - A[ 5]*t[11];
  s[ 7] +=
    A[ 0]*t[ 7] + A[ 1]*t[ 6] +
    A[ 2]*t[ 9] + A[ 3]*t[ 8] +
    A[ 4]*t[11] + A[ 5]*t[10];

  // s4 = s4 + U10*t3 + U11*t4 + U12*t5
  s[ 8] +=
    A[ 6]*t[ 6] - A[ 7]*t[ 7] +
    A[ 8]*t[ 8] - A[ 9]*t[ 9] +
    A[10]*t[10] - A[11]*t[11];
  s[ 9] +=
    A[ 6]*t[ 7] + A[ 7]*t[ 6] +
    A[ 8]*t[ 9] + A[ 9]*t[ 8] +
    A[10]*t[11] + A[11]*t[10];

  // s5 = s5 + U20*t3 + U21*t4 + U22*t5
  s[10] +=
    A[12]*t[ 6] - A[13]*t[ 7] +
    A[14]*t[ 8] - A[15]*t[ 9] +
    A[16]*t[10] - A[17]*t[11];
  s[11] +=
    A[12]*t[ 7] + A[13]*t[ 6] +
    A[14]*t[ 9] + A[15]*t[ 8] +
    A[16]*t[11] + A[17]*t[10];


  // s6 = s6 + U00*t6 + U01*t7 + U02*t8
  s[12] +=
    A[ 0]*t[12] - A[ 1]*t[13] +
    A[ 2]*t[14] - A[ 3]*t[15] +
    A[ 4]*t[16] - A[ 5]*t[17];
  s[13] +=
    A[ 0]*t[13] + A[ 1]*t[12] +
    A[ 2]*t[15] + A[ 3]*t[14] +
    A[ 4]*t[17] + A[ 5]*t[16];

  // s7 = s7 + U10*t6 + U11*t7 + U12*t8
  s[14] +=
    A[ 6]*t[12] - A[ 7]*t[13] +
    A[ 8]*t[14] - A[ 9]*t[15] +
    A[10]*t[16] - A[11]*t[17];
  s[15] +=
    A[ 6]*t[13] + A[ 7]*t[12] +
    A[ 8]*t[15] + A[ 9]*t[14] +
    A[10]*t[17] + A[11]*t[16];

  // s8 = s8 + U20*t6 + U21*t7 + U22*t8
  s[16] +=
    A[12]*t[12] - A[13]*t[13] +
    A[14]*t[14] - A[15]*t[15] +
    A[16]*t[16] - A[17]*t[17];
  s[17] +=
    A[12]*t[13] + A[13]*t[12] +
    A[14]*t[15] + A[15]*t[14] +
    A[16]*t[17] + A[17]*t[16];


  // s9 = s9 + U00*t9 + U01*t10 + U02*t11
  s[18] +=
    A[ 0]*t[18] - A[ 1]*t[19] +
    A[ 2]*t[20] - A[ 3]*t[21] +
    A[ 4]*t[22] - A[ 5]*t[23];
  s[19] +=
    A[ 0]*t[19] + A[ 1]*t[18] +
    A[ 2]*t[21] + A[ 3]*t[20] +
    A[ 4]*t[23] + A[ 5]*t[22];

  // s10 = s10 + U10*t9 + U11*t10 + U12*t11
  s[20] +=
    A[ 6]*t[18] - A[ 7]*t[19] +
    A[ 8]*t[20] - A[ 9]*t[21] +
    A[10]*t[22] - A[11]*t[23];
  s[21] +=
    A[ 6]*t[19] + A[ 7]*t[18] +
    A[ 8]*t[21] + A[ 9]*t[20] +
    A[10]*t[23] + A[11]*t[22];

  // s11 = s11 + U20*t9 + U21*t10 + U22*t11
  s[22] +=
    A[12]*t[18] - A[13]*t[19] +
    A[14]*t[20] - A[15]*t[21] +
    A[16]*t[22] - A[17]*t[23];
  s[23] +=
    A[12]*t[19] + A[13]*t[18] +
    A[14]*t[21] + A[15]*t[20] +
    A[16]*t[23] + A[17]*t[22];
}



// ********************



// s = s - A * t.

inline void fv_mi_eq_cm_ti_fv(double *s, const double *A, const double *t)
{
  // s0 = s0 - U00*t0 - U01*t1 - U02*t2
  s[ 0] -=
    (A[ 0]*t[ 0] - A[ 1]*t[ 1] +
     A[ 2]*t[ 2] - A[ 3]*t[ 3] +
     A[ 4]*t[ 4] - A[ 5]*t[ 5]);
  s[ 1] -=
    (A[ 0]*t[ 1] + A[ 1]*t[ 0] +
     A[ 2]*t[ 3] + A[ 3]*t[ 2] +
     A[ 4]*t[ 5] + A[ 5]*t[ 4]);

  // s1 = s1 - U10*t0 - U11*t1 - U12*t2
  s[ 2] -=
    (A[ 6]*t[ 0] - A[ 7]*t[ 1] +
     A[ 8]*t[ 2] - A[ 9]*t[ 3] +
     A[10]*t[ 4] - A[11]*t[ 5]);
  s[ 3] -=
    (A[ 6]*t[ 1] + A[ 7]*t[ 0] +
     A[ 8]*t[ 3] + A[ 9]*t[ 2] +
     A[10]*t[ 5] + A[11]*t[ 4]);

  // s2 = s2 - U20*t0 - U21*t1 - U22*t2
  s[ 4] -=
    (A[12]*t[ 0] - A[13]*t[ 1] +
     A[14]*t[ 2] - A[15]*t[ 3] +
     A[16]*t[ 4] - A[17]*t[ 5]);
  s[ 5] -=
    (A[12]*t[ 1] + A[13]*t[ 0] +
     A[14]*t[ 3] + A[15]*t[ 2] +
     A[16]*t[ 5] + A[17]*t[ 4]);


  // s3 = s3 - U00*t3 - U01*t4 - U02*t5
  s[ 6] -=
    (A[ 0]*t[ 6] - A[ 1]*t[ 7] +
     A[ 2]*t[ 8] - A[ 3]*t[ 9] +
     A[ 4]*t[10] - A[ 5]*t[11]);
  s[ 7] -=
    (A[ 0]*t[ 7] + A[ 1]*t[ 6] +
     A[ 2]*t[ 9] + A[ 3]*t[ 8] +
     A[ 4]*t[11] + A[ 5]*t[10]);

  // s4 = s4 - U10*t3 - U11*t4 - U12*t5
  s[ 8] -=
    (A[ 6]*t[ 6] - A[ 7]*t[ 7] +
     A[ 8]*t[ 8] - A[ 9]*t[ 9] +
     A[10]*t[10] - A[11]*t[11]);
  s[ 9] -=
    (A[ 6]*t[ 7] + A[ 7]*t[ 6] +
     A[ 8]*t[ 9] + A[ 9]*t[ 8] +
     A[10]*t[11] + A[11]*t[10]);

  // s5 = s5 - U20*t3 - U21*t4 - U22*t5
  s[10] -=
    (A[12]*t[ 6] - A[13]*t[ 7] +
     A[14]*t[ 8] - A[15]*t[ 9] +
     A[16]*t[10] - A[17]*t[11]);
  s[11] -=
    (A[12]*t[ 7] + A[13]*t[ 6] +
     A[14]*t[ 9] + A[15]*t[ 8] +
     A[16]*t[11] + A[17]*t[10]);


  // s6 = s6 - U00*t6 - U01*t7 - U02*t8
  s[12] -=
    (A[ 0]*t[12] - A[ 1]*t[13] +
     A[ 2]*t[14] - A[ 3]*t[15] +
     A[ 4]*t[16] - A[ 5]*t[17]);
  s[13] -=
    (A[ 0]*t[13] + A[ 1]*t[12] +
     A[ 2]*t[15] + A[ 3]*t[14] +
     A[ 4]*t[17] + A[ 5]*t[16]);

  // s7 = s7 - U10*t6 - U11*t7 - U12*t8
  s[14] -=
    (A[ 6]*t[12] - A[ 7]*t[13] +
     A[ 8]*t[14] - A[ 9]*t[15] +
     A[10]*t[16] - A[11]*t[17]);
  s[15] -=
    (A[ 6]*t[13] + A[ 7]*t[12] +
     A[ 8]*t[15] + A[ 9]*t[14] +
     A[10]*t[17] + A[11]*t[16]);

  // s8 = s8 - U20*t6 - U21*t7 - U22*t8
  s[16] -=
    (A[12]*t[12] - A[13]*t[13] +
     A[14]*t[14] - A[15]*t[15] +
     A[16]*t[16] - A[17]*t[17]);
  s[17] -=
    (A[12]*t[13] + A[13]*t[12] +
     A[14]*t[15] + A[15]*t[14] +
     A[16]*t[17] + A[17]*t[16]);


  // s9 = s9 - U00*t9 - U01*t10 - U02*t11
  s[18] -=
    (A[ 0]*t[18] - A[ 1]*t[19] +
     A[ 2]*t[20] - A[ 3]*t[21] +
     A[ 4]*t[22] - A[ 5]*t[23]);
  s[19] -=
    (A[ 0]*t[19] + A[ 1]*t[18] +
     A[ 2]*t[21] + A[ 3]*t[20] +
     A[ 4]*t[23] + A[ 5]*t[22]);

  // s10 = s10 - U10*t9 - U11*t10 - U12*t11
  s[20] -=
    (A[ 6]*t[18] - A[ 7]*t[19] +
     A[ 8]*t[20] - A[ 9]*t[21] +
     A[10]*t[22] - A[11]*t[23]);
  s[21] -=
    (A[ 6]*t[19] + A[ 7]*t[18] +
     A[ 8]*t[21] + A[ 9]*t[20] +
     A[10]*t[23] + A[11]*t[22]);

  // s11 = s11 - U20*t9 - U21*t10 - U22*t11
  s[22] -=
    (A[12]*t[18] - A[13]*t[19] +
     A[14]*t[20] - A[15]*t[21] +
     A[16]*t[22] - A[17]*t[23]);
  s[23] -=
    (A[12]*t[19] + A[13]*t[18] +
     A[14]*t[21] + A[15]*t[20] +
     A[16]*t[23] + A[17]*t[22]);
}



// ********************



// s = s + A^\dagger * t.

inline void fv_pl_eq_cm_dag_ti_fv(double *s, const double *A, const double *t)
{
  // s0 = s0 + U_00*t0 + U_10*t1 + U_20*t2
  s[ 0] +=
    A[ 0]*t[ 0] + A[ 1]*t[ 1] +
    A[ 6]*t[ 2] + A[ 7]*t[ 3] +
    A[12]*t[ 4] + A[13]*t[ 5];
  s[ 1] +=
    A[ 0]*t[ 1] - A[ 1]*t[ 0] +
    A[ 6]*t[ 3] - A[ 7]*t[ 2] +
    A[12]*t[ 5] - A[13]*t[ 4];

  // s1 = s1 + U_01*t0 + U_11*t1 + U_21*t2
  s[ 2] +=
    A[ 2]*t[ 0] + A[ 3]*t[ 1] +
    A[ 8]*t[ 2] + A[ 9]*t[ 3] +
    A[14]*t[ 4] + A[15]*t[ 5];
  s[ 3] +=
    A[ 2]*t[ 1] - A[ 3]*t[ 0] +
    A[ 8]*t[ 3] - A[ 9]*t[ 2] +
    A[14]*t[ 5] - A[15]*t[ 4];

  // s2 = s2 + U_02*t0 + U_12*t1 + U_22*t2
  s[ 4] +=
    A[ 4]*t[ 0] + A[ 5]*t[ 1] +
    A[10]*t[ 2] + A[11]*t[ 3] +
    A[16]*t[ 4] + A[17]*t[ 5];
  s[ 5] +=
    A[ 4]*t[ 1] - A[ 5]*t[ 0] +
    A[10]*t[ 3] - A[11]*t[ 2] +
    A[16]*t[ 5] - A[17]*t[ 4];


  // s3 = s3 + U_00*t3 + U_10*t4 + U_20*t5
  s[ 6] +=
    A[ 0]*t[ 6] + A[ 1]*t[ 7] +
    A[ 6]*t[ 8] + A[ 7]*t[ 9] +
    A[12]*t[10] + A[13]*t[11];
  s[ 7] +=
    A[ 0]*t[ 7] - A[ 1]*t[ 6] +
    A[ 6]*t[ 9] - A[ 7]*t[ 8] +
    A[12]*t[11] - A[13]*t[10];

  // s4 = s4 + U_01*t3 + U_11*t4 + U_21*t5
  s[ 8] +=
    A[ 2]*t[ 6] + A[ 3]*t[ 7] +
    A[ 8]*t[ 8] + A[ 9]*t[ 9] +
    A[14]*t[10] + A[15]*t[11];
  s[ 9] +=
    A[ 2]*t[ 7] - A[ 3]*t[ 6] +
    A[ 8]*t[ 9] - A[ 9]*t[ 8] +
    A[14]*t[11] - A[15]*t[10];

  // s5 = s5 + U_02*t3 + U_12*t4 + U_22*t5
  s[10] +=
    A[ 4]*t[ 6] + A[ 5]*t[ 7] +
    A[10]*t[ 8] + A[11]*t[ 9] +
    A[16]*t[10] + A[17]*t[11];
  s[11] +=
    A[ 4]*t[ 7] - A[ 5]*t[ 6] +
    A[10]*t[ 9] - A[11]*t[ 8] +
    A[16]*t[11] - A[17]*t[10];


  // s6 = s6 + U_00*t6 + U_10*t7 + U_20*t8
  s[12] +=
    A[ 0]*t[12] + A[ 1]*t[13] +
    A[ 6]*t[14] + A[ 7]*t[15] +
    A[12]*t[16] + A[13]*t[17];
  s[13] +=
    A[ 0]*t[13] - A[ 1]*t[12] +
    A[ 6]*t[15] - A[ 7]*t[14] +
    A[12]*t[17] - A[13]*t[16];

  // s7 = s7 + U_01*t6 + U_11*t7 + U_21*t8
  s[14] +=
    A[ 2]*t[12] + A[ 3]*t[13] +
    A[ 8]*t[14] + A[ 9]*t[15] +
    A[14]*t[16] + A[15]*t[17];
  s[15] +=
    A[ 2]*t[13] - A[ 3]*t[12] +
    A[ 8]*t[15] - A[ 9]*t[14] +
    A[14]*t[17] - A[15]*t[16];

  // s8 = s8 + U_02*t6 + U_12*t7 + U_22*t8
  s[16] +=
    A[ 4]*t[12] + A[ 5]*t[13] +
    A[10]*t[14] + A[11]*t[15] +
    A[16]*t[16] + A[17]*t[17];
  s[17] +=
    A[ 4]*t[13] - A[ 5]*t[12] +
    A[10]*t[15] - A[11]*t[14] +
    A[16]*t[17] - A[17]*t[16];


  // s9 = s9 + U_00*t9 + U_10*t10 + U_20*t11
  s[18] +=
    A[ 0]*t[18] + A[ 1]*t[19] +
    A[ 6]*t[20] + A[ 7]*t[21] +
    A[12]*t[22] + A[13]*t[23];
  s[19] +=
    A[ 0]*t[19] - A[ 1]*t[18] +
    A[ 6]*t[21] - A[ 7]*t[20] +
    A[12]*t[23] - A[13]*t[22];

  // s10 = s10 + U_01*t9 + U_11*t10 + U_21*t11
  s[20] +=
    A[ 2]*t[18] + A[ 3]*t[19] +
    A[ 8]*t[20] + A[ 9]*t[21] +
    A[14]*t[22] + A[15]*t[23];
  s[21] +=
    A[ 2]*t[19] - A[ 3]*t[18] +
    A[ 8]*t[21] - A[ 9]*t[20] +
    A[14]*t[23] - A[15]*t[22];

  // s11 = s11 + U_02*t9 + U_12*t10 + U_22*t11
  s[22] +=
    A[ 4]*t[18] + A[ 5]*t[19] +
    A[10]*t[20] + A[11]*t[21] +
    A[16]*t[22] + A[17]*t[23];
  s[23] +=
    A[ 4]*t[19] - A[ 5]*t[18] +
    A[10]*t[21] - A[11]*t[20] +
    A[16]*t[23] - A[17]*t[22];
}



// ********************



// s = s - A^\dagger * t.

inline void fv_mi_eq_cm_dag_ti_fv(double *s, const double *A, const double *t)
{
  // s0 = s0 - U_00*t0 - U_10*t1 - U_20*t2
  s[ 0] -=
    (A[ 0]*t[ 0] + A[ 1]*t[ 1] +
     A[ 6]*t[ 2] + A[ 7]*t[ 3] +
     A[12]*t[ 4] + A[13]*t[ 5]);
  s[ 1] -=
    (A[ 0]*t[ 1] - A[ 1]*t[ 0] +
     A[ 6]*t[ 3] - A[ 7]*t[ 2] +
     A[12]*t[ 5] - A[13]*t[ 4]);

  // s1 = s1 - U_01*t0 - U_11*t1 - U_21*t2
  s[ 2] -=
    (A[ 2]*t[ 0] + A[ 3]*t[ 1] +
     A[ 8]*t[ 2] + A[ 9]*t[ 3] +
     A[14]*t[ 4] + A[15]*t[ 5]);
  s[ 3] -=
    (A[ 2]*t[ 1] - A[ 3]*t[ 0] +
     A[ 8]*t[ 3] - A[ 9]*t[ 2] +
     A[14]*t[ 5] - A[15]*t[ 4]);

  // s2 = s2 - U_02*t0 - U_12*t1 - U_22*t2
  s[ 4] -=
    (A[ 4]*t[ 0] + A[ 5]*t[ 1] +
     A[10]*t[ 2] + A[11]*t[ 3] +
     A[16]*t[ 4] + A[17]*t[ 5]);
  s[ 5] -=
    (A[ 4]*t[ 1] - A[ 5]*t[ 0] +
     A[10]*t[ 3] - A[11]*t[ 2] +
     A[16]*t[ 5] - A[17]*t[ 4]);


  // s3 = s3 - U_00*t3 - U_10*t4 - U_20*t5
  s[ 6] -=
    (A[ 0]*t[ 6] + A[ 1]*t[ 7] +
     A[ 6]*t[ 8] + A[ 7]*t[ 9] +
     A[12]*t[10] + A[13]*t[11]);
  s[ 7] -=
    (A[ 0]*t[ 7] - A[ 1]*t[ 6] +
     A[ 6]*t[ 9] - A[ 7]*t[ 8] +
     A[12]*t[11] - A[13]*t[10]);

  // s4 = s4 - U_01*t3 - U_11*t4 - U_21*t5
  s[ 8] -=
    (A[ 2]*t[ 6] + A[ 3]*t[ 7] +
     A[ 8]*t[ 8] + A[ 9]*t[ 9] +
     A[14]*t[10] + A[15]*t[11]);
  s[ 9] -=
    (A[ 2]*t[ 7] - A[ 3]*t[ 6] +
     A[ 8]*t[ 9] - A[ 9]*t[ 8] +
     A[14]*t[11] - A[15]*t[10]);

  // s5 = s5 - U_02*t3 - U_12*t4 - U_22*t5
  s[10] -=
    (A[ 4]*t[ 6] + A[ 5]*t[ 7] +
     A[10]*t[ 8] + A[11]*t[ 9] +
     A[16]*t[10] + A[17]*t[11]);
  s[11] -=
    (A[ 4]*t[ 7] - A[ 5]*t[ 6] +
     A[10]*t[ 9] - A[11]*t[ 8] +
     A[16]*t[11] - A[17]*t[10]);


  // s6 = s6 - U_00*t6 - U_10*t7 - U_20*t8
  s[12] -=
    (A[ 0]*t[12] + A[ 1]*t[13] +
     A[ 6]*t[14] + A[ 7]*t[15] +
     A[12]*t[16] + A[13]*t[17]);
  s[13] -=
    (A[ 0]*t[13] - A[ 1]*t[12] +
     A[ 6]*t[15] - A[ 7]*t[14] +
     A[12]*t[17] - A[13]*t[16]);

  // s7 = s7 - U_01*t6 - U_11*t7 - U_21*t8
  s[14] -=
    (A[ 2]*t[12] + A[ 3]*t[13] +
     A[ 8]*t[14] + A[ 9]*t[15] +
     A[14]*t[16] + A[15]*t[17]);
  s[15] -=
    (A[ 2]*t[13] - A[ 3]*t[12] +
     A[ 8]*t[15] - A[ 9]*t[14] +
     A[14]*t[17] - A[15]*t[16]);

  // s8 = s8 - U_02*t6 - U_12*t7 - U_22*t8
  s[16] -=
    (A[ 4]*t[12] + A[ 5]*t[13] +
     A[10]*t[14] + A[11]*t[15] +
     A[16]*t[16] + A[17]*t[17]);
  s[17] -=
    (A[ 4]*t[13] - A[ 5]*t[12] +
     A[10]*t[15] - A[11]*t[14] +
     A[16]*t[17] - A[17]*t[16]);


  // s9 = s9 - U_00*t9 - U_10*t10 - U_20*t11
  s[18] -=
    (A[ 0]*t[18] + A[ 1]*t[19] +
     A[ 6]*t[20] + A[ 7]*t[21] +
     A[12]*t[22] + A[13]*t[23]);
  s[19] -=
    (A[ 0]*t[19] - A[ 1]*t[18] +
     A[ 6]*t[21] - A[ 7]*t[20] +
     A[12]*t[23] - A[13]*t[22]);

  // s10 = s10 - U_01*t9 - U_11*t10 - U_21*t11
  s[20] -=
    (A[ 2]*t[18] + A[ 3]*t[19] +
     A[ 8]*t[20] + A[ 9]*t[21] +
     A[14]*t[22] + A[15]*t[23]);
  s[21] -=
    (A[ 2]*t[19] - A[ 3]*t[18] +
     A[ 8]*t[21] - A[ 9]*t[20] +
     A[14]*t[23] - A[15]*t[22]);

  // s11 = s11 - U_02*t9 - U_12*t10 - U_22*t11
  s[22] -=
    (A[ 4]*t[18] + A[ 5]*t[19] +
     A[10]*t[20] + A[11]*t[21] +
     A[16]*t[22] + A[17]*t[23]);
  s[23] -=
    (A[ 4]*t[19] - A[ 5]*t[18] +
     A[10]*t[21] - A[11]*t[20] +
     A[16]*t[23] - A[17]*t[22]);
}



// ********************



// s = gamma * t.

inline void fv_eq_gamma_ti_fv(double *s, int gamma_index, const double *t)
{
  s[ 0] = t[gamma_permutation[gamma_index][ 0]] * gamma_sign[gamma_index][ 0];
  s[ 1] = t[gamma_permutation[gamma_index][ 1]] * gamma_sign[gamma_index][ 1];
  s[ 2] = t[gamma_permutation[gamma_index][ 2]] * gamma_sign[gamma_index][ 2];
  s[ 3] = t[gamma_permutation[gamma_index][ 3]] * gamma_sign[gamma_index][ 3];
  s[ 4] = t[gamma_permutation[gamma_index][ 4]] * gamma_sign[gamma_index][ 4];
  s[ 5] = t[gamma_permutation[gamma_index][ 5]] * gamma_sign[gamma_index][ 5];

  s[ 6] = t[gamma_permutation[gamma_index][ 6]] * gamma_sign[gamma_index][ 6];
  s[ 7] = t[gamma_permutation[gamma_index][ 7]] * gamma_sign[gamma_index][ 7];
  s[ 8] = t[gamma_permutation[gamma_index][ 8]] * gamma_sign[gamma_index][ 8];
  s[ 9] = t[gamma_permutation[gamma_index][ 9]] * gamma_sign[gamma_index][ 9];
  s[10] = t[gamma_permutation[gamma_index][10]] * gamma_sign[gamma_index][10];
  s[11] = t[gamma_permutation[gamma_index][11]] * gamma_sign[gamma_index][11];

  s[12] = t[gamma_permutation[gamma_index][12]] * gamma_sign[gamma_index][12];
  s[13] = t[gamma_permutation[gamma_index][13]] * gamma_sign[gamma_index][13];
  s[14] = t[gamma_permutation[gamma_index][14]] * gamma_sign[gamma_index][14];
  s[15] = t[gamma_permutation[gamma_index][15]] * gamma_sign[gamma_index][15];
  s[16] = t[gamma_permutation[gamma_index][16]] * gamma_sign[gamma_index][16];
  s[17] = t[gamma_permutation[gamma_index][17]] * gamma_sign[gamma_index][17];

  s[18] = t[gamma_permutation[gamma_index][18]] * gamma_sign[gamma_index][18];
  s[19] = t[gamma_permutation[gamma_index][19]] * gamma_sign[gamma_index][19];
  s[20] = t[gamma_permutation[gamma_index][20]] * gamma_sign[gamma_index][20];
  s[21] = t[gamma_permutation[gamma_index][21]] * gamma_sign[gamma_index][21];
  s[22] = t[gamma_permutation[gamma_index][22]] * gamma_sign[gamma_index][22];
  s[23] = t[gamma_permutation[gamma_index][23]] * gamma_sign[gamma_index][23];
}

inline void fv_eq_gamma_ti_fv(float *s, int gamma_index, const float *t)
{
  s[ 0] = t[gamma_permutation[gamma_index][ 0]] * gamma_sign[gamma_index][ 0];
  s[ 1] = t[gamma_permutation[gamma_index][ 1]] * gamma_sign[gamma_index][ 1];
  s[ 2] = t[gamma_permutation[gamma_index][ 2]] * gamma_sign[gamma_index][ 2];
  s[ 3] = t[gamma_permutation[gamma_index][ 3]] * gamma_sign[gamma_index][ 3];
  s[ 4] = t[gamma_permutation[gamma_index][ 4]] * gamma_sign[gamma_index][ 4];
  s[ 5] = t[gamma_permutation[gamma_index][ 5]] * gamma_sign[gamma_index][ 5];

  s[ 6] = t[gamma_permutation[gamma_index][ 6]] * gamma_sign[gamma_index][ 6];
  s[ 7] = t[gamma_permutation[gamma_index][ 7]] * gamma_sign[gamma_index][ 7];
  s[ 8] = t[gamma_permutation[gamma_index][ 8]] * gamma_sign[gamma_index][ 8];
  s[ 9] = t[gamma_permutation[gamma_index][ 9]] * gamma_sign[gamma_index][ 9];
  s[10] = t[gamma_permutation[gamma_index][10]] * gamma_sign[gamma_index][10];
  s[11] = t[gamma_permutation[gamma_index][11]] * gamma_sign[gamma_index][11];

  s[12] = t[gamma_permutation[gamma_index][12]] * gamma_sign[gamma_index][12];
  s[13] = t[gamma_permutation[gamma_index][13]] * gamma_sign[gamma_index][13];
  s[14] = t[gamma_permutation[gamma_index][14]] * gamma_sign[gamma_index][14];
  s[15] = t[gamma_permutation[gamma_index][15]] * gamma_sign[gamma_index][15];
  s[16] = t[gamma_permutation[gamma_index][16]] * gamma_sign[gamma_index][16];
  s[17] = t[gamma_permutation[gamma_index][17]] * gamma_sign[gamma_index][17];

  s[18] = t[gamma_permutation[gamma_index][18]] * gamma_sign[gamma_index][18];
  s[19] = t[gamma_permutation[gamma_index][19]] * gamma_sign[gamma_index][19];
  s[20] = t[gamma_permutation[gamma_index][20]] * gamma_sign[gamma_index][20];
  s[21] = t[gamma_permutation[gamma_index][21]] * gamma_sign[gamma_index][21];
  s[22] = t[gamma_permutation[gamma_index][22]] * gamma_sign[gamma_index][22];
  s[23] = t[gamma_permutation[gamma_index][23]] * gamma_sign[gamma_index][23];
}

// s = gamma * t.

inline void fv_eq_gamma_ti_fv(double *s, const double *t, int const per[], int const sig[])
{
  s[ 0] = t[per[ 0]] * sig[ 0];
  s[ 1] = t[per[ 1]] * sig[ 1];
  s[ 2] = t[per[ 2]] * sig[ 2];
  s[ 3] = t[per[ 3]] * sig[ 3];
  s[ 4] = t[per[ 4]] * sig[ 4];
  s[ 5] = t[per[ 5]] * sig[ 5];

  s[ 6] = t[per[ 6]] * sig[ 6];
  s[ 7] = t[per[ 7]] * sig[ 7];
  s[ 8] = t[per[ 8]] * sig[ 8];
  s[ 9] = t[per[ 9]] * sig[ 9];
  s[10] = t[per[10]] * sig[10];
  s[11] = t[per[11]] * sig[11];

  s[12] = t[per[12]] * sig[12];
  s[13] = t[per[13]] * sig[13];
  s[14] = t[per[14]] * sig[14];
  s[15] = t[per[15]] * sig[15];
  s[16] = t[per[16]] * sig[16];
  s[17] = t[per[17]] * sig[17];

  s[18] = t[per[18]] * sig[18];
  s[19] = t[per[19]] * sig[19];
  s[20] = t[per[20]] * sig[20];
  s[21] = t[per[21]] * sig[21];
  s[22] = t[per[22]] * sig[22];
  s[23] = t[per[23]] * sig[23];
}

inline void gamma_eq_gamma_ti_gamma(int per[], int sig[], 
				    int const per1[], int const per2[], int const sig1[], int const sig2[]) {
  int isimag1 = 0, isimag2 = 0;
  int isimag = 1;
  if(per1[0]>per1[1]) isimag1 = 1;
  if(per2[0]>per2[1]) isimag2 = 1;
  if(isimag1 && isimag2)  isimag = 0;
  else if(!isimag1 && !isimag2)  isimag = 0;
  int p1[4], p2[4], p[4];
  int s1[4], s2[4], s[4];
  for(int mu = 0; mu < 4; mu++) {
    p1[mu] = per1[6*mu]/6;
    p2[mu] = per2[6*mu]/6;
    if(isimag1) s1[mu] = -sig1[6*mu];
    else s1[mu] = sig1[6*mu];
    if(isimag2) s2[mu] = -sig2[6*mu];
    else s2[mu] = sig2[6*mu];
  }

  for(int mu = 0; mu < 4; mu++) {
    for(int nu = 0; nu < 4; nu++) {
      if(p1[mu] == p2[nu]){
	p[mu] = nu;
	s[mu] = s1[mu]*s2[p2[nu]];
      }
    }
  }
  if(!isimag) {
    for(int mu = 0; mu < 4; mu++) {
      for(int i = 0; i < 6; i++) {
	per[6*mu+i] = p[mu]*6 + i;
	sig[6*mu+i] = s[mu];
	if(isimag1 && isimag2) {
	  sig[6*mu+i] = -sig[6*mu+i];
	}
      }
    }
  }
  else {
    for(int mu = 0; mu < 4; mu++) {
      for(int i = 0; i < 6; i+=2) {
	per[6*mu+i] = p[mu]*6 + i + 1;
	per[6*mu+i+1] = p[mu]*6 + i;
	sig[6*mu+i] = -s[mu];
	sig[6*mu+i+1] = s[mu];
      }
    }
  }
}


inline void i_times_gamma(int per[], int sig[]) {
  int tmp;
  for(int i = 0; i < 24; i+=2){
    tmp = per[i];
    per[i] = per[i+1];
    per[i+1] = tmp;
    tmp = sig[i];
    sig[i] = - sig[i+1];
    sig[i+1] = tmp;
  }
  return;
}

inline void minus_i_times_gamma(int per[], int sig[]) {
  int tmp;
  for(int i = 0; i < 24; i+=2){
    tmp = per[i];
    per[i] = per[i+1];
    per[i+1] = tmp;
    tmp = sig[i];
    sig[i] = sig[i+1];
    sig[i+1] = - tmp;
  }
  return;
}

inline void minus_gamma(int per[], int sig[]) {
  for(int i = 0; i < 24; i++) {
    sig[i] *= -1;
  }
  return;
}

// in our gamma-matrix convention this should be 
// (-g0-g5)/sqrt(2) = -(g0+g5)/sqrt(2)
inline int rotate_etmc_ukqcd(double * const field, const int T, const int L) {
  double psi0[24], psi1[24];
  for(int i = 0; i < L*L*L*T; i++) {
    // multiply with gamma_0
    fv_eq_gamma_ti_fv(psi0, 0, &field[i*24]);
    // multiply with gamma_5
    fv_eq_gamma_ti_fv(psi1, 5, &field[i*24]);
    fv_eq_fv_pl_fv(&field[i*24], psi1, psi0);
    fv_ti_eq_re(&field[i*24], -1./sqrt(2.));
  }
  return(0);
}

inline int rotate_etmc_ukqcd(float * const field, const int T, const int L) {
  float psi0[24], psi1[24];
  for(int i = 0; i < L*L*L*T; i++) {
    // multiply with gamma_0
    fv_eq_gamma_ti_fv(psi0, 0, &field[i*24]);
    // multiply with gamma_5
    fv_eq_gamma_ti_fv(psi1, 5, &field[i*24]);
    fv_eq_fv_pl_fv(&field[i*24], psi1, psi0);
    fv_ti_eq_re(&field[i*24], -1./sqrt(2.));
  }
  return(0);
}

// ********************



#endif



// ********************
