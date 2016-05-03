#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdlib.h>
// #include <getopt.h>
#ifdef OMP
# include <omp.h>
#endif

#include "fields.hh"
#include "geometry.hh"
#include "linear_algebra.hh"
#include "contract_twopoint.hh"

using namespace std;

// Performs contraction of 
// <gamma_nu gamma_5 chi^dagger gamma_5 gamma_mu phi>
// for all t
//
// this routine is for gamma_5 diagonal only!
void contract_twopoint(complex * CF, const int gamma_source, const int gamma_sink,
		       double * chi[], double * phi[],
		       const unsigned int T, const unsigned int L, 
		       const unsigned int ts) {
  double re, im;
  int tt=0;
  int isimag = (gamma_permutation[gamma_source][0])%2;
  int p[4];

  for(int mu = 0; mu < 4; mu++) {
    p[mu] = gamma_permutation[gamma_source][6*mu]/6;
  }
  double psi0[24], psi1[24];
  complex tmp;
  
  for(unsigned int t = 0; t < T; t++) {
    tmp.re = 0.;
    tmp.im = 0.;
    re = 0.;
    im = 0.;
    for(unsigned int i = 0; i < L*L*L; i++) {
      unsigned long int v = gsi( get_index(ts+t, 0, 0, 0, T, L) + i );
      for(unsigned int mu = 0; mu < 4; mu++) {
	  // here comes the gamma_sink and gamma_5 at sink
	fv_eq_gamma_ti_fv(psi0, gamma_sink, &phi[ mu ][ v ]);
	fv_eq_gamma_ti_fv(psi1, 5, psi0);
	
	// spinor index p[mu] takes care of gamma_source at source
	//co_eq_fv_dag_ti_fv(&tmp, psi1, &chi[ p ][ v ]);
	co_eq_fv_dag_ti_fv(&tmp, &chi[ p[mu] ][ v ], psi1);
	
	  // sign at source (also gamma_5) and possible factor of i
	if( !isimag ) {
	  re += gamma_sign[5][6*mu]*gamma_sign[gamma_source][6*mu]*tmp.re;
	  im += gamma_sign[5][6*mu]*gamma_sign[gamma_source][6*mu]*tmp.im;
	}
	else {
	  re += gamma_sign[5][6*mu]*gamma_sign[gamma_source][6*mu]*tmp.im;
	  im -= gamma_sign[5][6*mu]*gamma_sign[gamma_source][6*mu]*tmp.re;
	}
      }
    }
    
    CF[tt].re += re;
    CF[tt].im += im;
      tt++;
  }
  return;
}

void contract_twopoint(complex * CF, const int gamma_source, const int gamma_sink,
		       double * chi[], double * phi[],
		       const unsigned int T, const unsigned int L, 
		       const unsigned int ts, 
		       int * const perm[], int * const sign[]) {
  double re, im;
  int tt=0;
  int isimag = (perm[gamma_source][0])%2;
  int p[4];
  
  for(unsigned int mu = 0; mu < 4; mu++) {
    p[mu] = perm[gamma_source][6*mu]/6;
  }
  complex tmp;
  double psi0[24], psi1[24];    
  for(unsigned int t = 0; t < T; t++) {
    tmp.re = 0.;
    tmp.im = 0.;
    re = 0.;
    im = 0.;
    for(unsigned int i = 0; i < L*L*L; i++) {
      unsigned long int v = gsi( get_index(ts+t, 0, 0, 0, T, L) + i );
      for(unsigned int mu = 0; mu < 4; mu++) {
	
	// here comes the gamma_sink and gamma_5 at sink
	fv_eq_gamma_ti_fv(psi0, &phi[ mu ][ v ], perm[gamma_sink], sign[gamma_sink]);
	fv_eq_gamma_ti_fv(psi1, psi0, perm[5], sign[5]);
	
	// spinor index p[mu] takes care of gamma_source at source
	//co_eq_fv_dag_ti_fv(&tmp, psi1, &chi[ p ][ v ]);
	co_eq_fv_dag_ti_fv(&tmp, &chi[ p[mu] ][ v ], psi1);
	
	// sign at source (also gamma_5) and possible factor of i
	if( !isimag ) {
	  re += sign[5][6*mu]*sign[gamma_source][6*mu]*tmp.re;
	  im += sign[5][6*mu]*sign[gamma_source][6*mu]*tmp.im;
	}
	else {
	  re += sign[5][6*mu]*sign[gamma_source][6*mu]*tmp.im;
	  im -= sign[5][6*mu]*sign[gamma_source][6*mu]*tmp.re;
	}
      }
    }
    CF[tt].re += re;
    CF[tt].im += im;
    tt++;
  }
  return;
}


// next version is for general gamma_5
// should compute every gamma matrix combination at source and sink 
// using gamma_5 hermiticity
// this version can also cope with 12 components instead of only four

void contract_twopoint(complex * CF, int const per_source[], int const sig_source[],
		       int per_sink[], int sig_sink[],
		       int per_5[], int sig_5[],
		       double * chi[], double * phi[],
		       const unsigned int T, const unsigned int L, 
		       const unsigned int ts, const unsigned int n_c) {
  int isimag = (per_source[0])%2;
  int p[4];
  int psource[24], psink[24], ssource[24], ssink[24];
#ifndef OMP
  complex tmp;
  double psi0[24];
#endif

  // compute gamma_source gamma_5
  gamma_eq_gamma_ti_gamma(psource, ssource, per_source, per_5, sig_source, sig_5);
  // and gamma_5 gamma_sink
  gamma_eq_gamma_ti_gamma(psink, ssink, per_5, per_sink, sig_5, sig_sink);

  isimag = psource[0]%2;
  for(unsigned int mu = 0; mu < 4; mu++) {
    p[mu] = psource[6*mu]/6;
  }
  if(n_c == 1) {
#ifdef OMP
# pragma omp parallel for
# endif
    for(unsigned int t = 0; t < T; t++) {
#ifdef OMP
      complex tmp;
      double psi0[24];
#endif
      for(unsigned int i = 0; i < L*L*L; i++) {
	unsigned long int v = gsi( get_index(ts+t, 0, 0, 0, T, L) + i );
	for(unsigned int mu = 0; mu < 4; mu++) {
	  
	  // here comes the gamma_sink * gamma_5 at sink
	  fv_eq_gamma_ti_fv(psi0, &phi[ mu ][ v ], psink, ssink);
	  
	  // spinor index p[mu] takes care of gamma_source at source
	  co_eq_fv_dag_ti_fv(&tmp, &chi[ p[mu] ][ v ], psi0);
	  
	  // sign at source and possible factor of i
	  if( !isimag ) {
	    CF[t].re += ssource[6*mu]*tmp.re;
	    CF[t].im += ssource[6*mu]*tmp.im;
	  }
	  else {
	    CF[t].re += ssource[6*mu]*tmp.im;
	    CF[t].im -= ssource[6*mu]*tmp.re;
	  }
	}
      }
    }
  }
  else {
#ifdef OMP
# pragma omp parallel for
# endif
    for(unsigned int t = 0; t < T; t++) {
#ifdef OMP
      complex tmp;
      double psi0[24];
#endif
      for(unsigned int i = 0; i < L*L*L; i++) {
	unsigned long int v = gsi( get_index(ts+t, 0, 0, 0, T, L) + i );
	// spin indices
	for(unsigned int mu = 0; mu < 4; mu++) {
	  // colour indices
	  for(unsigned int c = 0; c < n_c; c++) {
	    
	    // here comes the gamma_sink * gamma_5 at sink
	    fv_eq_gamma_ti_fv(psi0, &phi[ mu*n_c + c ][ v ], psink, ssink);
	    
	    // spinor index p[mu] takes care of gamma_source at source
	    co_eq_fv_dag_ti_fv(&tmp, &chi[ p[mu]*n_c + c ][ v ], psi0);
	    
	    // sign at source and possible factor of i
	    if( !isimag ) {
	      CF[t].re += ssource[6*mu]*tmp.re;
	      CF[t].im += ssource[6*mu]*tmp.im;
	    }
	    else {
	      CF[t].re += ssource[6*mu]*tmp.im;
	      CF[t].im -= ssource[6*mu]*tmp.re;
	    }
	  }
	}
      }
    }
  }
  return;
}
