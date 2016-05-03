#ifndef _CONTRACT_TWOPOINT_HH
#define _CONTRACT_TWOPOINT_HH

void contract_twopoint(complex * CF, const int gamma_source, const int gamma_sink,
		       double * chi[], double * phi[],
		       const unsigned int T, const unsigned int L, const unsigned int ts);

void contract_twopoint(complex * CF, const int gamma_source, const int gamma_sink,
		       double * chi[], double * phi[],
		       const unsigned int T, const unsigned int L, const int ts, 
		       int * const perm[], int * const sign[]);

	/** added by S.D.
	
			computes <\chi^\dagger (0) \phi (t)>, where
			\phi = \bar{u} Gamma_{sink}   d
			\chi = \bar{u} Gamma_{source} d
			
			this routine can cope with 12 sources (n_c=3)
			or 4 sources, respectively (n_c=1)
			
			parameter				type								description											example
			------------------------------------------------------------------------------------------
	  	CF 							complex[T*n_c*4])		resulting two-point function		-
	  	per_source[]		int[24]							gamma permutation								&(gamma_permutation[i][0]), i=0,...,15
	  	sig_source[]		int[24]							gamma sign											&(gamma_sign[i][0]), i=0,...,15
	  	per_sink[]			int[24]							gamma permutation								&(gamma_permutation[i][0]), i=0,...,15
	  	sig_sink[]			int[24]							gamma sign											&(gamma_sign[i][0]), i=0,...,15
	  	per_5[]					int[24]							gamma5 permutation							&(gamma_permutation[5][0])
	  	sig_5[]					int[24]							gamma5 sign											&(gamma_sign[5][0])
			chi[] 					double *[n_c*4]	 		spinor field (propagator)				-
			phi[] 					double *[n_c*4]	 		spinor field (propagator)				-
			T, L: 					int									lattice extent									-
			ts							int									time slice of source						-
			n_c							int									number of colors 								-
	*/
void contract_twopoint(complex * CF, int const per_source[], int const sig_source[],
		       int per_sink[], int sig_sink[],
		       int per_5[], int sig_5[],
		       double * chi[], double * phi[],
		       const unsigned int T, const unsigned int L, const unsigned int ts,
		       const unsigned int n_c = 1);

#endif
