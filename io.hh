/* $Id: io.h,v 1.2 2006/04/18 15:29:27 urbach Exp $ */
#ifndef _IO_H
#define _IO_H

#include "linear_algebra.hh"

void read_gwc(double *v, const int L, const int T, const char * filename, int header=0);

void read_cmi(double *v, const int L, const int T, const char * filename, double kappa);
void read_ukqcd(double *v, const int L, const int T, const char * filename, double kappa);

int read_lime_gauge_field_doubleprec(double *config, const char * filename, 
				     const int T, const int LX, const int LY, const int LZ);

int read_lime_gauge_field_doubleprec_timeslices(double *config, const char *filename, const int T, const int LX, const int LY, const int LZ, const int slice_i, const int slice_f);

int read_lime_gauge_field_singleprec(double *config, const char * filename, 
				     const int T, const int LX, const int LY, const int LZ);


// Read in a gauge field configuration (perhaps some old gwc format ...).
void read_gauge_field(double *gauge_field, char *filename, int T, int L);

// reading a scalar field
int read_scalar_field(char * filename, complex ** const sf, const unsigned int T, const unsigned int L);


#endif
