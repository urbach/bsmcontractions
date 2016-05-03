/* $Id: propagator_io.h,v 1.2 2007/11/24 14:37:24 urbach Exp $ */

#ifndef _PROPAGATOR_IO_HH
#define _PROPAGATOR_IO_HH

int write_source(double * const s, char * filename, 
		 const int append, const int prec,
		 const int T, const int LX, const int LY, const int LZ);

int write_propagator(double * const s, char * filename, 
		     const int append, const int prec,
		     const int T, const int LX, const int LY, const int LZ);

int write_double_propagator(double * const s, double * const r, 
			    char * filename, const int append, const int prec,
			    const int T, const int LX, const int LY, const int LZ);

int write_propagator_format(char * filename, const int prec, const int no_flavours,
			    const int T, const int LX, const int LY, const int LZ);
int write_propagator_type(const int type, char * filename);
int write_source_type(const int type, char * filename);
int get_propagator_type(char * filename);
int write_lime_spinor(double * const s, char * filename, const int append, const int prec,
		      const unsigned int T, const unsigned int LX, const unsigned int LY, const unsigned int LZ);
int read_lime_spinor(double * const s, char * filename, const int position, const int ts,
		     const unsigned int T, const unsigned int LX, const unsigned int LY, const unsigned int LZ);

int read_lime_spinor(float * const s, char * filename, const int position, const int ts,
		     const unsigned int T, const unsigned int LX, const unsigned int LY, const unsigned int LZ);


int read_lime_spinor_timeslice(double *const s, char *filename, const int position, const int ts, const unsigned int T, const unsigned int LX, const unsigned int LY, const unsigned int LZ);

#endif
