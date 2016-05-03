/* $Id: io.c,v 1.2 2006/04/18 15:29:27 urbach Exp $ */

/****************************************************
 * IO routines:
 *
 * read_lime_gauge_field_doubleprec
 *
 * read_lime_gauge_field_singleprec
 *
 * Autor: 
 *        Carsten Urbach <urbach@ifh.de>
 *
 ****************************************************/

/*
 * Note:
 * Required version of lime: >= 1.2.3
 * n_uint64_t is a lime defined type!!
 *
 */

#define _FILE_OFFSET_BITS 64

#include"lime.h" 
//#include<complex>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <errno.h>
#include "geometry.hh"
#include "io.hh"
#include "io_utils.hh"
#include "linear_algebra.hh"

//using std::complex;

#define MAXBUF 1048576

void read_gwc(double *v, const int L, const int T, const char * filename, int header) {
  FILE * ifs;
  double tmp[24], tmp2[24];
  int words_bigendian;
  double beta, m0, s;
  int l, t;

  words_bigendian = big_endian();
  ifs = fopen(filename, "r");
  if(ifs != NULL ) {
    if(header == 1) {
      fscanf(ifs,"%lf %lf %lf %d %d\n", &beta, &m0, &s, &l, &t);
      if((l != L)||(t != T)) {
	fprintf(stderr, "probably wrong lattice size for %s\n", filename);
	exit(-1);
      }
    }
    for(int x = 0; x < L; x++) {
      for(int y = 0; y < L; y++) {
	for(int z = 0; z < L; z++) {
	  for(int t = 0; t < T; t++) {
	    int ix = (t*L*L*L + x*L*L + y*L + z)*12;
	    // complex<double> phase(cos(t*3.141593/T), sin(t*3.141593/T));
	    fread(tmp, 24*sizeof(double), 1, ifs);
	    
	    if(!words_bigendian) {
	      byte_swap_assign(tmp2, tmp, 24);
	    }
	    
	    for(int i = 0; i < 12; i++) {
	      // v(ix+i) = phase*complex<double>(tmp[2*i], tmp[2*i+1]);
	      // v(ix+i) = complex<double>(tmp[2*i], tmp[2*i+1]);
	      
	      int index = (ix+i)*2;
	      if(words_bigendian) {
		v[index  ] = tmp[2*i  ];
		v[index+1] = tmp[2*i+1];
	      }
	      else {
		v[index  ] = tmp2[2*i  ];
		v[index+1] = tmp2[2*i+1];
	      }
	    }
	  }
	}
      }
    }
  }
  fclose(ifs);
  return;
}

void read_cmi(double *v, const int L, const int T, const char * filename, double kappa) {
  FILE * ifs;
  float tmp[24];
  ifs = fopen(filename, "r");

  double _2_kappa = 2.0 * kappa;

  for(int x = 0; x < L; x++) {
    for(int y = 0; y < L; y++) {
      for(int z = 0; z < L; z++) {
	for(int t = 0; t < T; t++) {
	  int ix = (t*L*L*L + x*L*L + y*L + z)*12;
// 	  complex<double> phase(cos(t*3.141593/T), sin(t*3.141593/T));
	  fread(tmp, 24*sizeof(float), 1, ifs);

          // !!!!!!!!!!
          // !!!!!!!!!!
          // !!!!!!!!!!
          // !!!!!!!!!!
          // !!!!!!!!!!

          // fprintf(stderr, "before byte swap: %f\n", tmp[0]);
          // byte_swap(tmp, 24*sizeof(float) / 4);
          // fprintf(stderr, "after byte swap: %f\n", tmp[0]);

          // !!!!!!!!!!
          // !!!!!!!!!!
          // !!!!!!!!!!
          // !!!!!!!!!!
          // !!!!!!!!!!

	  for(int i = 0; i < 12; i++) {
	    // v(ix+i) = phase*complex<double>(tmp[2*i], tmp[2*i+1]);
	    // v(ix+i) = complex<double>(tmp[2*i], tmp[2*i+1]);

	    int index = (ix+i)*2;
	    v[index  ] = _2_kappa * tmp[2*i  ];
	    v[index+1] = _2_kappa * tmp[2*i+1];
	  }
	}
      }
    }
  }
  fclose(ifs);
  return;
}

void read_ukqcd(double *v, const int L, const int T, const char * filename, double kappa) {
  FILE * ifs;
  float tmp[24];
  ifs = fopen(filename, "r");

  double _2_kappa = 2.0 * kappa;

  for(int t = 0; t < T; t++) {
    for(int z = 0; z < L; z++) {
      for(int y = 0; y < L; y++) {
	for(int x = 0; x < L; x++) {

	  int ix = (t*L*L*L + x*L*L + y*L + z)*12;
// 	  complex<double> phase(cos(t*3.141593/T), sin(t*3.141593/T));
	  fread(tmp, 24*sizeof(float), 1, ifs);

          // fprintf(stderr, "before byte swap: %f\n", tmp[0]);
          // byte_swap(tmp, 24*sizeof(float) / 4);
          // fprintf(stderr, "after byte swap: %f\n", tmp[0]);

	  for(int i = 0; i < 4; i++) {
	    for(int j = 0; j < 3; j++) {
	      // v(ix+i) = phase*complex<double>(tmp[2*i], tmp[2*i+1]);
	      // v(ix+i) = complex<double>(tmp[2*i], tmp[2*i+1]);

	      int index = 2*(ix + 3*i + j);
// 	      if(i == 0) {
// 		  index = 2*(ix + 9 + j);
// 	      }
// 	      else {
// 		  index = 2*(ix + (i-1)*3 + j);
// 	      }
	      v[index  ] = _2_kappa * tmp[2*(i*3+j)  ];
	      v[index+1] = _2_kappa * tmp[2*(i*3+j)+1];
	    }
	  }
	}
      }
    }
  }
  fclose(ifs);
  return;
}


// int read_lime_gauge_field_doubleprec(matrix< complex<double> > &config, const char * filename, const int T, const int LX, const int LY, const int LZ) {
int read_lime_gauge_field_doubleprec(double *config, const char * filename, 
				     const int T, const int LX, const int LY, const int LZ) {
  FILE * ifs;
  int t, x, y, z, status;
  n_uint64_t bytes;
  char * header_type;
  LimeReader * limereader;
  double tmp[72], tmp2[72];
  int words_bigendian;

  words_bigendian = big_endian();
  ifs = fopen(filename, "r");
  if(ifs == (FILE *)NULL) {
    fprintf(stderr, "Could not open file %s\n Aborting...\n", filename);
    exit(500);
  }
  limereader = limeCreateReader( ifs );
  if( limereader == (LimeReader *)NULL ) {
    fprintf(stderr, "Unable to open LimeReader\n");
    exit(500);
  }
  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
    if(status != LIME_SUCCESS ) {
      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
      status = LIME_EOF;
      break;
    }
    header_type = limeReaderType(limereader);
    if(strcmp("ildg-binary-data",header_type) == 0) break;
  }
  if(status == LIME_EOF) {
    fprintf(stderr, "no ildg-binary-data record found in file %s\n",filename);
    limeDestroyReader(limereader);
    fclose(ifs);
    exit(-2);
  }
  bytes = limeReaderBytes(limereader);
  if(bytes != (n_uint64_t)LX*LY*LZ*T*72*(n_uint64_t)sizeof(double)) {
    if(bytes != (n_uint64_t)LX*LY*LZ*T*72*(n_uint64_t)sizeof(double)/2) {
      fprintf(stderr, "Probably wrong lattice size or precision (bytes=%lu) in file %s expected %lu\n", 
	      (unsigned long)((n_uint64_t)bytes), filename, (unsigned long)((n_uint64_t)LX*LY*LZ*T*72*(n_uint64_t)sizeof(double)));
      fprintf(stderr, "Aborting...!\n");
      fflush( stdout );
      exit(501);
    }
    else {
      fclose(ifs);
      fprintf(stderr, "single precision read!\n");
      return( read_lime_gauge_field_singleprec(config, filename, T, LX, LY, LZ) );
    }
  }

  bytes = (n_uint64_t)72*sizeof(double);

  for(t = 0; t < T; t++) {
    for(z = 0; z < LZ; z++) {
      for(y = 0; y < LY; y++) {
	for(x = 0; x < LX; x++) {
	  n_uint64_t p = (((t*LX+x)*LY+y)*LZ+z)*(n_uint64_t)12;
	  if(!words_bigendian) {
	    status = limeReaderReadData(tmp, &bytes, limereader);
	    byte_swap_assign(tmp2, tmp, 72);
	  }
	  else {
	    status = limeReaderReadData(tmp2, &bytes, limereader);
	  }
	  n_uint64_t k =0;
	  //ILDG has mu-order: x,y,z,t
	  for(int mu = 1; mu < 4; mu++) {
	    for(int i = 0; i < 3; i++) {
	      for(int j = 0; j < 3; j++) {
		// config (p+mu*3+i, j) = complex<double> (tmp2[2*k], tmp2[2*k+1]);

		n_uint64_t index = ((p+mu*3+i) * 3 + j) * 2;
		config[index  ] = tmp2[2*k];
		config[index+1] = tmp2[2*k+1];

		k++;
	      }
	    }
	  }
 	  for(int i = 0; i < 3; i++) {
 	    for(int j = 0; j < 3; j++) {
 	      // config (p+i, j) = complex<double> (tmp2[2*k], tmp2[2*k+1]);

	      n_uint64_t index = ((p+i) * 3 + j) * 2;
	      config[index  ] = tmp2[2*k];
	      config[index+1] = tmp2[2*k+1];

 	      k++; 	    
	    }
 	  }
	  if(status < 0 && status != LIME_EOR) {
	    fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n", 
		    status, filename);
	    exit(500);
	  }
	}
      }
    }
  }
  limeDestroyReader(limereader);
  fclose(ifs);
  return(0);
}



// **********
// **********
// **********
// **********
// **********

int read_lime_gauge_field_doubleprec_timeslices(double *config, const char *filename, const int T, const int LX, const int LY, const int LZ, const int slice_i, const int slice_f) {
  FILE * ifs;
  int t, x, y, z, status;
  n_uint64_t bytes;
  char * header_type;
  LimeReader * limereader;
  double tmp[72], tmp2[72];
  int words_bigendian;

  words_bigendian = big_endian();
  ifs = fopen(filename, "r");
  if(ifs == (FILE *)NULL) {
    fprintf(stderr, "Could not open file %s\n Aborting...\n", filename);
    exit(500);
  }
  limereader = limeCreateReader( ifs );
  if( limereader == (LimeReader *)NULL ) {
    fprintf(stderr, "Unable to open LimeReader\n");
    exit(500);
  }
  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
    if(status != LIME_SUCCESS ) {
      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
      status = LIME_EOF;
      break;
    }
    header_type = limeReaderType(limereader);
    if(strcmp("ildg-binary-data",header_type) == 0) break;
  }
  if(status == LIME_EOF) {
    fprintf(stderr, "no ildg-binary-data record found in file %s\n",filename);
    limeDestroyReader(limereader);
    fclose(ifs);
    exit(-2);
  }
  bytes = limeReaderBytes(limereader);
  if(bytes != (n_uint64_t)LX*LY*LZ*T*72*(n_uint64_t)sizeof(double)) {
    if(bytes != (n_uint64_t)LX*LY*LZ*T*72*(n_uint64_t)sizeof(double)/2) {
      fprintf(stderr, "Probably wrong lattice size or precision (bytes=%lu) in file %s expected %lu\n", 
	      (n_uint64_t)bytes, filename, (n_uint64_t)LX*LY*LZ*T*72*(n_uint64_t)sizeof(double));
      fprintf(stderr, "Aborting...!\n");
      fflush( stdout );
      exit(501);
    }
    else {
      fclose(ifs);
      fprintf(stderr, "single precision read!\n");

      fprintf(stderr, "Not implemented!\n");
      exit(EXIT_FAILURE);

      return( read_lime_gauge_field_singleprec(config, filename, T, LX, LY, LZ) );
    }
  }

  bytes = (n_uint64_t)72*sizeof(double);

  for(t = 0; t < T; t++) {
    for(z = 0; z < LZ; z++) {
      for(y = 0; y < LY; y++) {
	for(x = 0; x < LX; x++) {
	  n_uint64_t p = ((((t-slice_i)*LX+x)*LY+y)*LZ+z)*(n_uint64_t)12;
	  if(!words_bigendian) {
	    status = limeReaderReadData(tmp, &bytes, limereader);
	    byte_swap_assign(tmp2, tmp, 72);
	  }
	  else {
	    status = limeReaderReadData(tmp2, &bytes, limereader);
	  }

	  if(t < slice_i || t > slice_f)
	    continue;

	  n_uint64_t k =0;
	  //ILDG has mu-order: x,y,z,t
	  for(int mu = 1; mu < 4; mu++) {
	    for(int i = 0; i < 3; i++) {
	      for(int j = 0; j < 3; j++) {
		// config (p+mu*3+i, j) = complex<double> (tmp2[2*k], tmp2[2*k+1]);

		n_uint64_t index = ((p+mu*3+i) * 3 + j) * 2;
		config[index  ] = tmp2[2*k];
		config[index+1] = tmp2[2*k+1];

		k++;
	      }
	    }
	  }
 	  for(int i = 0; i < 3; i++) {
 	    for(int j = 0; j < 3; j++) {
 	      // config (p+i, j) = complex<double> (tmp2[2*k], tmp2[2*k+1]);

	      n_uint64_t index = ((p+i) * 3 + j) * 2;
	      config[index  ] = tmp2[2*k];
	      config[index+1] = tmp2[2*k+1];

 	      k++; 	    
	    }
 	  }
	  if(status < 0 && status != LIME_EOR) {
	    fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n", 
		    status, filename);
	    exit(500);
	  }
	}
      }
    }
  }
  limeDestroyReader(limereader);
  fclose(ifs);
  return(0);
}

// **********
// **********
// **********
// **********
// **********

 
int read_lime_gauge_field_singleprec(double *config, const char * filename, 
				     const int T, const int LX, const int LY, const int LZ) {
  FILE * ifs;
  int t, x, y, z, status;
  n_uint64_t bytes;
  char * header_type;
  LimeReader * limereader;
  float tmp[72], tmp2[72];
  int words_bigendian;

  words_bigendian = big_endian();
  ifs = fopen(filename, "r");
  if(ifs == (FILE *)NULL) {
    fprintf(stderr, "Could not open file %s\n Aborting...\n", filename);
    exit(500);
  }
  limereader = limeCreateReader( ifs );
  if( limereader == (LimeReader *)NULL ) {
    fprintf(stderr, "Unable to open LimeReader\n");
    exit(500);
  }
  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
    if(status != LIME_SUCCESS ) {
      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
      status = LIME_EOF;
      break;
    }
    header_type = limeReaderType(limereader);
    if( strcmp("ildg-binary-data",header_type) == 0) break;
  }
  if(status == LIME_EOF) {
    fprintf(stderr, "no ildg-binary-data record found in file %s\n",filename);
    limeDestroyReader(limereader);
    fclose(ifs);
    exit(-2);
  }
  bytes = limeReaderBytes(limereader);
  if(bytes != (n_uint64_t)LX*LY*LZ*T*72*(n_uint64_t)sizeof(float)) {
    fprintf(stderr, "Probably wrong lattice size or precision (bytes=%d) in file %s\n", (int)bytes, filename);
    fprintf(stderr, "Aborting...!\n");
    fflush( stdout );
    exit(501);
  }

  bytes = (n_uint64_t)72*sizeof(float);
  for(t = 0; t < T; t++){
    for(z = 0; z < LZ; z++){
      for(y = 0; y < LY; y++){
	for(x = 0; x < LX; x++) {
	  n_uint64_t p = (((t*LX+x)*LY+y)*LZ+z)*(n_uint64_t)12;
	  if(!words_bigendian) {
	    status = limeReaderReadData(tmp, &bytes, limereader);
	    byte_swap_assign_singleprec(tmp2, tmp, 72);
	  }
	  else {
	    status = limeReaderReadData(tmp2, &bytes, limereader);
	  }
	  n_uint64_t k =0;
	  //ILDG has mu-order: x,y,z,t
	  for(int mu = 1; mu < 4; mu++) {
	    for(int i = 0; i < 3; i++) {
	      for(int j = 0; j < 3; j++) {
		//config (p+mu*3+i, j) = complex<double> (tmp2[2*k], tmp2[2*k+1]);
		n_uint64_t index = ((p+mu*3+i) * 3 + j) * 2;
		config[index  ] = (double)tmp2[2*k];
		config[index+1] = (double)tmp2[2*k+1];

		k++;
	      }
 	    }
	  }
 	  for(int i = 0; i < 3; i++) {
	    for(int j = 0; j < 3; j++) {
	      //  config (p+i, j) = complex<double> (tmp2[2*k], tmp2[2*k+1]);
	      n_uint64_t index = ((p+i) * 3 + j) * 2;
	      config[index  ] = tmp2[2*k];
	      config[index+1] = tmp2[2*k+1];

	      k++;
	    }
	  }

	  if(status < 0 && status != LIME_EOR) {
	    fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n", 
		    status, filename);
	    exit(500);
	  }
	} 
      }
    }
  }
  limeDestroyReader(limereader);
  fclose(ifs);
  return(0);
}




// **********



// Read in a gauge field configuration (perhaps some old gwc format ...).

void read_gauge_field(double *gauge_field, char *filename, int T, int L)
{
  FILE *fd;

  if((fd = fopen(filename, "r")) == NULL)
    {
      fprintf(stderr, "(1) Error: int read_gauge_field(...\n");
      fprintf(stderr, "  Could not open file '%s'!\n", filename);
      exit(EXIT_FAILURE);
    }

  double beta;
  int L_, T_;
  fscanf(fd, "%lf %d %d\n", &beta, &L_, &T_);

  if(T != T_ || L != L_)
    {
      fprintf(stderr, "(2) Error: int read_gauge_field(...\n");
      fprintf(stderr, "  Wrong T (%2d, %2d) or L (%2d, %2d)!\n", T, T_, L, L_);
      exit(EXIT_FAILURE);
    }

  int ix, iy, iz, it;

  for(ix = 0; ix < L; ix++)
    {
      for(iy = 0; iy < L; iy++)
	{
	  for(iz = 0; iz < L; iz++)
	    {
	      for(it = 0; it < T; it++)
		{
		    double *links =
		      gauge_field + ggi(get_index(it, ix, iy, iz, T, L), 0);

		    /*
		        // no conversion

			  if(fread(links, 4*18*sizeof(double), 1, fd) != 1)
			      {
			            fprintf(stderr, "(3) Error: int read_gauge_field(...\n");
				          fprintf(stderr, "  Could not read link!\n");
					        exit(EXIT_FAILURE);
						    }
		    */

		    // /*
		    // conversion: big <--> little endian

		    double links_[4*18];

		    if(fread(links_, 4*18*sizeof(double), 1, fd) != 1)
		      {
			fprintf(stderr, "(3) Error: int read_gauge_field(...\n");
			fprintf(stderr, "  Could not read link!\n");
			exit(EXIT_FAILURE);
		      }

		    byte_swap_assign(links, links_, 4*18*sizeof(double)/8);
		    // */
		}
	    }
	}
    }

  fclose(fd);
}



// **********


int read_scalar_field(char * filename, complex ** const sf, 
		      const unsigned int T, const unsigned int L) {

  const int scalar_precision_read_flag = 64;
  FILE *ptr;

  const unsigned int VOL = T*L*L*L;
  const unsigned int count = 4*VOL;
  int scalarreadsize = ( scalar_precision_read_flag==64 ? 2*sizeof(double) : 2*sizeof(float) );

  ptr = fopen(filename,"rb");  // r for read, b for binary
  
  // read into buffer
  void *buffer;
  if((buffer = malloc(count*scalarreadsize)) == NULL) {
    printf ("malloc errno : %d\n", errno);
    errno = 0;
    return(2);
  }

  if( count > fread(buffer, scalarreadsize, count, ptr) )
    return(-1);

  // copy to sf
  int k,l;
  for( unsigned int s = 0; s < 2; s++ ) {
    for( unsigned int i = 0; i < VOL; i++ ) {
      if(s == 0) {
	k=0; l=3;
      }
      else {
	k=2; l=1;
      }
      if ( scalar_precision_read_flag == 64 ) {
        sf[s][i].re = ((double*)buffer)[4*i+k];
        sf[s][i].im = ((double*)buffer)[4*i+l];
      }
      else {
        sf[s][i].re = ((float*)buffer)[4*i+k];
	sf[s][i].im = ((float*)buffer)[4*i+l];
      }
    }
  }

  return(0);
}
