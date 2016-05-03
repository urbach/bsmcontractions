/* $Id: propagator_io.c,v 1.3 2007/11/24 14:37:24 urbach Exp $ */

#define _FILE_OFFSET_BITS 64

#include"lime.h" 
#include<iostream>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<time.h>
#include<sys/time.h> 
#include<sys/types.h>
#include<cmath>
#include"lime.h" 
#include"geometry.hh"
#include"io_utils.hh"
#include"dml.hh"
#include"propagator_io.hh"

using namespace std;
/* write a one flavour propagator to file */

int write_source(double * const s, char * filename, 
		 const int append, const int prec,
		 const int T, const int LX, const int LY, const int LZ) {

  int err = 0;

  write_source_type(0, filename);
  write_propagator_format(filename, prec, 1, T, LX, LY, LZ);
  err = write_lime_spinor(s, filename, 1, prec, T, LX, LY, LZ);
  return(err);
}

int write_propagator(double * const s, char * filename, 
		     const int append, const int prec,
		     const int T, const int LX, const int LY, const int LZ) {
  int err = 0;

//   write_propagator_type(0, filename);
  write_propagator_format(filename, prec, 1, T, LX, LY, LZ);
  err = write_lime_spinor(s, filename, append, prec, T, LX, LY, LZ);
  return(err);
}

/* write two flavour operator to file */

int write_double_propagator(double * const s, double * const r, 
			    char * filename, const int append, const int prec,
			    const int T, const int LX, const int LY, const int LZ) {
  int err = 0;

  write_propagator_format(filename, prec, 2, T, LX, LY, LZ);
  err = write_lime_spinor(s, filename, append, prec, T, LX, LY, LZ);
  err += write_lime_spinor(r, filename, append, prec, T, LX, LY, LZ);
  return(err);
}

DML_Checksum write_binary_spinor_data(double * const s, LimeWriter * limewriter,
				      const int prec,
				      const unsigned int T, const unsigned int LX, 
				      const unsigned int LY, const unsigned int LZ) {
  
  int status=0;
  double tmp[24];
  float tmp2[24];
  n_uint64_t bytes;
  DML_Checksum ans;
  DML_SiteRank rank;
  int words_bigendian;
  words_bigendian = big_endian();

  DML_checksum_init(&ans);
  rank = (DML_SiteRank) 0;
	  
  if(prec == 32) bytes = 24*sizeof(float);
  else bytes = 24*sizeof(double);
  for(unsigned int t = 0; t < T; t++) {
    for(unsigned int z = 0; z < LZ; z++) {
      for(unsigned int y = 0; y < LY; y++) {
	for(unsigned int x = 0; x < LX; x++) {
	  n_uint64_t ix = (t*LX*LY*LZ + x*LY*LZ + y*LZ + z)*(n_uint64_t)12;
	  rank = (DML_SiteRank) (((t*LZ + z)*LY + y)*LX + x);
	  if(!words_bigendian) {
	    if(prec == 32) {
	      byte_swap_assign_double2single((float*)tmp2, &s[2*ix], 24); 
	      DML_checksum_accum(&ans,rank,(char *) tmp2, bytes);
	      status = limeWriteRecordData((void*)tmp2, &bytes, limewriter);
	    }
	    else {
	      byte_swap_assign(tmp, &s[2*ix], 24);
	      DML_checksum_accum(&ans,rank,(char *) tmp, bytes);
	      status = limeWriteRecordData((void*)tmp, &bytes, limewriter);
	    }
	  }
	  else {
	    if(prec == 32) {
	      double2single((float*)tmp2, &s[2*ix], 24); 
	      DML_checksum_accum(&ans,rank, (char *) tmp2, bytes);
	      status = limeWriteRecordData((void*)tmp2, &bytes, limewriter);
	    }
	    else {
	      status = limeWriteRecordData((void*) &s[2*ix], &bytes, limewriter);
	      DML_checksum_accum(&ans,rank,(char *) &s[2*ix], bytes);
	    }
	  }
	}
      }
    }
  }
  printf("The final checksum is %#lx %#lx\n", (unsigned long)ans.suma, (unsigned long)ans.sumb);  
  return(ans);
}

/* if -1 < ts < T read only timeslice ts */
/* the other entries in s are untouched  */

int read_binary_spinor_data(double * const s, LimeReader * limereader, 
			    const double prec, const int ts, DML_Checksum &ans,
			    unsigned long int get_index(const int, const int, const int, const int, const int, const int),
			    const unsigned int T, const unsigned int LX, 
			    const unsigned int LY, const unsigned int LZ) {
  int status=0;
  n_uint64_t bytes;
  double tmp[24];
  DML_SiteRank rank;
  float tmp2[24];
  int words_bigendian;
  words_bigendian = big_endian();

  DML_checksum_init(&ans);
  rank = (DML_SiteRank) 0;
  
  if(prec == 32) bytes = 24*sizeof(float);
  else bytes = 24*sizeof(double);
  for(unsigned int t = 0; t < T; t++){
    if(ts > -1 && (unsigned int)abs(ts) < T) {
      t = ts;
      limeReaderSeek(limereader,(n_uint64_t) 
		     (t*LZ*LY*LX)*bytes,
		     SEEK_SET);
    }
    for(unsigned int z = 0; z < LZ; z++){
      for(unsigned int y = 0; y < LY; y++){
	for(unsigned int x = 0; x < LX; x++){
	  n_uint64_t ix = get_index(t, x, y, z, T, LX)*(n_uint64_t)12;
	  rank = (DML_SiteRank) (((t*LZ + z)*LY + y)*LX + x);
	  if(prec == 32) {
	    status = limeReaderReadData(tmp2, &bytes, limereader);
	    DML_checksum_accum(&ans,rank,(char *) tmp2, bytes);	    
	  }
	  else {
	    status = limeReaderReadData(tmp, &bytes, limereader);
	    DML_checksum_accum(&ans,rank,(char *) tmp, bytes);
	  }
	  if(!words_bigendian) {
	    if(prec == 32) {
	      byte_swap_assign_single2double(&s[2*ix], (float*)tmp2, 24);
	    }
	    else {
	      byte_swap_assign(&s[2*ix], tmp, 24);
	    }
	  }
	  else {
	    if(prec == 32) {
	      single2double(&s[2*ix], (float*)tmp2, 24);
	    }
	    else memcpy(&s[2*ix], tmp, bytes);
	  }
	  if(status < 0 && status != LIME_EOR) {
	    return(-1);
	  }
	}
      }
    }
    if(ts > -1 && ts < (int)T) {
      t = T;
    }
  }
  printf("The final checksum is %#lx %#lx\n", (unsigned long)ans.suma, (unsigned long)ans.sumb);
  return(0);
}


int read_binary_spinor_data(float * const s, LimeReader * limereader, 
			    const double prec, const int ts, DML_Checksum &ans,
			    unsigned long int get_index(const int, const int, const int, const int, const int, const int),
			    const unsigned int T, const unsigned int LX, 
			    const unsigned int LY, const unsigned int LZ) {
  int status=0;
  n_uint64_t bytes;
  double tmp[24];
  DML_SiteRank rank;
  float tmp2[24];
  int words_bigendian;
  words_bigendian = big_endian();

  DML_checksum_init(&ans);
  rank = (DML_SiteRank) 0;
  
  if(prec == 32) bytes = 24*sizeof(float);
  else return(-1);
  for(unsigned int t = 0; t < T; t++){
    if(ts > -1 && (unsigned int)abs(ts) < T) {
      t = ts;
      limeReaderSeek(limereader,(n_uint64_t) 
		     (t*LZ*LY*LX)*bytes,
		     SEEK_SET);
    }
    for(unsigned int z = 0; z < LZ; z++){
      for(unsigned int y = 0; y < LY; y++){
	for(unsigned int x = 0; x < LX; x++){
	  n_uint64_t ix = get_index(t, x, y, z, T, LX)*(n_uint64_t)12;
	  rank = (DML_SiteRank) (((t*LZ + z)*LY + y)*LX + x);
	  status = limeReaderReadData(tmp2, &bytes, limereader);
	  DML_checksum_accum(&ans,rank,(char *) tmp2, bytes);	    
	  if(!words_bigendian) {
	    byte_swap_assign_singleprec(&s[2*ix], tmp2, 24);
	  }
	  else {
	    memcpy(&s[2*ix], tmp2, bytes);
	  }
	  if(status < 0 && status != LIME_EOR) {
	    return(-1);
	  }
	}
      }
    }
    if(ts > -1 && ts < (int)T) {
      t = T;
    }
  }
  printf("The final checksum is %#lx %#lx\n", (unsigned long)ans.suma, (unsigned long)ans.sumb);
  return(0);
}



// **********
// **********
// **********
// **********
// **********

int read_binary_spinor_data_timeslice(double * const s, LimeReader * limereader, 
			    const double prec, const int ts, DML_Checksum &ans,
			    unsigned long int get_index(const int, const int, const int, const int, const int, const int),
			    const unsigned int T, const unsigned int LX, 
			    const unsigned int LY, const unsigned int LZ) {
  int status=0;
  n_uint64_t bytes;
  double tmp[24];
  DML_SiteRank rank;
  float tmp2[24];
  int words_bigendian;
  words_bigendian = big_endian();

  DML_checksum_init(&ans);
  rank = (DML_SiteRank) 0;
  
  if(prec == 32) bytes = 24*sizeof(float);
  else bytes = 24*sizeof(double);
  for(unsigned int t = 0; t < T; t++){
    if(ts > -1 && abs(ts) < T) {
      t = ts;
      limeReaderSeek(limereader,(n_uint64_t) 
		     (t*LZ*LY*LX)*bytes,
		     SEEK_SET);
    }
    else
      {
	fprintf(stderr, "Error: int read_binary_spinor_data_timeslice(...\n");
	exit(EXIT_SUCCESS);
      }
    for(unsigned int z = 0; z < LZ; z++){
      for(unsigned int y = 0; y < LY; y++){
	for(unsigned int x = 0; x < LX; x++){



	  // n_uint64_t ix = get_index(t, x, y, z, T, LX)*(n_uint64_t)12;
	  n_uint64_t ix = get_index(0, x, y, z, T, LX)*(n_uint64_t)12;



	  rank = (DML_SiteRank) (((t*LZ + z)*LY + y)*LX + x);
	  if(prec == 32) {
	    status = limeReaderReadData(tmp2, &bytes, limereader);
	    DML_checksum_accum(&ans,rank,(char *) tmp2, bytes);	    
	  }
	  else {
	    status = limeReaderReadData(tmp, &bytes, limereader);
	    DML_checksum_accum(&ans,rank,(char *) tmp, bytes);
	  }
	  if(!words_bigendian) {
	    if(prec == 32) {
	      byte_swap_assign_single2double(&s[2*ix], (float*)tmp2, 24);
	    }
	    else {
	      byte_swap_assign(&s[2*ix], tmp, 24);
	    }
	  }
	  else {
	    if(prec == 32) {
	      single2double(&s[2*ix], (float*)tmp2, 24);
	    }
	    else memcpy(&s[2*ix], tmp, bytes);
	  }
	  if(status < 0 && status != LIME_EOR) {
	    return(-1);
	  }
	}
      }
    }
    if(ts > -1 && ts < T) {
      t = T;
    }
  }
  fprintf(stderr, "The final checksum is %#lx %#lx\n", ans.suma, ans.sumb);
  return(0);
}

// **********
// **********
// **********
// **********
// **********



int write_checksum(DML_Checksum * c, char * filename) {

  return(0);
}

int write_propagator_type(const int type, char * filename) {

  FILE * ofs = NULL;
  LimeWriter * limewriter = NULL;
  LimeRecordHeader * limeheader = NULL;
  int status = 0;
  int ME_flag=1, MB_flag=1;
  char message[500];
  n_uint64_t bytes;

  ofs = fopen(filename, "w");
  
  if(ofs == (FILE*)NULL) {
    fprintf(stderr, "Could not open file %s for writing!\n Aboring...\n", filename);
    exit(500);
  }
  limewriter = limeCreateWriter( ofs );
  if(limewriter == (LimeWriter*)NULL) {
    fprintf(stderr, "LIME error in file %s for writing!\n Aborting...\n", filename);
    exit(500);
  }
  
  if(type == 0) {
    sprintf(message,"DiracFermion_Sink");
    bytes = strlen( message );
  }
  else if (type == 1) {
    sprintf(message,"DiracFermion_Source_Sink_Pairs");
    bytes = strlen( message );
  }
  else if (type == 2) {
    sprintf(message,"DiracFermion_ScalarSource_TwelveSink");
    bytes = strlen( message );
  }
  else if (type == 3) {
    sprintf(message,"DiracFermion_ScalarSource_FourSink");
    bytes = strlen( message );
  }
  
  limeheader = limeCreateHeader(MB_flag, ME_flag, "etmc-propagator-type", bytes);
  status = limeWriteRecordHeader( limeheader, limewriter);
  if(status < 0 ) {
    fprintf(stderr, "LIME write header error %d\n", status);
    exit(500);
  }
  limeDestroyHeader( limeheader );
  limeWriteRecordData(message, &bytes, limewriter);
  
  limeDestroyWriter( limewriter );
  fclose(ofs);
  fflush(ofs);
  return(0);
}

int write_source_type(const int type, char * filename) {

  FILE * ofs = NULL;
  LimeWriter * limewriter = NULL;
  LimeRecordHeader * limeheader = NULL;
  int status = 0;
  int ME_flag=1, MB_flag=1;
  char message[500];
  n_uint64_t bytes;

  ofs = fopen(filename, "w");
  
  if(ofs == (FILE*)NULL) {
    fprintf(stderr, "Could not open file %s for writing!\n Aboring...\n", filename);
    exit(500);
  }
  limewriter = limeCreateWriter( ofs );
  if(limewriter == (LimeWriter*)NULL) {
    fprintf(stderr, "LIME error in file %s for writing!\n Aborting...\n", filename);
    exit(500);
  }
  
  sprintf(message,"DiracFermion_Source");
  bytes = strlen( message );
  
  limeheader = limeCreateHeader(MB_flag, ME_flag, "source-type", bytes);
  status = limeWriteRecordHeader( limeheader, limewriter);
  if(status < 0 ) {
    fprintf(stderr, "LIME write header error %d\n", status);
    exit(500);
  }
  limeDestroyHeader( limeheader );
  limeWriteRecordData(message, &bytes, limewriter);
  
  limeDestroyWriter( limewriter );
  fclose(ofs);
  fflush(ofs);
  return(0);
}

int write_propagator_format(char * filename, const int prec, const int no_flavours,
			    const int T, const int LX, const int LY, const int LZ) {
  FILE * ofs = NULL;
  LimeWriter * limewriter = NULL;
  LimeRecordHeader * limeheader = NULL;
  int status = 0;
  int ME_flag=0, MB_flag=1;
  char message[500];
  n_uint64_t bytes;
  /*   char * message; */


  ofs = fopen(filename, "a");
  
  if(ofs == (FILE*)NULL) {
    fprintf(stderr, "Could not open file %s for writing!\n Aborting...\n", filename);
    exit(500);
  }
  limewriter = limeCreateWriter( ofs );
  if(limewriter == (LimeWriter*)NULL) {
    fprintf(stderr, "LIME error in file %s for writing!\n Aborting...\n", filename);
    exit(500);
  }
  
  sprintf(message, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<etmcFormat>\n<field>diracFermion</field>\n<precision>%d</precision>\n<flavours>%d</flavours>\n<lx>%d</lx>\n<ly>%d</ly>\n<lz>%d</lz>\n<lt>%d</lt>\n</etmcFormat>", prec, no_flavours, LX, LY, LZ, T);
  bytes = strlen( message );
  limeheader = limeCreateHeader(MB_flag, ME_flag, "etmc-propagator-format", bytes);
  status = limeWriteRecordHeader( limeheader, limewriter);
  if(status < 0 ) {
    fprintf(stderr, "LIME write header error %d\n", status);
    exit(500);
  }
  limeDestroyHeader( limeheader );
  limeWriteRecordData(message, &bytes, limewriter);
  
  limeDestroyWriter( limewriter );
  fclose(ofs);
  fflush(ofs);
  
  return(0);
}


int write_lime_spinor(double * const s, char * filename, 
		      const int append, const int prec,
		      const unsigned int T, const unsigned int LX, const unsigned int LY, const unsigned int LZ) {

  FILE * ofs = NULL;
  LimeWriter * limewriter = NULL;
  LimeRecordHeader * limeheader = NULL;
  int status = 0;
  int ME_flag=0, MB_flag=0;
  n_uint64_t bytes;
  DML_Checksum checksum;


  if(append) {
    ofs = fopen(filename, "a");
  }
  else {
    ofs = fopen(filename, "w");
  }
  if(ofs == (FILE*)NULL) {
    fprintf(stderr, "Could not open file %s for writing!\n Aborting...\n", filename);
    exit(500);
  }
  limewriter = limeCreateWriter( ofs );
  if(limewriter == (LimeWriter*)NULL) {
    fprintf(stderr, "LIME error in file %s for writing!\n Aborting...\n", filename);
    exit(500);
  }
  
  bytes = LX*LY*LZ*T*(n_uint64_t)24*sizeof(double)*prec/64;
  MB_flag=0; ME_flag=1;
  limeheader = limeCreateHeader(MB_flag, ME_flag, "scidac-binary-data", bytes);
  status = limeWriteRecordHeader( limeheader, limewriter);
  if(status < 0 ) {
    fprintf(stderr, "LIME write header (scidac-binary-data) error %d\n", status);
    exit(500);
  }
  limeDestroyHeader( limeheader );
  
  checksum = write_binary_spinor_data(s, limewriter, prec, T, LX, LY, LZ);
  printf("Final check sum is (%#lx  %#lx)\n", (unsigned long)checksum.suma, (unsigned long)checksum.sumb);
  if(ferror(ofs)) {
    fprintf(stderr, "Warning! Error while writing to file %s \n", filename);
  }
  limeDestroyWriter( limewriter );
  fclose(ofs);
  fflush(ofs);
  return(0);
}

int get_propagator_type(char * filename) {
  FILE * ifs;
  int status=0, ret=-1;
  n_uint64_t bytes;
  char * tmp;
  LimeReader * limereader;
  
  if((ifs = fopen(filename, "r")) == (FILE*)NULL) {
    fprintf(stderr, "Error opening file %s\n", filename);
    return(ret);
  }
  
  limereader = limeCreateReader( ifs );
  if( limereader == (LimeReader *)NULL ) {
    fprintf(stderr, "Unable to open LimeReader\n");
    return(ret);
  }
  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
    if(status != LIME_SUCCESS ) {
      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
      status = LIME_EOF;
      break;
    }
    if(strcmp("etmc-propagator-type", limeReaderType(limereader)) == 0) break;
  }
  if(status == LIME_EOF) {
    fprintf(stderr, "no etmc-propagator-type record found in file %s\n",filename);
    limeDestroyReader(limereader);
    fclose(ifs);
    return(ret);
  }
  tmp = (char*) calloc(500, sizeof(char));
  bytes = limeReaderBytes(limereader);
  status = limeReaderReadData(tmp, &bytes, limereader);
  limeDestroyReader(limereader);
  fclose(ifs);
  if(strcmp("DiracFermion_Sink", tmp) == 0) ret = 0;
  else if(strcmp("DiracFermion_Source_Sink_Pairs", tmp) == 0) ret = 1;
  else if(strcmp("DiracFermion_ScalarSource_TwelveSink", tmp) == 0) ret = 2;
  else if(strcmp("DiracFermion_ScalarSource_FourSink", tmp) == 0) ret = 3;
  free(tmp);
  return(ret);
}

 
int read_lime_spinor(double * const s, char * filename, const int position,
		     const int ts,
		     const unsigned int T, const unsigned int LX, const unsigned int LY, const unsigned int LZ) {
  FILE * ifs;
  int status=0, getpos=-1;
  n_uint64_t bytes;
  char * header_type;
  LimeReader * limereader;
  n_uint64_t prec = 32;
  DML_Checksum checksum;
  
  if((ifs = fopen(filename, "r")) == (FILE*)NULL) {
    fprintf(stderr, "Error opening file %s\n", filename);
    return(-1);
  }

  limereader = limeCreateReader( ifs );
  if( limereader == (LimeReader *)NULL ) {
    fprintf(stderr, "Unable to open LimeReader\n");
    return(-1);
  }
  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
    if(status != LIME_SUCCESS ) {
      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
      status = LIME_EOF;
      break;
    }
    header_type = limeReaderType(limereader);

    if(strcmp("scidac-binary-data",header_type) == 0) getpos++;
    printf("... found record of type %s pos = %d!\n", header_type, getpos);
    if(getpos == position) break;
  }
  if(status == LIME_EOF) {
    fprintf(stderr, "no scidac-binary-data record found in file %s\n",filename);
    limeDestroyReader(limereader);
    fclose(ifs);
    return(-1);
  }
  bytes = limeReaderBytes(limereader);
  if(bytes == LX*LY*LZ*T*(uint64_t)(24*sizeof(double))) prec = 64;
  else if(bytes == LX*LY*LZ*T*(uint64_t)(24*sizeof(float))) prec = 32;
  else {
    fprintf(stderr, "wrong length in eospinor: bytes = %lu, not %lu. Aborting read!\n", 
	    (unsigned long)bytes, (unsigned long)(LX*LY*LZ*T*(uint64_t)(24*sizeof(double))));
    return(-1);
  }
  printf("# %lu Bit precision read\n", (unsigned long)prec);

  status = read_binary_spinor_data(s, limereader, prec, ts, checksum, &get_index, T, LX, LY, LZ);

  if(status < 0) {
    fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n", 
	    status, filename);
    exit(500);
  }

  limeDestroyReader(limereader);
  fclose(ifs);
  return(0);
}

int read_lime_spinor(float * const s, char * filename, const int position,
		     const int ts,
		     const unsigned int T, const unsigned int LX, const unsigned int LY, const unsigned int LZ) {
  FILE * ifs;
  int status=0, getpos=-1;
  n_uint64_t bytes;
  char * header_type;
  LimeReader * limereader;
  n_uint64_t prec = 32;
  DML_Checksum checksum;
  
  if((ifs = fopen(filename, "r")) == (FILE*)NULL) {
    fprintf(stderr, "Error opening file %s\n", filename);
    return(-1);
  }

  limereader = limeCreateReader( ifs );
  if( limereader == (LimeReader *)NULL ) {
    fprintf(stderr, "Unable to open LimeReader\n");
    return(-1);
  }
  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
    if(status != LIME_SUCCESS ) {
      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
      status = LIME_EOF;
      break;
    }
    header_type = limeReaderType(limereader);

    if(strcmp("scidac-binary-data",header_type) == 0) getpos++;
    printf("... found record of type %s pos = %d!\n", header_type, getpos);
    if(getpos == position) break;
  }
  if(status == LIME_EOF) {
    fprintf(stderr, "no scidac-binary-data record found in file %s\n",filename);
    limeDestroyReader(limereader);
    fclose(ifs);
    return(-1);
  }
  bytes = limeReaderBytes(limereader);
  if(bytes == LX*LY*LZ*T*(uint64_t)(24*sizeof(double))) prec = 64;
  else if(bytes == LX*LY*LZ*T*(uint64_t)(24*sizeof(float))) prec = 32;
  else {
    fprintf(stderr, "wrong length in eospinor: bytes = %lu, not %lu. Aborting read!\n", 
	    (unsigned long)bytes, (unsigned long)(LX*LY*LZ*T*(uint64_t)(24*sizeof(double))));
    return(-1);
  }
  printf("# %lu Bit precision read\n", (unsigned long)prec);

  status = read_binary_spinor_data(s, limereader, prec, ts, checksum, &get_index, T, LX, LY, LZ);

  if(status < 0) {
    fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n", 
	    status, filename);
    exit(500);
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

// difference to read_lime_spinor: data is stored at the beginning of s.

int read_lime_spinor_timeslice(double *const s, char *filename, const int position, const int ts, const unsigned int T, const unsigned int LX, const unsigned int LY, const unsigned int LZ) {
  FILE * ifs;
  int status=0, getpos=-1;
  n_uint64_t bytes;
  char * header_type;
  LimeReader * limereader;
  n_uint64_t prec = 32;
  DML_Checksum checksum;
  
  if((ifs = fopen(filename, "r")) == (FILE*)NULL) {
    fprintf(stderr, "Error opening file %s\n", filename);
    return(-1);
  }

  limereader = limeCreateReader( ifs );
  if( limereader == (LimeReader *)NULL ) {
    fprintf(stderr, "Unable to open LimeReader\n");
    return(-1);
  }
  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
    if(status != LIME_SUCCESS ) {
      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
      status = LIME_EOF;
      break;
    }
    header_type = limeReaderType(limereader);
    if(strcmp("scidac-binary-data",header_type) == 0) getpos++;
    if(getpos == position) break;
  }
  if(status == LIME_EOF) {
    fprintf(stderr, "no scidac-binary-data record found in file %s\n",filename);
    limeDestroyReader(limereader);
    fclose(ifs);
    return(-1);
  }
  bytes = limeReaderBytes(limereader);
  if(bytes == LX*LY*LZ*T*(uint64_t)(24*sizeof(double))) prec = 64;
  else if(bytes == LX*LY*LZ*T*(uint64_t)(24*sizeof(float))) prec = 32;
  else {
    fprintf(stderr, "wrong length in eospinor: bytes = %llu, not %llu. Aborting read!\n", 
	    bytes, LX*LY*LZ*T*(uint64_t)(24*sizeof(double)));
    return(-1);
  }
  fprintf(stderr, "# %llu Bit precision read\n", prec);

  status = read_binary_spinor_data_timeslice(s, limereader, prec, ts, checksum, &get_index, T, LX, LY, LZ);

  if(status < 0) {
    fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n", 
	    status, filename);
    exit(500);
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
