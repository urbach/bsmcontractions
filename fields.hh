// ********************



// fields.hh

// Author: Marc Wagner
// Date: October 2007



// ********************



#ifndef __FIELDS_HH__

#define __FIELDS_HH__



// ********************



#include <stdio.h>



// ********************



// Allocates and frees memory for a gauge field (lattice size T * L^3).

void Gauge_Field_Alloc(double **gauge_field, int T, int L);

void Mu_Fixed_Gauge_Field_Alloc(double **gauge_field, int T, int L);

void Timeslice_Gauge_Field_Alloc(double **gauge_field, int L);

void Gauge_Field_Free(double **gauge_field);



// Copies a gauge field (lattice size T * L^3).

void Gauge_Field_Copy(double *gauge_field_dst, double *gauge_field_src, int T, int L);

// Copies a timeslice of a gauge field (lattice size T * L^3).

void Timeslice_Gauge_Field_Copy(double *timeslice_gauge_field_dst, double *gauge_field_src, int T, int L, int timeslice);



// Generates a unit gauge field.

void Gauge_Field_Unity(double *gauge_field, int T, int L);



// Allocates and frees memory for a spinor field (lattice size T * L^3).

void Spinor_Field_Alloc(double **spinor_field, int T, int L);

void Timeslice_Spinor_Field_Alloc(double **spinor_field, int L);

void Spinor_Field_Free(double **spinor_field);



// Copies a timeslice of a spinor field (lattice size T * L^3).

void Timeslice_Spinor_Field_Copy(double *timeslice_spinor_field_dst, double *spinor_field_src, int T, int L, int timeslice);

void Timeslice_Spinor_Field_Copy_(double *spinor_field_dst, double *timeslice_spinor_field_src, int T, int L, int timeslice);

void Spinor_Field_Copy(double *timeslice_spinor_field_dst, double *spinor_field_src, int T, int L);



// Prints the non-zero entries of a spinor field to file.
void fprintf_Spinor_Field_Not_0(FILE *file, double *spinor_field, int T, int L, double epsilon = 0.000001);

// Prints all entries of a gauge field to file.
void fprintf_Gauge_Field(FILE *file, double *gauge_field, int T, int L);



// Generates a random gauge transformation (allocate memory with
// void Mu_Fixed_Gauge_Field_Alloc(...).
void Gauge_Trafo_Random(double *g, int T, int L);

void Gauge_Trafo_Apply_gauge_field(double *g, double *gauge_field, int T, int L);
void Gauge_Trafo_Apply_spinor_field(double *g, double *spinor_field, int T, int L);



// ********************



#endif



// ********************
