// ********************



// fields.cc

// Author: Marc Wagner
// Date: October 2007



// ********************



#include "fields.hh"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "geometry.hh"
#include "linear_algebra.hh"



// ********************



// Allocates and frees memory for a gauge field (lattice size T * L^3).

void Gauge_Field_Alloc(double **gauge_field, int T, int L)
{
  fprintf(stdout,
	  "void Gauge_Field_Alloc(...   -->   Trying to allocate %d M ...",
	  T*L*L*L * 4 * 18 * sizeof(double) / 1000000);

  if((*gauge_field = (double *)malloc(T*L*L*L * 4 * 18 * sizeof(double))) ==
     NULL)
    {
      fprintf(stderr, "\nError: void Gauge_Field_Alloc(...\n");
      exit(EXIT_FAILURE);
    }

  fprintf(stdout, " o.k.\n");
}

void Mu_Fixed_Gauge_Field_Alloc(double **gauge_field, int T, int L)
{
  fprintf(stdout,
    "void Mu_Fixed_Gauge_Field_Alloc(...   -->   Trying to allocate %d M ...",
    T*L*L*L * 18 * sizeof(double) / 1000000);

  if((*gauge_field = (double *)malloc(T*L*L*L * 18 * sizeof(double))) ==
     NULL)
    {
      fprintf(stderr, "\nError: void Mu_Fixed_Gauge_Field_Alloc(...\n");
      exit(EXIT_FAILURE);
    }

  fprintf(stdout, " o.k.\n");
}

void Timeslice_Gauge_Field_Alloc(double **gauge_field, int L)
{
  fprintf(stdout,
    "void Timeslice_Gauge_Field_Alloc(...   -->   Trying to allocate %d M ...",
    L*L*L * 4 * 18 * sizeof(double) / 1000000);

  if((*gauge_field = (double *)malloc(L*L*L * 4 * 18 * sizeof(double))) ==
     NULL)
    {
      fprintf(stderr, "\nError: void Gauge_Field_Alloc(...\n");
      exit(EXIT_FAILURE);
    }

  fprintf(stdout, " o.k.\n");
}

void Gauge_Field_Free(double **gauge_field)
{
  if(*gauge_field == NULL)
    {
      fprintf(stderr, "Error: void Gauge_Field_Free(...\n");
      exit(EXIT_FAILURE);
    }

  free(*gauge_field);
  *gauge_field = NULL;
}



// ********************



// Copies a gauge field (lattice size T * L^3).

void Gauge_Field_Copy(double *gauge_field_dst, double *gauge_field_src, int T, int L)
{
  memcpy(gauge_field_dst, gauge_field_src, T*L*L*L * 4 * 18 * sizeof(double));
}

// Copies a timeslice of a gauge field (lattice size T * L^3).

void Timeslice_Gauge_Field_Copy(double *timeslice_gauge_field_dst, double *gauge_field_src, int T, int L, int timeslice)
{
  memcpy(timeslice_gauge_field_dst,
	 gauge_field_src + ggi(get_index(timeslice, 0, 0, 0, T, L), 0),
	 L*L*L * 4 * 18 * sizeof(double));
}



// ********************



// Generates a unit gauge field.

void Gauge_Field_Unity(double *gauge_field, int T, int L)
{
  int i1;
  int it, ix, iy, iz;

  for(it = 0; it < T; it++)
    {
      for(ix = 0; ix < L; ix++)
	{
	  for(iy = 0; iy < L; iy++)
	    {
	      for(iz = 0; iz < L; iz++)
		{
		  for(i1 = 0; i1 < 4; i1++)
		    {
		      int index = ggi(get_index(it, ix, iy, iz, T, L), i1);

		      cm_eq_id(gauge_field + index);
		    }
		}
	    }
	}
    }  
}



// ********************



// Allocates and frees memory for a spinor field (lattice size T * L^3).

void Spinor_Field_Alloc(double **spinor_field, int T, int L)
{
  fprintf(stdout,
	  "void Spinor_Field_Alloc(...   -->   Trying to allocate %e M ...",
	  T*L*L*L * 24 * sizeof(double) / (double)1000000);

  if((*spinor_field = (double *)malloc(T*L*L*L * 24 * sizeof(double))) ==
     NULL)
    {
      fprintf(stderr, "\nError: void Spinor_Field_Alloc(...\n");
      exit(EXIT_FAILURE);
    }

  fprintf(stdout, " o.k.\n");
}

void Timeslice_Spinor_Field_Alloc(double **spinor_field, int L)
{
  fprintf(stdout,
    "void Timeslice_Spinor_Field_Alloc(...   -->   Trying to allocate %d M ...",
    L*L*L * 24 * sizeof(double) / 1000000);

  if((*spinor_field = (double *)malloc(L*L*L * 24 * sizeof(double))) ==
     NULL)
    {
      fprintf(stderr, "\nError: void Timeslice_Spinor_Field_Alloc(...\n");
      exit(EXIT_FAILURE);
    }

  fprintf(stdout, " o.k.\n");
}

void Spinor_Field_Free(double **spinor_field)
{
  if(*spinor_field == NULL)
    {
      fprintf(stderr, "Error: void Spinor_Field_Free(...\n");
      exit(EXIT_FAILURE);
    }

  free(*spinor_field);
  *spinor_field = NULL;
}



// ********************



// Copies a timeslice of a spinor field (lattice size T * L^3).

void Timeslice_Spinor_Field_Copy(double *timeslice_spinor_field_dst, double *spinor_field_src, int T, int L, int timeslice)
{
  memcpy(timeslice_spinor_field_dst,
	 spinor_field_src + gsi(get_index(timeslice, 0, 0, 0, T, L)),
	 L*L*L * 24 * sizeof(double));
}

void Timeslice_Spinor_Field_Copy_(double *spinor_field_dst, double *timeslice_spinor_field_src, int T, int L, int timeslice)
{
  memcpy(spinor_field_dst + gsi(get_index(timeslice, 0, 0, 0, T, L)),
	 timeslice_spinor_field_src,
	 L*L*L * 24 * sizeof(double));
}

void Spinor_Field_Copy(double *timeslice_spinor_field_dst, double *spinor_field_src, int T, int L)
{
  memcpy(timeslice_spinor_field_dst,
	 spinor_field_src,
	 L*L*L*T * 24 * sizeof(double));
}



// ********************



// Prints the non-zero entries of a spinor field to file.

void fprintf_Spinor_Field_Not_0(FILE *file, double *spinor_field, int T, int L, double epsilon)
{
  int i1;
  int it, ix, iy, iz;

  for(it = 0; it < T; it++)
    {
      for(ix = 0; ix < L; ix++)
	{
	  for(iy = 0; iy < L; iy++)
	    {
	      for(iz = 0; iz < L; iz++)
		{
		  int index = gsi(get_index(it, ix, iy, iz, T, L));

		  for(i1 = 0; i1 < 24; i1++)
		    {
		      if(fabs(spinor_field[index+ i1]) > epsilon)
			break;
		    }

		  if(i1 < 24)
		    {
		      fprintf(file, "t=%2d, x=%2d, y=%2d, z=%2d:\n",
			      it, ix, iy, iz);

		      fv_fprintf(spinor_field + index, file);
		    }
		}
	    }
	}
    }
}



// ********************



// Prints all entries of a gauge field to file. (added by S.D.)

void fprintf_Gauge_Field(FILE *file, double *gauge_field, int T, int L)
{
	int it, ix, iy, iz, mu;

	for(it = 0; it < T; it++)
	{
		for(ix = 0; ix < L; ix++)
		{
			for(iy = 0; iy < L; iy++)
			{
				for(iz = 0; iz < L; iz++)
				{
					for(mu = 0; mu < 4; mu++)
					{
						int index = ggi(get_index(it, ix, iy, iz, T, L),mu);
						fprintf(file, "t=%2d, x=%2d, y=%2d, z=%2d, mu=%d:\n",
						        it, ix, iy, iz, mu);
						cm_fprintf(gauge_field + index, file);
					}
				}
			}
		}
	}
}



// ********************



// Generates a random gauge transformation (allocate memory with
// void Mu_Fixed_Gauge_Field_Alloc(...).

double DRand(double min, double max)
{
  return min + (max-min) * ( (rand() + 0.5) / (RAND_MAX + 1.0) );
}

void Gauge_Trafo_Random(double *g, int T, int L)
{
  int i1;
  int it, ix, iy, iz;

  for(it = 0; it < T; it++)
    {
      for(ix = 0; ix < L; ix++)
	{
	  for(iy = 0; iy < L; iy++)
	    {
	      for(iz = 0; iz < L; iz++)
		{
		  unsigned long int index = get_index(it, ix, iy, iz, T, L) * 18;
		  double *g_ = g + index;

		  for(i1 = 0; i1 < 18; i1++)
		    {
		      g_[i1] = DRand(-1.0, +1.0);
		    }

		  cm_proj(g_);
		}
	    }
	}
    }
}



// ********************



void Gauge_Trafo_Apply_gauge_field(double *g, double *gauge_field, int T, int L)
{
  // U --> g U g^\dagger

  int it, ix, iy, iz, mu;
  double SU3_1[18];

  for(it = 0; it < T; it++)
    {
      for(ix = 0; ix < L; ix++)
	{
	  for(iy = 0; iy < L; iy++)
	    {
	      for(iz = 0; iz < L; iz++)
		{
		  for(mu = 0; mu < 4; mu++)
		    {
		      unsigned long int index_gf =
			ggi(get_index(it, ix, iy, iz, T, L), mu);

		      unsigned long int index_g_l =
			get_index(it, ix, iy, iz, T, L) * 18;

		      unsigned long int index_g_r;

		      if(mu == 0)
			index_g_r = get_index(it+1, ix, iy, iz, T, L) * 18;
		      if(mu == 1)
			index_g_r = get_index(it, ix+1, iy, iz, T, L) * 18;
		      if(mu == 2)
			index_g_r = get_index(it, ix, iy+1, iz, T, L) * 18;
		      if(mu == 3)
			index_g_r = get_index(it, ix, iy, iz+1, T, L) * 18;

  cm_eq_cm_ti_cm(SU3_1, g + index_g_l, gauge_field + index_gf);
  cm_eq_cm_ti_cm_dag(gauge_field + index_gf, SU3_1, g + index_g_r);
		    }
		}
	    }
	}
    }
}

void Gauge_Trafo_Apply_spinor_field(double *g, double *spinor_field, int T, int L)
{
  // s --> g s

  int it, ix, iy, iz;
  double spinor1[24];

  for(it = 0; it < T; it++)
    {
      for(ix = 0; ix < L; ix++)
	{
	  for(iy = 0; iy < L; iy++)
	    {
	      for(iz = 0; iz < L; iz++)
		{
		  unsigned long int index_s =
		    gsi(get_index(it, ix, iy, iz, T, L));

		  unsigned long int index_g =
		    get_index(it, ix, iy, iz, T, L) * 18;

  fv_eq_cm_ti_fv(spinor1, g + index_g, spinor_field + index_s);
  fv_eq_fv(spinor_field + index_s, spinor1);
		}
	    }
	}
    }
}



// ********************
