// ********************



// geometry.hh

// Author: Carsten Urbach
// Date: ???

// Author: Marc Wagner
// Date: October 2007



// ********************



#ifndef __GEOMETRY_HH__

#define __GEOMETRY_HH__



// ********************



inline unsigned long int get_index(const int t, const int x, const int y, const int z, const int T, const int L)
{
  unsigned long int tt = (t+T)%T;
  unsigned long int xx = (x+L)%L;
  unsigned long int yy = (y+L)%L;
  unsigned long int zz = (z+L)%L;

  return ((tt*L+xx)*L+yy)*L+zz;
}

inline unsigned long int get_index_timeslice(const int x, const int y, const int z, const int T, const int L)
{
  unsigned long int xx = (x+L)%L;
  unsigned long int yy = (y+L)%L;
  unsigned long int zz = (z+L)%L;

  return (xx*L+yy)*L+zz;
}

inline unsigned long int get_index_timeslice_t(const int t, const int x, const int y, const int z, const int T, const int L)
{
  unsigned long int xx = (x+L)%L;
  unsigned long int yy = (y+L)%L;
  unsigned long int zz = (z+L)%L;

  return (xx*L+yy)*L+zz;
}



inline unsigned long int ggi(const unsigned long int ix, const int mu)
{
  return (4*ix+mu) * (unsigned long int)18;
}



inline unsigned long gsi(const unsigned long ix)
{
  return ix*12 * 2;
}



// ********************



#endif



// ********************
