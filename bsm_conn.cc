// Code to perform heavy-heavy (1+1) contractions
// 
// AUTHOR: Konstantin Ottnad, Carsten Urbach
//
//
// DATE:   20101215
//
//============================== INCLUDE SECTION ======================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <getopt.h>
#include <errno.h>
#include "fields.hh"
#include "geometry.hh"
#include "io.hh"
#include "propagator_io.hh"
#include "linear_algebra.hh"
//#include "Q_phi.hh"
#include "contract_twopoint.hh"
//#include "fuzz.hh"
//#include "smearing_techniques.hh"


//============================== GLOBAL DECLARATIONS AND STUFF ========================

using namespace std;

int verbose = 0;
int no_fields = 24; 
const int no_scalars = 4;                                       // number of (pseudo-) scalar gamma combinations
const int no_vectors = 0;                                       // number of (axial-) vector gamma combinations
const int n_s = 4;                                              // number of spin components
const int num_level = 2; // it's just a dummy variable right now ... 

// further variables are defined local to main()



//============================== FUNCTIONS ============================================ 



// This function prints program info and options and terminates the program
void usage()
{
  cout << endl << "This program performs heavy-heavy neutral connected contractions" << endl << endl << "Usage: bsm_conn [options]" << endl << endl;
  cout << "Options:" << endl << endl << "-L (required) set the spatial extension of the lattice" << endl;
  cout << "-N (required) no of the gauge configuration" << endl;
  cout << "-T (required) set the temporal extension of the lattice" << endl;
  cout << "-c (optional) set n_c [default=3]" << endl;
  cout << "-s (optional) read propagators from separate files [default=off]" << endl;
  cout << "-t (optional) timeslice [default=0]" << endl;
  cout << "-n (optional) number of scalar configs per gauge [default=1]" << endl;
  cout << "-v (optional) level of verbosity; values: 0 (none) [default], 1 (some), 2 (all)" << endl << endl;
  exit(0);
}



// This function processes the commandline params and performs some error / consistency checks on them
int read_params(int argc, char **argv, int *N, int *L, int *T, int *n_c, int *splitted, int *timeslice, int *nscalar)
{
  int i = 0;

  while ((i = getopt(argc,argv,"h?fsN:L:T:c:l:n:t:v:")) !=-1) {
    switch (i) {
    case 'N':
      *N = atoi(optarg); 
      break;
    case 'L':
      *L = atoi(optarg); 
      break;
    case 'T':
      *T = atoi(optarg); 
      break;
    case 'c':
      *n_c = atoi(optarg);
      break;
    case 's':
      *splitted = 1;
      break;
    case 't':
      *timeslice = atoi(optarg); 
      break;
    case 'v':
      verbose = atoi(optarg); 
      break; 
    case 'n':
      *nscalar = atoi(optarg);
      break;
    case 'h':
    case '?':
    default:
      usage(); 
    }
  }
  
  if (*N<0)
  { 
    cout << endl << "FAIL: the number of the gauge configuration >N< must be specified and has to be larger or equal to zero" << endl;
    return -1;
  }

  if (*L<=0)
  {
    cout << endl << "FAIL: L must be specified and has to be positive" << endl;
    return -1;
  }
 
  if (*T<=0) 
  {
    cout << endl << "FAIL: T must specified and has to be larger than zero" << endl;
    return -1;
  }

  if (*n_c<=0)
  {
    cout << endl << "FAIL: n_c must be positive" << endl;
    return -1;
  }

  if ((*timeslice<0)||(*timeslice>=*T))
  {
    cout << endl << "FAIL: t must be in the range [0,T=" << *T << ")" << endl;
    return -1;
  }

  verbose = (verbose>2) ? 2 : verbose; // check for (and if necessary set) the correct verbosity value 
  verbose = (verbose<0) ? 0 : verbose;

  return 0;
}



// This function reads a gauge configuration from a file 'conf.gggg' where gggg denotes the number of the gauge configuration given by N
int read_gauge_field(double **gauge_field, const int N, const int L, const int T)
{
  stringstream filename;

  filename << "conf.";
  filename.width(4);
  filename.fill('0');
  filename << N << ends;
  if (verbose>=1)
  {
    cout << "Reading gauge configuration no " << N << " ... ";
  }
  Gauge_Field_Alloc(gauge_field, T, L);  // get memory for the gauge field
  if (read_lime_gauge_field_doubleprec((*gauge_field), filename.str().c_str(), T, L, L, L)!=0) return -1; // read it

  if (verbose>=1)  
  {
    cout << "done!" << endl;
  }
  return 0;
}



//  This function reads the heavy sources from a file 'source.gggg.tt.ii.hinverted'
//  with 'gggg' being the number of the cgauge configuration, 'tt' the timeslice and 'ii' the sample number
int read_heavy_field(double *heavy_field, const int N, const int L, const int T, const int timeslice, const int ii, const int pos)
{
  stringstream filename;
  
  filename.str(""); // reset the filename string object 
  filename.width(1);
  filename << "propagator.";
  filename.width(4);
  filename.fill('0');
  filename << N;
  filename.width(1);
  filename << ".";
  filename.width(2);
  filename << timeslice;
  filename.width(1);
  filename << ".";
  filename.width(2);
  filename << ii;
  filename.width(1);
  filename << ".inverted" << ends;

  if (verbose>=1)  {
    cout << "Reading heavy propagators from file " << filename.str().c_str() << " ... "; 
  } 

  if (read_lime_spinor(heavy_field,(char *)filename.str().c_str(), pos, -1, T, L, L, L) != 0) return -1; // read it
  
  if (verbose>=1) {
    cout << "done!" << endl;
  }
  return 0;
}


// AUXILIARY FUNCTIONS

inline int get_pos1(int sector)
// returns the position of the first propagator in a source. file corresponding to the sector index
{
  int i = (sector / 4);
  if(i > 1) return(i + 2);
  return i;
}

inline int get_pos2(int sector)
// return the position of the second propagator in a source. file corresponding to the sector index
{ 
  int i = (sector % 4) + 2;
  if(i > 3) return(i + 2);
  return i;
} 

// END OF AUXILIARY FUNCTIONS



//============================== MAIN PROGRAM =========================================

int main(int argc, char *argv[])
{
  int N = -1, L = -1, T = -1;
  int splitted = 0, timeslice = 0, n_c = 3;
 
  int * permutation[16]; // required for the generation of the gamma matrix combinations
  int * sign[16];
  const int sfnr = 2;

  double *gauge_field, *auxiliary_field;

  complex * sca = NULL;
  complex * sf[sfnr];
  // the number of scalar configs to loop over with the same gauge.
  int nscalar = 1;

  // INITIALIZATION

  if (read_params(argc, argv, &N, &L, &T, &n_c, &splitted, &timeslice, &nscalar)!=0) 
  {
    return -1;
  }
  no_fields *= n_c;

  double *heavy_field[no_fields]; // need the correct no_fields fot this one (in case of fuzzing)  
  
  if((void*)(sca = (complex*)calloc(sfnr*T*L*L*L, sizeof(complex))) == NULL) {
    printf ("scalar field malloc errno : %d\n",errno); 
    errno = 0;
    return(1);
  }
  for(int i = 0; i < sfnr; i++) {
    sf[i] = &sca[i*T*L*L*L];
  }
  //read_scalar_field("hallo", sf, T, L);

  // END OF INITALIZATION



  // GAMMA MATRICES
  //
  // Index mapping conventions are as follows:
  //
  // index --> gamma combination
  // 
  // 0  --> g_0
  // 1  --> g_1
  // 2  --> g_2
  // 3  --> g_3
  // 4  --> Id
  // 5  --> g_5
  // 6  --> g_0 g_5
  // 7  --> g_5 g_1
  // 8  --> g_5 g_2
  // 9  --> g_5 g_3
  // 10 --> g_0 g_1
  // 11 --> g_0 g_2
  // 12 --> g_0 g_3
  // 13 --> g_1 g_2
  // 14 --> g_1 g_3
  // 15 --> g_2 g_3
   
  for(int i = 0; i < 16; i++) // there are always 16 of them... 
  {
    permutation[i] = (int *)malloc(24*sizeof(int)); // memory allocation
  }

  for(int i = 0; i < 16; i++) 
  {
    sign[i] = (int *)malloc(24*sizeof(int));
  }

  for(int i = 0; i < 6; i++) 
  {
    for(int j = 0; j < 24; j++) 
    {
      permutation[i][j] = gamma_permutation[i][j];
      sign[i][j] = gamma_sign[i][j];
    }
  }

  // these are the remaining products (10), same def as in neutral.cc / heavylight.cc
  gamma_eq_gamma_ti_gamma(permutation[6], sign[6], permutation[0], permutation[5], sign[0], sign[5]);    // g0*g5 
  gamma_eq_gamma_ti_gamma(permutation[7], sign[7], permutation[5], permutation[1], sign[5], sign[1]);    // g5*g1
  gamma_eq_gamma_ti_gamma(permutation[8], sign[8], permutation[5], permutation[2], sign[5], sign[2]);    // g5*g2
  gamma_eq_gamma_ti_gamma(permutation[9], sign[9], permutation[5], permutation[3], sign[5], sign[3]);    // g5*g3
  gamma_eq_gamma_ti_gamma(permutation[10], sign[10], permutation[0], permutation[1], sign[0], sign[1]);  // g0*g1
  gamma_eq_gamma_ti_gamma(permutation[11], sign[11], permutation[0], permutation[2], sign[0], sign[2]);  // g0*g2
  gamma_eq_gamma_ti_gamma(permutation[12], sign[12], permutation[0], permutation[3], sign[0], sign[3]);  // g0*g3
  gamma_eq_gamma_ti_gamma(permutation[13], sign[13], permutation[1], permutation[2], sign[1], sign[2]);  // g1*g2
  gamma_eq_gamma_ti_gamma(permutation[14], sign[14], permutation[1], permutation[3], sign[1], sign[3]);  // g1*g3
  gamma_eq_gamma_ti_gamma(permutation[15], sign[15], permutation[2], permutation[3], sign[2], sign[3]);  // g2*g3


  // these are the gamma combinations which are required for the neutral two point contractions
  // the first elements are (pseudo-) scalars followed by (axial-) vectors. These two groups must NOT be mixed!
  // NOTE: do not forget to update 'no_scalars', 'no_vectors' etc.  whenever changing the size of this array.
  const unsigned int neutral_gamma_index[][2] = {{5,5},{5,4},{4,5},{4,4}};


  // these arrays contain info about wether the resulting correlator for a gamma combination is imaginary and needs to be corrected (i.e. Xchg real and imag part of the result) 
  int neutral_is_imaginary[] = {0, 0, 0, 0};  // #total == no_scalars + no_vectors - "i"-vectors are probably wrong


  // these arrays contain the relative signs for the addition of the three 'i'-components of (axial-) vector gamma combinations
  double neutral_vsign[] = {1., 1., 1., 1.};  // #total == 3 * no_vectors


  // END OF GAMMA MATRICES


  // INPUT + SMEARING AND FUZZING

  Gauge_Field_Alloc(&gauge_field, T, L);  // get memory for the gauge field
  if (read_gauge_field(&gauge_field, N, L, T)!=0)
  {
    return -1;
  }

  for (int i = 0; i < no_fields; i++) 
    Spinor_Field_Alloc(&(heavy_field[i]),T,L); // get memory 

  complex C[1 * 16 * (no_scalars + no_vectors) * T];  // memory for the contraction results
  // Need: 1 for local-local only right now
  //       16 for all quark combinations (u,d,u',d')
  //       (no_scalars + no_vectors) for all source-sink gamma pairs (mind the summation over vector components, if applicable!)
  //       T for T :)
  for (int t = 0; t < 1*16*(no_scalars + no_vectors)*T; t++) {
    C[t].re = 0.;
    C[t].im = 0.;
  }


  for (int sector = 0; sector < 16; sector++)
  {
    cout << "Sector " << sector << endl;
    // this loop runs over all 16 relevant (see below) connected quark combinations, i.e.
    // lets call the two flavour components u and d
    //
    // sector_id | sector | pos1 | pos2 in source. file
    //
    //  uuuu         0        0      0
    //  uuud         1        0      1
    //  uudu         2        0      2
    //  uudd         3        0      3
    //  uduu         4        1      0
    //  udud         5        1      1
    //  uddu         6        1      2
    //  uddd         7        1      3
    //  duuu         8        2      0
    //  duud         9        2      1
    //  dudu         10       2      2
    //  dudd         11       2      3
    //  dduu         12       3      0
    //  ddud         13       3      1
    //  dddu         14       3      2
    //  dddd         15       3      3
    //
    // the order in the propagator file is
    // source-sink
    //   uu
    //   ud
    //  (uu)^+
    //  (ud)^+
    //   du
    //   dd
    //  (du)^+
    //  (dd)^+
    //
    // Note that these rules describe the corresponding loops clockwise:
    //
    //         propagator1
    //       ________\______
    //      /        /	   \
    //     /		    \
    //  quark3             quark4
    //    |
    // Gamma_sink        Gamma_source
    //    |                   |
    //  quark2             quark1  <- start here
    //     \                 /
    //      \________/______/
    //		     \
    //        propagator2
    

    int ll = 0; // index to fix the position of local / smearing combinations in C[] 

    for (int j = 0; j < 1; j++) {
      // loop over all possible local-local, local-smeared, smeared-local, smeared-smeared combinations (or the fuzzed instead of the smeared version thereof)
      
      if (j==0) { // local - local
	
        ll = 0;
        for (int ii = 0; ii < n_s*n_c ; ii++) {
	  if ((read_heavy_field(heavy_field[ii], N, L, T, timeslice, ii, get_pos1(sector)) != 0) ||            // propagator 1
	      (read_heavy_field(heavy_field[ii + n_s*n_c], N, L, T, timeslice, ii, get_pos2(sector)) != 0)) {  // propagator 2
	    return -1;
	  }
	}
      }
      
      // END OF INPUT 
      
      // CONTRACTIONS
      // now we perform the contractions, looping over all gamma combinations
      // note that the flavor-off-diagonal elements sectors aquire an additional minus sign due to the gamma_5 trick, which is NOT being taken care of in the code !!! (It should be taken care of "analytically")

      unsigned int result_index = sector * 1 * (no_scalars + no_vectors) * T + ll * (no_scalars + no_vectors) * T;

      // (pseudo-) scalars first...
      for (int gamma = 0; gamma < no_scalars; gamma++) {
        contract_twopoint(&C[result_index], 
			  permutation[neutral_gamma_index[gamma][0]], sign[neutral_gamma_index[gamma][0]], 
			  permutation[neutral_gamma_index[gamma][1]], sign[neutral_gamma_index[gamma][1]], 
			  permutation[4], sign[4], // no gamma5 trick here, so this is not g5 but 1!
			  heavy_field, &heavy_field[n_s*n_c], T, L, timeslice, n_c);  
	// this is the most general version of the contract_twopoint function...
        result_index+=T;
      }

      // END OF CONTRACTIONS
    } // end of j-loop
  }   // end of sector-loop



  // OUTPUT SECTION
  for (int sector = 0; sector < 16; sector++) {
    // write a separate file for each quark-combination
    
    ofstream ofs;
    ostringstream filename;
    unsigned int result_index = 0;
    const string sector_id[] = {"uuuu","uuud","uudu","uudd","uduu","udud","uddu","uddd","duuu","duud","dudu","dudd","dduu","ddud","dddu","dddd"};
    filename.str("");  // reset the filename string object
    filename << "outprcvn."; // base filename
    filename << sector_id[sector] << "."; // add the sector identifier 
    filename.width(2);
    filename.fill('0');
    filename << timeslice; // the timeslice index
    filename.width(1);
    filename << ".";
    filename.width(4);
    filename << N << ends; // number of the gauge configuration
    filename.width(1);
    ofs.open(filename.str().c_str());
    //    ofs << " " << N << " " << timeslice << " " << L << " " << L << " " << L << " " << T << " " << scientific << kappa << " " << sector_id[sector] << endl; // write some info in the first line of the file
    for (int gamma = 0; gamma < (no_scalars + no_vectors); gamma++) {
      for (int ll = 0; ll<1; ll++) {
        for (int t = 0; t<T/2+1; t++) {
          result_index = 1 * sector * (no_scalars + no_vectors) * T + gamma * T + ll * (no_scalars + no_vectors) * T;
          
          ofs.width(3);
          ofs.fill(' ');
          ofs << gamma + 1;  // gamma combination index -- mind the "+1" convention for connected data
          ofs.width(3);
          ofs << 2 * ll + 1;  // local / smearing / fuzzing combination
          ofs.width(4);
          ofs << t;
          ofs << " ";
          ofs.width(15);


	  ofs << scientific << C[result_index + t].re / L / L / L << " ";
	  ofs.width(15);
	  ofs << scientific << C[result_index + t].im / L / L / L << " ";
	  if ((t==0) || (t==T/2)) {
	    ofs.width(15);
	    ofs << scientific << 0.0 << " " << 0.0 << endl;
	  }
          else {
	    ofs.width(15);
	    ofs << scientific << C[result_index + T - t].re / L / L / L << " ";
	    ofs.width(15);
	    ofs << scientific << C[result_index + T - t].im / L / L / L << endl;
	  }
        }
      }
    }
    
    ofs.close();
  }
  // END OF OUTPUT
  
  // CLEAN UP
  Gauge_Field_Free(&gauge_field);
  for (int i = 0; i < no_fields; i++) Spinor_Field_Free(&(heavy_field[i]));

  for (int i=0; i<16; i++) {
    delete [] permutation[i];
    delete [] sign[i];
  }

  return 0;
}
// END OF MAIN PROGRAM



