#include <stdio.h>
#include <math.h>
#include <string.h>

#include <vector>
#include <algorithm>
using namespace std;
#include "symbol.h"

#include <stdarg.h>
extern "C" int vsscanf( const char*, const char*, va_list );

struct EigenState
{
  double energy;
  vector<double> vcoeff;

  void set( char* ptr ){
    double coeff;
    int offset;

    sscanf( ptr, "EIGEN %lf%n", &energy, &offset );
    ptr += offset;

    while( 1 == sscanf( ptr, "%lf%n", &coeff, &offset ) ){
      ptr += offset;
      vcoeff.push_back(coeff);
    }
  }

  void output( void ){
    printf("EIGEN %+f", energy );
    for( int i=0; i<(int)vcoeff.size(); i++ ){
      printf(" %+f", vcoeff[i] );
    }
    printf("\n");
  }

  void arrangeSign( const EigenState& es ) {
    double max_v=0.0;
    int    max_i=0;
    for( int i=0; i<(int)es.vcoeff.size(); i++ ){
      if( fabs(max_v)<fabs(es.vcoeff[i]) ){
	max_v=es.vcoeff[i];
	max_i=i;
      }
    }

    if( max_v * vcoeff[max_i] < 0.0 ){
      for( int i=0; i<(int)vcoeff.size(); i++ ){
	vcoeff[i] *= -1.0;
      }
    }
  }

  static bool checkDegeneracy( const EigenState& es1, const EigenState& es2 ){
    // Disabled this function by MIZUHO-IR. 2012/08/22.
    // Degenerated bands can not be identified correctly by the present method.
    return false;
    //    return fabs(es1.energy-es2.energy)<1.0e-4;
  }

  static void mixDegeneracy( vector<EigenState>& ves, vector<int>& vindex ){
    if( vindex.empty() ) return;

    double sum;
    sum = 0.0;
    for( int m=0; m<(int)vindex.size(); m++ ){
      int n = vindex[m];
      sum += ves[n].energy;
    }
    sum /= vindex.size();
    for( int m=0; m<(int)vindex.size(); m++ ){
      int n = vindex[m];
      ves[n].energy = sum;
    }

    for( int i=0; i<(int)ves.front().vcoeff.size(); i++ ){
      sum = 0.0;
      for( int m=0; m<(int)vindex.size(); m++ ){
	int n = vindex[m];
	sum += ves[n].vcoeff[i] * ves[n].vcoeff[i];
      }
      sum = sqrt(sum/vindex.size());
      for( int m=0; m<(int)vindex.size(); m++ ){
	int n = vindex[m];
	ves[n].vcoeff[i] = sum;
      }
    }
  }

  static double compare( const EigenState& es1, const EigenState& es2 ){
    double norm1=0.0;
    double norm2=0.0;
    double prod12=0.0;

    for( int i=0; i<(int)es1.vcoeff.size(); i++ ){
      norm1  += es1.vcoeff[i] * es1.vcoeff[i];
      norm2  += es2.vcoeff[i] * es2.vcoeff[i];
      prod12 += es1.vcoeff[i] * es2.vcoeff[i];
    }

    return prod12*prod12/(norm1*norm2);
  }
};

namespace ELSES
{
  double total;
  vector<EigenState> ves;

  int checkNumberOfBasis( char* elem ){
    static vector<char*> velem;
    static vector<int  > vbasis;

    for( int i=0; i<(int)velem.size(); i++ ){
      if( strcmp( velem[i], elem ) == 0 ) return vbasis[i];
    }

    char fname[32];
    sprintf(fname,"%s.xml",elem);
    FILE* fptr = fopen(fname,"rb");
    if( fptr == NULL ){
      velem.push_back(strdup(elem));
      vbasis.push_back(0);
      return 0;
    }

    char word[64];
    int  basis=0;
    while( 1 == fscanf(fptr,"%s",word) ){
      if( strcmp(word,"<principal_quantum_number>") == 0 ){
	break;
      }
    }
    while( 1 == fscanf(fptr,"%s",word) ){
      if( strcmp(word,"</principal_quantum_number>") == 0 ){
	break;
      }
      basis++;
    }

    // begin modify by MIZUHO-IR 2012/08/22
    // fixed a bug that miss-understanded the number of kinds of basis
    switch( basis ){
    case 1 : // 2S or 1*2S
      basis = 1; // S only
      break;
    case 2 : // 2S 3*2P
      basis = 2; // S and P
      break;
    case 3 : // 3S 3*3P 5*3D
      basis = 3; // S P D
      break;
    case 4 : // 2S 2P 2P 2P
      basis = 2; // S and P
      break;
    case 5 : // 3S 3P 3P 3P 5*3D
      basis = 3; // S P D
      break;
    case 7 : // 3S 3*3P 3D 3D 3D 3D 3D
      basis = 3; // S P D
      break;
    case 9 : // 3S 3P 3P 3P 3D 3D 3D 3D 3D
      basis = 3; // S P D
      break;
    }
    // end modify by MIZUHO-IR 2012/08/22
    //    fprintf(stderr,"DEBUG %s %d\n", elem, basis );

    velem.push_back(strdup(elem));
    vbasis.push_back(basis);

    return basis;
  }

  bool getData(  const char* file_mol,
		 const char* file_energy,
		const char* file_vcoeff,
		const char* file_total ){
    char buf[1024], term1[64], term2[64];
    double energy, occup;

    double energy_homo;
    int num_eigen;
    {
      FILE* fptr_energy = fopen( file_energy, "r" );
      if( fptr_energy == NULL ){
	fprintf( stderr,"Error: can not open file %s\n", file_energy );
	return false;
      }

      num_eigen = 0;
      while( fgets(buf,sizeof(buf),fptr_energy) ){
	if( 2 == sscanf(buf,"%*d %lf %lf", &energy, &occup ) ){
	  energy *= 1.0/27.211384588;

	  num_eigen++;

	  ves.push_back( EigenState() );
	  ves.back().energy = energy;
	  if( occup > 0.1 ){
	    energy_homo = energy;
	  }
	}
      }
      fclose(fptr_energy);
    }

    for( int i=0; i<(int)ves.size(); i++ ){
      ves[i].vcoeff.resize(num_eigen);
    }


    vector<bool> voutershell;

    {
      FILE* fptr_mol = fopen( file_mol, "r" );
      if( fptr_mol == NULL ){
	fprintf( stderr,"Error: can not open file %s\n", file_mol );
	return false;
      }

      char elem[4];
      while( fgets(buf,sizeof(buf),fptr_mol) ){
	if( 1 == sscanf( buf, "<atom element=\"%[^\"]", elem ) ){
	  int basis = checkNumberOfBasis(elem);
	  if( basis == 0 ) basis = 2; // supporse SP as default 

	  // skip outer bases
	  if( false );
	  else if( !strcmp( elem, "H" ) ||
		   !strcmp( elem, "He" ) ){
	    voutershell.push_back(false); // 1S
	  }
	  else if( !strcmp( elem, "Li" ) ||
		   !strcmp( elem, "Be" ) ||
		   !strcmp( elem, "B"  ) ||
		   !strcmp( elem, "C"  ) ||
		   !strcmp( elem, "N"  ) ||
		   !strcmp( elem, "O"  ) ||
		   !strcmp( elem, "F"  ) ||
		   !strcmp( elem, "Ne" ) ){
	    if( basis >= 1 ){
	      voutershell.push_back(false); // 2S
	    }
	    if( basis >= 2 ){
	      voutershell.push_back(false); // 2Px
	      voutershell.push_back(false); // 2Py
	      voutershell.push_back(false); // 2Pz
	    }
	    if( basis >= 3 ){
	      voutershell.push_back(true); // 3D
	      voutershell.push_back(true); // 3D
	      voutershell.push_back(true); // 3D
	      voutershell.push_back(true); // 3D
	      voutershell.push_back(true); // 3D
	    }
	  }
	  else if( !strcmp( elem, "Na" ) ||
		   !strcmp( elem, "Mg" ) ||
		   !strcmp( elem, "Al" ) ||
		   !strcmp( elem, "Si" ) ||
		   !strcmp( elem, "P"  ) ||
		   !strcmp( elem, "S"  ) ||
		   !strcmp( elem, "Cl" ) ||
		   !strcmp( elem, "Ar" ) ){
	    if( basis >= 1 ){
	      voutershell.push_back(false); // 3S
	    }
	    if( basis >= 2 ){
	      voutershell.push_back(false); // 3Px
	      voutershell.push_back(false); // 3Py
	      voutershell.push_back(false); // 3Pz
	    }
	    if( basis >= 3 ){
	      voutershell.push_back(true); // 3D
	      voutershell.push_back(true); // 3D
	      voutershell.push_back(true); // 3D
	      voutershell.push_back(true); // 3D
	      voutershell.push_back(true); // 3D
	    }
	  }
	  else if( !strcmp( elem, "K"  ) ||
		   !strcmp( elem, "Ca" ) ||
		   !strcmp( elem, "Sc" ) ||
		   !strcmp( elem, "Ti" ) ||
		   !strcmp( elem, "V"  ) ||
		   !strcmp( elem, "Cr" ) ||
		   !strcmp( elem, "Mn" ) ||
		   !strcmp( elem, "Fe" ) ||
		   !strcmp( elem, "Co" ) ||
		   !strcmp( elem, "Ni" ) ||
		   !strcmp( elem, "Cu" ) ||
		   !strcmp( elem, "Zn" ) ||
		   !strcmp( elem, "Ga" ) ||
		   !strcmp( elem, "Ge" ) ||
		   !strcmp( elem, "As" ) ||
		   !strcmp( elem, "Se" ) ||
		   !strcmp( elem, "Br" ) ||
		   !strcmp( elem, "Kr" ) ){
	    if( basis >= 1 ){
	      voutershell.push_back(false); // 3S
	    }
	    if( basis >= 2 ){
	      voutershell.push_back(false); // 3Px
	      voutershell.push_back(false); // 3Py
	      voutershell.push_back(false); // 3Pz
	    }
	    if( basis >= 3 ){
	      voutershell.push_back(false); // 3D
	      voutershell.push_back(false); // 3D
	      voutershell.push_back(false); // 3D
	      voutershell.push_back(false); // 3D
	      voutershell.push_back(false); // 3D
	    }
	  }
	}
      }

      fclose(fptr_mol);
    }

    int skip_num=0;

    /*
      // do not skip, use higher unoccupied states for comparison with Gaussian
    for( int i=0; i<(int)voutershell.size(); i++ ){
      if( voutershell[i] ) skip_num++;
    }
    */

    FILE* fptr_vcoeff = fopen( file_vcoeff, "r" );
    if( fptr_vcoeff == NULL ){
      fprintf( stderr,"Error: can not open file %s\n", file_vcoeff );
      return false;
    }

    int base_number;
    double coeff;
    
    while( fgets(buf,sizeof(buf),fptr_vcoeff) ){
      if( buf[19] == '.' ) break;
    }

    for( int i=0; i+skip_num<(int)ves.size(); i++ ){
      for( int j=0; j<num_eigen; j++ ){
	strncpy( term1, buf+ 8,  3 );
	strncpy( term2, buf+16, 24 );

	if( 1 == sscanf( term1, "%d",  &base_number ) &&
	    1 == sscanf( term2, "%lf", &coeff ) ){
	  ves[i].vcoeff[j] = coeff;
	}
	fgets(buf,sizeof(buf),fptr_vcoeff);
      }
    }


    fclose(fptr_vcoeff);

    FILE* fptr_total = fopen( file_total, "r" );
    if( fptr_total == NULL ){
      fprintf( stderr,"Error: can not open file %s\n", file_total );
      return false;
    }

    total = 0.0;
    while( fgets(buf,sizeof(buf),fptr_total) ){
      if( 1 == sscanf( buf, " EBD+ECSC+ECC Energy : ETOT [au]:"
		       "             %lf", &total ) ){
	break;
      }
    }
    fclose(fptr_total);

    if( total == 0.0 ){
      return false;
    }

    // output data
    printf("TOTAL %+f\n", total );
    printf("NUM   %d\n", (int)ves.size()-skip_num );
    for( int i=0; i+skip_num<(int)ves.size(); i++ ){
      printf("EIGEN %+f", ves[i].energy );
      for( int j=0; j<(int)ves[i].vcoeff.size(); j++ ){
	if( voutershell[j] ) continue;
	printf(" %+f", ves[i].vcoeff[j]  );
      }
      printf("\n");
    }
    printf("HOMO  %+f\n", energy_homo );

    return true;
  }

}

namespace GAUSS
{
  double total;
  vector<EigenState> ves;

  int mscanf( char* buf, size_t is, size_t ie, const char* const format, ... ){
    char work[1024];
    strncpy( work, buf+is, ie-is+1 );
    work[ie-is+1] = '\0';

    va_list argptr;
    va_start( argptr, format );

    return vsscanf( work, format, argptr );
  }

  // Comparative list of STO-3G bases for ELSES
  bool isComparableBaseSTO3G( bool& inner, const char* elemname, const char* basename ){
    unsigned int elem_id = ElementSymbol(elemname);

    inner = false;

    if( false );
    else if( ElementSymbol::H  <= elem_id && elem_id <= ElementSymbol::He ){
      if( !strcmp( basename, "1S" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Li <= elem_id && elem_id <= ElementSymbol::Ne ){
      if( !strcmp( basename, "1S" ) ){
	inner = true;
      }
      if( !strcmp( basename, "2S" ) || !strcmp( basename, "2P" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Na <= elem_id && elem_id <= ElementSymbol::Ar ){
      if( !strcmp( basename, "1S" ) ||
	  !strcmp( basename, "2S" ) || !strcmp( basename, "2P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "3S" ) || !strcmp( basename, "3P" ) ){
	return true;
      }
    }
    else if( ElementSymbol::K  <= elem_id && elem_id <= ElementSymbol::Ca ){
      if( !strcmp( basename, "1S" ) ||
	  !strcmp( basename, "2S" ) || !strcmp( basename, "2P" ) ||
	  !strcmp( basename, "3S" ) || !strcmp( basename, "3P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "4S" ) || !strcmp( basename, "4P" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Sc <= elem_id && elem_id <= ElementSymbol::Kr ){
      if( !strcmp( basename, "1S" ) ||
	  !strcmp( basename, "2S" ) || !strcmp( basename, "2P" ) ||
	  !strcmp( basename, "3S" ) || !strcmp( basename, "3P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "4S" ) || !strcmp( basename, "4P" ) ||
	  !strcmp( basename, "5D" ) ){
	return true;
      }
    }
    return false;
  }

  // Comparative list of 3-21G bases for ELSES
  bool isComparableBase321G( bool& inner, const char* elemname, const char* basename ){
    unsigned int elem_id = ElementSymbol(elemname);

    inner = false;

    if( false );
    else if( ElementSymbol::H  <= elem_id && elem_id <= ElementSymbol::He ){
      if( !strcmp( basename, "1S" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Li <= elem_id && elem_id <= ElementSymbol::Ne ){
      if( !strcmp( basename, "1S" ) ){
	inner = true;
      }
      if( !strcmp( basename, "2S" ) || !strcmp( basename, "2P" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Na <= elem_id && elem_id <= ElementSymbol::Ar ){
      if( !strcmp( basename, "1S" ) ||
	  !strcmp( basename, "2S" ) || !strcmp( basename, "2P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "3S" ) || !strcmp( basename, "3P" ) ){
	return true;
      }
    }
    else if( ElementSymbol::K  <= elem_id && elem_id <= ElementSymbol::Ca ){
      if( !strcmp( basename, "1S" ) ||
	  !strcmp( basename, "2S" ) || !strcmp( basename, "2P" ) ||
	  !strcmp( basename, "3S" ) || !strcmp( basename, "3P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "4S" ) || !strcmp( basename, "4P" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Sc <= elem_id && elem_id <= ElementSymbol::Zn ){
      if( !strcmp( basename, "1S" ) ||
	  !strcmp( basename, "2S" ) || !strcmp( basename, "2P" ) ||
	  !strcmp( basename, "3S" ) || !strcmp( basename, "3P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "4S" ) || !strcmp( basename, "4P" ) ||
	  !strcmp( basename, "6D" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Ga <= elem_id && elem_id <= ElementSymbol::Kr ){
      if( !strcmp( basename, "1S" ) ||
	  !strcmp( basename, "2S" ) || !strcmp( basename, "2P" ) ||
	  !strcmp( basename, "3S" ) || !strcmp( basename, "3P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "4S" ) || !strcmp( basename, "4P" ) ){
	return true;
      }
    }
    return false;
  }

  // Comparative list of 3-21G* bases for ELSES
  bool isComparableBase321Gd( bool& inner, const char* elemname, const char* basename ){
    unsigned int elem_id = ElementSymbol(elemname);

    inner = false;

    if( false );
    else if( ElementSymbol::H  <= elem_id && elem_id <= ElementSymbol::He ){
      if( !strcmp( basename, "1S" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Li <= elem_id && elem_id <= ElementSymbol::Ne ){
      if( !strcmp( basename, "1S" ) ){
	inner = true;
      }
      if( !strcmp( basename, "2S" ) || !strcmp( basename, "2P" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Na <= elem_id && elem_id <= ElementSymbol::Ar ){
      if( !strcmp( basename, "1S" ) ||
	  !strcmp( basename, "2S" ) || !strcmp( basename, "2P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "3S" ) || !strcmp( basename, "3P" ) ){
	return true;
      }
    }
    else if( ElementSymbol::K  <= elem_id && elem_id <= ElementSymbol::Ca ){
      if( !strcmp( basename, "1S" ) ||
	  !strcmp( basename, "2S" ) || !strcmp( basename, "2P" ) ||
	  !strcmp( basename, "3S" ) || !strcmp( basename, "3P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "4S" ) || !strcmp( basename, "4P" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Sc <= elem_id && elem_id <= ElementSymbol::Zn ){
      if( !strcmp( basename, "1S" ) ||
	  !strcmp( basename, "2S" ) || !strcmp( basename, "2P" ) ||
	  !strcmp( basename, "3S" ) || !strcmp( basename, "3P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "4S" ) || !strcmp( basename, "4P" ) ||
	  !strcmp( basename, "6D" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Ga <= elem_id && elem_id <= ElementSymbol::Kr ){
      if( !strcmp( basename, "1S" ) ||
	  !strcmp( basename, "2S" ) || !strcmp( basename, "2P" ) ||
	  !strcmp( basename, "3S" ) || !strcmp( basename, "3P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "4S" ) || !strcmp( basename, "4P" ) ){
	return true;
      }
    }
    return false;
  }

  // Comparative list of 6-311G bases for ELSES
  bool isComparableBase6311G( bool& inner, const char* elemname, const char* basename ){
    unsigned int elem_id = ElementSymbol(elemname);

    inner = false;

    if( false );
    else if( ElementSymbol::H  <= elem_id && elem_id <= ElementSymbol::He ){
      if( !strcmp( basename, "1S" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Li <= elem_id && elem_id <= ElementSymbol::Ne ){
      if( !strcmp( basename, "1S" ) ){
	inner = true;
      }
      if( !strcmp( basename, "2S" ) || !strcmp( basename, "2P" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Na <= elem_id && elem_id <= ElementSymbol::Ar ){
      if( !strcmp( basename, "1S" ) || !strcmp( basename, "2S" ) ||
	  !strcmp( basename, "7P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "3S" ) || !strcmp( basename, "8P" ) ){
	return true;
      }
    }
    else if( ElementSymbol::K  <= elem_id && elem_id <= ElementSymbol::Ca ){
      if( !strcmp( basename, "1S" ) || !strcmp( basename, "2S" ) ||
	  !strcmp( basename, "3S" ) ||
	  !strcmp( basename, "9P" ) || !strcmp( basename, "10P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "4S" ) || !strcmp( basename, "11P" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Sc <= elem_id && elem_id <= ElementSymbol::Zn ){
      if( !strcmp( basename, "1S" ) || !strcmp( basename, "2S" ) ||
	  !strcmp( basename, "3S" ) ||
	  !strcmp( basename, "10P" ) || !strcmp( basename, "11P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "4S" ) || !strcmp( basename, "12P" ) ||
	  !strcmp( basename, "15D" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Ga <= elem_id && elem_id <= ElementSymbol::Kr ){
      if( !strcmp( basename, "1S" ) || !strcmp( basename, "2S" ) ||
	  !strcmp( basename, "3S" ) ||
	  !strcmp( basename, "9P" ) || !strcmp( basename, "10P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "4S" ) || !strcmp( basename, "11P" ) ){
	return true;
      }
    }
    return false;
  }

  // Comparative list of 6-311G* bases for ELSES
  bool isComparableBase6311Gd( bool& inner, const char* elemname, const char* basename ){
    unsigned int elem_id = ElementSymbol(elemname);

    inner = false;

    if( false );
    else if( ElementSymbol::H  <= elem_id && elem_id <= ElementSymbol::He ){
      if( !strcmp( basename, "1S" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Li <= elem_id && elem_id <= ElementSymbol::Ne ){
      if( !strcmp( basename, "1S" ) ){
	inner = true;
      }
      if( !strcmp( basename, "2S" ) || !strcmp( basename, "2P" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Na <= elem_id && elem_id <= ElementSymbol::Ar ){
      if( !strcmp( basename, "1S" ) || !strcmp( basename, "2S" ) ||
	  !strcmp( basename, "7P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "3S" ) || !strcmp( basename, "8P" ) ){
	return true;
      }
    }
    else if( ElementSymbol::K  <= elem_id && elem_id <= ElementSymbol::Ca ){
      if( !strcmp( basename, "1S" ) || !strcmp( basename, "2S" ) ||
	  !strcmp( basename, "3S" ) ||
	  !strcmp( basename, "9P" ) || !strcmp( basename, "10P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "4S" ) || !strcmp( basename, "11P" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Sc <= elem_id && elem_id <= ElementSymbol::Zn ){
      if( !strcmp( basename, "1S" ) || !strcmp( basename, "2S" ) ||
	  !strcmp( basename, "3S" ) ||
	  !strcmp( basename, "10P" ) || !strcmp( basename, "11P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "4S" ) || !strcmp( basename, "12P" ) ||
	  !strcmp( basename, "15D" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Ga <= elem_id && elem_id <= ElementSymbol::Kr ){
      if( !strcmp( basename, "1S" ) || !strcmp( basename, "2S" ) ||
	  !strcmp( basename, "3S" ) ||
	  !strcmp( basename, "9P" ) || !strcmp( basename, "10P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "4S" ) || !strcmp( basename, "11P" ) ){
	return true;
      }
    }
    return false;
  }


  // Comparative list of 6-311G** bases for ELSES
  bool isComparableBase6311Gdp( bool& inner, const char* elemname, const char* basename ){
    unsigned int elem_id = ElementSymbol(elemname);

    inner = false;

    if( false );
    else if( ElementSymbol::H  <= elem_id && elem_id <= ElementSymbol::He ){
      if( !strcmp( basename, "1S" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Li <= elem_id && elem_id <= ElementSymbol::Ne ){
      if( !strcmp( basename, "1S" ) ){
	inner = true;
      }
      if( !strcmp( basename, "2S" ) || !strcmp( basename, "2P" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Na <= elem_id && elem_id <= ElementSymbol::Ar ){
      if( !strcmp( basename, "1S" ) || !strcmp( basename, "2S" ) ||
	  !strcmp( basename, "7P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "3S" ) || !strcmp( basename, "8P" ) ){
	return true;
      }
    }
    else if( ElementSymbol::K  <= elem_id && elem_id <= ElementSymbol::Ca ){
      if( !strcmp( basename, "1S" ) || !strcmp( basename, "2S" ) ||
	  !strcmp( basename, "3S" ) ||
	  !strcmp( basename, "9P" ) || !strcmp( basename, "10P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "4S" ) || !strcmp( basename, "11P" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Sc <= elem_id && elem_id <= ElementSymbol::Zn ){
      if( !strcmp( basename, "1S" ) || !strcmp( basename, "2S" ) ||
	  !strcmp( basename, "3S" ) ||
	  !strcmp( basename, "10P" ) || !strcmp( basename, "11P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "4S" ) || !strcmp( basename, "12P" ) ||
	  !strcmp( basename, "15D" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Ga <= elem_id && elem_id <= ElementSymbol::Kr ){
      if( !strcmp( basename, "1S" ) || !strcmp( basename, "2S" ) ||
	  !strcmp( basename, "3S" ) ||
	  !strcmp( basename, "9P" ) || !strcmp( basename, "10P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "4S" ) || !strcmp( basename, "11P" ) ){
	return true;
      }
    }
    return false;
  }

  // Comparative list of 6-311+G** bases for ELSES
  bool isComparableBase6311pGdp( bool& inner, const char* elemname, const char* basename ){
    unsigned int elem_id = ElementSymbol(elemname);

    inner = false;

    if( false );
    else if( ElementSymbol::H  <= elem_id && elem_id <= ElementSymbol::He ){
      if( !strcmp( basename, "1S" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Li <= elem_id && elem_id <= ElementSymbol::Ne ){
      if( !strcmp( basename, "1S" ) ){
	inner = true;
      }
      if( !strcmp( basename, "2S" ) || !strcmp( basename, "2P" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Na <= elem_id && elem_id <= ElementSymbol::Ar ){
      if( !strcmp( basename, "1S" ) || !strcmp( basename, "2S" ) ||
	  !strcmp( basename, "7P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "3S" ) || !strcmp( basename, "8P" ) ){
	return true;
      }
    }
    else if( ElementSymbol::K  <= elem_id && elem_id <= ElementSymbol::Ca ){
      if( !strcmp( basename, "1S" ) || !strcmp( basename, "2S" ) ||
	  !strcmp( basename, "3S" ) ||
	  !strcmp( basename, "9P" ) || !strcmp( basename, "10P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "4S" ) || !strcmp( basename, "11P" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Sc <= elem_id && elem_id <= ElementSymbol::Zn ){
      if( !strcmp( basename, "1S" ) || !strcmp( basename, "2S" ) ||
	  !strcmp( basename, "3S" ) ||
	  !strcmp( basename, "11P" ) || !strcmp( basename, "12P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "4S" ) || !strcmp( basename, "13P" ) ||
	  !strcmp( basename, "18D" ) ){
	return true;
      }
    }
    else if( ElementSymbol::Ga <= elem_id && elem_id <= ElementSymbol::Kr ){
      if( !strcmp( basename, "1S" ) || !strcmp( basename, "2S" ) ||
	  !strcmp( basename, "3S" ) ||
	  !strcmp( basename, "9P" ) || !strcmp( basename, "10P" ) ){
	inner = true;
      }
      if( !strcmp( basename, "4S" ) || !strcmp( basename, "11P" ) ){
	return true;
      }
    }
    return false;
  }

  bool getData( const char* file_output, const char* basetype ){
    FILE* fptr = fopen( file_output, "r" );
    if( fptr == NULL ){
      fprintf( stderr,"Error: can not open file %s\n", file_output );
      return false;
    }

    char buf[1024], term[16];
    double eigen[5];
    char   occup[5][2];
    double coeff[5];
    int base_number;
    char elem_name[8];
    char base_name[8];


    vector<double> veigen;
    vector<char>   voccup;
    vector< vector<double> > vvcoeff;
    vector<char*>   velemname;
    vector<char*>   vbasename;

    total = 0.0;
    bool flag=false;
    char elem_prev[3]="  ";
    while( fgets(buf,sizeof(buf),fptr) ){
      // start of molecular orbital coefficients
      if( !strncmp(buf,"     Alpha Molecular Orbital Coefficients:", 42 ) ||
	  !strncmp(buf,"     Molecular Orbital Coefficients:", 36 ) ){
	flag  = true;
	continue;
      }

      // end of molecular orbital coefficients
      if( !strncmp(buf,"     Beta Molecular Orbital Coefficients:",  41 ) ||
	  !strncmp(buf,"     Density Matrix:", 20 ) ){
	flag  = false;
	continue;
      }

      // inside molecular orbital coefficients
      if( flag ){
	int status;

	if( false );
	// get eigenvalues
	else if( (status=sscanf( buf,
				 "     Eigenvalues -- %lf %lf %lf %lf %lf",
				 &eigen[0], &eigen[1],&eigen[2], &eigen[3], &eigen[4] ) )>0 ){
	  for( int s=0; s<status; s++ ){
	    veigen.push_back(eigen[s]);
	  }
	}

	// get occupasions
	else if( (status=sscanf( buf,
				 " %[OV] %[OV] %[OV] %[OV] %[OV]",
				 occup[0], occup[1], occup[2], occup[3], occup[4] ) )>0 ){
	  for( int s=0; s<status; s++ ){
	    voccup.push_back(occup[s][0]);
	  }
	}


	// get coefficients
	// example of this line
	//   1 1   P  1S          0.00000   0.00000   0.00000   0.00000   0.00000
	//   2        2S          0.00000   0.00000   0.00000   0.00000   0.00000
	//01234567890123456789012345678901234567890123456789012345678901234567890123456789
	else if( 1==mscanf( buf,0,3,"%d", &base_number ) &&
		 1==mscanf( buf,11,13,"%s", base_name ) &&
		 (status=mscanf( buf,20,72,"%lf %lf %lf %lf %lf",
				 &coeff[0], &coeff[1], &coeff[2],
				 &coeff[3], &coeff[4] ) )>0 ){
	  if( base_number > (int)vvcoeff.size() ){
	    vvcoeff.resize(base_number);
	    velemname.resize(base_number);
	    vbasename.resize(base_number);
	  }

	  if( 1==mscanf( buf, 9,10, "%s", elem_name ) ){
	    // backup elem for the next elem
	    strcpy( elem_prev, elem_name );
	  }else{
	    // in case of empty elem, use prev elem name
	    strcpy( elem_name, elem_prev );
	  }

	  velemname[base_number-1] = strdup(elem_name);
	  vbasename[base_number-1] = strdup(base_name);

	  for( int s=0; s<status; s++ ){
	    vvcoeff[base_number-1].push_back(coeff[s]);
	  }
	}

	// ignore indeces of molecular orbital
	else{
	}
      }

      // outside molecular orbital coefficients
      else{
	if( !strncmp( buf, " SCF Done:", 10 ) &&
	    1 == sscanf( buf+24, "%lf", &total ) ){
	}
      }
    }

    fclose(fptr);

    for( int j=0; j<(int)vvcoeff.size(); j++ ){
      if( (int)vvcoeff[j].size() != base_number ){
	fprintf( stderr,"Error: eigen functions broken in file %s.\n", file_output );
	return false;
      }
    }

    // Check each bases if it is comparable for ELSES base.
    vector<bool> vcompbase(vbasename.size());

    int skip_inner_levels=0;
    int skip_outer_levels=0;
    bool inner;

    for( int j=0; j<(int)vbasename.size(); j++ ){
      if( false );
      else if( !strcmp( basetype, "STO-3G" ) ){
	vcompbase[j] = isComparableBaseSTO3G( inner, velemname[j], vbasename[j] );
      }
      else if( !strcmp( basetype, "3-21G" ) ){
	vcompbase[j] = isComparableBase321G( inner, velemname[j], vbasename[j] );
      }
      else if( !strcmp( basetype, "3-21G*" ) ){
	vcompbase[j] = isComparableBase321Gd( inner, velemname[j], vbasename[j] );
      }
      else if( !strcmp( basetype, "6-311G" ) ){
	vcompbase[j] = isComparableBase6311G( inner, velemname[j], vbasename[j] );
      }
      else if( !strcmp( basetype, "6-311G*" ) ){
	vcompbase[j] = isComparableBase6311Gd( inner, velemname[j], vbasename[j] );
      }
      else if( !strcmp( basetype, "6-311G**" ) ){
	vcompbase[j] = isComparableBase6311Gdp( inner, velemname[j], vbasename[j] );
      }
      else if( !strcmp( basetype, "6-311+G**" ) ){
	vcompbase[j] = isComparableBase6311pGdp( inner, velemname[j], vbasename[j] );
      }
      else{
	fprintf( stderr,"Error: unknown base type %s\n", basetype );
	return false;
      }
      if( !vcompbase[j] ){
	if( inner ){
	  skip_inner_levels++;
	}
	else{
	  skip_outer_levels++;
	}
      }
    }

    // for debug
    printf("# the basetype is %s\n", basetype );
    printf("# bases comparable or not for ELSES\n");
    for( int j=0; j<(int)vbasename.size(); j++ ){
      printf("# %2d %2s %2s %s\n", j, velemname[j], vbasename[j],
	     vcompbase[j] ? "comparable" : "not comparable" );
    }
    printf("# eigen energies comparable or not for ELSES\n");
    for( int i=0; i<(int)veigen.size(); i++ ){
      if( i<skip_inner_levels || (int)veigen.size()-skip_outer_levels<=i ){
	printf("# %2d %+f not comarable\n", i, veigen[i] );
      }
      else{
	printf("# %2d %+f comarable\n", i, veigen[i] );
      }
    }

    double energy_homo;
    for( int i=0; i<(int)veigen.size(); i++ ){
      if( i<skip_inner_levels || (int)veigen.size()-skip_outer_levels<=i ) continue;

      EigenState es;
      es.energy = veigen[i];
      if( voccup[i] == 'O' ){
	energy_homo = veigen[i];
      }

      for( int j=0; j<(int)vvcoeff.size(); j++ ){
	if( !vcompbase[j] ) continue;

	es.vcoeff.push_back( vvcoeff[j][i] );
      }
      
      ves.push_back(es);
    }

    if( total == 0.0 ){
      return false;
    }

    // output data
    printf("TOTAL %+f\n", total );
    printf("NUM   %d\n", (int)ves.size() );
    for( int i=0; i<(int)ves.size(); i++ ){
      printf("EIGEN %+f", ves[i].energy );
      for( int j=0; j<(int)ves[i].vcoeff.size(); j++ ){
	printf(" %+f", ves[i].vcoeff[j]  );
      }
      printf("\n");
    }
    printf("HOMO  %+f\n", energy_homo );

    return true;
  }
}

namespace ELSES_SORT
{
  double total;
  double energy_homo;
  vector<EigenState> ves_elses;
  vector<EigenState> ves_gauss;

  bool getData( const char* file_elses, const char* file_gauss ){
    char buf[1024];

    FILE* fptr_elses = fopen(file_elses, "r");
    FILE* fptr_gauss = fopen(file_gauss, "r");

    if( fptr_elses == NULL ){
      fprintf( stderr,"Error: can not open file %s\n", file_elses );
      return 1;
    }
    if( fptr_gauss == NULL ){
      fprintf( stderr,"Error: can not open file %s\n", file_gauss );
      return 1;
    }


    while( fgets(buf,sizeof(buf),fptr_elses ) ){
      if( 1 == sscanf( buf, "TOTAL %lf", &total ) ){
      }
      else if( 1 == sscanf( buf, "HOMO %lf", &energy_homo ) ){
      }
      else if( !strncmp( buf, "EIGEN ", 6 ) ){
	EigenState es;
	es.set(buf);
	ves_elses.push_back(es);
      }
    }

    while( fgets(buf,sizeof(buf),fptr_gauss ) ){
      if( strncmp( buf, "EIGEN ", 6 ) ) continue;

      EigenState es;
      es.set(buf);
      ves_gauss.push_back(es);
    }

    fclose(fptr_elses);
    fclose(fptr_gauss);

    {
      vector<int> vindex;

      vindex.clear();
      for( int n=0; n+1<(int)ves_gauss.size(); n++ ){
	if( EigenState::checkDegeneracy( ves_gauss[n], ves_gauss[n+1]) ){
	  if( vindex.empty() ){
	    vindex.push_back(n);
	  }
	  vindex.push_back(n+1);
	}
	else if( !vindex.empty() ){
	  EigenState::mixDegeneracy( ves_gauss, vindex );
	  vindex.clear();
	}
      }
      if( !vindex.empty() ){
	EigenState::mixDegeneracy( ves_gauss, vindex );
	vindex.clear();
      }

      vindex.clear();
      for( int n=0; n+1<(int)ves_elses.size(); n++ ){
	if( EigenState::checkDegeneracy( ves_elses[n], ves_elses[n+1]) ){
	  if( vindex.empty() ){
	    vindex.push_back(n);
	  }
	  vindex.push_back(n+1);
	}
	else if( !vindex.empty() ){
	  EigenState::mixDegeneracy( ves_elses, vindex );
	  vindex.clear();
	}
      }
      if( !vindex.empty() ){
	EigenState::mixDegeneracy( ves_elses, vindex );
	vindex.clear();
      }
    }


    int mo_index[(int)ves_gauss.size()];
    for(int i=0;i<(int)ves_gauss.size();i++){
      mo_index[i]=i;
    }

    printf("TOTAL %+f\n", total );
    printf("NUM   %d\n", (int)ves_gauss.size() );
    for( int j=0; j<(int)ves_gauss.size(); j++ ){
      double match_max=0.0;
      int    match_i  =j;
      for( int i=j; i<(int)ves_elses.size(); i++ ){
	double v = EigenState::compare( ves_elses[i], ves_gauss[j] );
	if( match_max<v ){
	  match_max=v;
	  match_i  =i;
	}
      }
      std::swap( ves_elses[j], ves_elses[match_i] );
      std::swap( mo_index[j], mo_index[match_i]);

      ves_elses[j].arrangeSign( ves_gauss[j] );
      ves_elses[j].output();

      printf("MO match: gaussian %d: ELSES %d: production %f \n",
	     j, mo_index[j], match_max);

    }
    printf("HOMO  %f\n", energy_homo );
  }
}


main( int argc, char* argv[] )
{
  if( argc <= 1 ){
    printf("%s elses ...\n", argv[0] );
    printf("%s gaussian ...\n", argv[0] );
    printf("%s elses-sort ...\n", argv[0] );
    return 1;
  }

  // in case of ELSES
  if( !strcmp( argv[1], "elses" ) ){
    if( argc != 6 ){
      printf("%s elses molecule.xml output_eigen_levels.txt output_wavefunction.txt Output.txt\n", argv[0] );
      return 1;
    }
    if( !ELSES::getData(argv[2],argv[3],argv[4],argv[5]) ){
      return 1;
    }

  }

  // in case of Gaussian
  else if( !strcmp( argv[1], "gaussian" ) ){
    if( argc != 4 ){
      printf("%s gaussian molecule.out basetype\n", argv[0] );
      return 1;
    }
    if( !GAUSS::getData(argv[2],argv[3]) ){
      return 1;
    }
  }


  // in case of Gaussian
  else if( !strcmp( argv[1], "elses-sort" ) ){
    if( argc != 4 ){
      printf("%s elses-sort elses_data.dat gaussian_data.dat\n", argv[0] );
      return 1;
    }
    if( !ELSES_SORT::getData(argv[2],argv[3]) ){
      return 1;
    }
  }



  else{
    printf("%s elses molecule.xml output_eigen_levels.txt output_wavefunction.txt Output.txt\n", argv[0] );
    printf("%s gaussian molecule.out\n", argv[0] );
    printf("%s sort elses.dat gaussian.dat\n", argv[0] );
    return 1;
  }

}

