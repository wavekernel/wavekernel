#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/wait.h>
#include <vector>
#include <algorithm>
#include <functional>
using namespace std;

enum ParamID {
  PARAM_Ds,
  PARAM_Dp,
  PARAM_Dd,
  PARAM_Z1s,
  PARAM_Z1p,
  PARAM_Z1d,
  PARAM_Z2s,
  PARAM_Z2p,
  PARAM_Z2d,
  PARAM_CAs,
  PARAM_CAp,
  PARAM_CAd,
  PARAM_Rs,
  PARAM_Rp,
  PARAM_Rd,
  PARAM_CRs,
  PARAM_CRp,
  PARAM_CRd,
};

struct Param {
  ParamID id;
  const char* shortname;
  const char* longname;
  bool  set;

  Param( ParamID id, const char* shortname, const char* longname ){
    this->id = id;
    this->shortname = shortname;
    this->longname  = longname;
    set = false;
    min = 0.0;
    max = 0.0;
  }
  double min, max, init;

  bool setValue( char* value ){
    if( 3 == sscanf(value,"%lf %lf %lf", &min, &max, &init ) ){
      set = true;
    }
    else if( 2 == sscanf(value,"%lf %lf", &min, &max ) ){
      init = (min+max)*0.5;
      set = true;
    }
    else if( 1 == sscanf(value,"%lf", &min ) ){
      max = min;
      init = min;
      set = true;
    }
    else{
      fprintf(stderr,"Error: variable %s is broekn.\n", longname );
      set = false;
    }
    return set;
  }
};


Param param[] = {
  Param( PARAM_Ds,  "Ds",  "param_diag_s" ),
  Param( PARAM_Dp,  "Dp",  "param_diag_p" ),
  Param( PARAM_Dd,  "Dd",  "param_diag_d" ),
  Param( PARAM_Z1s, "Z1s", "param_zeta_s" ),
  Param( PARAM_Z1p, "Z1p", "param_zeta_p" ),
  Param( PARAM_Z1d, "Z1d", "param_zeta_d" ),
  Param( PARAM_Z2s, "Z2s", "param_zeta2_s" ),
  Param( PARAM_Z2p, "Z2p", "param_zeta2_p" ),
  Param( PARAM_Z2d, "Z2d", "param_zeta2_d" ),
  Param( PARAM_CAs, "CAs", "param_carg_s" ),
  Param( PARAM_CAp, "CAp", "param_carg_p" ),
  Param( PARAM_CAd, "CAd", "param_carg_d" ),
  Param( PARAM_Rs,  "Rs",  "param_repulsive_s" ),
  Param( PARAM_Rp,  "Rp",  "param_repulsive_p" ),
  Param( PARAM_Rd,  "Rd",  "param_repulsive_d" ),
  Param( PARAM_CRs, "CRs", "param_cratio_s" ),
  Param( PARAM_CRp, "CRp", "param_cratio_p" ),
  Param( PARAM_CRd, "CRd", "param_cratio_d" )
};

const int DIM = 15;

double weight_eigenerror=1.0, weight_totalerror=10.0;
double criterion_error = 1.0e-7;

struct Parameter
{
  double p[DIM];
  double err,err_ebd,err_etot;

  Parameter( void ){
    for( int d=0; d<DIM; d++ ){
      p[d] = 0.0;
    }
    err = 0.0;
    err_ebd=0.0;
    err_etot=0.0;
  }

  bool operator < ( const Parameter& P ) const {
    return err > P.err;
  }

  Parameter operator + ( const Parameter& P ) const {
    Parameter P2;
    for( int d=0; d<DIM; d++ ){
      P2.p[d] = p[d] + P.p[d];
    }
    return P2;
  }

  Parameter operator - ( const Parameter& P ) const {
    Parameter P2;
    for( int d=0; d<DIM; d++ ){
      P2.p[d] = p[d] - P.p[d];
    }
    return P2;
  }

  Parameter operator * ( const double& f ) const {
    Parameter P2;
    for( int d=0; d<DIM; d++ ){
      P2.p[d] = p[d] * f;
    }
    return P2;
  }

  bool checkRange( void ){
    for( int d=0; d<DIM; d++ ){
      if( p[d] < param[d].min || param[d].max < p[d] ){
	return false;
      }
    }
    return true;
  }

  void evaluate( void ){
    const double err_large = 10000.0;

    if( !checkRange() ){
      printf("# Info: parameter out of range.\n");
      err = err_large; // set very large value to the err
      return;
    }

    if( fork() == 0 ){
      char* arg[1+DIM+1];
      
      arg[0] = "./run.sh";
      for( int d=0; d<DIM; d++ ){
	char tmp[64];
	sprintf( tmp, "%f", p[d]);
	arg[1+d] = strdup(tmp);
      }
      arg[1+DIM] = NULL;

      execv(arg[0],arg);
      perror("execl");
      exit(1);
    }
    int status;
    wait(&status);

    if( status != 0 ){
      printf("# Info: parameter not good for elses.\n");
      err = err_large; // set very large value to the err
      return;
    }

    char fname[1024];
    char buf[128];

    snprintf( fname, sizeof(fname), "energy_difference.dat");
    FILE* fptr = fopen( fname,"r");
    if( fptr == NULL ){
      fprintf(stderr,"# Error: not found file %s\n", fname );
      exit(1);
    }
    // read the err
    err = 0.0;
    err_ebd=0.0;
    err_etot=0.0;

    double e;
    fgets(buf,sizeof(buf),fptr);
    if( 1 != sscanf(buf,"%lf", &e ) ){
      fprintf(stderr,"# Error: can not read data from file %s\n", fname );
      exit(1);
    }
    err += weight_eigenerror*e;
    err_ebd += weight_eigenerror*e;

    fgets(buf,sizeof(buf),fptr);
    if( 1 != sscanf(buf,"%lf", &e ) ){
      fprintf(stderr,"# Error: can not read data from file %s\n", fname );
      exit(1);
    }
    err += weight_totalerror*e;
    err_etot += weight_totalerror*e;

    fclose(fptr);

    print();
  }

  void print( void ){
    for( int d=0; d<DIM; d++ ){
      printf(" %+7.3f", p[d] );
    }
    printf("  %+15.8f %+15.8f %+15.8f\n", err, err_ebd, err_etot );
  }
};


bool isComment( char* buf )
{
  char dummy[1024];

  return buf[0] == '#' || sscanf(buf,"%s", dummy ) == -1;
}

bool getValue( char* buf, const char* name, char* value )
{
  char format[1024];

  sprintf( format, " %s%%*[ =] (%%[^)]", name );

  if( 1 == sscanf(buf,format, value ) ) return true;

  sprintf( format, " %s%%*[ =] %%s", name );
  if( 1 == sscanf(buf,format, value ) ) return true;

  return false;
}

int getEnv( void )
{
  char buf[1024];
  char tag[256];
  char name[1024],value[1024];

  FILE* fptr = fopen("../param.cfg","rb");

  while( fgets(buf,sizeof(buf),fptr) ){
    if( isComment(buf) ) continue;
    if( 1 != sscanf(buf,"%s",tag) ) continue;

    if( strcmp(tag,"config{") == 0 ){
      while( fgets(buf,sizeof(buf),fptr) ){
	if( isComment(buf) ) continue;

	if( 1 != sscanf(buf,"%s",tag) ) continue;

	if( strcmp(tag,"}config") == 0 ){
	  break;
	}
	if( strcmp(tag,"material{") == 0 ){
	  while( fgets(buf,sizeof(buf),fptr) ){
	    if( isComment(buf) ) continue;

	    if( 1 == sscanf(buf,"%s",tag) &&
		strcmp(tag,"}material") == 0 ){
	      break;
	    }

	    if( false );
	    else if( getValue(buf,"element", value ) ){
	      //	      printf("element=\"%s\"\n", value );
	    }
	    else if( getValue(buf,"molecule", value ) ){
	      //	      printf("molecule=\"%s\"\n", value );
	    }
	    else if( getValue(buf,"spin", value ) ){
	      //	      printf("spin=\"%s\"\n", value );
	    }
	    else{
	      fprintf(stderr,"broken line %s", buf);
	    }
	  }
	  continue;
	} // material

	if( strcmp(tag,"conformation{") == 0 ){
	  while( fgets(buf,sizeof(buf),fptr) ){
	    if( isComment(buf) ) continue;

	    if( 1 == sscanf(buf,"%s",tag) &&
		strcmp(tag,"}conformation") == 0 ){
	      break;
	    }

	    if( false );
	    else if( getValue(buf,"bond_optimized", value ) ){
	      //	      printf("bond_optimized=\"%s\"\n", value );
	    }
	    else if( getValue(buf,"angle_optimized", value ) ){
	      //	      printf("angle_optimized=\"%s\"\n", value );
	    }
	    else if( getValue(buf,"bonds_conformation", value ) ){
	      //	      printf("bonds_conformation=\"%s\"\n", value );
	    }
	    else if( getValue(buf,"angles_conformation", value ) ){
	      //	      printf("angles_conformation=\"%s\"\n", value );
	    }
	    else{
	      fprintf(stderr,"broken line %s", buf);
	    }
	  }
	  continue;
	} // conformation

	if( strcmp(tag,"difference{") == 0 ){
	  while( fgets(buf,sizeof(buf),fptr) ){
	    if( isComment(buf) ) continue;

	    if( 1 == sscanf(buf,"%s",tag) &&
		strcmp(tag,"}difference") == 0 ){
	      break;
	    }

	    if( false );
	    else if( getValue(buf,"weight_eigen", value ) ){
	      if( 1 != sscanf( value,"%lf", &weight_eigenerror ) ){
		fprintf(stderr,"Error: variable %s is broekn.\n",
			"weight_eigen" );
		fclose(fptr);
		return 1;
	      }
	      //	      printf("weight_eigen  = [%s]\n", value );
	    }
	    else if( getValue(buf,"weight_total", value ) ){
	      if( 1 != sscanf( value,"%lf", &weight_totalerror ) ){
		fprintf(stderr,"Error: variable %s is broekn.\n",
			"weight_total" );
		fclose(fptr);
		return 1;
	      }
	      //	      printf("weight_total  = [%s]\n", value );
	    }
	    else if( getValue(buf,"convergence", value ) ){
	      if( 1 != sscanf( value,"%lf", &criterion_error ) ){
		fprintf(stderr,"Error: variable %s is broekn.\n",
			"convergence" );
		fclose(fptr);
		return 1;
	      }
	      //	      printf("convergence   = [%s]\n", value );
	    }
	  }
	  continue;
	} // difference

	fprintf(stderr,"broken line %s", buf);
      }
      continue;
    }	// config

    if( strcmp(tag,"param{") == 0 ){
      while( fgets(buf,sizeof(buf),fptr) ){
	if( isComment(buf) ) continue;

	if( 1 != sscanf(buf,"%s",tag) ) continue;

	if( strcmp(tag,"}param") == 0 ){
	  break;
	}
	if( strcmp(tag,"S{") == 0 ){
	  while( fgets(buf,sizeof(buf),fptr) ){
	    if( isComment(buf) ) continue;

	    if( 1 == sscanf(buf,"%s",tag) &&
		strcmp(tag,"}S") == 0 ){
	      break;
	    }

	    if( false );
	    else if( getValue(buf,"n", value ) ){
	      //	      printf("n  = [%s]\n", value );
	    }
	    else if( getValue(buf,"diag", value ) ){
	      if( !param[PARAM_Ds].setValue(value) ){
		fclose(fptr);
		return 1;
	      }
	      //	      printf("diag  = [%s]\n", value );
	    }
	    else if( getValue(buf,"zeta1", value ) ){
	      if( !param[PARAM_Z1s].setValue(value) ){
		fclose(fptr);
		return 1;
	      }
	      //	      printf("zeta1  = [%s]\n", value );
	    }
	    else if( getValue(buf,"zeta2", value ) ){
	      if( !param[PARAM_Z2s].setValue(value) ){
		fclose(fptr);
		return 1;
	      }
	      //	      printf("zeta2  = [%s]\n", value );
	    }
	    else if( getValue(buf,"carg", value ) ){
	      if( !param[PARAM_CAs].setValue(value) ){
		fclose(fptr);
		return 1;
	      }
	      //	      printf("carg  = [%s]\n", value );
	    }
	    else if( getValue(buf,"cratio", value ) ){
	      if( !param[PARAM_CRs].setValue(value) ){
		fclose(fptr);
		return 1;
	      }
	      //	      printf("cratio = [%s]\n", value );
	    }
	    else if( getValue(buf,"repulsive", value ) ){
	      if( !param[PARAM_Rs].setValue(value) ){
		fclose(fptr);
		return 1;
	      }
	      //	      printf("repulsive = [%s]\n", value );
	    }
	    else{
	      fprintf(stderr,"broken line %s", buf);
	    }
	  }
	  continue;
	} // S

	if( strcmp(tag,"P{") == 0 ){
	  while( fgets(buf,sizeof(buf),fptr) ){
	    if( isComment(buf) ) continue;

	    if( 1 == sscanf(buf,"%s",tag) &&
		strcmp(tag,"}P") == 0 ){
	      break;
	    }

	    if( false );
	    else if( getValue(buf,"n", value ) ){
	      //	      printf("n  = [%s]\n", value );
	    }
	    else if( getValue(buf,"diag", value ) ){
	      if( !param[PARAM_Dp].setValue(value) ){
		fclose(fptr);
		return 1;
	      }
	      //	      printf("diag  = [%s]\n", value );
	    }
	    else if( getValue(buf,"zeta1", value ) ){
	      if( !param[PARAM_Z1p].setValue(value) ){
		fclose(fptr);
		return 1;
	      }
	      //	      printf("zeta1  = [%s]\n", value );
	    }
	    else if( getValue(buf,"zeta2", value ) ){
	      if( !param[PARAM_Z2p].setValue(value) ){
		fclose(fptr);
		return 1;
	      }
	      //	      printf("zeta2  = [%s]\n", value );
	    }
	    else if( getValue(buf,"carg", value ) ){
	      if( !param[PARAM_CAp].setValue(value) ){
		fclose(fptr);
		return 1;
	      }
	      //	      printf("carg  = [%s]\n", value );
	    }
	    else if( getValue(buf,"cratio", value ) ){
	      if( !param[PARAM_CRp].setValue(value) ){
		fclose(fptr);
		return 1;
	      }
	      //      	      printf("cratio = [%s]\n", value );
	    }
	    else if( getValue(buf,"repulsive", value ) ){
	      if( !param[PARAM_Rp].setValue(value) ){
		fclose(fptr);
		return 1;
	      }
	      //      	      printf("repulsive = [%s]\n", value );
	    }
	    else{
	      fprintf(stderr,"broken line %s", buf);
	    }
	  }
	  continue;
	} // P

	if( strcmp(tag,"D{") == 0 ){
	  while( fgets(buf,sizeof(buf),fptr) ){
	    if( isComment(buf) ) continue;

	    if( 1 == sscanf(buf,"%s",tag) &&
		strcmp(tag,"}D") == 0 ){
	      break;
	    }

	    if( false );
	    else if( getValue(buf,"n", value ) ){
	      //	      printf("n  = [%s]\n", value );
	    }
	    else if( getValue(buf,"diag", value ) ){
	      if( !param[PARAM_Dd].setValue(value) ){
		fclose(fptr);
		return 1;
	      }
	      //	      printf("diag  = [%s]\n", value );
	    }
	    else if( getValue(buf,"zeta1", value ) ){
	      if( !param[PARAM_Z1d].setValue(value) ){
		fclose(fptr);
		return 1;
	      }
	      //	      printf("zeta1  = [%s]\n", value );
	    }
	    else if( getValue(buf,"zeta2", value ) ){
	      if( !param[PARAM_Z2d].setValue(value) ){
		fclose(fptr);
		return 1;
	      }
	      //	      printf("zeta2  = [%s]\n", value );
	    }
	    else if( getValue(buf,"carg", value ) ){
	      if( !param[PARAM_CAd].setValue(value) ){
		fclose(fptr);
		return 1;
	      }
	      //	      printf("carg  = [%s]\n", value );
	    }
	    else if( getValue(buf,"cratio", value ) ){
	      if( !param[PARAM_CRd].setValue(value) ){
		fclose(fptr);
		return 1;
	      }
	      //	      printf("cratio = [%s]\n", value );
	    }
	    else if( getValue(buf,"repulsive", value ) ){
	      if( !param[PARAM_Rd].setValue(value) ){
		fclose(fptr);
		return 1;
	      }
	      //	      printf("repulsive = [%s]\n", value );
	    }
	    else{
	      fprintf(stderr,"broken line %s", buf);
	    }
	  }
	  continue;
	} // D

	fprintf(stderr,"broken line %s", buf);
      }
      continue;
    } // param

    fprintf(stderr,"broken line %s", buf);
  } //

  fclose(fptr);

#if 0
  FILE* fptr = fopen("../config.sh","rb");

  while( fgets(buf,sizeof(buf),fptr) ){
    if( buf[0] == '#' ) continue;
    if( 1 != sscanf(buf,"export %[^=]",name) &&
	1 != sscanf(buf," %[^ =]",name) ) continue;

    bool match=false;
    for( int p=PARAM_Ds; p<=PARAM_CRd; p++ ){
      if( strcmp( name, param[p].longname ) == 0 ){
	if( 2 == sscanf(buf,"%*[^=] = \" %lf %lf\"",
			&param[p].min, &param[p].max ) ){
	  match=true;
	  param[p].set = true;
	  break;
	}
	else{
	  fprintf(stderr,"Error: environment variable %s is broekn.\n",
	      param[p].longname );
	  fclose(fptr);
	  return 1;
	}
      }
    }
    if(!match){
      if( false );
      else if( strcmp( name, "weight_eigenerror" ) == 0 ){
	if( 1 == sscanf(buf,"%*[^=]=\"%lf\"", &weight_eigenerror ) ){
	  match=true;
	}
	else{
	  fprintf(stderr,"Error: environment variable %s is broekn.\n",
		  "weight_eigenerror" );
	  fclose(fptr);
	  return 1;
	}
      }
      else if( strcmp( name, "weight_totalerror" ) == 0 ){
	if( 1 == sscanf(buf,"%*[^=]=\"%lf\"", &weight_totalerror ) ){
	  match=true;
	}
	else{
	  fprintf(stderr,"Error: environment variable %s is broekn.\n",
		  "weight_totalerror" );
	  fclose(fptr);
	  return 1;
	}
      }
      else if( strcmp( name, "criterion_error" ) == 0 ){
	if( 1 == sscanf(buf,"%*[^=]=\"%lf\"", &criterion_error ) ){
	  match=true;
	}
	else{
	  fprintf(stderr,"Error: environment variable %s is broekn.\n",
		  "criterion_error" );
	  fclose(fptr);
	  return 1;
	}
      }
    }
  }
  fclose(fptr);
#endif


  if( param[PARAM_CRs].set && !param[PARAM_CAs].set ){
    param[PARAM_CAs].min = 180.0/M_PI * atan(param[PARAM_CRs].min);
    param[PARAM_CAs].max = 180.0/M_PI * atan(param[PARAM_CRs].max);
    param[PARAM_CAs].set = true;
  }
  if( param[PARAM_CRp].set && !param[PARAM_CAp].set ){
    param[PARAM_CAp].min = 180.0/M_PI * atan(param[PARAM_CRp].min);
    param[PARAM_CAp].max = 180.0/M_PI * atan(param[PARAM_CRp].max);
    param[PARAM_CAp].set = true;
  }
  if( param[PARAM_CRd].set && !param[PARAM_CAd].set ){
    param[PARAM_CAd].min = 180.0/M_PI * atan(param[PARAM_CRd].min);
    param[PARAM_CAd].max = 180.0/M_PI * atan(param[PARAM_CRd].max);
    param[PARAM_CAd].set = true;
  }

  for( int p=PARAM_Ds; p<=PARAM_Rd; p++ ){
    if( !param[p].set ){
      fprintf(stderr,"Error: variable %s is not defined.\n",
	      param[p].longname );
    }
    else{
      //      fprintf(stderr,"%s %f %f\n", param[p].longname, param[p].min, param[p].max );
    }
  }
  //  fprintf(stderr,"%s %f\n", "weight_eigenerror", weight_eigenerror );
  //  fprintf(stderr,"%s %f\n", "weight_totalerror", weight_totalerror );
  //  fprintf(stderr,"%s %e\n", "criterion_error", criterion_error );


  /*
  for( int p=0; p<DIM; p++ ){
    char* envname  = param_str[p].longname;
    char* envvalue = getenv( envname );
    if( envvalue == NULL ){
      if( strcmp( envname, "param_carg_s") == 0 ){
	envname  = "param_cratio_s";
	envvalue = getenv( envname );
      }
      if( strcmp( envname, "param_carg_p") == 0 ){
	envname  = "param_cratio_p";
	envvalue = getenv( envname );
      }
      if( strcmp( envname, "param_carg_d") == 0 ){
	envname  = "param_cratio_d";
	envvalue = getenv( envname );
      }

      if( envvalue == NULL ){
	fprintf(stderr,"Error: environment variable %s is not defined.\n", param_str[p].longname );
	return 1;
      }
    }

    if( 2 != sscanf(envvalue,"%lf %lf",&prange[p].min,&prange[p].max ) ){
      fprintf(stderr,"Error: environment variable %s is broekn [%s].\n", envname, envvalue );
      return 1;
    }

    if( strcmp( envname, "param_cratio_s") == 0 ||
	strcmp( envname, "param_cratio_p") == 0 ||
	strcmp( envname, "param_cratio_d") == 0 ){
      prange[p].min = 180.0/M_PI * atan(prange[p].min);
      prange[p].max = 180.0/M_PI * atan(prange[p].max);
    }
  }

  {
    char* env = getenv("weight_eigenerror");
    if( env == NULL || 1 != sscanf(env,"%lf",&weight_eigenerror ) ){
      fprintf(stderr,"Error: environment variable weight_eigenerror is not defined.\n");
      return 1;
    }
  }
  {
    char* env = getenv("weight_totalerror");
    if( env == NULL || 1 != sscanf(env,"%lf",&weight_totalerror ) ){
      fprintf(stderr,"Error: environment variable weight_totalerror is not defined.\n");
      return 1;
    }
  }
  {
    char* env = getenv("criterion_error");
    if( env == NULL || 1 != sscanf(env,"%lf",&criterion_error ) ){
      fprintf(stderr,"Error: environment variable criterion_error is not defined.\n");
      return 1;
    }
  }
  */


  return 0;
}

int init( vector<Parameter>& vP )
{

  // the 0th trial
  {
    Parameter Pt;
    for( int d=0; d<DIM; d++ ){
      Pt.p[d] = param[d].init;
      //      printf("Pt.p[%d]: %f \n",d,Pt.p[d]);
    }
    Pt.evaluate();
    vP.push_back(Pt);
  }

  // other trials
  for( int m=0; m<DIM; m++ ){
    // skip if the range is not available
    if( param[m].min == param[m].max ) continue;

    Parameter Pt;
    for( int d=0; d<DIM; d++ ){
      if( d==m ){
	Pt.p[d] = param[d].min;
      }
      else{
	Pt.p[d] = param[d].max;
      }
      //      printf("Pt.p[%d]: %f \n",d,Pt.p[d]);
    }
    Pt.evaluate();
    vP.push_back(Pt);
  }

  // return the size of vP
  return (int)vP.size();
}

void sort( vector<Parameter>& vP )
{
  std::sort( vP.begin(), vP.end() );
}

bool converged( const vector<Parameter>& vP, const double criterion_error )
{
  // if the best and the worst is almost the same, it has converged.
  return fabs(vP.back().err - vP.front().err) < criterion_error;
}

double calcVolume( const vector<Parameter>& vP )
{
  double volume = 1.0;
  double min,max;
  for( int d=0; d<DIM; d++ ){
    min=vP.front().p[d];
    max=vP.front().p[d];
    for( int m=0; m<(int)vP.size(); m++ ){
      if( min > vP[m].p[d] ) min = vP[m].p[d];
      if( max < vP[m].p[d] ) max = vP[m].p[d];
    }

    if( min == max ) continue;
    volume *= max-min;
  }
  return volume;
}

int main()
{
  vector<Parameter> vP;

  if( getEnv() ){
    fprintf( stderr, "Error: failed to start optimizing parameters.\n" );
    return 1;
  }

  printf("#");
  for( int p=0; p<DIM; p++ ){
    printf(" %s", param[p].shortname );
  }
  //  printf(", sum of errors\n");
  printf(", eval_f_all eval_f_ebd eval_f_etot\n");

  const int M = init(vP);

  do{
    // sort, then vP[0] is the worst, vP[M-1] is the best.
    sort(vP);

    // calc center except the worst
    Parameter Pg;
    for( int m=1; m<M; m++ ){
      Pg = Pg + vP[m];
    }
    Pg = Pg * (1.0/(M-1));

    // calc the first trial
    Parameter Pt1 = vP[0] + (Pg - vP[0])*2.0;
    Pt1.evaluate();

    // if the first trial is better than the second worst
    if( Pt1.err < vP[1].err ){
      // replace the worst by the first trial
      vP[0] = Pt1;
    }
    // if the first trial is worse than the second worst
    else{
      // calc the second trial
      Parameter Pt2 = vP[0] + (Pg - vP[0])*0.5;
      Pt2.evaluate();

      // if the second trial is better than the worst
      if( Pt2.err < vP[0].err ){
	// replace the worst by the second trial
	vP[0] = Pt2;
      }
      // if the second trial is worse than the worst
      else{
	// shrink except the best
	for( int m=0; m+1<M; m++ ){
	  vP[m] = (vP[m] + vP[M-1])*0.5;
	  vP[m].evaluate();
	}
      }
    }

    //    printf("# Vol = %e\n", calcVolume(vP) );
  } while( !converged(vP,criterion_error) );

  printf("# Info: converged\n");

  return 0;
}
