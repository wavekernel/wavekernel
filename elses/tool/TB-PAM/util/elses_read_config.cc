#include <stdio.h>
#include <string.h>

bool isComment( char* buf )
{
  char dummy[1024];

  return buf[0] == '#' || sscanf(buf,"%s", dummy ) == -1;
}

bool getValue( const char* buf, const char* name, char* value )
{
  char format[1024];

  sprintf( format, " %s%%*[ =] (%%[^)]", name );

  if( 1 == sscanf(buf,format, value ) ) return true;

  sprintf( format, " %s%%*[ =] %%s", name );
  if( 1 == sscanf(buf,format, value ) ) return true;

  return false;
}

int main( int argc, char* argv[] )
{
  char buf[1024];
  char tag[256];
  char value[256];

  char* fname = argv[1];
  FILE* fptr = fopen( fname, "rb");
  if( fptr == NULL ){
    return 1;
  }

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
	      printf("element=\"%s\"\n", value );
	    }
	    else if( getValue(buf,"molecule", value ) ){
	      printf("molecule=\"%s\"\n", value );
	    }
	    else if( getValue(buf,"spin", value ) ){
	      printf("spin=\"%s\"\n", value );
	    }
	    else if( getValue(buf,"basetype", value ) ){
	      printf("basetype=\"%s\"\n", value );
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
	      printf("bond_optimized=\"%s\"\n", value );
	    }
	    else if( getValue(buf,"angle_optimized", value ) ){
	      printf("angle_optimized=\"%s\"\n", value );
	    }
	    else if( getValue(buf,"bonds_conformation", value ) ){
	      printf("bonds_conformation=\"%s\"\n", value );
	    }
	    else if( getValue(buf,"angles_conformation", value ) ){
	      printf("angles_conformation=\"%s\"\n", value );
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
	    else if( getValue(buf,"nbands", value ) ){
	      	      printf("nbands=\"%s\"\n", value );
	    }
	    else if( getValue(buf,"weight_eigen", value ) ){
	      //	      printf("weight_eigen  = [%s]\n", value );
	    }
	    else if( getValue(buf,"weight_total", value ) ){
	      //	      printf("weight_total  = [%s]\n", value );
	    }
	    else if( getValue(buf,"convergence", value ) ){
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
	      printf("doublezeta_n_s=\"%s\"\n", value );
	    }
	    else if( getValue(buf,"diag", value ) ){
	      //	      printf("diag  = [%s]\n", value );
	    }
	    else if( getValue(buf,"zeta1", value ) ){
	      //	      printf("zeta1  = [%s]\n", value );
	    }
	    else if( getValue(buf,"zeta2", value ) ){
	      //	      printf("zeta2  = [%s]\n", value );
	    }
	    else if( getValue(buf,"carg", value ) ){
	      //	      printf("carg  = [%s]\n", value );
	    }
	    else if( getValue(buf,"cratio", value ) ){
	      //	      printf("cratio = [%s]\n", value );
	    }
	    else if( getValue(buf,"repulsive", value ) ){
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
	      printf("doublezeta_n_p=\"%s\"\n", value );
	    }
	    else if( getValue(buf,"diag", value ) ){
	      //	      printf("diag  = [%s]\n", value );
	    }
	    else if( getValue(buf,"zeta1", value ) ){
	      //	      printf("zeta1  = [%s]\n", value );
	    }
	    else if( getValue(buf,"zeta2", value ) ){
	      //	      printf("zeta2  = [%s]\n", value );
	    }
	    else if( getValue(buf,"carg", value ) ){
	      //	      printf("carg  = [%s]\n", value );
	    }
	    else if( getValue(buf,"cratio", value ) ){
	      //	      printf("cratio = [%s]\n", value );
	    }
	    else if( getValue(buf,"repulsive", value ) ){
	      //	      printf("repulsive = [%s]\n", value );
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
	      printf("doublezeta_n_d=\"%s\"\n", value );
	    }
	    else if( getValue(buf,"diag", value ) ){
	      //	      printf("diag  = [%s]\n", value );
	    }
	    else if( getValue(buf,"zeta1", value ) ){
	      //	      printf("zeta1  = [%s]\n", value );
	    }
	    else if( getValue(buf,"zeta2", value ) ){
	      //	      printf("zeta2  = [%s]\n", value );
	    }
	    else if( getValue(buf,"carg", value ) ){
	      //	      printf("carg  = [%s]\n", value );
	    }
	    else if( getValue(buf,"cratio", value ) ){
	      //	      printf("cratio = [%s]\n", value );
	    }
	    else if( getValue(buf,"repulsive", value ) ){
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

  return 0;
}
