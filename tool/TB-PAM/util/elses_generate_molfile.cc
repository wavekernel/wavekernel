#include <stdio.h>
#include <string.h>
#include <math.h>

namespace ELSES
{
  void outputHeader( const int ion, const bool spin ){
    printf( "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" );
    printf( "<structure name=\"test\">\n" );

    const double Lcell = 10.00; // ang
    printf( "<unitcell>\n"
	    "  <vector unit=\"ang\">  %11.8f  0.00  0.00 </vector>\n"
	    "  <vector unit=\"ang\">   0.00 %11.8f  0.00 </vector>\n"
	    "  <vector unit=\"ang\">   0.00  0.00 %11.8f </vector>\n"
	    "</unitcell>\n", Lcell, Lcell, Lcell );
  }
  void outputFooter( void ){
    printf( "</structure>\n" );
  }
  void outputAtom( const char* elem, const double x, const double y, const double z ){
    printf( "<atom element=\"%s\">\n"
	    "  <position unit=\"ang\"> %11.8f %11.8f %11.8f </position>\n"
	    "</atom>\n", elem, x, y, z );
  }

  bool output( const char* mol, const char* basetype, const double bond, const double angle, const bool spin ){
    outputHeader(0,spin);
    if( false ){
    }
    else if( !strcmp( mol, "H2O" ) ){
      const double x = bond*cos(angle*0.5*M_PI/180.0);
      const double y = bond*sin(angle*0.5*M_PI/180.0);

      outputAtom( "H",   x,  +y, 0.0 );
      outputAtom( "H",   x,  -y, 0.0 );
      outputAtom( "O", 0.0, 0.0, 0.0 );
    }
    else if( !strcmp( mol, "CO2" ) ){
      const double x = bond;

      outputAtom( "C", 0.0, 0.0, 0.0 );
      outputAtom( "O",  +x, 0.0, 0.0 );
      outputAtom( "O",  -x, 0.0, 0.0 );
    }
    else if( !strcmp( mol, "NH3" ) ){
      outputAtom( "N", 0.0, 0.0, 0.0 );
      outputAtom( "H", bond, 0.0, 0.0 );
      outputAtom( "H",(-0.5)*bond, (+0.5*sqrt(3.0))*bond, 0.0 );
      outputAtom( "H",(-0.5)*bond, (-0.5*sqrt(3.0))*bond, 0.0 );
    }
    else if( !strcmp( mol, "CH4" ) ){
      const double Lb = sqrt(3.0)/3.0*bond;
      const double Lc = sqrt(6.0)/3.0*bond;

      outputAtom( "C", 0.0, 0.0, 0.0 );
      outputAtom( "H", -Lc, 0.0, -Lb );
      outputAtom( "H", +Lc, 0.0, -Lb );
      outputAtom( "H", 0.0, -Lc, +Lb );
      outputAtom( "H", 0.0, +Lc, +Lb );
    }
    else if( !strcmp( mol, "PS4_3-" ) ){
      const double Lb = sqrt(3.0)/3.0*bond;
      const double Lc = sqrt(6.0)/3.0*bond;

      outputAtom( "P", 0.0, 0.0, 0.0 );
      outputAtom( "S", -Lc, 0.0, -Lb );
      outputAtom( "S", +Lc, 0.0, -Lb );
      outputAtom( "S", 0.0, -Lc, +Lb );
      outputAtom( "S", 0.0, +Lc, +Lb );
    }
    else if( !strcmp( mol, "NO2" ) ){
      const double x = bond*cos(angle*0.5*M_PI/180.0);
      const double y = bond*sin(angle*0.5*M_PI/180.0);

      outputAtom( "N", 0.0, 0.0, 0.0 );
      outputAtom( "O",   x,  +y, 0.0 );
      outputAtom( "O",   x,  -y, 0.0 );
    }
    else if( !strcmp( mol, "SO2" ) ){
      const double x = bond*cos(angle*0.5*M_PI/180.0);
      const double y = bond*sin(angle*0.5*M_PI/180.0);

      outputAtom( "S", 0.0, 0.0, 0.0 );
      outputAtom( "O",   x,  +y, 0.0 );
      outputAtom( "O",   x,  -y, 0.0 );
    }
    else if( !strcmp( mol, "CF4" ) ){
      const double Lb = sqrt(3.0)/3.0*bond;
      const double Lc = sqrt(6.0)/3.0*bond;

      outputAtom( "C", 0.0, 0.0, 0.0 );
      outputAtom( "F", -Lc, 0.0, -Lb );
      outputAtom( "F", +Lc, 0.0, -Lb );
      outputAtom( "F", 0.0, -Lc, +Lb );
      outputAtom( "F", 0.0, +Lc, +Lb );
    }
    else if( !strcmp( mol, "C6H6" ) ){
      const double rC=bond;
      const double rCH=angle;
      const double rH=rC+rCH;
      int n;
      for(n=1;n<7;n++){
	outputAtom( "C", rC*cos(2*M_PI*n/6.0), rC*sin(2*M_PI*n/6.0), 0.0 );
	outputAtom( "H", rH*cos(2*M_PI*n/6.0), rH*sin(2*M_PI*n/6.0), 0.0 );
      }
    }
    else if( !strcmp( mol, "C2H6" ) ){
      const double rCC=bond;
      const double diangle=2*M_PI*angle/360.0;
      const double rCH=1.09;
      const double phi=2.0*M_PI*109.0/360.0;
      const double theta=2*M_PI/3.0;
      int n;

      outputAtom("C", 0.0, 0.0, 0.0);
      outputAtom("C", 0.0, 0.0, rCC);

      for(n=0;n<3;n++){
	outputAtom("H", rCH*(cos(M_PI/2.0-phi))*cos(n*theta),
		   rCH*(cos(M_PI/2.0-phi))*sin(n*theta),rCH*sin(M_PI/2.0-phi));
      }

      for(n=0;n<3;n++){
	outputAtom("H", rCH*(cos(M_PI/2.0-phi))*cos(n*theta+diangle),
		   rCH*(cos(M_PI/2.0-phi))*sin(n*theta+diangle),rCC-rCH*sin(M_PI/2.0-phi));
      }
    }
    else{
      fprintf(stderr,"Error: uknown molecule %s\n", mol );
      return false;
    }
    outputFooter();

    return true;
  }
}

namespace GAUSS
{
  void outputHeader( const char* basetype, const int ion, const bool spin ){
    if( spin ){
      printf( "#p Sp, ROB3LYP/%s formcheck pop=full punch=archive scf=Tight nosymm\n"
	      "\n"
	      "comment\n"
	      "\n"
	      "%d 2\n", basetype, ion );
    }
    else{
      printf( "#p Sp, B3LYP/%s formcheck pop=full punch=archive scf=Tight nosymm\n"
	      "\n"
	      "comment\n"
	      "\n"
	      "%d 1\n", basetype, ion );
    }
  }
  void outputFooter( void )
  {
    printf("\n\n");
  }
  void outputAtom( const char* elem, const double x, const double y, const double z )
  {
    printf(" %2s %+11.8f %+11.8f %+11.8f\n", elem, x, y, z );
  }

  bool output( const char* mol, const char* basetype, const double bond, const double angle, const bool spin ){
    if( false ){
    }
    else if( !strcmp( mol, "H2O" ) ){
      const double x = bond*cos(angle*0.5*M_PI/180.0);
      const double y = bond*sin(angle*0.5*M_PI/180.0);

      outputHeader(basetype,0,spin);
      outputAtom( "H",   x,  +y, 0.0 );
      outputAtom( "H",   x,  -y, 0.0 );
      outputAtom( "O", 0.0, 0.0, 0.0 );
    }
    else if( !strcmp( mol, "CO2" ) ){
      const double x = bond;

      outputHeader(basetype,0,spin);
      outputAtom( "C", 0.0, 0.0, 0.0 );
      outputAtom( "O",  +x, 0.0, 0.0 );
      outputAtom( "O",  -x, 0.0, 0.0 );
    }
    else if( !strcmp( mol, "NH3" ) ){
      outputHeader(basetype,0,spin);
      outputAtom( "N", 0.0, 0.0, 0.0 );
      outputAtom( "H", bond, 0.0, 0.0 );
      outputAtom( "H",(-0.5)*bond, (+0.5*sqrt(3.0))*bond, 0.0 );
      outputAtom( "H",(-0.5)*bond, (-0.5*sqrt(3.0))*bond, 0.0 );
    }
    else if( !strcmp( mol, "CH4" ) ){
      const double Lb = sqrt(3.0)/3.0*bond;
      const double Lc = sqrt(6.0)/3.0*bond;

      outputHeader(basetype,0,spin);
      outputAtom( "C", 0.0, 0.0, 0.0 );
      outputAtom( "H", -Lc, 0.0, -Lb );
      outputAtom( "H", +Lc, 0.0, -Lb );
      outputAtom( "H", 0.0, -Lc, +Lb );
      outputAtom( "H", 0.0, +Lc, +Lb );
    }
    else if( !strcmp( mol, "PS4_3-" ) ){
      const double Lb = sqrt(3.0)/3.0*bond;
      const double Lc = sqrt(6.0)/3.0*bond;

      outputHeader(basetype,-3,spin);
      outputAtom( "P", 0.0, 0.0, 0.0 );
      outputAtom( "S", -Lc, 0.0, -Lb );
      outputAtom( "S", +Lc, 0.0, -Lb );
      outputAtom( "S", 0.0, -Lc, +Lb );
      outputAtom( "S", 0.0, +Lc, +Lb );
    }
    else if( !strcmp( mol, "NO2" ) ){
      const double x = bond*cos(angle*0.5*M_PI/180.0);
      const double y = bond*sin(angle*0.5*M_PI/180.0);

      outputHeader(basetype,0,spin);
      outputAtom( "N", 0.0, 0.0, 0.0 );
      outputAtom( "O",   x,  +y, 0.0 );
      outputAtom( "O",   x,  -y, 0.0 );
    }
    else if( !strcmp( mol, "SO2" ) ){
      const double x = bond*cos(angle*0.5*M_PI/180.0);
      const double y = bond*sin(angle*0.5*M_PI/180.0);

      outputHeader(basetype,0,spin);
      outputAtom( "S", 0.0, 0.0, 0.0 );
      outputAtom( "O",   x,  +y, 0.0 );
      outputAtom( "O",   x,  -y, 0.0 );
    }
    else if( !strcmp( mol, "CF4" ) ){
      const double Lb = sqrt(3.0)/3.0*bond;
      const double Lc = sqrt(6.0)/3.0*bond;

      outputHeader(basetype,0,spin);
      outputAtom( "C", 0.0, 0.0, 0.0 );
      outputAtom( "F", -Lc, 0.0, -Lb );
      outputAtom( "F", +Lc, 0.0, -Lb );
      outputAtom( "F", 0.0, -Lc, +Lb );
      outputAtom( "F", 0.0, +Lc, +Lb );
    }
    else if( !strcmp( mol, "C6H6" ) ){
      const double rC=bond;
      const double rCH=angle;
      const double rH=rC+rCH;
      int n;

      outputHeader(basetype,0,spin);
      for(n=1;n<7;n++){
	outputAtom( "C", rC*cos(2*M_PI*n/6.0), rC*sin(2*M_PI*n/6.0), 0.0 );
	outputAtom( "H", rH*cos(2*M_PI*n/6.0), rH*sin(2*M_PI*n/6.0), 0.0 );
      }
    }
    else if( !strcmp( mol, "C2H6" ) ){
      const double rCC=bond;
      const double diangle=2*M_PI*angle/360.0;
      const double rCH=1.09;
      const double phi=2.0*M_PI*109.0/360.0;
      const double theta=2*M_PI/3.0;
      int n;

      outputHeader(basetype,0,spin);

      outputAtom("C", 0.0, 0.0, 0.0);
      outputAtom("C", 0.0, 0.0, rCC);

      for(n=0;n<3;n++){
	outputAtom("H", rCH*(cos(M_PI/2.0-phi))*cos(n*theta),
		   rCH*(cos(M_PI/2.0-phi))*sin(n*theta),rCH*sin(M_PI/2.0-phi));
      }

      for(n=0;n<3;n++){
	outputAtom("H", rCH*(cos(M_PI/2.0-phi))*cos(n*theta+diangle),
		   rCH*(cos(M_PI/2.0-phi))*sin(n*theta+diangle),rCC-rCH*sin(M_PI/2.0-phi));
      }
    }
    else{
      fprintf(stderr,"Error: uknown molecule %s\n", mol );
      return false;
    }

    outputFooter();

    return true;
  }
}

int main( const int argc, const char* argv[] )
{
  double bond, angle;
  bool spin;

  if( argc != 7 ){
    printf("%s molecule_name [elses|gaussian] basetype bond angle [nospin|spin]\n", argv[0] );
    printf(" Basetype is only effective in gaussian.\n");
    printf(" acceptable: STO-3G,3-21G,3-21G*,6-311G,6311G*,6-311G**,6-311+G**.\n");
    return 1;
  }

  if( 1 != sscanf( argv[4], "%lf", &bond ) ){
    printf("bond(arg3) is not given\n");
    return 1;
  }
  if( 1 != sscanf( argv[5], "%lf", &angle ) ){
    printf("angle(arg4) is not given\n");
    return 1;
  }

  if( strcmp( argv[6], "spin" ) == 0 ){
    spin = true;
  }
  else if( strcmp( argv[6], "nospin" ) == 0 ){
    spin = false;
  }
  else{
    printf("spin(arg6) is not given\n");
    return 1;
  }


  if( false ){
  }
  else if( !strcmp( argv[2], "elses" ) ){
    if( !ELSES::output( argv[1], argv[3], bond, angle, spin ) ){
      return 1;
    }
  }
  else if( !strcmp( argv[2], "gaussian" ) ){
    if( !GAUSS::output( argv[1], argv[3], bond, angle, spin ) ){
      return 1;
    }
  }
  else{
    printf("arg2 is not given\n");
    return 1;
  }

  return 0;
}
