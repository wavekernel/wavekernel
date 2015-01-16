#include <stdio.h>
#include <math.h>

main( int argc, char* argv[] )
{
  if( argc != 5 && argc != 6 ){
    fprintf(stderr,"# usage typeA: give angle of (c1,c2)\n");
    fprintf(stderr,"  %s zeta1 zeta2 n cang\n", argv[0] );
    fprintf(stderr,"# usage typeB: give (c1,c2) directoly\n");
    fprintf(stderr,"  %s zeta1 zeta2 n c1 c2\n", argv[0] );
    return 1;
  }

  double z1, z2, cang, c1, c2;
  int n;

  if( 1 != sscanf( argv[1], "%lf", &z1 ) ){ // zeta1
    printf("%f %f\n", 0.0, 0.0 );
    return 0;
  }
  if( 1 != sscanf( argv[2], "%lf", &z2 ) ){ // zeta1
    printf("%f %f\n", 0.0, 0.0 );
    return 0;
  }

  if( z2 == 0.0 ){
    printf("%f %f\n", 0.0, 0.0 );
    return 0;
  }

  if( 1 != sscanf( argv[3], "%d",  &n ) ){ // n
    printf("%f %f\n", 0.0, 0.0 );
    return 0;
  }

  double am = (z1+z2)/2;
  double gm = sqrt(z1*z2);
  double gamma = pow( gm/am, 2*n+1 );

  if( argc == 5 && 1 == sscanf( argv[4], "%lf", &cang ) ){
    double coef = 1.0/sqrt(1.0+gamma*sin(2*cang*M_PI/180.0));
    c1 = coef * cos(cang*M_PI/180.0); 
    c2 = coef * sin(cang*M_PI/180.0);

    //  double f  = c1*c1 + c2*c2 + 2*gamma*c1*c2;

    printf("%f %f\n", c1, c2 );
    //  printf("%f\n", f );
  }
  else if( argc == 6 && 2 == sscanf( argv[4], "%lf", &c1, &c2 ) ){
    if( c1 == 0.0 ){
      c2 = 1.0;
    }
    else{
      double alpha = c2/c1;
      double coef = 1.0/sqrt(1.0 + alpha*alpha + 2*gamma*alpha);

      c1 = coef * 1.0;
      c2 = coef * alpha;
    }
    //  double f  = c1*c1 + c2*c2 + 2*gamma*c1*c2;
    printf("%f %f\n", c1, c2 );
    //  printf("%f\n", f );
  }
  else{
    printf("%f %f\n", 0.0, 0.0 );
    return 0;
  }

  return 0;
}
