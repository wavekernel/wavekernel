#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main( const int argc, const char* argv[] )
{
  if( argc != 4 ){
    printf("%s fileA fileB nbands\n", argv[0] );
    return 1;
  }

  FILE* fptr1 = fopen(argv[1], "r");
  FILE* fptr2 = fopen(argv[2], "r");

  int nbands;  
  nbands=atoi(argv[3]);
  //  printf("%d \n",nbands);

  if( fptr1 == NULL ){
    printf("file %s not found\n", argv[1]);
    return 1;
  }
  if( fptr2 == NULL ){
    printf("file %s not found\n", argv[2]);
    return 1;
  }

  double value1, value2;
  int status1, status2;

  int n=0; // number of configuration
  double diff=0.0;
  char buf1[1024], buf2[1024];
  while( fgets(buf1,sizeof(buf1),fptr1) ){
    fgets(buf2,sizeof(buf2),fptr2);

    if( buf1[0] == '#' || buf1[0] == '\n' || buf1[0] == '\r' ) continue;
    // modified by MIZUHO-IR 2013/10/22

    char* ptr1 = buf1;
    char* ptr2 = buf2;
    int offset1;
    int offset2;
    int i_level=0;  // index of level
    double value_test;

    while( 1 == sscanf( ptr1,"%lf%n",&value1,&offset1) ){
    //while(i_level < nbands){
    //  if( 1 == sscanf( ptr1,"%lf%n",&value1,&offset1)) exit;
      ptr1 += offset1;
      sscanf( ptr2,"%lf%n",&value2,&offset2);
      ptr2 += offset2;

      //     printf("value1:% f  value2:% f\n", value1,value2);

      //      value_test= fabs(value1 - value2);
      //      printf("fabs(value1-value2)=%f\n",value_test);

      //      if( value1 != value2 ){
      if( fabs(value1 - value2) > 1e-8 && i_level < nbands){
	diff += (value1-value2)*(value1-value2);
	//printf("value1:% f  lavlu2:% f diff:% f \n",value1,value2,diff);
	//	n++; // removed by MIZUHO-IR 2013/10/22
	i_level+=1;
      }
    }
    n++; // add by MIZUHO-IR 2013/10/22
  }
  diff = sqrt(diff/n);

  //printf("Number of configuration: %d \n",n);
  printf("%f\n", diff );

  return 0;
}
