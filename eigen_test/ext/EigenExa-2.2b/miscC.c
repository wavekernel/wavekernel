#include <stdio.h>
#include <stdlib.h>

#if __IBM_REGISTER_VARS
#	define	GET_DELTA	get_delta
#	define	REDIST1		dc_redist1
#	define	REDIST1C	dc_redist1c
#	define	REDIST2		dc_redist2
#	define	REDIST2C	dc_redist2c
#else
#	define	GET_DELTA	get_delta_
#	define	REDIST1		dc_redist1_
#	define	REDIST1C	dc_redist1c_
#	define	REDIST2		dc_redist2_
#	define	REDIST2C	dc_redist2c_
#endif


void
GET_DELTA(void *a1, void *a2, int *i)
{
	long	b1 = (long)a1;
	long	b2 = (long)a2;
	long	delta;

	delta = b1 - b2;
	*i = (int)delta;
}


void REDIST1(
	int *n, int *NB,
	double *z, double *work, int *lddz, double *iwork, int *liwork );
void REDIST2(
	int *n, int *NB,
	double *work, int *lddz, double *z, int *ldz, double *iwork, int *liwork );


void REDIST1C(
	int *n, int *NB,
	double *z, double *work, int *lddz, int *iwork, int *liwork )
{
	double	*dwork = (double *)iwork;
	REDIST1( n, NB, z, work, lddz, dwork, liwork );
}

void REDIST2C(
	int *n, int *NB,
	double *work, int *lddz, double *z, int *ldz, int *iwork, int *liwork )
{
	double	*dwork = (double *)iwork;
	REDIST2( n, NB, work, lddz, z, ldz, dwork, liwork );

}

#if __IBM_REGISTER_VARS
void flush( int *unit ) {
	fflush( stdout );
	flush_( unit );
	sleep( 0 );
}
#endif

