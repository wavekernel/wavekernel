/TRD-BLK /{TRD=$3}
/D\&C/{DC=$2;if(NF>3)DC=$3}
/TRDBAK/{BK=$3}
/Total/{print TRD", "DC", "BK", "$2}
