#! /usr/bin/perl
# -*- mode: Perl -*-
#================================================================
# ELSES version 0.02
# Copyright (C) ELSES. 2007-2010 all rights reserved
#================================================================

use strict;
use warnings;
use diagnostics;
use POSIX;
use XML::LibXML;

my $bandCalcInput="ForBandCalc.txt";
my $tolerance = 1.0e-12;
my $twigState = "";
my (@unitVecs, @transVecs, @positions, @elements, @classes, @motions, @trans_unitCoeff );

if( $#ARGV < 0 || $#ARGV > 2 ) {
    die <<EOT;
    
  Usage:
    mkSupercell.pl < -u > unitcell.xml supercell.xml
      Result:
	1. read an unitcell size from <unitcell> tag in unitcell.xml
	2. read a supercell size from <unitcell> tag in supercell.xml
	3. check compatibility between the unitcell and the supercell
	4. write supercell structure to supercell_new.xml
	(in order to avoid overwriting, _new is added to the basename of supercell.xml)
    5. Add "-u" option to skip pretty-format part.
        If you don't have twig-perl module, try "-u".
        If you feel difficulty to handle large system, try "-u".

EOT
    }
  if( $ARGV[0] ne "-u" ){
	$twigState="necessary";
  }else{
	shift @ARGV;
  }
my($unitCellFile,$outputFile)=@ARGV;
# foreach(@ARGV){
#   print "$_\n";
# }
# die "debugging\n";

my $parser = XML::LibXML->new();

foreach my $file ( $unitCellFile, $outputFile ) {
    if( $file !~ /\.xml$/ ) {
	die "Given filename $file should end with extention \".xml\"\n";
    }
}

sub toAU($$){
    my $Bohr = 0.529177e0;  #my $Bohr = 5.291772108e-1;
    my ($unit,$vr) = @_; # $vr = (reference of vector)
    if($unit eq 'angstrom'){
	  foreach(@$vr) {
		# print "$_\n"; # for debug
	    $_ /= $Bohr;
	  }
	  # change unit -> a.u. # Not implemented yet
    }elsif($unit ne 'a.u.'){
	  die "\"unit\" attribute must be 'a.u.' or 'angstrom'.\n";
    }
}

sub getUnitVecs($$){ # read unit vector & convert to Atomic Unit
    my $parser = XML::LibXML->new();
    my ($file,$vr) = @_; # $vr = (reference of vector)
    my $doc    = $parser->parse_file($file)
	or die "can't parse \"$file\"\n";
    #my $root   = $doc->documentElement(); #DOM
    #print $doc->toString(1);
    
    my @unit_vec=
	$doc->getElementsByTagName('unitcell')->
	shift->getElementsByTagName('vector');
    
    if($#unit_vec != 2){
	die "Wrong number of unit_vecs\n";
    }
    
    my (@vectors, $unit);
    for ( my $i=0 ; $i < 3 ; $i++ ) {
	$unit=$unit_vec[$i]->getAttribute('unit');
	$unit_vec[$i]->string_value =~ /^\s+(\S+)\s+(\S+)\s+(\S+)/ ;
	$unit_vec[$i] = [$1,$2,$3];
	&toAU($unit,$unit_vec[$i]);
	@$vr=@unit_vec;
    }
}

sub outProd ($) {
    my ($vr) = @_ ; # $vr = (reference of vector)
    return [
	    $vr->[0][1] * $vr->[1][2] - $vr->[0][2] * $vr->[1][1],
	    $vr->[0][2] * $vr->[1][0] - $vr->[0][0] * $vr->[1][2],
	    $vr->[0][0] * $vr->[1][1] - $vr->[0][1] * $vr->[1][0]
	    ] ;
}

sub inProd ($$){
    my ($vr1, $vr2) = @_ ; # $vr = (reference of vector)
    my $sum = 0;
    for(my $i=0; $i < 3; $i++){
	$sum += $vr1->[$i] * $vr2->[$i];
    }
    return $sum;
}

sub det3 ($) {
    my ($vr) =  @_ ; # $vr = (reference of vector)
    return &inProd(&outProd($vr),$vr->[2]);
}

sub norm2Vec($){
    my ($vr) =  @_ ; # $vr = (reference of vector)
    my $norm2 = 0.0e0;
    for(my $i=0 ; $i < 3 ; $i++){
	$norm2 += $vr->[$i] * $vr->[$i];
    }
    return $norm2;
}

sub normVec($){
    my ($vr) =  @_ ; # $vr = (reference of vector)
    return sqrt(&norm2Vec($vr));
}

sub linearCombination($$){
    my ($cr,$vr)=@_;
 # $cr = (reference of coefficients),  $vr = (reference of vectors or matrix)
    my @tv;
    for(my $i=0; $i < 3; $i++) {
	$tv[$i]=$cr->[0]*$vr->[0][$i]+$cr->[1]*$vr->[1][$i]+$cr->[2]*$vr->[2][$i];
    }
    @$cr=@tv;
}

sub solveCoeff($$){
    my($vecRef,$unitVecsRef)=@_;
    my $volUnit=&det3($unitVecsRef);
    my @ansVec;
    for(my $j=0; $j < 3; $j++) {
	my @tmpVecs = (map [@$_] , @$unitVecsRef ); # deep copy 2-D array
	@{$tmpVecs[$j]}=@$vecRef;
	my $volTmp = &det3(\@tmpVecs);
	$ansVec[$j] = $volTmp/$volUnit;
    }
    @$vecRef = @ansVec;
}

sub printVec($){
    my ($vr)=@_;
    printf "( %10e, %10e, %10e )\n", @$vr;
}

sub print3Vecs($){
    my ($vr)=@_;
    for(my $i = 0 ; $i < 3 ; $i++){
	&printVec($vr->[$i])
	}
}

sub nint($){ # nint means nearest integer
    my ($realNum) = @_ ;
    my $int  = floor($realNum);
    my $frac = $realNum - $int;
    return (( $frac > 0.5 ) ? $int + 1 : $int ) ;
}

sub gcd{     # greatest common diviser
    my ($a,$b) = @_;
    #print "$a, $b\n";
    $a=abs($a);
    $b=abs($b);
    ($a,$b) = ( $a > $b ) ? ($b,$a) : ($a,$b) ;
    if ($a == 0) {
	return $b;
    }
    my $m = $b % $a;
    if ($m == 0){
	return $a;
    }else{
	return &gcd($a,$m)
	}
}

sub gcd3 {
    my $m = &gcd(@_);
    return &gcd($m,$_[2]);
}

sub lcm {    # least common multiple
    my($a,$b) = @_;
    return abs($a*$b)/&gcd($a,$b);
}

sub lcm3 {
    my ($a,$b,$c) = @_;
    return &lcm(&lcm($a,$b),$c);
}

sub equivalent($$){
    my($vr1,$vr2)=@_; # vr? : reference of vector
    for(my $i=0; $i < 3; $i++){
	if(abs($vr2->[$i] - $vr1->[$i] - &nint($vr2->[$i] - $vr1->[$i])) > $tolerance ){
	    return (0!=0);
	}
    }
    return (0==0);
}

sub onRectangularLattice($) {
    my($vr)=@_; # vr : reference of vector
    return &equivalent($vr,[0,0,0]) ;
}

&getUnitVecs($unitCellFile,\@unitVecs);
&getUnitVecs($outputFile,\@transVecs);

my $volUnit  = &det3(\@unitVecs);
my $volTrans = &det3(\@transVecs);
my $numUnit = &nint($volTrans/$volUnit);
if ( $numUnit == 0 ) {
    die "The volume of the supercell is less than that of the unitcell.";
}

print "translation vectors of unitcell (primitive translation vectors)\n";
&print3Vecs(\@unitVecs);
print "translation vectors of supercell\n";
&print3Vecs(\@transVecs);

print "lengths of translation vectors of supercell (in atomic unit)\n";
printf " %20.17e, %20.17e, %20.17e \n",
  &normVec($transVecs[0]),
  &normVec($transVecs[1]),
  &normVec($transVecs[2]);

printf "VOL(supercell)=%10e, VOL(unitcell)=%10e\n",$volTrans,$volUnit ;
if ( abs($volTrans/$volUnit/$numUnit-1.0e0 ) > $tolerance ){
    die "Super cell does not consistent with unit cell.\n";
}else{
    printf "VOL(supercell)/VOL(unitcell)=%10e\n",
    &nint($volTrans/$volUnit) ;
}

# check if two lattices are commensurate or not
for(my $i=0; $i < 3; $i++){
    my @vec = @{$transVecs[$i]};
    &solveCoeff(\@vec,\@unitVecs);
    
    for(my $j=0; $j < 3; $j++){
	my $inum = &nint($vec[$j]);
	if ($inum == 0) {
	    if( abs($vec[$j]) > $tolerance ){
		printf "\n %10e \n", $vec[$j];
		die "Incommensurate lattice (1)\n";
	    }
	}elsif(abs($vec[$j]/$inum-1.0e0) > $tolerance){
	    printf "\n %10e \n", $vec[$j];
	    die "Incommensurate lattice (2)\n";
	}
	$trans_unitCoeff[$i][$j]=$inum;
    }
    printf
	"(Translation vector %2d) = (%2d, %2d, %2d) in the unit of primitive translation vectors\n",
	$i+1,  $trans_unitCoeff[$i][0],$trans_unitCoeff[$i][1],$trans_unitCoeff[$i][2];
}

my $doc = $parser->parse_file($unitCellFile)
    or die "can't parse \"$unitCellFile\"\n";

print "Structure in unitcell\n";
my ($nodeList, $node, $unit, $atomNode);
foreach $atomNode ($doc->getElementsByTagName('atom')) {
    $nodeList = $atomNode->getElementsByTagName('position');
    if( $nodeList->size >1 ) {
	die "Too many posisions in ",$atomNode->getAttribute('element'),"\n";
    }
    $node = $nodeList->get_node(0);
    $unit = $node->getAttribute('unit');
    my @tv;
    $node->string_value =~ /^\s*(\S+)\s+(\S+)\s+(\S+)/;
    eval {
	local $SIG{__WARN__} =
	    sub {
		print STDERR "atom position\n";
		print STDERR "$1,$2,$3\n";
		die "Numeric format error\n";
	    };
	my $string = sprintf "Read %10e,%10e,%10e\n", $1, $2, $3;
	# checking numeric format
    };
    if($@) { die $@; }
    @tv=($1,$2,$3);
    if($unit eq 'angstrom') {
	&toAU($unit,\@tv);
	$unit='a.u.';
    }
    if($unit eq 'a.u.') {
	&solveCoeff(\@tv,\@unitVecs);
    }
    push @positions, \@tv;
    push @elements, $atomNode->getAttribute('element');
    push @classes,  $atomNode->getAttribute('class');
    push @motions,  $atomNode->getAttribute('motion');
    printf "Position%6d = ( %10e, %10e, %10e ) : %6s\n",
    $#positions+1, @{$positions[$#positions]}, $elements[$#positions];
}

# analyze the arrangement of unitcells
my @partition;
for(my $i=0; $i < 3; $i++){
    $partition[$i]=&gcd3(@{$trans_unitCoeff[$i]});
}
print "We assume supercell is rectangular.\n";
printf
    "The minimal rectangular cell consists of %4d unitcells.\n",
    $numUnit/($partition[0] * $partition[1] * $partition[2]);
printf
    "Given supercell consists of %4d x %4d x %4d minimal rectangular cells\n",
    @partition;

my @rectVecs;
for(my $i=0; $i < 3; $i++){
    for(my $j=0; $j < 3; $j++){
	$trans_unitCoeff[$i][$j] /= $partition[$i];
	$rectVecs[$i][$j] = $transVecs[$i][$j]/$partition[$i];
    }
}

sub elementWiseRecipro($){
    my($vr)=@_; # $vr : reference of vector
    my(@ansVec);
    for(my $i=0; $i<3; $i++){
	if($vr->[$i]==0){
	    $ansVec[$i] = 1;
	}else{
	    $ansVec[$i] = &nint(1/($vr->[$i]));
	}
    }
    return \@ansVec; # return the reference of a resultant vector
}

my @rectTrans;
my @unit_rectCoeff;
for(my $i=0; $i < 3; $i++){
    $unit_rectCoeff[$i]=[];
    @{$unit_rectCoeff[$i]} = @{$unitVecs[$i]};
    &solveCoeff($unit_rectCoeff[$i],\@rectVecs);
    my $ansVecRef = elementWiseRecipro($unit_rectCoeff[$i]);
#     &printVec($ansVecRef);
#     printf
# 	"How many unit-translaion matches rectangular translation = %d\n",
# 	&lcm3(@$ansVecRef);
    $rectTrans[$i]=&lcm3(@$ansVecRef);
}

my @BravaisLattice;
for(my $i =0; $i < $rectTrans[0]; $i++){
    for(my $j =0; $j < $rectTrans[1]; $j++){
      outer: for(my $k =0; $k < $rectTrans[2]; $k++){
	  my @vec = ($i,$j,$k);
	  &linearCombination(\@vec,\@unit_rectCoeff);
	  for(my $n=0; $n < $#BravaisLattice+1; $n++){
	      if( &equivalent(\@vec,$BravaisLattice[$n]) ){
		  if(&norm2Vec(\@vec) < &norm2Vec($BravaisLattice[$n])){
		      $BravaisLattice[$n]=\@vec;
		      # shorter distance to the origin is preferable.
		  }
		  next outer;
	      }
	  }
	  push @BravaisLattice, \@vec;
      }
    }
}
printf
    "(# of Bravais lattice points) = %d\n",
    $#BravaisLattice+1;
print
    "List of Bravais lattice points (in the unit of minimal rectangular cell)\n";
foreach(@BravaisLattice){
    &printVec($_);
}

# Output
$doc = $parser->parse_file($outputFile)
    or die "can't parse \"$outputFile\"\n";

my $structureNode = $doc->getDocumentElement(); # (root node) = (<structure />)
if($structureNode->nodeName() ne 'structure' ){
    die 'root node does not <structure/>\n';
}
# print $structureNode->nodePath(); print "\n";
foreach my $atomNode ($doc->getElementsByTagName('atom')){
    $structureNode->removeChild($atomNode);
}

sub putUnitCell{ # internal co-ordinate in the unit of @rectVecs
    my @disp=@_;
    for(my $n=0; $n < $#positions+1; $n++){
	# construct new node
	my $newAtomNode = XML::LibXML::Element->new('atom');
	$newAtomNode->setAttribute('element',$elements[$n]);
	$newAtomNode->setAttribute('class',$classes[$n]);
	$newAtomNode->setAttribute('motion',$motions[$n]);
	my $newPositionNode = XML::LibXML::Element->new('position');
	$newPositionNode->setAttribute('unit','internal');
	my @vec=@{$positions[$n]};           # unit : @unitVecs
	&linearCombination(\@vec,\@unitVecs);
	&solveCoeff(\@vec,\@rectVecs);        # unit : @rectVecs
	for( my $i = 0; $i < 3 ; $i++){
	    $vec[$i] /= $partition[$i];      # unit : Given supercell
	    $vec[$i] += $disp[$i];
	}
	$newPositionNode->addChild($doc->createTextNode(sprintf " %22.15e %22.15e %22.15e ",@vec));
	$newAtomNode->addChild($newPositionNode);
	$structureNode->addChild($newAtomNode);
    }
}

print "List of all unit-cell origins (in the unit of given supercell)\n";
for(my $i =0; $i < $partition[0]; $i++){
    for(my $j =0; $j < $partition[1]; $j++){
	for(my $k =0; $k < $partition[2]; $k++){
	    for(my $n=0; $n < $#BravaisLattice+1; $n++){
		my @vec=($i+$BravaisLattice[$n][0],
			 $j+$BravaisLattice[$n][1],
			 $k+$BravaisLattice[$n][2]);
		for(my $m=0; $m<3; $m++){
		    $vec[$m] /= $partition[$m];
		}
		&printVec(\@vec); # internal co-ordinate
		putUnitCell(@vec);
	    }
	}
    }
}

#
#undef @unitVecs, @transVecs, @positions, @elements, @classes, @motions, @trans_unitCoeff ;


# Formatting & writing to "$outputFile"
$outputFile =~ s/\.xml/_new.xml/;
open(OUT,">$outputFile");

eval 'local $SIG{__WARN__} = sub { print STDERR "try -u option\n"; die ; };use XML::Twig;';
if($@){
    print "? error Message: $@ \n";
    $twigState="No Twig";
}else{
	if($twigState eq "necessary"){
	    my $twig = new XML::Twig;
	    $twig->set_indent(" "x2);
	    $twig->parse($doc->toString());
	    $twig->set_pretty_print("indented");
	    print OUT $twig->sprint();
	    $twigState="done";
	}
}
if($twigState ne "done"){
    print OUT $doc->toString();
}
close(OUT);

# output for Band calculation
open(OUT,">$bandCalcInput");

print OUT "# translation vectors of Bravais lattice (a.u.)\n";
for(my $i=0; $i<3; $i++){
    printf OUT "%23.15e %23.15e %23.15e\n",
    @{$rectVecs[$i]};
}

print OUT "# number of positions in the unitcell\n";
printf OUT "%4d\n", $#positions+1;

print OUT "# positions in the unitcell (internal co-ordinate)\n";
for(my $n=0; $n < $#positions+1; $n++){
    my @vec=@{$positions[$n]};           # unit : @unitVecs
    &linearCombination(\@vec,\@unitVecs);
    &solveCoeff(\@vec,\@rectVecs);        # unit : @rectVecs
    printf OUT "%23.15e %23.15e %23.15e\n", @vec;
}

print OUT "# number of unitcell\n";
printf OUT "%4d\n", &nint($volTrans/$volUnit) ;

print OUT "# translation vectors of the origin of unitcell\n";
print OUT "# (internal co-ordinate, with shifted center)\n";
my(@orgShift)=@partition;
my($odd)=(0==1);
foreach(@orgShift){
    if($_ % 2 != 0) {
	$odd = (1==1);
	$_ /= 2.0;
    }else{
	$_ /= 2;
    }
}
if($odd) {
    print  "** CAUTION ** : sizes of supercell includes Odd number.\n";
    printf "   The origin is shifted to (%9.1e, %9.1e, %9.1e).\n", 
	$orgShift[0]-&nint($orgShift[0]),
	$orgShift[1]-&nint($orgShift[1]),
	$orgShift[2]-&nint($orgShift[2]);
}
my($counter,$whereIsOrigin)=(0,0);
for(my $i =0; $i < $partition[0]; $i++){
    for(my $j =0; $j < $partition[1]; $j++){
	for(my $k =0; $k < $partition[2]; $k++){
	    for(my $n=0; $n < $#BravaisLattice+1; $n++){
		$counter++;
		if($i==$orgShift[0] && $j==$orgShift[1] && $k==$orgShift[2] && $n==0){
		    $whereIsOrigin=$counter;
		}
		my @vec=($i+$BravaisLattice[$n][0]-$orgShift[0],
			 $j+$BravaisLattice[$n][1]-$orgShift[1],
			 $k+$BravaisLattice[$n][2]-$orgShift[2]);
		printf OUT "%23.15e %23.15e %23.15e\n", @vec;
	    }
	}
    }
}
if(!$odd){
    print  OUT  "# k-th unitcell includes origin.\n";
    print  OUT  "# k is wrtten in the following line.\n";
    printf OUT  "%4d\n",$whereIsOrigin;
}
