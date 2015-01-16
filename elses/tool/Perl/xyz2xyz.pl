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

my($xyz_file,$xml_file,$output);

if( $#ARGV < 1 ) {
  die <<EOT;

Purpose: Wrapping around atom positions into periodic cell.

Usage: 
  xyz2xyz.pl structure.xml positions.xyz

Output:
  positionsn_new.xyz

EOT

}elsif( $#ARGV == 1 ){
  ($xml_file,$xyz_file)= @ARGV;
  $output = $xyz_file;
  $output =~ s/\.xyz/_new.xyz/;
}

if( $xml_file !~ /\.xml$/ ) {
  die "The 1st argument must be \"*.xml\"";
}
if( $xyz_file !~ /\.xyz$/ ) {
  die "The 2nd argument must be \"*.xyz\"";
}

# read_xyz()
open(IN,$xyz_file);
open(OUT,">$output");

sub read_xyz($){
  my ($ref_pos)=@_;
  my @pos;
  eval{
	local $SIG{__WARN__} = sub { die $_[0]; };
	my $line=<IN>;
	$line =~ /(\S+)/;
	my $number_of_positions=$1;
	$line = <IN>; # skip single comment line
	for( my $i=0; $i < $number_of_positions; $i++ ){
	  $line = <IN>;
	  $line =~ s/^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)//;
	  push (@pos,[$1,$2,$3,$4]);
	}
  };
  if(! $@){
	@{$ref_pos}=@pos;
  }
  return $@;
}

sub write_xyz($){
  my ($ref_pos)=@_;
  my $number_of_positions=$#$ref_pos+1;
  print OUT "$number_of_positions\n";
  print OUT "Created by xyz2xyz.pl using $xml_file \& $xyz_file\n";
  foreach (@{$ref_pos}){
	my ($element, $x, $y, $z) = @$_;
	my $output=sprintf("%s   %22.15e %22.15e %22.15e\n",$element,$x,$y,$z);
	print OUT $output;
  }
}

sub toAngstrom($$){
    my $Bohr = 0.529177e0;  #my $Bohr = 5.291772108e-1;
    my ($unit,$vr) = @_; # $vr = (reference of vector)
    if($unit eq 'a.u.'){
	  foreach(@$vr) {
		# print "$_\n"; # for debug
	    $_ *= $Bohr;
	  }
	  # change unit -> angstrom
    }elsif($unit ne 'angstrom'){
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
	&toAngstrom($unit,$unit_vec[$i]);
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

my @unitVecs;

sub wraparoundPositions($){
  my ($refPos)=@_ ; # $refPos means reference to @positions
  my @vec;
  foreach(@$refPos){
      my($a,$x,$y,$z)=@$_;
      @vec=($x,$y,$z);
      &solveCoeff(\@vec,\@unitVecs);
      foreach(@vec){
	  $_ = $_ - floor($_);
      }
      &linearCombination(\@vec,\@unitVecs);
      @$_=($a,$vec[0],$vec[1],$vec[2]);
  }
}

# read unit vectors
&getUnitVecs($xml_file,\@unitVecs);

my @positions;
while(! &read_xyz(\@positions)){
  &wraparoundPositions(\@positions);
  &write_xyz(\@positions);
}
close(IN);
close(OUT)
