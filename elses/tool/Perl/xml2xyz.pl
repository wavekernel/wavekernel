#! /usr/bin/perl
# -*- mode: Perl -*-
#================================================================
# ELSES version 0.02
# Copyright (C) ELSES. 2007-2010 all rights reserved
#================================================================

use strict;
use warnings;
use XML::LibXML;

my $Bohr = 0.529177e0;  #my $Bohr = 5.291772108e-1; # length unit in xyz= "Angstrom"

if( $#ARGV < 0 ){
  die "Usage: xml2xyz.pl filename.xml\n";
}

my $filename = $ARGV[0];
my $parser = XML::LibXML->new();
my $doc    = $parser->parse_file($filename)
  or die "can't parse \"$filename\"\n";
#my $root   = $doc->documentElement(); #DOM
#print $doc->toString(1);

$filename =~ s/\.xml/.xyz/ or die 'Extension must be ".xml"';

my @vectors=
  $doc->getElementsByTagName('unitcell')->
  shift->getElementsByTagName('vector');

if($#vectors != 2){
  die "Wrong number of vectors\n";
}

for ( my $i=0 ; $i < 3 ; $i++ ) {
  my $vector = shift(@vectors);
  $vector->string_value =~ /^\s+(\S+)\s+(\S+)\s+(\S+)/ ;
  if($vector->getAttribute('unit') eq 'a.u.'){
    push(@vectors,[$1,$2,$3]);
  }elsif($vector->getAttribute('unit') eq 'angstrom'){
    push(@vectors,[$1/$Bohr,$2/$Bohr,$3/$Bohr]);
  }else{
     die "\"unit\" attribute must be 'internal' or 'angstrom'.\n";
  }
}

# for ( my $i=0 ; $i < 3 ; $i++ ) {
# 	print "$vectors[$i][0], $vectors[$i][1], $vectors[$i][2] \n";
# }

my @xyzLines = ();
my @original = ();
my @atomList = $doc->getElementsByTagName('atom');


foreach my $atomNode (@atomList) {
  my  $nodeList = $atomNode->getElementsByTagName('position');
  if( $nodeList->size >1 ) { 	#      debugging  >=1, original >1
	die "Too many positions in ",$atomNode->getAttribute('element'),"\n";
  }

  my ($node, $line, $unit, $x, $y, $z);

  $node = $nodeList->get_node(0);
  $line = $node->string_value;
  $unit = $node->getAttribute('unit');

  push @original, $line;
  $line =~ /^\s*(\S+)\s+(\S+)\s+(\S+)/;
#   my $format = "";
#   for (my $i=0; $i<12; $i++){
# 	$format .= " %.2e";
#   }

  eval {
	local $SIG{__WARN__} =
	  sub {
		print STDERR "unit vectors\n";
		print STDERR "$vectors[0][0],$vectors[0][1],$vectors[0][2]\n";
		print STDERR "$vectors[1][0],$vectors[1][1],$vectors[1][2]\n";
		print STDERR "$vectors[2][0],$vectors[2][1],$vectors[2][2]\n";
		print STDERR "atom position\n";
		print STDERR "$1,$2,$3\n";
		die "Numeric format error\n";
	  };
	if($unit eq 'internal') {
	  $x=($1*$vectors[0][0]+$2*$vectors[1][0]+$3*$vectors[2][0])*$Bohr;
	  $y=($1*$vectors[0][1]+$2*$vectors[1][1]+$3*$vectors[2][1])*$Bohr;
	  $z=($1*$vectors[0][2]+$2*$vectors[1][2]+$3*$vectors[2][2])*$Bohr;
	}elsif($unit eq 'angstrom'){
	  $x=$1*1.0e0;
	  $y=$2*1.0e0;
	  $z=$3*1.0e0;
	}else{
	  die "\"unit\" attribute must be 'internal' or 'angstrom'.\n";
	}
  };
  if($@) { die $@; }
  push @xyzLines, [ $atomNode->getAttribute('element'),
					$x, $y, $z ];
}

# Checking distances ...
for(my $i = 0 ; $i <= $#xyzLines ; $i++) {
  for(my $j = $i +1 ; $j <= $#xyzLines ; $j++) {
	my ($x,$y,$z) =
	  (${$xyzLines[$i]}[1]-${$xyzLines[$j]}[1],
	   ${$xyzLines[$i]}[2]-${$xyzLines[$j]}[2],
	   ${$xyzLines[$i]}[3]-${$xyzLines[$j]}[3]);
    my $d = sqrt($x*$x+$y*$y+$z*$z);
    if( $d < 1e-12 ) {
	  print "The same position $i (",
		join(',',@{$xyzLines[$i]}),
		  ") appeared twice at $j.\n";
	  print "Line: ",$original[$i],"\n";
	}
  }
}

open(OUT,">$filename");
print OUT $#xyzLines +1,"\n";
print OUT "# converted from \"$ARGV[0]\"\n";
foreach (@xyzLines){
  my ($element, $x, $y, $z) = @$_;
  my $output=sprintf("%s   %22.15e %22.15e %22.15e\n",$element,$x,$y,$z);
  print OUT $output;
}
