#! /usr/bin/perl
# -*- mode: Perl -*-

use strict;
use warnings;
use diagnostics;
use XML::LibXML;

my $Bohr = 0.529177e0;  #my $Bohr = 5.291772108e-1;
# length unit in xyz= "Angstrom"

my($xyz_file,$xml_file,$slice_number);

if( $#ARGV < 0 ) {
  die <<EOT;

Usage:
  xyz2xml.pl positions.xyz template.xml [slice #]
  If slice # is omitted, the last slice is cut.

EOT
}elsif( $#ARGV == 1 ){
  ($xyz_file,$xml_file)= @ARGV;
  undef($slice_number);
}elsif( $#ARGV == 2 ){
  ($xyz_file,$xml_file,$slice_number)= @ARGV;
}

if( $xyz_file !~ /\.xyz$/ ) {
  die "The 1st argument must be \"*.xyz\"";
}
if( $xml_file !~ /\.xml$/ ) {
  die "The 2nd argument must be \"*.xml\"";
}

# read_xyz()
open(IN,$xyz_file);
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

my @positions;
my $slice_count=0;
while(! &read_xyz(\@positions)){
  if(defined($slice_number)){
	$slice_count++;
	if($slice_count==$slice_number){
	  print "The ${slice_number}-th slice is cut out\n";
	  last;
	}
  }
}
close(IN);

my $parser = XML::LibXML->new();
my $doc    = $parser->parse_file($xml_file)
  or die "can't parse \"$xml_file\"\n";
#my $root   = $doc->documentElement(); #DOM
#print $doc->toString(1);

my @unit_vec=
  $doc->getElementsByTagName('unitcell')->
  shift->getElementsByTagName('vector');

if($#unit_vec != 2){
  die "Wrong number of unit_vecs\n";
}

my @vectors;
for ( my $i=0 ; $i < 3 ; $i++ ) {
  $unit_vec[$i]->string_value =~ /^\s+(\S+)\s+(\S+)\s+(\S+)/ ;
  if($unit_vec[$i]->getAttribute('unit') eq 'a.u.'){
	$unit_vec[$i] = [$1,$2,$3];
  }elsif($unit_vec[$i]->getAttribute('unit') eq 'angstrom'){
	# construct new node
	my $newVectorNode = XML::LibXML::Element->new('vector');
	$newVectorNode->setAttribute('unit','a.u.');
	$newVectorNode->appendTextNode(sprintf(' %22.15e %22.15e %22.15e',$1/$Bohr,$2/$Bohr,$3/$Bohr));
	my $parentNode=$unit_vec[$i]->parentNode;
	# replace node
	$parentNode->replaceChild($newVectorNode,$unit_vec[$i]);
	# Now $unit_vec[$i] can be disposed.
	$unit_vec[$i] = [$1/$Bohr,$2/$Bohr,$3/$Bohr];
	# change unit -> a.u. # Not implemented yet
	
  }else{
	die "\"unit\" attribute must be 'a.u.' or 'angstrom'.\n";
  }
  @vectors=@unit_vec;
}

sub det3 ($) {
  my ($vr) =  @_ ; # $vr = (reference of vector)
#   print "passed vectors\n";
#   foreach (@$vr) {
# 	foreach (@$_) {
#  	  print $_,",";
#  	}
#  	print "\n";
#   }
  return  $vr->[0][0] * $vr->[1][1] * $vr->[2][2] - $vr->[0][0] * $vr->[1][2] * $vr->[2][1]
	- $vr->[0][1] * $vr->[1][0] * $vr->[2][2] + $vr->[0][1] * $vr->[1][2] * $vr->[2][0]
	+ $vr->[0][2] * $vr->[1][0] * $vr->[2][1] - $vr->[0][2] * $vr->[1][1] * $vr->[2][0];
}
		  
my $vol = det3(\@vectors);
#print " Volume=$vol\n";

my $atom_count=0;

my ($nodeList, $node, $unit, $atomNode, $newPositionNode);

foreach $atomNode ($doc->getElementsByTagName('atom')) {
  $nodeList = $atomNode->getElementsByTagName('position');
  if( $nodeList->size >1 ) {
        die "Too many positions in ",$atomNode->getAttribute('element'),"\n";
  }

  $vol = det3(\@vectors);
  #print "xVolume=$vol\n";

  $node = $nodeList->get_node(0);
  $unit = $node->getAttribute('unit');
  
  my $x = $positions[$atom_count][1] / $Bohr;
  my $y = $positions[$atom_count][2] / $Bohr;
  my $z = $positions[$atom_count][3] / $Bohr;

  $newPositionNode = XML::LibXML::Element->new('position');
  $newPositionNode->setAttribute('unit','a.u.');
  $newPositionNode->appendTextNode( sprintf(' %22.15e %22.15e %22.15e', $x, $y, $z));

  $atomNode->replaceChild($newPositionNode,$node);

  $atom_count++;
  
}
# data consistency check
if($atom_count != $#positions+1) {
  die "The number of positions does not match\n.";
}

# Output "*_new.xml" file
my $new_file=$xml_file;
$new_file =~ s/\.xml/_new.xml/;
$doc->toFile($new_file);
