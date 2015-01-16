#!/usr/bin/perl
# -*- mode: Perl -*-
opendir(DIR,".");

#$search_str ="^[c,C](\\s+implicit none)";
#$replace_str ="      implicit none";
$search_str  ="lengthScale";
$replace_str ="myLength";

foreach(readdir(DIR))
{
  
	next unless /(.+)\.f90$/;
	my $bak="$_~";
	my $src=$_;

# Read source
	open(IN,"$src")  or die "Can't open $src\n";
	my $f="0";
	undef @buff;
	while(<IN>){
	# String Substitution
		$f="1" if s/$search_str/$replace_str/g;
		push @buff,$_;
	}
	close(IN,"$src");

	# 
	next if $f == "0";

	unlink $bak;
	rename $src,$bak;
	print "renamed: $src -> $bak\n";
    # Write destination
	open(OUT,">$src")or die "Can't open $src\n";
	foreach(@buff)
	{
		print OUT;
	}
	close(OUT);
}
