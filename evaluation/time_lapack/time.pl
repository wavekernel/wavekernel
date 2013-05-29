#!/usr/bin/perl
use warnings;
use strict;

my @num_threads = (1, 2, 3, 4);
my $size = $ARGV[0];

my @times_all = ();
foreach my $p (@num_threads) {
    my $filename = "p${p}_n$size.log";
    open IN, $filename or die;
    my %times = ();
    while (<IN>) {
        if (/(\w+) TIME\s+=\s+([\d\.]+)/) {
            $times{$1} = $2;
        }
    }
    while (my ($name, $time) = each %times) {
#        print $name, $time, "\n";
    }
    push @times_all, \%times;
    close IN;
}
print "# p\tDSYTRD\tDSTEDC\tDORMTR\tTOTAL\n";
foreach my $p (@num_threads) {
    my %times = %{shift @times_all};
    print "$p\t$times{DSYTRD}\t$times{DSTEDC}\t$times{DORMTR}\t$times{TOTAL}\n";
}
