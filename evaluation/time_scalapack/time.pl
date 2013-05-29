#!/usr/bin/perl
use warnings;
use strict;

my @output_files = `ls out/o*`;

my @info_all = ();
foreach my $file (@output_files) {
    open IN, $file or die;
    my %info = ();
    my $is_elapse_time_region = 0;
    while (<IN>) {
        if (/procs:\s*\d+\s*x\s*\d+\s*\(\s*(\d+)\s*\)/) {
            $info{'procs'} = $1;
        }
        $is_elapse_time_region = 1 if (/Elapse time/);
        if ($is_elapse_time_region && /(\w+)\s*:\s*([\d\.]+)/) {
            $info{"t_$1"} = $2;
        }
    }
    while (my ($name, $value) = each %info) {
        #print $name, $value, "\n";
    }
    push @info_all, \%info;
    close IN;
}

print "n_procs\tPDSYTRD\tPDSTEDC\tPDORMTR\tTOTAL\n";

@info_all = sort {$a->{'procs'} <=> $b->{'procs'}} @info_all;
foreach my $i (@info_all) {
    my $p = $i->{'procs'};
    print "$p\t$i->{t_pdsytrd}\t$i->{t_pdstedc}\t$i->{t_pdormtr}\t$i->{t_total}\n";
}
