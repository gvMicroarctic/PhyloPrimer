#! /usr/bin/perl -w

use strict;

#connected to tree_check_pp.pl : check uploaded gene similarity

#usage: ./genesCheck_pp.pl -file ./phyloprimer/analysesPhyloprimer/

my %args = @ARGV;

my $fileInfo = $args{-file};

open(FILEIN, "<$fileInfo") or die "Couldn't open the file: $!"; #input file with aligned sequences

#create the hash: each key is the base number: the values of the keys are 1/5 arrays (one array per possible character)

my $num = 0;
my $acc;
my %all;
my $num_seq = 0;
my @seq1;

while (defined(my $input = <FILEIN>)) {
    chomp($input);
    if ($input =~ "^>") {
        $num_seq++;
        if (defined $acc) {
            foreach my $i (@seq1) {
                $num++;
                ###create hash with num - array with letter for that position
                if ($i eq 'A') {
                    $all{$num}{'A'}++;
                } elsif ($i eq 'C') {
                    $all{$num}{'C'}++;
                } elsif ($i eq 'T') {
                    $all{$num}{'T'}++;
                } elsif ($i eq 'G') {
                    $all{$num}{'G'}++;
                } elsif (($i eq '-') or ($i eq '.')) {
                    $all{$num}{'GAP'}++;
                } elsif ($i eq 'R') {
                    $all{$num}{'A'} += 0.5;
                    $all{$num}{'G'} += 0.5;
                } elsif ($i eq 'Y') {
                    $all{$num}{'C'} += 0.5;
                    $all{$num}{'T'} += 0.5;
                } elsif ($i eq 'S') {
                    $all{$num}{'G'} += 0.5;
                    $all{$num}{'C'} += 0.5;
                } elsif ($i eq 'W') {
                    $all{$num}{'A'} += 0.5;
                    $all{$num}{'T'} += 0.5;
                } elsif ($i eq 'K') {
                    $all{$num}{'G'} += 0.5;
                    $all{$num}{'T'} += 0.5;
                } elsif ($i eq 'M') {
                    $all{$num}{'A'} += 0.5;
                    $all{$num}{'C'} += 0.5;
                } elsif ($i eq 'B') {
                    $all{$num}{'C'} += 0.333;
                    $all{$num}{'G'} += 0.333;
                    $all{$num}{'T'} += 0.333;
                } elsif ($i eq 'D') {
                    $all{$num}{'A'} += 0.333;
                    $all{$num}{'G'} += 0.333;
                    $all{$num}{'T'} += 0.333;
                } elsif ($i eq 'H') {
                    $all{$num}{'A'} += 0.333;
                    $all{$num}{'C'} += 0.333;
                    $all{$num}{'T'} += 0.333;
                } elsif ($i eq 'V') {
                    $all{$num}{'A'} += 0.333;
                    $all{$num}{'C'} += 0.333;
                    $all{$num}{'G'} += 0.333;
                } elsif ($i eq 'N') {
                    $all{$num}{'C'} += 0.25;
                    $all{$num}{'A'} += 0.25;
                    $all{$num}{'T'} += 0.25;
                    $all{$num}{'G'} += 0.25;
                }
            }
            undef @seq1;
            $num=0;
        }
        ($acc) = $input =~ />(\S+)/;
    } else {
        $input = uc $input; #all to upper case
        my @seq = split //, $input;
        push @seq1, @seq;
    }
}
# insert also the last sequence in the hash
foreach my $i (@seq1) {
    $num++;
    ###create hash with num - array with letter for that position
    if ($i eq 'A') {
        $all{$num}{'A'}++;
    } elsif ($i eq 'C') {
        $all{$num}{'C'}++;
    } elsif ($i eq 'T') {
        $all{$num}{'T'}++;
    } elsif ($i eq 'G') {
        $all{$num}{'G'}++;
    } elsif (($i eq '-') or ($i eq '.')) {
        $all{$num}{'GAP'}++;
    } elsif ($i eq 'R') {
        $all{$num}{'A'} += 0.5;
        $all{$num}{'G'} += 0.5;
    } elsif ($i eq 'Y') {
        $all{$num}{'C'} += 0.5;
        $all{$num}{'T'} += 0.5;
    } elsif ($i eq 'S') {
        $all{$num}{'G'} += 0.5;
        $all{$num}{'C'} += 0.5;
    } elsif ($i eq 'W') {
        $all{$num}{'A'} += 0.5;
        $all{$num}{'T'} += 0.5;
    } elsif ($i eq 'K') {
        $all{$num}{'G'} += 0.5;
        $all{$num}{'T'} += 0.5;
    } elsif ($i eq 'M') {
        $all{$num}{'A'} += 0.5;
        $all{$num}{'C'} += 0.5;
    } elsif ($i eq 'B') {
        $all{$num}{'C'} += 0.33333;
        $all{$num}{'G'} += 0.33333;
        $all{$num}{'T'} += 0.33333;
    } elsif ($i eq 'D') {
        $all{$num}{'A'} += 0.33333;
        $all{$num}{'G'} += 0.33333;
        $all{$num}{'T'} += 0.33333;
    } elsif ($i eq 'H') {
        $all{$num}{'A'} += 0.33333;
        $all{$num}{'C'} += 0.33333;
        $all{$num}{'T'} += 0.33333;
    } elsif ($i eq 'V') {
        $all{$num}{'A'} += 0.33333;
        $all{$num}{'C'} += 0.33333;
        $all{$num}{'G'} += 0.33333;
    } elsif ($i eq 'N') {
        $all{$num}{'C'} += 0.25;
        $all{$num}{'A'} += 0.25;
        $all{$num}{'T'} += 0.25;
        $all{$num}{'G'} += 0.25;
    }
}

#assign to each column a letter (IUPAC nomenclature)

my $thr = 0.20;
my $thr1 = 0.15;
my $final_consensus = '';

my $finalGap = 0;
my $finalMatch = 0;
my $finalMis = 0;
my $finalAll = 0;

for my $num (1..(keys %all)) { #for each position
    my $a = 0;
    my $c = 0;
    my $t = 0;
    my $g = 0;
    my $gap = 0;
    
    if (defined($all{$num}{'A'})) {
        $a = ($all{$num}{'A'}/$num_seq);
    }
    if (defined($all{$num}{'C'})) {
        $c = ($all{$num}{'C'}/$num_seq);
    }
    if (defined($all{$num}{'T'})) {
        $t = ($all{$num}{'T'}/$num_seq);
    }
    if (defined($all{$num}{'G'})) {
        $g = ($all{$num}{'G'}/$num_seq);
    }
    if (defined($all{$num}{'GAP'})) {
        $gap = ($all{$num}{'GAP'}/$num_seq);
    }
    $finalAll++;
    if ($gap != 0) {
        $finalGap++;
    } elsif (($a == 1) or ($c == 1) or ($g == 1) or ($t == 1)) { ## A
            $finalMatch++;
    } else { ## C
            $finalMis++;
    }
}

my $finalGap_diff = sprintf("%.1f",($finalGap/$finalAll *100));
my $finalMatch_diff = sprintf("%.1f",($finalMatch/$finalAll *100));
my $finalMis_diff = sprintf("%.1f",($finalMis/$finalAll *100));

my $finalPerc = $finalGap_diff . "-" . $finalMatch_diff . "-" . $finalMis_diff;

print "$finalPerc";

