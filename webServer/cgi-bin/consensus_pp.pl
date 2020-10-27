#! /usr/bin/perl -w

use strict;

#connected to tree_consensus_pp.pl and sequence_consensus_pp.pl : calculate consensus sequence from DNA aligned sequences

#usage: ./consensus.pl file

open(FILEIN, "<$ARGV[0]") or die "Couldn't open the file: $!"; #input file with aligned sequences

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
    if ($gap >=0.50) {
        $final_consensus .= "-";
    }   else {
        if (($a >=$thr) && ($c <$thr) && ($g <$thr) && ($t <$thr)) { ## A
            $final_consensus .= "A";
        } elsif (($c >=$thr) && ($a <$thr) && ($g <$thr) && ($t <$thr)) { ## C
            $final_consensus .= "C";
        } elsif (($t >=$thr) && ($a <$thr) && ($c <$thr) && ($g <$thr)) { ## G
            $final_consensus = $final_consensus . "T";
        } elsif (($g >=$thr) && ($a <$thr) && ($c <$thr) && ($t <$thr)) { ## T
            $final_consensus = $final_consensus . "G";
        } elsif (($a >=$thr) && ($t >=$thr) && ($c <$thr) && ($g <$thr)) { ## A and T
            $final_consensus = $final_consensus . "W";
        } elsif (($c >=$thr) && ($g >=$thr) && ($a <$thr) && ($t <$thr)) { ## C and G
            $final_consensus = $final_consensus . "S";
        } elsif (($a >=$thr) && ($g >=$thr) && ($c <$thr) && ($t <$thr)) { ## A and G
            $final_consensus = $final_consensus . "R";
        } elsif (($c >=$thr) && ($t >=$thr) && ($a <$thr) && ($g <$thr)) { ## C and T
            $final_consensus = $final_consensus . "Y";
        } elsif (($g >=$thr) && ($t >=$thr) && ($a <$thr) && ($c <$thr)) { ## G and T
            $final_consensus = $final_consensus . "K";
        } elsif (($a >=$thr) && ($c >=$thr) && ($g <$thr) && ($t <$thr)) { ## A and C
            $final_consensus = $final_consensus . "M";
        } elsif (($a >=$thr) && ($c >=$thr) && ($g >=$thr) && ($t <$thr)) { ## A and C and G
            $final_consensus = $final_consensus . "V";
        } elsif (($a >=$thr) && ($c >=$thr) && ($t >=$thr) && ($g <$thr)) { ## A and C and T
            $final_consensus = $final_consensus . "H";
        } elsif (($c >=$thr) && ($g >=$thr) && ($t >=$thr) && ($a <$thr)) { ## C and G and T
            $final_consensus = $final_consensus . "B";
        } elsif (($a >=$thr) && ($g >=$thr) && ($t >=$thr) && ($c <$thr)) { ## A and G and T
            $final_consensus = $final_consensus . "D";
        } elsif (($a >=$thr1) && ($c >=$thr1) && ($g >=$thr1) && ($t >=$thr1)) { ## A and C and G and T
            $final_consensus = $final_consensus . "N";
        }
    }
}

print "$final_consensus";
