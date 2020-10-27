#!/usr/bin/perl

use strict;

#connected to Primer_design_pp.pl, Probe_design_pp.pl, Oligo_design_pp.pl, Primer_check_pp.pl, Probe_check_pp.pl and Oligo_check_pp.pl : self-dimer calculation

#Usage: ./Self_dimer_checking_web.pl -primer -mg -mon -oligo -dntp -t

#set path to folders:
my $path_cgi = 'path_to_cgi_folder';

my %args = @ARGV;

my $folder = $args{-folder};
my $file_primer = $args{-primer};
my $mg_tot = $args{-mg}; #Mg ion concentration
my $monovalent = $args{-mon}; #monovalent ion concentration
my $C = $args{-oligo}; #oligo concentration
my $dNTP_tot = $args{-dntp}; #dNTP concentration
my $temperature_celsius = $args{-t}; #temperature

#all structures
open(my $file, ">$folder/tmp/SecondaryStructure_${file_primer}") or die; #file with all the necessary parameters
open(my $file1, ">$folder/tmp/SecondaryStructure_${file_primer}.delta") or die; #file with all the necessary parameters

my %corr;
$corr{'A'} = 'T';
$corr{'T'} = 'A';
$corr{'G'} = 'C';
$corr{'C'} = 'G';


$C = 0.000002; #oligo concentration
$monovalent = 0.005; #monovalent ion concentration
$mg_tot = 0.0015; #Mg ion concentration
$dNTP_tot = 0.00001; #dNTP concentration


#degenarate oligo check: create all the oligo possible alternatives
#wildcards
my %wildcard;
$wildcard{'R'}{'A'} = '';
$wildcard{'R'}{'G'} = '';

$wildcard{'Y'}{'C'} = '';
$wildcard{'Y'}{'T'} = '';

$wildcard{'S'}{'G'} = '';
$wildcard{'S'}{'C'} = '';

$wildcard{'W'}{'A'} = '';
$wildcard{'W'}{'T'} = '';

$wildcard{'K'}{'G'} = '';
$wildcard{'K'}{'T'} = '';

$wildcard{'M'}{'A'} = '';
$wildcard{'M'}{'C'} = '';

$wildcard{'B'}{'C'} = '';
$wildcard{'B'}{'G'} = '';
$wildcard{'B'}{'T'} = '';

$wildcard{'D'}{'A'} = '';
$wildcard{'D'}{'G'} = '';
$wildcard{'D'}{'T'} = '';

$wildcard{'H'}{'A'} = '';
$wildcard{'H'}{'C'} = '';
$wildcard{'H'}{'T'} = '';

$wildcard{'V'}{'A'} = '';
$wildcard{'V'}{'C'} = '';
$wildcard{'V'}{'G'} = '';

$wildcard{'N'}{'C'} = '';
$wildcard{'N'}{'A'} = '';
$wildcard{'N'}{'T'} = '';
$wildcard{'N'}{'G'} = '';

#open file with all oligos
open(IN, "<$folder/tmp/$file_primer") or die;

while (my $input = <IN>) {
    chomp($input);
    my $primer_input = $input;
    print $file1 "$primer_input\t";
    my $len = length($primer_input);
    my $primer_inputRev = reverse($primer_input);

    #check if degenerancy Seq
    my ($primerSeqRef) = degenerateAlt($primer_input);
    my @primerSeq = @{$primerSeqRef};

    my $pair = 0;
    my %done;
    my %allDeg;

    foreach my $seq (@primerSeq) {
        if ($pair >= 1) {
            undef %done;
            foreach my $pair1 (keys %allDeg) {
                if (!(defined($done{$allDeg{$pair1}{'R'}}))) {
                    $pair++;
                    $allDeg{$pair}{'F'} = $seq;
                    $allDeg{$pair}{'R'} = $allDeg{$pair1}{'R'};
                    $done{$allDeg{$pair1}{'R'}} = '';
                }
            }
        }
        $pair++;
        my $rev = reverse($seq); #reverse of the oligo
        $allDeg{$pair}{'F'} = $seq;
        $allDeg{$pair}{'R'} = $rev;
    }

    #foreach my $pair (keys %allDeg) {
    #    print "all\t$pair$allDeg{$pair}{'F'}\t$allDeg{$pair}{'R'}\n";
    #}


    my %final_max;
    my %done;
    my $passPrint = 0;

    my $degPos = 0;
    if ($primer_input =~ /[RYSWKMBDHVN]/) {
        $degPos = 1;
    }

    print $file "\n>\t$primer_input\n";

    #foreach my $seq (@primerSeq) {
    #    my $rev = reverse($seq); #reverse of the oligo

    foreach my $pair (keys %allDeg) {
        my $seq = $allDeg{$pair}{'F'};
        my $rev = $allDeg{$pair}{'R'};
        
        my @seq1 = split /(?=[A-Z])/i, $seq;
        my @rev1 = split /(?=[A-Z])/i, $rev;
        #sequence without shift
        my $count=0;
        my $consecutive = 0;
        my $consecutive_gap = 0;
        my $gap = 0;
        my $num = 0;
        my $inside = 0;
        my $more = 0;
        my $time = 0;
        my %small;
        my %max;
        my @alignment;
        my $start;
        my $step;
        foreach my $position (0..($len-1)) {
            if ($seq1[$position] eq $corr{$rev1[$position]}) {
                my $real = $position + 1;
                $count++;
                push @alignment, '|';
                if ($inside == 0) {
                    $start = $position;
                }
                if ($gap == 0) {
                    $consecutive++;
                    $inside = 1;
                } elsif ($gap == 1) {
                    $consecutive++;
                    $consecutive_gap++;
                    $inside = 1;
                    $gap = 0;
                }
            } else {
                if ($inside == 1) {
                    $gap++;
                    if ($gap >=2) {
                        $num++;
                        $step = $consecutive + $consecutive_gap;
                        $small{$num} = $start . "_" . $step . "_" . $time;
                        if ($step > 1) {
                            $more = 1;
                        }
                        $gap = 0;
                        $inside = 0;
                        $consecutive = 0;
                        $consecutive_gap = 0;
                        $inside = 0;
                    }
                }
                push @alignment, ' ';
            }
        }
        if ($inside == 1) {
            $num++;
            $step = $consecutive + $consecutive_gap;
            $small{$num} = $start . "_" . $step . "_" . $time;
            if ($step > 1) {
                $more = 1;
            }
        }
        
        if (($num > 0) && ($more == 1)) {
            undef %max;
            foreach $num (sort keys %small) {
                my @data = split(/_/, $small{$num}); #start, step, shift
                my $check1 = substr $seq, $data[0], $data[1]; #portion alignment F
                my $check2 = substr $rev, $data[0], $data[1]; #portion alignment R
                if (length($check1) > 1) {
                    my $dG = `$path_cgi/DeltaG_calculation_pp.pl -primerF $check1 -primerR $check2 -type self -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius`;
                    
                    chomp($dG);
                    $dG = sprintf("%.2f", $dG);
    #                print "dG\t$check1 and $check2\t$dG\n";
                    $max{$dG}{$small{$num}} = $small{$num} . '-' . $check1 . '-' . $check2 . '-' . $dG;
                }
            }
        }
        my $max = 0;
        my $max_dG;
        my $infoMax;
        my $deg = 0;
        foreach my $dG (sort {$a <=> $b} keys %max) {
            if ($max == 0) {
                foreach my $info (sort {$a <=> $b} keys %{$max{$dG}}) {
                    if ($max == 0) {
                        $infoMax = $info;
                        $max_dG = $dG;
                        my $pos = 0;
                        my @data = split(/_/, $info); #start, step, shift
                        if ($degPos == 1) {
                            $infoMax = $max{$dG}{$info};
                        }
                        foreach my $match (@alignment) {
                            if ($match eq "|") {
                                if (($pos >= $data[0]) && ($pos <= ($data[0] + $data[1]))) {
                                } else {
                                    $alignment[$pos] = "¦";
                                }
                            }
                            $pos++;
                        }
                        $max++;
                    } else {
                        last;
                    }
                }
            } else {
                last;
            }
        }
        if ($max_dG < 0) {
            if (!(defined($done{$infoMax}))) { #check that it has not already be saved
                $passPrint++;
                push @{$final_max{$max_dG}{$passPrint}{'seq1'}}, @seq1;
                push @{$final_max{$max_dG}{$passPrint}{'al'}}, @alignment;
                push @{$final_max{$max_dG}{$passPrint}{'seq2'}}, @rev1;
                $done{$infoMax} = "";
            }
        }
        
        #Shift to the left
        my @rev2 = @rev1;
        my @seq2 = @seq1;
        foreach my $time (0..($len-2)) {
            shift(@rev2); #shift array to the left
            my @alignment;
            my @spaces;
            my $count=0;
            my $consecutive = 0;
            my $consecutive_gap = 0;
            my $gap = 0;
            my $num = 0;
            my $inside = 0;
            my $more = 0;
            undef %small;
            foreach my $position (0..($len-1)) {
                if ($seq2[$position] eq $corr{$rev2[$position]}) {
                    my $real = $position + 1;
                    $count++;
                    push @alignment, '|';
                    if ($inside == 0) {
                        $start = $position;
                    }
                    if ($gap == 0) {
                        $consecutive++;
                        $inside = 1;
                    } elsif ($gap == 1) {
                        $consecutive++;
                        $consecutive_gap++;
                        $inside = 1;
                        $gap = 0;
                    }
                } else {
                    if ($inside == 1) {
                        $gap++;
                        if ($gap >=2) {
                            $num++;
                            $step = $consecutive + $consecutive_gap;
                            $small{$num} = $start . "_" . $step . "_" . $time;
                            if ($step > 1) {
                                $more = 1;
                            }
                            $gap = 0;
                            $inside = 0;
                            $consecutive = 0;
                            $consecutive_gap = 0;
                            $inside = 0;
                        }
                    }
                    push @alignment, ' ';
                }
            }
            if ($inside == 1) {
                $num++;
                $step = $consecutive + $consecutive_gap;
                $small{$num} = $start . "_" . $step . "_" . $time;
                if ($step > 1) {
                    $more = 1;
                }
            }
            if (($num > 0) && ($more == 1)) {
                my $real = $time + 1;
                foreach my $time1 (1 .. $real) {
                    push @spaces, ' ';
                }
                my @seq3 = @seq2;
                unshift @seq3, @spaces; #shift array to the right
                unshift @alignment, @spaces; #shift array to the right
                undef %max;
                foreach my $num (sort keys %small) {
                    my ($check1, $check2) = checkCalc($seq,$rev,"left", $small{$num});
                    if (length($check1) > 1) {
                        my $dG = `$path_cgi/DeltaG_calculation_pp.pl -primerF $check1 -primerR $check2 -type self  -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius`;
                        chomp($dG);
                        $dG = sprintf("%.2f", $dG);
    #                    print "dG\t$check1 and $check2\t$dG\n";
                        $max{$dG}{$small{$num}} = $small{$num} . '-' . $check1 . '-' . $check2 . '-' . $dG;
                    }
                }
                my $max = 0;
                my $max_dG;
                my $infoMax;
                my $deg = 0;
                foreach my $dG (sort {$a <=> $b} keys %max) {
                    if ($max == 0) {
                        foreach my $info (sort {$a <=> $b} keys %{$max{$dG}}) {
                            if ($max == 0) {
                                $infoMax = $info;
                                $max_dG = $dG;
                                my $pos = 0;
                                my @data = split(/_/, $info); #start, step, shift
                                if ($degPos == 1) {
                                    $infoMax = $max{$dG}{$info};
                                }
                                foreach my $match (@alignment) {
                                    if ($match eq "|") {
                                        my $one = $data[0] + $data[2] + 1;
                                        my $two = $data[0] + $data[1] + $data[2] + 1;
                                        if (($pos >= ($data[0] + $data[2] + 1)) && ($pos <= ($data[0] + $data[1] + $data[2] + 1))) {
                                        } else {
                                            $alignment[$pos] = "¦";
                                        }
                                    }
                                    $pos++;
                                }
                                $max++;
                            } else {
                                last;
                            }
                        }
                    } else {
                        last;
                    }
                }
                if ($max_dG < 0) {
                    if (!(defined($done{$infoMax}))) { #check that it has not already be saved
                        $passPrint++;
                        push @{$final_max{$max_dG}{$passPrint}{'seq1'}}, @seq3;
                        push @{$final_max{$max_dG}{$passPrint}{'al'}}, @alignment;
                        push @{$final_max{$max_dG}{$passPrint}{'seq2'}}, @rev1;
                        $done{$infoMax} = "";
                    }
                }
            }
        }
        #Shift to the right
        my @rev2 = @rev1;
        my @seq2 = @seq1;
        foreach my $time (0..($len-2)) {
            unshift @rev2, ' '; #shift array to the right
            my @alignment;
            my @spaces;
            my $count=0;
            my $consecutive = 0;
            my $consecutive_gap = 0;
            my $gap = 0;
            my $num = 0;
            my $inside = 0;
            my $more = 0;
            undef %small;
            foreach my $position (0..($len-1)) {
                if ($seq2[$position] eq $corr{$rev2[$position]}) {
                    my $real = $position + 1;
                    $count++;
                    push @alignment, '|';
                    if ($inside == 0) {
                        $start = $position;
                    }
                    if ($gap == 0) {
                        $consecutive++;
                        $inside = 1;
                    } elsif ($gap == 1) {
                        $consecutive++;
                        $consecutive_gap++;
                        $inside = 1;
                        $gap = 0;
                    }
                } else {
                    if ($inside == 1) {
                        $gap++;
                        if ($gap >=2) {
                            $num++;
                            $step = $consecutive + $consecutive_gap;
                            $small{$num} = $start . "_" . $step . "_" . $time;
                            if ($step > 1) {
                                $more = 1;
                            }
                            $gap = 0;
                            $inside = 0;
                            $consecutive = 0;
                            $consecutive_gap = 0;
                            $inside = 0;
                        }
                    }
                    push @alignment, ' ';
                }
            }
            if ($inside == 1) {
                $num++;
                $step = $consecutive + $consecutive_gap;
                $small{$num} = $start . "_" . $step . "_" . $time;
                if ($step > 1) {
                    $more = 1;
                }
            }
            if (($num > 0) && ($more == 1)) {
                undef %max;
                foreach my $num (sort keys %small) {
                    my ($check1, $check2) = checkCalc($seq,$rev,"right", $small{$num});
                    if (length($check1) > 1) {
                        my $dG = `$path_cgi/DeltaG_calculation_pp.pl -primerF $check1 -primerR $check2 -type self  -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius`;
                        $dG = sprintf("%.2f", $dG);
                        chomp($dG);
    #                    print "dG\t$check1 and $check2\t$dG\n";
                        $max{$dG}{$small{$num}} = $small{$num} . '-' . $check1 . '-' . $check2 . '-' . $dG;
                    }
                }
                my $max = 0;
                my $max_dG;
                my $infoMax;
                my $deg = 0;
                foreach my $dG (sort {$a <=> $b} keys %max) {
                    if ($max == 0) {
                        foreach my $info (sort {$a <=> $b} keys %{$max{$dG}}) {
                            if ($max == 0) {
                                $infoMax = $info;
                                $max_dG = $dG;
                                my $pos = 0;
                                my @data = split(/_/, $info); #start, step, shift
                                if ($degPos == 1) {
                                    $infoMax = $max{$dG}{$info};
                                }
                                foreach my $match (@alignment) {
                                    if ($match eq "|") {
                                        if (($pos >= $data[0]) && ($pos <= ($data[0] + $data[1]))) {
                                        } else {
                                            $alignment[$pos] = "¦";
                                        }
                                    }
                                    $pos++;
                                }
                                $max++;
                            }
                        }
                    } else {
                        last;
                    }
                }
                if ($max_dG < 0) {
                    if (!(defined($done{$infoMax}))) { #check that it has not already be saved
                        $passPrint++;
                        push @{$final_max{$max_dG}{$passPrint}{'seq1'}}, @seq1;
                        push @{$final_max{$max_dG}{$passPrint}{'al'}}, @alignment;
                        push @{$final_max{$max_dG}{$passPrint}{'seq2'}}, @rev2;
                        $done{$infoMax} = "";
                    }
                }
            }
        }
    }

    my $max = 0;
    my $dG_def = 0;
    my $inside = 0;
    foreach my $dG (sort {$a <=> $b} keys %final_max) {
        $inside = 1;
        foreach my $passPrint (keys %{$final_max{$dG}}) {
            if ($max == 0) {
                $dG_def = $dG;
                print $file1 "$dG_def\n";
                $max++;
            }
            print $file "\n1:\t@{$final_max{$dG}{$passPrint}{'seq1'}}\n";
            print $file "A:\t@{$final_max{$dG}{$passPrint}{'al'}}\n";
            print $file "2:\t@{$final_max{$dG}{$passPrint}{'seq2'}}\n";
            print $file "dG:\t$dG\n";
        }
    }
    print $file "\n@\t$dG_def\n";
    
    if ($inside == 0) {
        print $file1 "\n";
    }
    
}

`rm $folder/tmp/$file_primer`;

#subroutines
sub degenerateAlt { #if $primer_input has degenerate bases I need to retrieve all the possible alternatives
    my ($primer) = $_[0];
    
    my %degenerate;
    my %degenerateNew;
    my $count = 0;
    my $inside = 0;
    my @all;
    my $primer0;
    
    if ($primer =~ /[RYSWKMBDHVN]/) {
        my @each = split(//, $primer);
        foreach my $e (@each) {
            if (defined($wildcard{$e})) {
                undef %degenerateNew;
                foreach my $w (keys %{$wildcard{$e}}) {
                    $count++;
                    if ($inside == 0) {
                        $degenerate{$count} = $primer0 . $w;
                        $degenerateNew{$count} = $degenerate{$count};
                    } else {
                        foreach my $c (keys %degenerate) {
                            $count++;
                            $degenerateNew{$count} = $degenerate{$c} . $w;
                        }
                    }
                }
                undef %degenerate;
                foreach my $c (keys %degenerateNew) {
                    $degenerate{$c} = $degenerateNew{$c};
                }
                undef %degenerateNew;
                $inside++;
            } else {
                if ($inside == 0) {
                    $primer0 .= $e;
                } elsif ($inside == 1) {
                    foreach my $c (keys %degenerate) {
                        $degenerate{$c} = $degenerate{$c} . $e;
                    }
                    undef %degenerateNew;
                } else {
                    foreach my $c (keys %degenerate) {
                        $degenerateNew{$c} = $degenerate{$c} . $e;
                    }
                    undef %degenerate;
                    foreach my $c (keys %degenerateNew) {
                        $degenerate{$c} = $degenerateNew{$c};
                    }
                    undef %degenerateNew;
                }
            }
        }
    } else {
        $degenerate{'1'} = $primer;
    }
    foreach my $c (keys %degenerate) {
        push @all, $degenerate{$c};
    }
    my $allRef = \@all;
    return($allRef);
}

sub checkCalc { #retrieve alignment portion - also considering of dangling ends if necessary
    my ($seq) = $_[0];
    my ($rev) = $_[1];
    my ($direction) = $_[2];
    my @data = split(/_/, $_[3]); #start, step, shift
    
    my $check1_start = $data[0];
    my $check1_end = $check1_start + $data[1];
    my $check2_start;
    if ($direction eq "left") {
        $check2_start = $data[2] + $data[0] + 1;
    } elsif ($direction eq "right") {
        $check2_start = $data[0] - $data[2] - 1;
    } elsif ($direction eq "no") {
        $check2_start = $data[0];
    }
    my $check2_end = $check2_start + $data[1];
    
    my $check1 = substr $seq, $check1_start, $data[1]; #portion alignment F
    my $check2 = substr $rev, $check2_start, $data[1]; #portion alignment R
    
    my $s1;
    if ($check1_start == 0) {
        $s1 = "";
    } else {
        $s1 = substr($seq, ($check1_start-1), 1); #base before seq
    }
    my $e1 = substr($seq, $check1_end, 1); #base after seq
    my $s2;
    if ($check2_start == 0) {
        $s2 = "";
    } else {
        $s2 = substr($rev, ($check2_start-1), 1); #base before rev
    }
    my $e2 = substr($rev, $check2_end, 1); #base after rev
    
    ###Dangling ends
    if (($s1 eq "") && ($s2 ne "")) {
        $check1 = "X" . $check1;
        $check2 = $s2 . $check2;
    } elsif (($s1 ne "") && ($s2 eq "")) {
        $check1 = $s1 . $check1;
        $check2 = "X" . $check2;
    }
    if (($e1 eq "") && ($e2 ne "")) {
        $check1 = $check1 . "X";
        $check2 = $check2 . $e2;
    } elsif (($e1 ne "") && ($e2 eq "")) {
        $check1 = $check1 . $e1;
        $check2 = $check2 . "X";
    }
    return ($check1,$check2);
}
