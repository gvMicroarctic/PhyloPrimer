#!/usr/bin/perl

use strict;

#connected to Primer_design_pp.pl, Probe_design_pp.pl, Oligo_design_pp.pl, Primer_check_pp.pl, Probe_check_pp.pl and Oligo_check_pp.pl : hairpin calculation

#Usage: ./Hairpin_checking_web_pp.pl -primer -mg -mon -oligo -dntp -t

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

    my $degPos = 0;
    if ($primer_input =~ /[RYSWKMBDHVN]/) {
        $degPos = 1;
    }

    #check if degenerancy Seq
    my ($primerSeqRef) = degenerateAlt($primer_input);
    my @primerSeq = @{$primerSeqRef};

    my %final_max;
    my %done;
    my $passPrint = 0;

    print $file "\n>\t$primer_input\n";

    foreach my $seq (@primerSeq) {
        #even loops
        foreach my $shift (1..($len-1)) {
            my $trim = $len - $shift;
            my $up = substr $seq, 0, $trim; #portion alignment F
            my $down = substr $seq, $trim, $shift; #portion alignment F
            $down = reverse($down);
            my $down1;
            my @down1;
            my $up1;
            my @up1;
            if (length($up) < length($down)) {
                my $diff = abs(length($down) - length($up));
                my $add = "X" x $diff;
                $up1 =  $add . $up;
                @up1 = split /(?=[A-Z])/i, $up1;
                $down1 = $down;
                @down1 = split /(?=[A-Z])/i, $down1;
            } else {
                my $diff = abs(length($down) - length($up));
                my $add = "X" x $diff;
                $down1 = $add . $down;
                @down1 = split /(?=[A-Z])/i, $down1;
                $up1 = $up;
                @up1 = split /(?=[A-Z])/i, $up1;
            }

            my $count=0;
            my $consecutive = 0;
            my $consecutive_gap = 0;
            my $gap = 0;
            my $num = 0;
            my $inside = 0;
            my $more = 0;
            my $time = 0;
            my %small = ();
            my @alignment = ();
            my @spaces = ();
            my %max = ();
            my $start;
            my $step;
            foreach my $position (0..(length($down1)-1)) {
                if ($up1[$position] eq $corr{$down1[$position]}) {
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
                            $small{$num} = $start . "_" . $step;
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
                $small{$num} = $start . "_" . $step;
                if ($step > 1) {
                    $more = 1;
                }
            }
            my $deg = 0;
            if (($num > 0) && ($more == 1)) {
                foreach my $num (sort keys %small) {
                    my @data = split(/_/, $small{$num}); #start, step, shift
                    my $check1 = substr $up1, $data[0], $data[1]; #portion alignment F
                    my $check2 = substr $down1, $data[0], $data[1]; #portion alignment R
                    if (length($check1) > 1) {
                        #Dangling ends
                        if ($data[0] != 0) {
                            my $s1 = substr($up1, ($data[0]-1), 1); #base before seq
                            my $s2 = substr($down1, ($data[0]-1), 1); #base before rev
                            if (($s1 eq "X") && ($s2 ne "X")) {
                                $check1 = "X" . $check1;
                                $check2 = $s2 . $check2;
                            } elsif (($s1 ne "X") && ($s2 eq "X")) {
                                $check1 = $s1 . $check1;
                                $check2 = "X" . $check2;
                            }
                        }
                        #Loop
                        my $loop1 = substr $up1, ($data[0]+$data[1]); #loop base F
                        my $loop2 = substr $down1, ($data[0]+$data[1]); #loop base R
                        my $loop3 = reverse($loop2);
                        my $loop = $loop1 . $loop3;
                        my $loop_len = length($loop);

                        if ($loop_len > 2) {
                            my $dG = `$path_cgi/DeltaG_calculation_hairpin_pp.pl -primerF $check1 -primerR $check2 -loop $loop -loopLen $loop_len -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot`;
                            chomp($dG);
                            $dG = sprintf("%.2f", $dG);
    #                        print "beg\t$check1\t$check2\t$loop\t$loop_len\t$dG\n";
                            $max{$dG}{$small{$num}} = $small{$num} . '-' . $check1 . '-' . $check2 . '-' . $loop . '-' . $dG . '-' . $shift;
                        }
                    }
                }
            }
            my $max = 0;
            my $dG_def;
            my $pos;
            my $infoMax;
            foreach my $dG (sort {$a <=> $b} keys %max) {
                if ($max == 0) {

                    foreach my $info (sort {$a <=> $b} keys %{$max{$dG}}) {
                        chomp($dG);
                        if ($max == 0) {
                            $dG_def = $dG;
                            $pos = 0;
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
            $pos = 0;
            foreach my $u (@up1) {
                if ($u eq 'X') {
                    $up1[$pos] = " ";
                }
                $pos++;
            }
            $pos = 0;
            foreach my $d (@down1) {
                if ($d eq 'X') {
                    $down1[$pos] = " ";
                }
                $pos++;
            }
            if ($dG_def < 0) {
                if (!(defined($done{$infoMax}))) {
                    $passPrint++;
                    push @up1, "\\";
                    push @down1, "\/";
                    push @{$final_max{$dG_def}{$passPrint}{'seq1'}}, @up1;
                    push @{$final_max{$dG_def}{$passPrint}{'al'}}, @alignment;
                    push @{$final_max{$dG_def}{$passPrint}{'seq2'}}, @down1;
                    $done{$infoMax} = "";
                }
            }
        }

        #uneven loops
        foreach my $shift (1..($len-1)) {
            my $trim = $len - $shift;
            my $up = substr $seq, 0, $trim; #portion alignment F
            my $middle = substr $seq, $trim, 1; #portion alignment F
            my $down = substr $seq, ($trim+1), $shift; #portion alignment F
            $down = reverse($down);
            my $down1;
            my @down1;
            my $up1;
            my @up1;
            if (length($up) < length($down)) {
                my $diff = abs(length($down) - length($up));
                my $add = "X" x $diff;
                $up1 =  $add . $up;
                @up1 = split /(?=[A-Z])/i, $up1;
                $down1 = $down;
                @down1 = split /(?=[A-Z])/i, $down1;
            } else {
                my $diff = abs(length($down) - length($up));
                my $add = "X" x $diff;
                $down1 = $add . $down;
                @down1 = split /(?=[A-Z])/i, $down1;
                $up1 = $up;
                @up1 = split /(?=[A-Z])/i, $up1;
            }

            my $count=0;
            my $consecutive = 0;
            my $consecutive_gap = 0;
            my $gap = 0;
            my $num = 0;
            my $inside = 0;
            my $more = 0;
            my $time = 0;
            my %small = ();
            my @alignment = ();
            my @spaces = ();
            my %max = ();
            my $start;
            my $step;
            foreach my $position (0..(length($down1)-1)) {
                if ($up1[$position] eq $corr{$down1[$position]}) {
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
                            $small{$num} = $start . "_" . $step;
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
                $small{$num} = $start . "_" . $step;
                if ($step > 1) {
                    $more = 1;
                }
            }
            my $deg = 0;
            if (($num > 0) && ($more == 1)) {
                foreach my $num (sort keys %small) {
                    my @data = split(/_/, $small{$num}); #start, step, shift
                    my $check1 = substr $up1, $data[0], $data[1]; #portion alignment F
                    my $check2 = substr $down1, $data[0], $data[1]; #portion alignment R
                    if (length($check1) > 1) {
                        #Dangling ends
                        if ($data[0] != 0) {
                            my $s1 = substr($up1, ($data[0]-1), 1); #base before seq
                            my $s2 = substr($down1, ($data[0]-1), 1); #base before rev
                            if (($s1 eq "X") && ($s2 ne "X")) {
                                $check1 = "X" . $check1;
                                $check2 = $s2 . $check2;
                            } elsif (($s1 ne "X") && ($s2 eq "X")) {
                                $check1 = $s1 . $check1;
                                $check2 = "X" . $check2;
                            }
                        }
                        #Loop
                        my $loop1 = substr $up1, ($data[0]+$data[1]); #loop base F
                        my $loop2 = substr $down1, ($data[0]+$data[1]); #loop base R
                        my $loop3 = reverse($loop2);
                        my $loop = $loop1 . $middle . $loop3;
                        my $loop_len = length($loop);

                        if ($loop_len > 2) {
                            my $dG = `$path_cgi/DeltaG_calculation_hairpin_pp.pl -primerF $check1 -primerR $check2 -loop $loop -loopLen $loop_len -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot`;
                            chomp($dG);
                            $dG = sprintf("%.2f", $dG);
                            $max{$dG}{$small{$num}} = $small{$num} . '-' . $check1 . '-' . $check2 . '-' . $loop . '-' . $dG . '-' . $shift;
                        }
                    }
                }
            }
            my $max = 0;
            my $dG_def;
            my $pos;
            my $infoMax;
            foreach my $dG (sort {$a <=> $b} keys %max) {
                if ($max == 0) {

                    foreach my $info (sort {$a <=> $b} keys %{$max{$dG}}) {
                        chomp($dG);
                        if ($max == 0) {
                            $dG_def = $dG;
                            $pos = 0;
                            if ($degPos == 1) {
                                $infoMax = $max{$dG}{$info};
                            }
                            my @data = split(/_/, $info); #start, step, shift
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
                }
            }
            $pos = 0;
            foreach my $u (@up1) {
                if ($u eq 'X') {
                    $up1[$pos] = " ";
                }
                $pos++;
            }
            $pos = 0;
            foreach my $d (@down1) {
                if ($d eq 'X') {
                    $down1[$pos] = " ";
                }
                $pos++;
            }
            if ($dG_def < 0) {
                if (!(defined($done{$infoMax}))) {
                    $passPrint++;
                    push @up1, "\\";
                    push @alignment, " ";
                    push @alignment, $middle;
                    push @down1, "\/";
                    push @{$final_max{$dG_def}{$passPrint}{'seq1'}}, @up1;
                    push @{$final_max{$dG_def}{$passPrint}{'al'}}, @alignment;
                    push @{$final_max{$dG_def}{$passPrint}{'seq2'}}, @down1;
                    $done{$infoMax} = "";
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
            print $file "\n@{$final_max{$dG}{$passPrint}{'seq1'}}\n";
            print $file "@{$final_max{$dG}{$passPrint}{'al'}}\n";
            print $file "@{$final_max{$dG}{$passPrint}{'seq2'}}\n";
            print $file "dG:\t$dG\tkcal/mol\n";
        }
    }
    print $file "\n@\t$dG_def\tkcal/mol\n";

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
