#!/usr/bin/perl

use strict;

#connected to Primer_design_pp.pl, Probe_design_pp.pl, Oligo_design_pp.pl, Primer_check_pp.pl, Probe_check_pp.pl and Oligo_check_pp.pl : cross-dimer calculation

#Usage: ./Cross_dimer_checking_self_web.pl -primer_input -mg -mon -oligo -dntp -t

#set path to folders:
my $path_cgi = 'path_to_cgi_folder';

my %args = @ARGV;

my $folder = $args{-folder}; #folder
my $primer_input = $args{-primer}; #primer file
my $mg_tot = $args{-mg}; #Mg ion concentration
my $monovalent = $args{-mon}; #monovalent ion concentration
my $C = $args{-oligo}; #oligo concentration
my $dNTP_tot = $args{-dntp}; #dNTP concentration
my $temperature_celsius = $args{-t}; #temperature


my %corr;
$corr{'A'} = 'T';
$corr{'T'} = 'A';
$corr{'G'} = 'C';
$corr{'C'} = 'G';

open(IN, "<$folder/tmp/$primer_input");
my ($index) = $primer_input =~ /CrossSelf_(.*).tmp/;
open(my $file, ">$folder/tmp/SecondaryStructure_crossSelf_${index}.tmp") or die; #file with all the necessary parameters

while(defined(my $input = <IN>)) {
    chomp($input);
    my @pair = split(/\t/, $input);
    my $countPair = 0;
    
    my %final_max;
    my $passPrint = 0;
    
    print $file "\n>\t$pair[0]\n"; #print the original oligo

    shift(@pair); #remove first element array
    foreach my $p (@pair) { #each pair
        $countPair++;

        my ($for,$rev) = split(/,/, $p);

        $rev = reverse($rev);

                my %done;
    
                my $len_min;
                my $len_max;
                if (length($for) < length($rev)) {
                    $len_min = length($for);
                    $len_max = length($rev);
                } else {
                    $len_max = length($for);
                    $len_min = length($rev);
                }
                my @for1 = split /(?=[A-Z])/i, $for;
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
                foreach my $position (0..($len_min-1)) {
                    if ($for1[$position] eq $corr{$rev1[$position]}) {
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
                        my $check1_start = $data[0];
                        my $check1_end = $check1_start + $data[1];
                        my $check1 = substr $for, $check1_start, $data[1]; #portion alignment F
                        if (length($check1) > 1) {
                            my $check2_start = $data[0];
                            my $check2_end = $check2_start + $data[1];
                            my $check2 = substr $rev, $check2_start, $data[1]; #portion alignment R
                            
                            my $s1;
                            if ($check1_start == 0) {
                                $s1 = "";
                            } else {
                                $s1 = substr($for, ($check1_start-1), 1); #base before seq
                            }
                            my $e1 = substr($for, $check1_end, 1); #base after seq
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
                            my $dG = `$path_cgi/DeltaG_calculation_pp.pl -primerF $check1 -primerR $check2 -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius`;
                            chomp($dG);
                            $dG = sprintf("%.2f", $dG);
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
                                $max_dG = $dG;
                                my $pos = 0;
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
                        push @{$final_max{$max_dG}{$passPrint}{'seq1'}}, @for1;
                        push @{$final_max{$max_dG}{$passPrint}{'al'}}, @alignment;
                        push @{$final_max{$max_dG}{$passPrint}{'seq2'}}, @rev1;
                        $done{$infoMax} = "";
                    }
                }
                
                #Shift to the left
                my @rev2 = @rev1;
                my @for2 = @for1;
                foreach my $time (0..($len_max-2)) {
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
                    foreach my $position (0..($len_min-1)) {
                        if ($for2[$position] eq $corr{$rev2[$position]}) {
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
                        my @for3 = @for2;
                        unshift @for3, @spaces; #shift array to the right
                        unshift @alignment, @spaces; #shift array to the right
                        undef %max;
                        foreach my $num (sort keys %small) {
                            my @data = split(/_/, $small{$num}); #start, step, shift
                            my $check1_start = $data[0];
                            my $check1_end = $check1_start + $data[1];
                            my $check1 = substr $for, $check1_start, $data[1]; #portion alignment F
                            
                            if (length($check1) > 1) {
                                my $check2_start = $data[2] + $data[0] + 1;
                                my $check2_end = $check2_start + $data[1];
                                my $check2 = substr $rev, $check2_start, $data[1]; #portion alignment R
                                
                                my $s1;
                                if ($check1_start == 0) {
                                    $s1 = "";
                                } else {
                                    $s1 = substr($for, ($check1_start-1), 1); #base before seq
                                }
                                my $e1 = substr($for, $check1_end, 1); #base after seq
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
                                my $dG = `$path_cgi/DeltaG_calculation_pp.pl -primerF $check1 -primerR $check2 -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius`;
                                chomp($dG);
                                $dG = sprintf("%.2f", $dG);
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
                                        $max_dG = $dG;
                                        my $pos = 0;
                                        my @data = split(/_/, $info); #start, step, shift
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
                                push @{$final_max{$max_dG}{$passPrint}{'seq1'}}, @for3;
                                push @{$final_max{$max_dG}{$passPrint}{'al'}}, @alignment;
                                push @{$final_max{$max_dG}{$passPrint}{'seq2'}}, @rev1;
                                $done{$infoMax} = "";
                            }
                        }
                    }
                }
                
                #Shift to the right
                my @rev2 = @rev1;
                my @for2 = @for1;
                foreach my $time (0..($len_max-2)) {
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
                    foreach my $position (0..($len_min-1)) {
                        if ($for2[$position] eq $corr{$rev2[$position]}) {
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
                            my @data = split(/_/, $small{$num}); #start, step, shift
                            my $check1_start = $data[0];
                            my $check1_end = $check1_start + $data[1];
                            my $check1 = substr $for, $check1_start, $data[1]; #portion alignment
                            if (length($check1) > 1) {
                                my $check2_start = $data[0] - $data[2] - 1;
                                my $check2_end = $check2_start + $data[1];
                                my $check2 = substr $rev, $check2_start, $data[1]; #portion alignment R
                                
                                my $s1;
                                if ($check1_start == 0) {
                                    $s1 = "";
                                } else {
                                    $s1 = substr($for, ($check1_start-1), 1); #base before seq
                                }
                                my $e1 = substr($for, $check1_end, 1); #base after seq
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
                                my $dG = `$path_cgi/DeltaG_calculation_pp.pl -primerF $check1 -primerR $check2 -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius`;
                                chomp($dG);
                                $dG = sprintf("%.2f", $dG);
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
                                        $max_dG = $dG;
                                        my $pos = 0;
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
                                push @{$final_max{$max_dG}{$passPrint}{'seq1'}}, @for1;
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
    my %already; #no same structure twice
    foreach my $dG (sort {$a <=> $b} keys %final_max) {
        foreach my $passPrint (keys %{$final_max{$dG}}) {
            if ($max == 0) {
                $dG_def = $dG;
                $max++;
            }
            my $count1 = 0;
            foreach my $c (@{$final_max{$dG}{$passPrint}{'al'}}) {
                if ($c eq "|") {
                    $count1++;
                }
            }
            my $count2 = 0;
            foreach my $c (@{$final_max{$dG}{$passPrint}{'al'}}) {
                if ($c eq "¦") {
                    $count2++;
                }
            }
            my $al = $dG . scalar(@{$final_max{$dG}{$passPrint}{'al'}}) . $count1 . $count2;
            if (!(defined($already{$al}))) {
                print $file "\n1:\t@{$final_max{$dG}{$passPrint}{'seq1'}}\n";
                print $file "A:\t@{$final_max{$dG}{$passPrint}{'al'}}\n";
                print $file "2:\t@{$final_max{$dG}{$passPrint}{'seq2'}}\n";
                print $file "dG:\t$dG\n";
                $already{$al} = '';
            }
        }
    }
    print $file "\n@\t$dG_def\n";
    undef %already;
}
close(IN);

`rm $folder/tmp/$primer_input 2> /dev/null`;
