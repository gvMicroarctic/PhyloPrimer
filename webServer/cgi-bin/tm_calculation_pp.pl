#!/usr/bin/perl

#connected to Primer_design_pp.pl, Probe_design_pp.pl, Oligo_design_pp.pl, Primer_check_pp.pl, Probe_check_pp.pl and Oligo_check_pp.pl : melting temperature calculation

#Usage: ./tm_calculation -primer -mg -mon -oligo -dntp

my %args = @ARGV;

$primer_input = $args{-primer};
$sense = $args{-sense}; #F or R
$dangling = $args{-dang}; #dangling base
$type = $args{-type}; #primer or probe
$mg_tot = $args{-mg}; #Mg ion concentration
$monovalent = $args{-mon}; #monovalent ion concentration
$C = $args{-oligo}; #oligo concentration
$dNTP_tot = $args{-dntp}; #dNTP concentration

#fixed parameters
$R = 1.9872; #cal/K x mol
$x = 4;
$Ka = 30000;

#ΔH° kcal/mol
$dH{AA_TT} = -7.6;
$dH{TT_AA} = -7.6;

$dH{AT_TA} = -7.2;

$dH{TA_AT} = -7.2;

$dH{CA_GT} = -8.5;
$dH{TG_AC} = -8.5;

$dH{GT_CA} = -8.4;
$dH{AC_TG} = -8.4;

$dH{CT_GA} = -7.8;
$dH{AG_TC} = -7.8;

$dH{GA_CT} = -8.2;
$dH{TC_AG} = -8.2;

$dH{CG_GC} = -10.6;

$dH{GC_CG} = -9.8;

$dH{GG_CC} = -8.0;
$dH{CC_GG} = -8.0;

#any case
$dH{initial} = 0.2;

#only if A_T or T_A
$dH{terminal} = 2.2;

#ΔS° cal/k·mol
$dS{AA_TT} = -21.3;
$dS{TT_AA} = -21.3;

$dS{AT_TA} = -20.4;

$dS{TA_AT} = -21.3;

$dS{CA_GT} = -22.7;
$dS{TG_AC} = -22.7;

$dS{GT_CA} = -22.4;
$dS{AC_TG} = -22.4;

$dS{CT_GA} = -21.0;
$dS{AG_TC} = -21.0;

$dS{GA_CT} = -22.2;
$dS{TC_AG} = -22.2;

$dS{CG_GC} = -27.2;

$dS{GC_CG} = -24.4;

$dS{GG_CC} = -19.9;
$dS{CC_GG} = -19.9;

#any case
$dS{initial} = -5.7;

#only if A_T or T_A
$dS{terminal} = 6.9;

#dangling end
#ΔH° kcal/mol
$dH{AA_XT} = 0.2;
$dH{TX_AA} = 0.2;

$dH{AC_XG} = -6.3;
$dH{GX_CA} = -6.3;

$dH{AG_XC} = -3.7;
$dH{CX_GA} = -3.7;

$dH{AT_XA} = -2.9;
$dH{AX_TA} = -2.9;

$dH{CA_XT} = 0.6;
$dH{TX_AC} = 0.6;

$dH{CC_XG} = -4.4;
$dH{GX_CC} = -4.4;

$dH{CG_XC} = -4;
$dH{CX_GC} = -4;

$dH{CT_XA} = -4.1;
$dH{AX_TC} = -4.1;

$dH{GA_XT} = -1.1;
$dH{TX_AG} = -1.1;

$dH{GC_XG} = -5.1;
$dH{GX_CG} = -5.1;

$dH{GG_XC} = -3.9;
$dH{CX_GG} = -3.9;

$dH{GT_XA} = -4.2;
$dH{AX_TG} = -4.2;

$dH{TA_XT} = -6.9;
$dH{TX_AT} = -6.9;

$dH{TC_XG} = -4;
$dH{GX_CT} = -4;

$dH{TG_XC} = -4.9;
$dH{CX_GT} = -4.9;

$dH{TT_XA} = -0.2;
$dH{AX_TT} = -0.2;

$dH{XA_AT} = -0.7;
$dH{TA_AX} = -0.7;

$dH{XC_AG} = -2.1;
$dH{GA_CX} = -2.1;

$dH{XG_AC} = -5.9;
$dH{CA_GX} = -5.9;

$dH{XT_AA} = -0.5;
$dH{AA_TX} = -0.5;

$dH{XA_CT} = 4.4;
$dH{TC_AX} = 4.4;

$dH{XC_CG} = -0.2;
$dH{GC_CX} = -0.2;

$dH{XG_CC} = -2.6;
$dH{CC_GX} = -2.6;

$dH{XT_CA} = 4.7;
$dH{AC_TX} = 4.7;

$dH{XA_GT} = -1.6;
$dH{TG_AX} = -1.6;

$dH{XC_GG} = -3.9;
$dH{GG_CX} = -3.9;

$dH{XG_GC} = -3.2;
$dH{CG_GX} = -3.2;

$dH{XT_GA} = -4.1;
$dH{AG_TX} = -4.1;

$dH{XA_TT} = 2.9;
$dH{TT_AX} = 2.9;

$dH{XC_TG} = -4.4;
$dH{GT_CX} = -4.4;

$dH{XG_TC} = -5.2;
$dH{CT_GX} = -5.2;

$dH{XT_TA} = -3.8;
$dH{AT_TX} = -3.8;

#ΔS° cal/k·mol
$dS{AA_XT} = 2.3;
$dS{TX_AA} = 2.3;

$dS{AC_XG} = -17.1;
$dS{GX_CA} = -17.1;

$dS{AG_XC} = -10;
$dS{CX_GA} = -10;

$dS{AT_XA} = -7.6;
$dS{AX_TA} = -7.6;

$dS{CA_XT} = 3.3;
$dS{TX_AC} = 3.3;

$dS{CC_XG} = -12.6;
$dS{GX_CC} = -12.6;

$dS{CG_XC} = -11.9;
$dS{CX_GC} = -11.9;

$dS{CT_XA} = -13;
$dS{AX_TC} = -13;

$dS{GA_XT} = -1.6;
$dS{TX_AG} = -1.6;

$dS{GC_XG} = -14;
$dS{GX_CG} = -14;

$dS{GG_XC} = -10.9;
$dS{CX_GG} = -10.9;

$dS{GT_XA} = -15;
$dS{AX_TG} = -15;

$dS{TA_XT} = -20;
$dS{TX_AT} = -20;

$dS{TC_XG} = -10.9;
$dS{GX_CT} = -10.9;

$dS{TG_XC} = -13.8;
$dS{CX_GT} = -13.8;

$dS{TT_XA} = -0.5;
$dS{AX_TT} = -0.5;

$dS{XA_AT} = -0.8;
$dS{TA_AX} = -0.8;

$dS{XC_AG} = -3.9;
$dS{GA_CX} = -3.9;

$dS{XG_AC} = -16.5;
$dS{CA_GX} = -16.5;

$dS{XT_AA} = -1.1;
$dS{AA_TX} = -1.1;

$dS{XA_CT} = 14.9;
$dS{TC_AX} = 14.9;

$dS{XC_CG} = -0.1;
$dS{GC_CX} = -0.1;

$dS{XG_CC} = -7.4;
$dS{CC_GX} = -7.4;

$dS{XT_CA} = 14.2;
$dS{AC_TX} = 14.2;

$dS{XA_GT} = -3.6;
$dS{TG_AX} = -3.6;

$dS{XC_GG} = -11.2;
$dS{GG_CX} = -11.2;

$dS{XG_GC} = -10.4;
$dS{CG_GX} = -10.4;

$dS{XT_GA} = -13.1;
$dS{AG_TX} = -13.1;

$dS{XA_TT} = 10.4;
$dS{TT_AX} = 10.4;

$dS{XC_TG} = -13.1;
$dS{GT_CX} = -13.1;

$dS{XG_TC} = -15;
$dS{CT_GX} = -15;

$dS{XT_TA} = -12.6;
$dS{AT_TX} = -12.6;

my %wildcard;
#wildcards
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

my $len = length($primer_input);
$len_dg = $len;


#if $primer_input has degenerate bases I need to retrieve all the possible alternatives
if ($primer_input =~ /[RYSWKMBDHVN]/) {
    my @each = split(//, $primer_input);
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
    $degenerate{'1'} = $primer_input;
}

$count = 0;

foreach my $c (keys %degenerate) {
    $primer_input = $degenerate{$c};
#    print "primer $primer_input\n";
    if ($type eq 'primer') { #print "Primer mode: $primer1 and $primer2\n";
        if ($dangling ne '') {
            if ($sense eq "F") { #F
                $primer1 = $primer_input; #forward primer
                $primer2 = $primer_input; #template
                $primer2 =~ tr/CATGYKBDRMVH/GTACRMVHYKBD/;
                $primer1 .= "X";
                if ($dangling =~ /[RYSWKMBDHVN]/) {
                    foreach my $w (keys %{$wildcard{$dangling}}) {
                        $count++;
                        $template{$count} = $primer2 . $w;
                    }
                } else {
                    $count++;
                    $template{$count} = $primer2 . $dangling;
                }
                $primer2 = "";
                $len_dg += 1;
            } else { #R
                $primer2 = reverse($primer_input); #reverse primer
                $primer1 = $primer2;
                $primer1 =~ tr/CATGYKBDRMVH/GTACRMVHYKBD/; #template
                if ($dangling =~ /[RYSWKMBDHVN]/) {
                    foreach my $w (keys %{$wildcard{$dangling}}) {
                        $count++;
                        $template{$count} = $w . $primer1;
                    }
                } else {
                    $count++;
                    $template{$count} = $dangling . $primer1;
                }
                $primer1 = "";
                $primer2 = "X" . $primer2;
                $len_dg += 1;
            }
        } else {
            $count++;
            if ($sense eq "F") { #F
                $primer1 = $primer_input; #forward primer
                $primer2 = $primer_input; #template
                $primer2 =~ tr/CATGYKBDRMVH/GTACRMVHYKBD/;
                $template{$count} = $primer2;
                $primer2 = "";
            } else { #R
                $primer2 = reverse($primer_input); #reverse primer
                $primer1 = $primer2;
                $primer1 =~ tr/CATGYKBDRMVH/GTACRMVHYKBD/; #template
                $template{$count} = $primer1;
                $primer1 = "";
            }
        }
    } elsif ($type eq 'probe') { #print "Probe mode: $primer1 and $primer2\n";
        if ($dangling ne '') {
            @dan = split(/,/, $dangling);
            $count++;
            if ($sense eq "F") { #F
                $primer1 = $primer_input; #forward primer
                $primer2 = $primer_input; #template
                $primer2 =~ tr/CATGYKBDRMVH/GTACRMVHYKBD/;
                $primer1 = "X" . $primer_input . "X";
                
                if ($dan[0] =~ /[RYSWKMBDHVN]/) {
                    foreach my $w (keys %{$wildcard{$dan[0]}}) {
                        $count++;
                        $template{$count} = $w . $primer2 . $dan[1];
                    }
                    
                } elsif ($dan[1] =~ /[RYSWKMBDHVN]/) {
                    foreach my $w (keys %{$wildcard{$dan[1]}}) {
                        $count++;
                        $template{$count} = $dan[0] . $primer2 . $w;
                    }
                } else {
                    $count++;
                    $template{$count} = $dan[0] . $primer2 . $dan[1];
                }
                $primer2 = "";
                $len_dg += 2;
            } else {
                $primer2 = reverse($primer_input); #reverse primer
                $primer1 = $primer2;
                $primer1 =~ tr/CATGYKBDRMVH/GTACRMVHYKBD/; #template
                if ($dan[0] =~ /[RYSWKMBDHVN]/) {
                    foreach my $w (keys %{$wildcard{$dan[0]}}) {
                        $count++;
                        $template{$count} = $w . $primer1 . $dan[1];
                    }
                } elsif ($dan[1] =~ /[RYSWKMBDHVN]/) {
                    foreach my $w (keys %{$wildcard{$dan[1]}}) {
                        $count++;
                        $template{$count} = $dan[0] . $primer1 . $w;
                    }
                } else {
                    $count++;
                    $template{$count} = $dan[0] . $primer1 . $dan[1];
                }
                $primer1 = "";
                $primer2 = "X" . $primer2 . "X";
                $len_dg += 2;
            }
        } else {
            if ($sense eq "F") {
                $primer1 = $primer_input;
                $primer2 = $primer_input;
                $primer2 =~ tr/CATGYKBDRMVH/GTACRMVHYKBD/;
                $template{$count} = $primer2;
                $primer2 = "";
            } else {
                $primer1 = $primer_input;
                $primer1 =~ tr/CATGYKBDRMVH/GTACRMVHYKBD/;
                $primer2 = $primer_input;
                $template{$count} = $primer1;
                $primer1 = "";
            }
        }
    }
    
    #loop for dangling with degenerate base
    $inside = 0;
    foreach my $t (keys %template) {
        if ($primer1 eq "") {
            $primer1 = $template{$t};
            $primer2 = $primer2;
            $inside = 1;
        } else {
            $primer1 = $primer1;
            $primer2 = $template{$t};
        }
        
        my $offset = 0;
        my $sum_dH = 0;
        my $sum_dS = 0;
        while ($offset <= ($len_dg - 2)) {
            $sub_1 = substr($primer1, $offset, 2);
            $sub_2 = substr($primer2, $offset, 2);
            $pair = $sub_1 . "_" . $sub_2;
            $sum_dH += $dH{$pair};
            $sum_dS += $dS{$pair};
            $offset++;
        }
        
        $sum_dH += $dH{initial};
        $sum_dS += $dS{initial};
        
        $terminal1 = substr($primer1,0,1);
        $terminal2 = substr($primer1,-1);
        
        $terminal3 = substr($primer2,0,1);
        $terminal4 = substr($primer2,-1);
        
        if ((($terminal1 eq "T") or ($terminal1 eq "A")) && ($terminal1 ne "X") && ($terminal3 ne "X")) {
            $sum_dH += $dH{terminal};
            $sum_dS += $dS{terminal};
        }
        
        if ((($terminal2 eq "T") or ($terminal2 eq "A")) && ($terminal2 ne "X") && ($terminal4 ne "X")) {
            $sum_dH += $dH{terminal};
            $sum_dS += $dS{terminal};
        }
        
        #Calculate GC%
        $gc_base=0;
        while ($primer_input =~ /[GC]/g) {
            $gc_base++;
        }
        $gc = $gc_base/$len;
        $T = $sum_dH*1000/($sum_dS + ($R*log($C/$x))); # in kelvin
        
        $med = $R*log($C/$x);

        #mg concentration
        if ($dNTP_tot < (0.8 * $mg_tot)) {
            $mg = $mg_tot - $dNTP_tot;
        } else {
            $mg = (-($Ka*$dNTP_tot -$Ka*$mg_tot + 1) + sqrt(($Ka*$dNTP_tot -$Ka*$mg_tot + 1)**2 + 4 * $Ka * $mg_tot))/(2*$Ka);
            #print "here1\t$mg\n";
        }

        $a = 3.92 * 10**-5;
        $b = -9.11 * 10**-6;
        $c = 6.26 * 10**-5;
        $d = 1.42 * 10**-5;
        $e = -4.82 * 10**-4;
        $f = 5.25 * 10**-4;
        $g = 8.31 * 10**-5;
        
        #Nbp - the duplex length is calculated excluding the dangling ends! It fives an extimate of phospore bonds so counts of bases part of the inner duplex!!
        
        if ($monovalent == 0) {
            $Tb = (1/$T) + $a + ($b*log($mg)) + $gc*($c+$d*log($mg)) + (1/(2*($len - 1)))*($e+($f*log($mg))+$g*((log($mg)**2)));
        } else {
            $ratio = sqrt($mg)/$monovalent;
#            print "ratio\t$ratio\n";
            if ($ratio < 0.22) {
                $Tb = (1/$T) + ((4.29*$gc-3.95) * 10**-5 * log($monovalent)) + (9.4 * 10**-6 * (log($monovalent)**2));
            } elsif ($ratio < 6) {
                $Tb = (1/$T) + $a + ($b*log($mg)) + $gc*($c+$d*log($mg)) + (1/(2*($len - 1)))*($e+($f*log($mg))+$g*((log($mg)**2)));
            } else {
                $a = 3.92 * 10**-5 * (0.843-0.352*sqrt($monovalent)*log($monovalent));
                $d = 1.42 * 10**-5 * (1.279-0.00403*log($monovalent)-0.00803*(log($monovalent)**2));
                $g = 8.31 * 10**-5 * (0.486-0.258*log($monovalent)+0.00525*(log($monovalent)**3));
                $Tb = (1/$T) + $a + ($b*log($mg)) + $gc*($c+$d*log($mg)) + (1/(2*($len - 1)))*($e+($f*log($mg))+$g*((log($mg)**2)));
            }
        }

        $Tm = (1/$Tb) - 273.15; #convert to celsius
        $TmSum += $Tm;
        $TmCount++;
        if ($inside == 1) {
            $primer1 = "";
        } else {
            $primer2 = "";
        }
    }

}

$Tm = $TmSum/$TmCount;
print "$Tm";
