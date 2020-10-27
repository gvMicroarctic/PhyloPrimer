#!/usr/bin/perl

#connected to Cross_dimer_checking_web_pp.pl, Cross_dimer_checking_self_web_pp.pl, Cross_dimer_checking_probe_web_pp.pl and Self_dimer_checking_web_pp.pl : dG calculation for self- and cross- dimers

#Usage: ./DeltaG_calculation.pl -primerF -primerR -type -mg -mon -oligo -dntp -t

my %args = @ARGV;

my $seq1 = $args{-primerF};
my $seq2 = $args{-primerR};
my $type = $args{-type};
my $mg_tot = $args{-mg}; #Mg ion concentration
my $monovalent = $args{-mon}; #monovalent ion concentration
my $C = $args{-oligo}; #oligo concentration
my $dNTP_tot = $args{-dntp}; #dNTP concentration
my $temperature_celsius = $args{-t}; #temperature

#give option if self-dimer or between same primer
#A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics

#standard parameters
$R = 1.9872; #cal/K x mol
$x = 4;
$Ka = 30000;

$temperature = $temperature_celsius + 273,15;

#NN
#ΔH° kcal/mol
$dH{AA_TT} = -7.6;
$dH{AT_TA} = -7.2;
$dH{CA_GT} = -8.5;
$dH{GT_CA} = -8.4;
$dH{CT_GA} = -7.8;
$dH{GA_CT} = -8.2;
$dH{CG_GC} = -10.6;
$dH{GC_CG} = -9.8;
$dH{GG_CC} = -8.0;
#any case
$dH{initial} = 0.2;
#only if A_T or T_A
$dH{terminal} = 2.2;
#symmetry correction
$dH{symmetry} = 0;

#ΔS° cal/k·mol
$dS{AA_TT} = -21.3;
$dS{AT_TA} = -20.4;
$dS{TA_AT} = -21.3;
$dS{CA_GT} = -22.7;
$dS{GT_CA} = -22.4;
$dS{CT_GA} = -21.0;
$dS{GA_CT} = -22.2;
$dS{CG_GC} = -27.2;
$dS{GC_CG} = -24.4;
$dS{GG_CC} = -19.9;
#any case
$dS{initial} = -5.7;
#only if A_T or T_A
$dS{terminal} = 6.9;
#symmetry correction
$dS{symmetry} = -1.4;

##mismatches A.A and G.G and C.C and T.T
#ΔH° kcal/mol
$dH{AA_TA} = 1.2;
$dH{CA_GA} = -0.9;
$dH{GA_CA} = -2.9;
$dH{TA_AA} = 4.7;
$dH{AC_TC} = 0;
$dH{CC_GC} = -1.5;
$dH{GC_CC} = 3.6;
$dH{TC_AC} = 6.1;
$dH{AG_TG} = -3.1;
$dH{CG_GG} = -4.9;
$dH{GG_CG} = -6;
$dH{TG_AG} = 1.6;
$dH{AT_TT} = -2.7;
$dH{CT_GT} = -5;
$dH{GT_CT} = -2.2;
$dH{TT_AT} = 0.2;
#ΔS° cal/k·mol
$dS{AA_TA} = 1.7;
$dS{CA_GA} = -4.2;
$dS{GA_CA} = -9.8;
$dS{TA_AA} = 12.9;
$dS{AC_TC} = -4.4;
$dS{CC_GC} = -7.2;
$dS{GC_CC} = 8.9;
$dS{TC_AC} = 16.4;
$dS{AG_TG} = -9.5;
$dS{CG_GG} = -15.3;
$dS{GG_CG} = -15.8;
$dS{TG_AG} = 3.6;
$dS{AT_TT} = -10.8;
$dS{CT_GT} = -15.8;
$dS{GT_CT} = -8.4;
$dS{TT_AT} = -1.5;

#mismatches G.T
#ΔH° kcal/mol
$dH{AG_TT} = 1;
$dH{AT_TG} = -2.5;
$dH{CG_GT} = -4.1;
$dH{CT_GG} = -2.8;
$dH{GG_CT} = 3.3;
$dH{GT_CG} = -4.4;
$dH{TG_AT} = -0.1;
$dH{TT_AG} = -1.3;
#ΔS° cal/k·mol
$dS{AG_TT} = 0.9;
$dS{AT_TG} = -8.3;
$dS{CG_GT} = -11.7;
$dS{CT_GG} = -8;
$dS{GG_CT} = 10.4;
$dS{GT_CG} = -12.3;
$dS{TG_AT} = -1.7;
$dS{TT_AG} = -5.3;

#mismatches A . C
#ΔH° kcal/mol
$dH{AA_TC} = 2.3;
$dH{AC_TA} = 5.3;
$dH{CA_GC} = 1.9;
$dH{CC_GA} = 0.6;
$dH{GA_CC} = 5.2;
$dH{GC_CA} = -0.7;
$dH{TA_AC} = 3.4;
$dH{TC_AA} = 7.6;
#ΔS° cal/k·mol
$dS{AA_TC} = 4.6;
$dS{AC_TA} = 14.6;
$dS{CA_GC} = 3.7;
$dS{CC_GA} = -0.6;
$dS{GA_CC} = 14.2;
$dS{GC_CA} = -3.8;
$dS{TA_AC} = 8;
$dS{TC_AA} = 20.2;

#mismatches G . A
#ΔH° kcal/mol
$dH{AA_TG} = -0.6;
$dH{AG_TA} = -0.7;
$dH{CA_GG} = -0.7;
$dH{CG_GA} = -4;
$dH{GA_CG} = -0.6;
$dH{GG_CA} = 0.5;
$dH{TA_AG} = 0.7;
$dH{TG_AA} = 3;
#ΔS° cal/k·mol
$dS{AA_TG} = -2.3;
$dS{AG_TA} = -2.3;
$dS{CA_GG} = -2.3;
$dS{CG_GA} = -13.2;
$dS{GA_CG} = -1;
$dS{GG_CA} = 3.2;
$dS{TA_AG} = 0.7;
$dS{TG_AA} = 7.4;

#mismatches C . T
#ΔH° kcal/mol
$dH{AC_TT} = 0.7;
$dH{AT_TC} = -1.2;
$dH{CC_GT} = -0.8;
$dH{CT_GC} = -1.5;
$dH{GC_CT} = 2.3;
$dH{GT_CC} = 5.2;
$dH{TC_AT} = 1.2;
$dH{TT_AC} = 1;
#ΔS° cal/k·mol
$dS{AC_TT} = 0.2;
$dS{AT_TC} = -6.2;
$dS{CC_GT} = -4.5;
$dS{CT_GC} = -6.1;
$dS{GC_CT} = 5.4;
$dS{GT_CC} = 13.5;
$dS{TC_AT} = 0.7;
$dS{TT_AC} = 0.7;

#dangling ends
#TABLE2
#ΔH° kcal/mol
$dH{AA_XT} = 0.2;
$dH{AC_XG} = -6.3;
$dH{AG_XC} = -3.7;
$dH{AT_XA} = -2.9;
$dH{CA_XT} = 0.6;
$dH{CC_XG} = -4.4;
$dH{CG_XC} = -4;
$dH{CT_XA} = -4.1;
$dH{GA_XT} = -1.1;
$dH{GC_XG} = -5.1;
$dH{GG_XC} = -3.9;
$dH{GT_XA} = -4.2;
$dH{TA_XT} = -6.9;
$dH{TC_XG} = -4;
$dH{TG_XC} = -4.9;
$dH{TT_XA} = -0.2;
$dH{XA_AT} = -0.7;
$dH{XC_AG} = -2.1;
$dH{XG_AC} = -5.9;
$dH{XT_AA} = -0.5;
$dH{XA_CT} = 4.4;
$dH{XC_CG} = -0.2;
$dH{XG_CC} = -2.6;
$dH{XT_CA} = 4.7;
$dH{XA_GT} = -1.6;
$dH{XC_GG} = -3.9;
$dH{XG_GC} = -3.2;
$dH{XT_GA} = -4.1;
$dH{XA_TT} = 2.9;
$dH{XC_TG} = -4.4;
$dH{XG_TC} = -5.2;
$dH{XT_TA} = -3.8;
#ΔS° cal/k·mol
$dS{AA_XT} = 2.3;
$dS{AC_XG} = -17.1;
$dS{AG_XC} = -10;
$dS{AT_XA} = -7.6;
$dS{CA_XT} = 3.3;
$dS{CC_XG} = -12.6;
$dS{CG_XC} = -11.9;
$dS{CT_XA} = -13;
$dS{GA_XT} = -1.6;
$dS{GC_XG} = -14;
$dS{GG_XC} = -10.9;
$dS{GT_XA} = -15;
$dS{TA_XT} = -20;
$dS{TC_XG} = -10.9;
$dS{TG_XC} = -13.8;
$dS{TT_XA} = -0.5;
$dS{XA_AT} = -0.8;
$dS{XC_AG} = -3.9;
$dS{XG_AC} = -16.5;
$dS{XT_AA} = -1.1;
$dS{XA_CT} = 14.9;
$dS{XC_CG} = -0.1;
$dS{XG_CC} = -7.4;
$dS{XT_CA} = 14.2;
$dS{XA_GT} = -3.6;
$dS{XC_GG} = -11.2;
$dS{XG_GC} = -10.4;
$dS{XT_GA} = -13.1;
$dS{XA_TT} = 10.4;
$dS{XC_TG} = -13.1;
$dS{XG_TC} = -15;
$dS{XT_TA} = -12.6;

$gc_base=0;
while ($seq1 =~ /[GC]/g) {
    $gc_base++;
}
$len = length($seq1);
$gc = $gc_base/$len;

if ($dNTP_tot < (0.8 * $mg_tot)) {
    $mg = $mg_tot - $dNTP_tot;
} else {
    $mg = (-($Ka*$dNTP_tot -$Ka*$mg_tot + 1) + sqrt(($Ka*$dNTP_tot -$Ka*$mg_tot + 1)**2 + 4 * $Ka * $mg_tot))/(2*$Ka);
}

# set which correction is needed
$a = 3.92 * 10**-5;
$b = -9.11 * 10**-6;
$c = 6.26 * 10**-5;
$d = 1.42 * 10**-5;
$e = -4.82 * 10**-4;
$f = 5.25 * 10**-4;
$g = 8.31 * 10**-5;

if ($monovalent == 0) {
    $D = $a + ($b*log($mg)) + $gc*($c+$d*log($mg)) + (1/(2*($len - 1)))*($e+($f*log($mg))+$g*((log($mg)**2)));
} else {
    $ratio = sqrt($mg)/$monovalent;
#    print "ratio\t$ratio\n";
    if ($ratio < 0.22) {
        $D = ((4.29*$gc-3.95) * 10**-5 * log($monovalent)) + (9.4 * 10**-6 * (log($monovalent)**2));
    } elsif ($ratio < 6) {
        $D = $a + ($b*log($mg)) + $gc*($c+$d*log($mg)) + (1/(2*($len - 1)))*($e+($f*log($mg))+$g*((log($mg)**2)));
    } else {
        $a = 3.92 * 10**-5 * (0.843-0.352*sqrt($monovalent)*log($monovalent));
        $d = 1.42 * 10**-5 * (1.279-0.00403*log($monovalent)-0.00803*(log($monovalent)**2));
        $g = 8.31 * 10**-5 * (0.486-0.258*log($monovalent)+0.00525*(log($monovalent)**3));
        $D = $a + ($b*log($mg)) + $gc*($c+$d*log($mg)) + (1/(2*($len - 1)))*($e+($f*log($mg))+$g*((log($mg)**2)));
    }
}

#symmetry
if ($type eq 'self') { #symmetry only if self-dimers
    $sum_dS += $dS{symmetry};
    $sum_dH += $dH{symmetry};
}

#initial bases
$sum_dS += $dS{initial};
$sum_dH += $dH{initial};

#terminal bases corrected if A_T ot T_A
$start_1 =  substr($seq1, 0, 1);
$start_2 =  substr($seq2, 0, 1);
$pair = $start_1 . "_" . $start_2;
if (($pair eq 'A_T') or ($pair eq 'T_A')) {
    $sum_dS += $dS{terminal};
    $sum_dH += $dH{terminal};
}

$end_1 =  substr($seq1, -1);
$end_2 =  substr($seq2, -1);
$pair = $end_1 . "_" . $end_2;
if (($pair eq 'A_T') or ($pair eq 'T_A')) {
    $sum_dS += $dS{terminal};
    $sum_dH += $dH{terminal};
}

$offset = 0;
while ($offset <= ($len - 2)) {
    $sub_1 = substr($seq1, $offset, 2);
    $sub_2 = substr($seq2, $offset, 2);
    $pair = $sub_1 . "_" . $sub_2;
    if (!(defined($dS{$pair}))) {
        $pair = reverse($pair);
    }
    $sum_dS += $dS{$pair};
    $sum_dH += $dH{$pair};
    $sing = $dH{$pair} - $temperature*($dS_corrected/1000);
    $offset++;
}

$dS_corrected = $sum_dS + ($sum_dH*$D);

$dG = $sum_dH - $temperature*($dS_corrected/1000);

print "$dG\n";
