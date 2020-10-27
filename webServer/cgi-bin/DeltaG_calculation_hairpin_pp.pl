#!/usr/bin/perl

#connected to Hairpin_checking_web_pp.pl : dG calculation for hairpins

#Usage: ./DeltaG_calculation_hairpin.pl -primerF -primerR -loop -loopLen -mg -mon -oligo -dntp -t
my %args = @ARGV;

my $seq1 = $args{-primerF};
my $seq2 = $args{-primerR};
my $loop = $args{-loop};
my $loopLen = $args{-loopLen};
my $mg_tot = $args{-mg}; #Mg ion concentration
my $monovalent = $args{-mon}; #monovalent ion concentration
my $C = $args{-oligo}; #oligo concentration
my $dNTP_tot = $args{-dntp}; #dNTP concentration
my $temperature_celsius = $args{-t}; #temperature

#give option if self-dimer or between same primer

#A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics

$C = 0.000002;
$monovalent = 0.005;
$mg_tot = 0.0015;
$dNTP_tot = 0.00001;
$temperature_celsius = 25;

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

#triloops
#ΔH° kcal/mol
$dH{AGAAT} = -1.5;
$dH{AGCAT} = -1.5;
$dH{AGGAT} = -1.5;
$dH{AGTAT} = -1.5;
$dH{CGAAG} = -2;
$dH{CGCAG} = -2;
$dH{CGGAG} = -2;
$dH{CGTAG} = -2;
$dH{GGAAC} = -2;
$dH{GGCAC} = -2;
$dH{GGGAC} = -2;
$dH{GGTAC} = -2;
$dH{TGAAA} = -1.5;
$dH{TGCAA} = -1.5;
$dH{TGGAA} = -1.5;
$dH{TGTAA} = -1.5;
#ΔS° cal/k·mol
$dS{AGAAT} = 0;
$dS{AGCAT} = 0;
$dS{AGGAT} = 0;
$dS{AGTAT} = 0;
$dS{CGAAG} = 0;
$dS{CGCAG} = 0;
$dS{CGGAG} = 0;
$dS{CGTAG} = 0;
$dS{GGAAC} = 0;
$dS{GGCAC} = 0;
$dS{GGGAC} = 0;
$dS{GGTAC} = 0;
$dS{TGAAA} = 0;
$dS{TGCAA} = 0;
$dS{TGGAA} = 0;
$dS{TGTAA} = 0;

#tetraloops
#ΔH° kcal/mol
# I have matker with '#' all the dH that differ from dG, so that dS needs to be calculated
$dH{AAAAAT} = 0.5; #
$dH{AAAACT} = 0.7; #
$dH{AAACAT} = 1; #
$dH{ACTTGT} = 0; #
$dH{AGAAAT} = -1.1; #
$dH{AGAGAT} = -1.1; #
$dH{AGATAT} = -1.5; #
$dH{AGCAAT} = -1.6; #
$dH{AGCGAT} = -1.1; #
$dH{AGCTTT} = 0.2; #
$dH{AGGAAT} = -1.1; #
$dH{AGGGAT} = -1.1; #
$dH{AGGGGT} = 0.5; #
$dH{AGTAAT} = -1.6; #
$dH{AGTGAT} = -1.1; #
$dH{AGTTCT} = 0.8; #
$dH{ATTCGT} = -0.2; #
$dH{ATTTGT} = 0; #
$dH{ATTTTT} = -0.5; #
$dH{CAAAAG} = 0.5; #
$dH{CAAACG} = 0.7;
$dH{CAACAG} = 1;
$dH{CAACCG} = 0;
$dH{CCTTGG} = 0; #
$dH{CGAAAG} = -1.1;
$dH{CGAGAG} = -1.1;
$dH{CGATAG} = -1.5;
$dH{CGCAAG} = -1.6;
$dH{CGCGAG} = -1.1;
$dH{CGCTTG} = 0.2;
$dH{CGGAAG} = -1.1;
$dH{CGGGAG} = -1;
$dH{CGGGGG} = 0.5; #
$dH{CGTAAG} = -1.6;
$dH{CGTGAG} = -1.1;
$dH{CGTTCG} = 0.8;
$dH{CTTCGG} = -0.2;
$dH{CTTTGG} = 0;
$dH{CTTTTG} = -0.5;
$dH{GAAAAC} = 0.5; #
$dH{GAAACC} = 0.7;
$dH{GAACAC} = 1;
$dH{GCTTGC} = 0; #
$dH{GGAAAC} = -1.1;
$dH{GGAGAC} = -1.1;
$dH{GGATAC} = -1.6;
$dH{GGCAAC} = -1.6;
$dH{GGCGAC} = -1.1;
$dH{GGCTTC} = 0.2;
$dH{GGGAAC} = -1.1;
$dH{GGGGAC} = -1.1;
$dH{GGGGGC} = 0.5;#
$dH{GGTAAC} = -1.6;
$dH{GGTGAC} = -1.1;
$dH{GGTTCC} = 0.8;
$dH{GTTCGC} = -0.2;
$dH{GTTTGC} = 0;
$dH{GTTTTC} = -0.5;
$dH{GAAAAT} = 0.5;#
$dH{GAAACT} = 1;
$dH{GAACAT} = 1;
$dH{GCTTGT} = 0; #
$dH{GGAAAT} = -1.1;
$dH{GGAGAT} = -1.1;
$dH{GGATAT} = -1.6;
$dH{GGCAAT} = -1.6;
$dH{GGCGAT} = -1.1;
$dH{GGCTTT} = -0.1;
$dH{GGGAAT} = -1.1;
$dH{GGGGAT} = -1.1;
$dH{GGGGGT} = 0.5;#
$dH{GGTAAT} = -1.6;
$dH{GGTGAT} = -1.1;
$dH{GTATAT} = -0.5;
$dH{GTTCGT} = -0.4;
$dH{GTTTGT} = -0.4;
$dH{GTTTTT} = -0.5;
$dH{TAAAAA} = 0.5;#
$dH{TAAACA} = 0.7;#
$dH{TAACAA} = 1;#
$dH{TCTTGA} = 0;#
$dH{TGAAAA} = -1.1;#
$dH{TGAGAA} = -1.1;#
$dH{TGATAA} = -1.6;#
$dH{TGCAAA} = -1.6;#
$dH{TGCGAA} = -1.1;#
$dH{TGCTTA} = 0.2;#
$dH{TGGAAA} = -1.1;#
$dH{TGGGAA} = -1.1;#
$dH{TGGGGA} = 0.5;#
$dH{TGTAAA} = -1.6;#
$dH{TGTGAA} = -1.1;#
$dH{TGTTCA} = 0.8;#
$dH{TTTCGA} = -0.2;#
$dH{TTTTGA} = 0;#
$dH{TTTTTA} = -0.5;#
$dH{TAAAAG} = 0.5;#
$dH{TAAACG} = 1;#
$dH{TAACAG} = 1;#
$dH{TCTTGG} = 0;#
$dH{TGAAAG} = -1;#
$dH{TGAGAG} = -1;#
$dH{TGATAG} = -1.5;#
$dH{TGCAAG} = -1.5;#
$dH{TGCGAG} = -1;#
$dH{TGCTTG} = -0.1;#
$dH{TGGAAG} = -1;#
$dH{TGGGAG} = -1;#
$dH{TGGGGG} = 0.5;#
$dH{TGTAAG} = -1.5;#
$dH{TGTGAG} = -1;#
$dH{TTTCGG} = -0.4;#
$dH{TTTTAG} = -1;#
$dH{TTTTGG} = -0.4;#
$dH{TTTTTG} = -0.5;#
#ΔS° cal/mol
$dS{TGAAAA} = 1.6121231662099;
$dS{TTTTTG} = 1.6121231662099;
$dS{ATTTTT} = 1.6121231662099;
$dS{TGCGAA} = 1.6121231662099;
$dS{TTTTTA} = 1.6121231662099;
$dS{AGCAAT} = 1.6121231662099;
$dS{TAACAG} = 1.6121231662099;
$dS{TGTAAA} = 1.6121231662099;
$dS{TGCAAA} = 1.6121231662099;
$dS{TCTTGA} = 4.19152023214574;
$dS{AAAAAT} = -0.644849266483959;
$dS{AAAACT} = 1.6121231662099;
$dS{TGTAAG} = 1.6121231662099;
$dS{TGGAAA} = 1.6121231662099;
$dS{AGATAT} = 1.6121231662099;
$dS{AGCTTT} = 1.6121231662099;
$dS{TGAAAG} = 1.6121231662099;
$dS{TAAAAG} = -1.6121231662099;
$dS{TAACAA} = 1.6121231662099;
$dS{TGCGAG} = 1.6121231662099;
$dS{AGTAAT} = 1.6121231662099;
$dS{CCTTGG} = 2.57939706593584;
$dS{AGTTCT} = 1.6121231662099;
$dS{TGGGAA} = 1.6121231662099;
$dS{TTTCGG} = 1.6121231662099;
$dS{AGCGAT} = 1.6121231662099;
$dS{TAAACG} = 1.6121231662099;
$dS{TGCAAG} = 1.6121231662099;
$dS{TAAAAA} = 0.32242463324198;
$dS{CAAAAG} = -1.28969853296792;
$dS{GCTTGT} = 1.6121231662099;
$dS{TCTTGG} = 3.2242463324198;
$dS{AGTGAT} = 1.6121231662099;
$dS{TGGAAG} = 1.6121231662099;
$dS{TGGGGG} = 0.644849266483959;
$dS{TGCTTA} = 1.6121231662099;
$dS{TGTGAG} = 1.6121231662099;
$dS{TGATAA} = 1.6121231662099;
$dS{AGGGGT} = 0.644849266483959;
$dS{TAAACA} = 1.6121231662099;
$dS{GAAAAT} = -3.2242463324198;
$dS{GGGGGC} = -0.967273899725939;
$dS{TTTTGG} = 1.6121231662099;
$dS{GGGGGT} = -0.967273899725939;
$dS{TGCTTG} = 1.6121231662099;
$dS{AGAAAT} = 1.6121231662099;
$dS{AGGGAT} = 1.6121231662099;
$dS{CGGGGG} = -0.967273899725939;
$dS{TTTTAG} = 1.6121231662099;
$dS{TGAGAA} = 1.6121231662099;
$dS{TTTCGA} = 1.6121231662099;
$dS{ACTTGT} = 4.19152023214574;
$dS{ATTCGT} = 1.6121231662099;
$dS{AGAGAT} = 1.6121231662099;
$dS{GAAAAC} = -3.2242463324198;
$dS{AGGAAT} = 1.6121231662099;
$dS{TGATAG} = 1.6121231662099;
$dS{GCTTGC} = 2.57939706593584;
$dS{ATTTGT} = 1.6121231662099;
$dS{TGTTCA} = 1.6121231662099;
$dS{TTTTGA} = 1.6121231662099;
$dS{TGGGGA} = 0.644849266483959;
$dS{AAACAT} = 1.6121231662099;
$dS{TGAGAG} = 1.6121231662099;
$dS{TGTGAA} = 1.6121231662099;
$dS{TGGGAG} = 1.6121231662099;

#Length dependence for hairpins
#ΔH° kcal/mol
$dH{3} = 0;
$dH{4} = 0;
$dH{5} = 0;
$dH{6} = 0;
$dH{7} = 0;
$dH{8} = 0;
$dH{9} = 0;
$dH{10} = 0;
$dH{11} = 0; ###approx
$dH{12} = 0;
$dH{13} = 0; ###approx
$dH{14} = 0;
$dH{15} = 0; ###approx
$dH{16} = 0;
$dH{17} = 0; ###approx
$dH{18} = 0;
$dH{19} = 0; ###approx
$dH{20} = 0;
#ΔS° cal/k·mol
$dS{3} = -11.3;
$dS{4} = -11.3;
$dS{5} = -10.6;
$dS{6} = -12.9;
$dS{7} = -13.5;
$dS{8} = -13.9;
$dS{9} = -14.5;
$dS{10} = -14.8;
$dS{11} = -15.5; ###approx
$dS{12} = -16.1;
$dS{13} = -16.4; ###approx
$dS{14} = -16.4;
$dS{15} = -16.8; ###approx
$dS{16} = -17.1;
$dS{17} = -17.4; ###approx
$dS{18} = -17.7;
$dS{19} = -18.1; ###approx
$dS{20} = -18.4;

#hairpin terminal mismatches
#ΔH° kcal/mol
$dH{t_AA_TA} = -3.2;
$dH{t_AA_TC} = -0.9;
$dH{t_AA_TG} = -2.3;
#$dH{t_AA_TT} = -5;
$dH{t_AC_TA} = -2.2;
$dH{t_AC_TC} = -0.5;
#$dH{t_AC_TG} = -6;
$dH{t_AC_TT} = -1.2;
$dH{t_AG_TA} = -2.7;
#$dH{t_AG_TC} = -6;
$dH{t_AG_TG} = -1.3;
$dH{t_AG_TT} = -2.9;
#$dH{t_AT_TA} = -5;
$dH{t_AT_TC} = -2.8;
$dH{t_AT_TG} = -3.5;
$dH{t_AT_TT} = -2.4;
$dH{t_CA_GA} = -2.8;
$dH{t_CA_GC} = -2;
$dH{t_CA_GG} = -3;
#$dH{t_CA_GT} = -6;
$dH{t_CC_GA} = -2.4;
$dH{t_CC_GC} = -1.4;
#$dH{t_CC_GG} = -7;
$dH{t_CC_GT} = -2.4;
$dH{t_CG_GA} = -5.1;
#$dH{t_CG_GC} = -7;
$dH{t_CG_GG} = -2.9;
$dH{t_CG_GT} = -2.9;
#$dH{t_CT_GA} = -6;
$dH{t_CT_GC} = -3.1;
$dH{t_CT_GG} = -5.9;
$dH{t_CT_GT} = -5.3;
$dH{t_GA_CA} = -6;
$dH{t_GA_CC} = -4;
$dH{t_GA_CG} = -3.4;
#$dH{t_GA_CT} = -6;
$dH{t_GC_CA} = -2.8;
$dH{t_GC_CC} = -2.5;
#$dH{t_GC_CG} = -7;
$dH{t_GC_CT} = -3.3;
$dH{t_GG_CA} = -4;
#$dH{t_GG_CC} = -7;
$dH{t_GG_CG} = -5.3;
$dH{t_GG_CT} = -3.7;
#$dH{t_GT_CA} = -6;
$dH{t_GT_CC} = -2.5;
$dH{t_GT_CG} = -4.5;
$dH{t_GT_CT} = -6.1;
$dH{t_GA_TA} = 0;
$dH{t_GA_TC} = 0;
$dH{t_GA_TG} = 0;
$dH{t_GA_TT} = -3.7;
$dH{t_GC_TA} = 0;
$dH{t_GC_TC} = 0;
$dH{t_GC_TG} = -4.5;
$dH{t_GC_TT} = 0;
$dH{t_GG_TA} = 0;
$dH{t_GG_TC} = -5.9;
$dH{t_GG_TG} = 0;
$dH{t_GG_TT} = -2;
$dH{t_GT_TA} = -3.5;
$dH{t_GT_TC} = 0;
$dH{t_GT_TG} = -2;
$dH{t_GT_TT} = 0;
$dH{t_TA_AA} = -2.7;
$dH{t_TA_AC} = -2.6;
$dH{t_TA_AG} = -2.4;
#$dH{t_TA_AT} = -5;
$dH{t_TC_AA} = -2.6;
$dH{t_TC_AC} = -0.5;
#$dH{t_TC_AG} = -6;
$dH{t_TC_AT} = -2.7;
$dH{t_TG_AA} = -1.9;
#$dH{t_TG_AC} = -6;
$dH{t_TG_AG} = -1.5;
$dH{t_TG_AT} = -2.3;
#$dH{t_TT_AA} = -5;
$dH{t_TT_AC} = -1.4;
$dH{t_TT_AG} = -3.7;
$dH{t_TT_AT} = -2.3;
$dH{t_TA_GA} = 0;
$dH{t_TA_GC} = 0;
$dH{t_TA_GG} = 0;
$dH{t_TA_GT} = -2.3;
$dH{t_TC_GA} = 0;
$dH{t_TC_GC} = 0;
$dH{t_TC_GG} = -3.7;
$dH{t_TC_GT} = 0;
$dH{t_TG_GA} = 0;
$dH{t_TG_GC} = -2.9;
$dH{t_TG_GG} = 0;
$dH{t_TG_GT} = -2;
$dH{t_TT_GA} = -2.9;
$dH{t_TT_GC} = 0;
$dH{t_TT_GG} = -2;
$dH{t_TT_GT} = 0;
#ΔS° cal/k·mol
$dS{t_GC_CC} = -6.12606803159761;
$dS{t_TA_GG} = 1.6121231662099;
$dS{t_TA_AA} = -6.77091729808157;
$dS{t_CA_GC} = -3.86909559890376;
$dS{t_CG_GG} = -6.44849266483959;
$dS{t_GT_TA} = -9.67273899725939;
$dS{t_TT_GT} = 0.644849266483959;
$dS{t_GG_CT} = -9.35031436401741;
$dS{t_TG_GA} = 1.6121231662099;
$dS{t_CC_GT} = -5.48121876511366;
$dS{t_GT_CC} = -6.12606803159761;
$dS{t_AG_TT} = -7.73819119780751;
$dS{t_AC_TT} = -2.90182169917782;
$dS{t_AG_TA} = -6.77091729808157;
$dS{t_GT_TC} = 0.644849266483959;
$dS{t_GC_TC} = 0.644849266483959;
$dS{t_GA_TG} = 1.6121231662099;
$dS{t_AA_TG} = -5.80364339835563;
$dS{t_CC_GC} = -2.90182169917782;
$dS{t_TG_AT} = -5.80364339835563;
$dS{t_GC_TG} = -11.6072867967113;
$dS{t_GT_TT} = 0.644849266483959;
$dS{t_TC_AT} = -7.09334193132355;
$dS{t_GG_TA} = 1.6121231662099;
$dS{t_CA_GG} = -6.77091729808157;
$dS{t_GA_TT} = -9.99516363050137;
$dS{t_CG_GT} = -6.12606803159761;
$dS{t_GC_CA} = -5.80364339835563;
$dS{t_TG_GC} = -6.12606803159761;
$dS{t_TA_GA} = 1.6121231662099;
$dS{t_GA_TC} = 0.644849266483959;
$dS{t_TT_GG} = -4.8363694986297;
$dS{t_TT_AT} = -6.44849266483959;
$dS{t_GG_CA} = -9.67273899725939;
$dS{t_AT_TC} = -8.06061583104949;
$dS{t_TC_GT} = 0.644849266483959;
$dS{t_TA_GC} = 0.644849266483959;
$dS{t_TT_AG} = -9.99516363050137;
$dS{t_GT_CG} = -11.6072867967113;
$dS{t_AG_TG} = -2.90182169917782;
$dS{t_CC_GA} = -5.15879413187167;
$dS{t_TC_AC} = -0.967273899725939;
$dS{t_CT_GC} = -8.06061583104949;
$dS{t_AA_TA} = -8.06061583104949;
$dS{t_GA_CA} = -16.121231662099;
$dS{t_TT_AC} = -3.54667096566178;
$dS{t_TC_AA} = -6.77091729808157;
$dS{t_GG_TT} = -4.8363694986297;
$dS{t_GG_TG} = 1.6121231662099;
$dS{t_TG_AA} = -4.19152023214574;
$dS{t_AT_TG} = -9.67273899725939;
$dS{t_GA_CC} = -10.6400128969853;
$dS{t_GA_CG} = -8.38304046429147;
$dS{t_AA_TC} = -1.93454779945188;
$dS{t_CT_GG} = -16.121231662099;
$dS{t_GC_TA} = 0.644849266483959;
$dS{t_TG_GG} = 1.6121231662099;
$dS{t_TC_GC} = 0.644849266483959;
$dS{t_GC_CT} = -8.38304046429147;
$dS{t_TA_AG} = -6.12606803159761;
$dS{t_GT_TG} = -4.8363694986297;
$dS{t_CA_GA} = -5.80364339835563;
$dS{t_CT_GT} = -14.1866838626471;
$dS{t_TG_GT} = -4.8363694986297;
$dS{t_TT_GA} = -7.73819119780751;
$dS{t_GA_TA} = 1.6121231662099;
$dS{t_AT_TT} = -6.44849266483959;
$dS{t_GG_CG} = -13.8642592294051;
$dS{t_AC_TA} = -5.15879413187168;
$dS{t_CG_GA} = -13.2194099629212;
$dS{t_GC_TT} = 0.644849266483959;
$dS{t_TC_GG} = -9.35031436401741;
$dS{t_TT_GC} = 0.644849266483959;
$dS{t_TA_GT} = -5.80364339835563;
$dS{t_TG_AG} = -3.54667096566178;
$dS{t_GT_CT} = -16.7660809285829;
$dS{t_AC_TC} = -0.967273899725939;
$dS{t_GG_TC} = -16.121231662099;
$dS{t_TA_AC} = -7.09334193132355;
$dS{t_TC_GA} = 0.644849266483959;

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

#initial bases
$sum_dS += $dS{initial};
$sum_dH += $dH{initial};

#no simmetry correction because it is not a simmetrical formation!!

if ($loopLen == 3) {
    
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
    
    #H and S values from DNA duplex
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
        $offset++;
    }
    
    #loop len
    $sum_dH += $dH{$loopLen};
    $sum_dS += $dS{$loopLen};
    
    #triloop bonus
    $sum_dH += $dH{$loop};
    $sum_dS += $dS{$loop};

} elsif ($loopLen == 4) {
    #terminal bases corrected if A_T ot T_A
    $start_1 =  substr($seq1, 0, 1);
    $start_2 =  substr($seq2, 0, 1);
    $pair = $start_1 . "_" . $start_2;
    if (($pair eq 'A_T') or ($pair eq 'T_A')) {
        $dS_corrected = $dS{$terminal}+($dH{$terminal}*$D);
        $sum_dS += $dS_corrected;
        $sum_dH += $dH{terminal};
    }
    #H and S values from DNA duplex
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
        $offset++;
    }
    
    #loop len
    $sum_dH += $dH{$loopLen};
    $sum_dS += $dS{$loopLen};
    
    #tetralopp bonus
    $sum_dH += $dH{$loop};
    $sum_dS += $dS{$loop};
    
    #terminal mismatch
    $seq1_last = substr($seq1, -1);
    $seq2_last = substr($seq2, -1);
    $loop_first = substr($loop, 1, 1);
    $loop_last = substr($loop, -2, 1);
    $mis = "t_" . $seq1_last . $loop_first . "_" . $seq2_last . $loop_last;
    $sum_dH += $dH{$mis};
    $sum_dS += $dS{$mis};
} elsif ($loopLen > 4) {
    #terminal bases corrected if A_T ot T_A
    $start_1 =  substr($seq1, 0, 1);
    $start_2 =  substr($seq2, 0, 1);
    $pair = $start_1 . "_" . $start_2;
    if (($pair eq 'A_T') or ($pair eq 'T_A')) {
        $dS_corrected = $dS{$terminal}+($dH{$terminal}*$D);
        $sum_dS += $dS_corrected;
        $sum_dH += $dH{terminal};
    }
    
    #H and S values from DNA duplex
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
        $offset++;
    }
    
    #loop len
    $sum_dH += $dH{$loopLen};
    $sum_dS += $dS{$loopLen};
    
    #terminal mismatch
    $seq1_last = substr($seq1, -1);
    $seq2_last = substr($seq2, -1);
    $loop_first = substr($loop, 1, 1);
    $loop_last = substr($loop, -2, 1);
    $mis = "t_" . $seq1_last . $loop_first . "_" . $seq2_last . $loop_last;
    $sum_dH += $dH{$mis};
    $sum_dS += $dS{$mis};
}

$dS_corrected = $sum_dS + ($sum_dH*$D);

$dG = $sum_dH - $temperature*($dS_corrected/1000);

print "$dG\n";

