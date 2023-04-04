#!/usr/bin/perl

use strict;
use DBD::mysql;
use POSIX;

#connected to backgroundDesignCheck.pl : oligo design pipeline

#set path to folders:
my $path_cgi = 'path_to_cgi_folder';
my $path_html = 'path_to_html_folder';

#set mySQL parameters
my $dsn = "mysql_database";
my $user_name = "mysql_user";
my $password = "mysql_password";

#Usage: ./Oligo_design_pp.pl -folder $path_html/analysesPhyloprimer/AkuVricWqWoYjmWDMStmTa -file AkuVricW.info

##NOTES:
#Tm would be the same for oligo and reverse if I was not considering dangling ends
#Indicated position of oligo and reverse must always be the position at 5' end so F = position and R = position+length

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

my %args = @ARGV;

my $fileInfo = $args{-file};
my $folder = $args{-folder};
my $nameFile = substr($fileInfo, 0, 8);

my $indexType = substr($folder, -2,1);

open(IN, "<$folder/$fileInfo") or die; #file with all the necessary parameters

#save input into array
my %all;
while (defined(my $input = <IN>)) {
    chomp($input);
    my @data = split(/\t/, $input);
    $all{$data[0]} = $data[1];
}

#make temporary directory
`mkdir $folder/tmp 2> /dev/null`;

#make directory that will be compressed at the end for the user
`mkdir $folder/results 2> /dev/null`;
`mkdir $folder/inputs 2> /dev/null`;


`mv $folder/$fileInfo $folder/tmp/$fileInfo`; #move the input parameter file read by PhyloPrimer to the tmp directory

my $tail = 2;

#wildcards - define minredundant and maxredundant
my $minredundant;
my $maxredundant;
my $minwild;
my $maxwild;
my @allowwild = ();

#define max redundancy
if ($all{'NUMBERWILD_MAX'} == 0) {
    $minredundant = 0;
    $maxredundant = 0;
} else {
    my %r;
    $r{'R'} = 2;
    $r{'Y'} = 2;
    $r{'S'} = 2;
    $r{'W'} = 2;
    $r{'K'} = 2;
    $r{'M'} = 2;
    $r{'B'} = 3;
    $r{'D'} = 3;
    $r{'H'} = 3;
    $r{'V'} = 3;
    $r{'N'} = 4;
    my @which = split(/,/, $all{'WHICHWILD'});
    $minwild = $which[0];
    foreach my $w (@which) {
        if ($w == 2) {
            push @allowwild, ('R', 'Y', 'S', 'W', 'K', 'M');
        }
        if ($w == 3) {
            push @allowwild, ('B', 'D', 'H', 'V');
        }
        if ($w == 4) {
            push @allowwild, 'N';
        }
        $maxwild = $w;
    }
    $minredundant = ($all{'NUMBERWILD_MIN'} * $minwild);
    $maxredundant = ($all{'NUMBERWILD_MAX'} * $maxwild);
}

#long runs and repeats
my $longrun = $all{'MAXRUN'}; ###include this
my $repeat = ($all{'MAXREPEAT'}/2); ###include this

my $a = 'A' x $longrun;
my $t = 'T' x $longrun;
my $c = 'C' x $longrun;
my $g = 'G' x $longrun;

my $ta = 'TA' x $repeat;
my $tg = 'TG' x $repeat;
my $tc = 'TC' x $repeat;
my $at = 'AT' x $repeat;
my $ag = 'AG' x $repeat;
my $ac = 'AC' x $repeat;
my $gt = 'GT' x $repeat;
my $ga = 'GA' x $repeat;
my $gc = 'GC' x $repeat;
my $ct = 'CT' x $repeat;
my $ca = 'CA' x $repeat;
my $cg = 'CG' x $repeat;


#GC clamp
my $gcclamp_min = ($all{'GCCLAMP_MIN'}/5)*100;
my $gcclamp_max = ($all{'GCCLAMP_MAX'}/5)*100;

#highlight - determine areas for primers design
my %area;
my $lenConsensus = length($all{'CONSENSUS'});

if ($all{'CONSERVED'} eq "yes") {
    $area{'F'} = $all{'HIGHF'};
} elsif ($all{'CONSERVED'} eq "no") {
    $area{'ALL'} = $all{'CONSENSUS'}; #if none of the areas was specified by the user
}

my %selectedTaxa;

my $nameFileSpecies = $nameFile . ".taxonomy";
if (-e $nameFileSpecies) {
    open(IN, "<$folder/$nameFileSpecies") or die;
    while(defined(my $input = <IN>)) {
        chomp($input);
        my @info = split(/\t/, $input);
        $selectedTaxa{$info[0]}{$info[1]} = '';
    }
    close(IN);
}

my @mismatch = (0,1,2,3,4,5,6); #number of mismatches - max allowed is 3 per oligo (6 per oligo pair

##parameters for dG and Tm calculation
my $monovalent =  $all{'MON_DG'}/1000;
my $mg_tot = $all{'MG_DG'}/1000;
my $C = $all{'OLIGO_DG'}/1000000;
my $dNTP_tot = $all{'DNTP_DG'}/1000;
my $temperature_celsius = $all{'T_DG'};

##diffences in position between positive and negative consensus
#check there are no gaps at the beginning or end of the consensus - these gaps are not suitable regions for the oligo design because could be the result of fragmented input sequences
my $diffPos = $all{'DIFFERENT_POS'};
my @colGap = split(/\|/, $diffPos);
my @gaps = split(/,/, $colGap[2]); #gaps - blue letters

my $limMin = 10;
my $limMax = $lenConsensus - 10;
my %gaps;
my %avoidGap;
foreach my $g (@gaps) {
    $gaps{$g} = '';
    if (($g <= $limMin) or ($g >= $limMax)) {
        $avoidGap{$g} = "";
    }
}

my $limMinNew = $limMin;
if (defined($avoidGap{$limMin})) {
    my $done = 0;
    while($done == 0) {
        $limMinNew++;
        if (!(defined($gaps{$limMinNew}))) {
            $done = 1;
        } else {
            $avoidGap{$limMinNew} = '';
        }
    }
}

my $limMaxNew = $limMax;
if (defined($avoidGap{$limMax})) {
    my $done = 0;
    while($done == 0) {
        $limMaxNew--;
        if (!(defined($gaps{$limMaxNew}))) {
            $done = 1;
        } else {
            $avoidGap{$limMaxNew} = '';
        }
    }
}

#Create hash which consensus differences will be taken in consideration in the scoring system
my %difference;
$diffPos = $all{'DIFFERENT_POS'};
$diffPos =~ s/\|//g;
my @difference = split(/,/, $diffPos);
foreach my $d (@difference) {
    if (!(defined($avoidGap{$d}))) {
        $difference{$d} = '';
    }
}

my %oligoAll;
my $pos;

#check for primers that are present more than once in different positions on the consensus
my %duplicate;
my %exclude;

my $countEX = 0;

OUTER:
foreach my $plength ($all{'LEN_MIN'} .. $all{'LEN_MAX'}) { #so then I check for duplicates across all consensus
    my $seq = $all{'CONSENSUS'};
    my $l = length($seq);
    my $i = 0;
    #Go through the consensus sequence with a sliding window of $plength, one base at a time
    while ($i < ($l - $plength)) {
        my $seq = substr($seq, $i, $plength);
        my $rev = reverse($seq); #reverse sequence here is intended as the sequence where a reverse primer would anneal to the sense strand.

        #ATGCCGTATGGTGTA
        #TACGGCATACCACAT

        #seq to discard is ATG and rev to discard is GTA (which is the reverse if seq). GTA is a portion to which the reverse primer TAC would anneal.

        my ($pRef) = degenerateAlt($seq);
        my @primerDeg = @{$pRef}; #all degenearte oligos
        foreach my $pD (@primerDeg) {
          if (defined($duplicate{$pD})) { #if this primer already existed
              if (!(defined($exclude{$seq}))) {
                $countEX++;
              }
              if (!(defined($exclude{$pD}))) {
                $countEX++;
              }
              #print "$countEX\n";
              if ($countEX > 100000) { #exit this loop if more than 100000 repetitions, otherwise %exclude takes too much RAM. PhyloPrimer gets many repetitions when many degenerate bases in consensus
                last OUTER;
              }
              $exclude{$seq} = '';
              $exclude{$pD} = '';
          }
          my $pD_rev = reverse($pD); #reverse seq
          if (defined($duplicate{$pD_rev})) {
              if (!(defined($exclude{$rev}))) {
                $countEX++;
              }
              if (!(defined($exclude{$pD_rev}))) {
                $countEX++;
              }
              #print "$countEX\n";
              if ($countEX > 100000) { #exit this loop if more than 100000 repetitions, otherwise %exclude takes too much RAM. PhyloPrimer gets many repetitions when many degenerate bases in consensus
                last OUTER;
              }
              $exclude{$rev} = '';
              $exclude{$pD_rev} = '';
          }
          $duplicate{$pD} = '';
          $duplicate{$pD_rev} = '';
        }
        $duplicate{$seq} = '';
        $duplicate{$rev} = '';
        $i++;
    }
    undef %duplicate; #empty duplicate every different length
}

#scan the consensus areas for each length and a shift of 1
my $presence0 = 0; #how many oligos can be potentially designed from the consensus
my $presence1 = 0;
my $diffPres = 0;
my @diffMessage;

foreach my $type (sort keys %area) {
    foreach my $plength ($all{'LEN_MIN'} .. $all{'LEN_MAX'}) {
        my $seq = $all{'CONSENSUS'};
        my $l;

        if ($type eq 'F') {
            $pos = $all{'HIGHF_START'} -1;
            $l = $all{'HIGHF_START'} + length($area{$type}) + 1;
        } else {
            $pos = $limMinNew; #start from one and finish one before to be sure!
            $l = $limMaxNew; #start from one and finish one before to be sure!
        }

        my $i = $pos;

        #Go through the consensus sequence with a sliding window of $plength, one base at a time
        while ($i < ($l - $plength)) {
            $presence0++;
            my $primer = substr($seq, $i, $plength);

            if (!(defined($exclude{$primer}))) {

                #Get the redundancy level for the primer so we can check it doesn't exceed the $maxredundant threshold below
                my ($redundancy,$numberwild) = redundancy($primer); #0 is no wildcards

                #Get the final starting coordinate for this primer, based on the input sequence position (5' end)
                my $posReal = $pos + 1;
                my $final = "${posReal}-$primer";

                #CHECK 1 - no ambiguous bases in the first or last $tail bases of the primer, and redundancy check
                if(($redundancy <= $maxredundant) && ($primer =~ /[CATG]{$tail}$/ || $primer =~ /^[CATG]{$tail}/) && ($numberwild <= $all{'NUMBERWILD_MAX'})) {

                    #CHECK 2 - long runs
                    if (($final !~ /$a/) && ($final !~ /$t/) && ($final !~ /$c/) && ($final !~ /$g/)) {

                        #CHECK 3 - repeats
                        if (($final !~ /$ta/) && ($final !~ /$tg/) && ($final !~ /$tc/)) {
                            if (($final !~ /$at/) && ($final !~ /$ag/) && ($final !~ /$ac/)) {
                                if (($final !~ /$gt/) && ($final !~ /$ga/) && ($final !~ /$gc/)) {
                                    if (($final !~ /$ct/) && ($final !~ /$ca/) && ($final !~ /$cg/)) {
                                        $oligoAll{$final} = ""; #primer = ""
                                        $presence1++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            $i++;
            $pos++;
        }
    }
}

undef %exclude;
if ($presence0 == 0) {
    $diffPres = 0;
} else {
    $diffPres = $presence1/$presence0;
}

my $message;
if (($diffPres < 0.3) and ($countEX > 100000)) {
    push @diffMessage, '<li>Try to increase the redundancy level and the homopolymer length.</li>';
    push @diffMessage, '<li>The consensus sequence reported a high level of degenerate bases; try to change the sequence selection.</li>';
    $message = "Content-Type: text/html; charset=ISO-8859-1\n\n<html><body>Hi,<br>PhyloPrimer did find any suitable oligos with the selected parameters. Try to widen the search criteria for redundancy and homopolymers. Also, the consensus sequence reported a high level of degenerate bases: try to change the sequence selection.<br><br><br>All the best,<br>Gilda</body>";
} elsif ($diffPres < 0.3) {
    push @diffMessage, '<li>Try to increase the redundancy level and the homopolymer length</li>';
    $message = "Content-Type: text/html; charset=ISO-8859-1\n\n<html><body>Hi,<br>PhyloPrimer did find any suitable oligos with the selected parameters. Try to widen the search criteria for redundancy and homopolymers.<br><br><br>All the best,<br>Gilda</body>";
}

#if no suitable oligos were found during the first screening
if (($presence0 == 0) or ($presence1 == 0)) {
    `mv ${folder}/${nameFile}* ${folder}/inputs/`;
    #send email
    my $to = $all{'EMAIL'};
    my $from = 'gv16363@bristol.ac.uk';
    my $subject = 'PhyloPrimer results - ' . $all{'PROJECT'};
    
    open(MAIL, "|/usr/sbin/sendmail -t");

    # Email Header
    print MAIL "To: $to\n";
    print MAIL "From: $from\n";
    print MAIL "Subject: $subject\n";
    # Email Body
    print MAIL $message;

    close(MAIL);
    exit;
}

my %oligo;
my $presence2=0;
if (($all{'STRAND'} eq 'anti') or ($all{'STRAND'} eq 'both')) {
    foreach my $primer (keys %oligoAll) {

        my ($pos, $seq) = split(/-/, $primer);
        my $gc = gplusc($seq);

        #CHECK 4 - GC content
        if ($gc >= $all{'GC_CONTENT_MIN'} && $gc <= $all{'GC_CONTENT_MAX'}) {
            my $plength = length($seq);
            my $tm;
            if ($all{'TYPE'} eq 'primer') {
                my $dang = substr($all{'CONSENSUS'}, ($pos+$plength-1), 1);
                $tm = `$path_cgi/tm_calculation_pp.pl -primer $seq -type primer -sense F -mg $mg_tot -dang $dang -mon $monovalent -oligo $C -dntp $dNTP_tot`; #Tm for oligo
            } else {
                my $dangS = substr($all{'CONSENSUS'}, ($pos-2), 1);
                my $dangE = substr($all{'CONSENSUS'}, ($pos+$plength-1), 1);
                $tm = `$path_cgi/tm_calculation_pp.pl -primer $seq -type probe -sense F -mg $mg_tot -dang ${dangS},${dangE} -mon $monovalent -oligo $C -dntp $dNTP_tot`; #Tm for oligo
            }
            chomp($tm);
            $tm = sprintf("%.2f", $tm);

            #CHECK 5 - melting temperature
            if (($tm >= $all{'TM_MIN'}) && ($tm <= $all{'TM_MAX'})) {
                my $tail = substr($seq, -5);
                my $tail_gc = gplusc($tail);

                #CHECK 6 - GC clamp
                if (($tail_gc >= $gcclamp_min) && ($tail_gc <= $gcclamp_max)) {
                    my $redundancy = redundancy($primer);
                    my $bad = 0;
                    if ($redundancy > 0) {
                        my $tail = substr($seq, -($all{'END3WILD'}));
                        if ($tail =~ /[RYSWKMBDHVN]/) {
                            $bad = 1;
                        }
                        my $head = substr($seq, 0 ,$all{'END5WILD'});
                        if ($head =~ /[RYSWKMBDHVN]/) {
                            $bad = 1;
                        }
                    }
                    if ($bad == 0) {
                        $presence2++;
                        #oligo information
                        $oligo{$seq}{'POS'} = $pos;
                        $oligo{$seq}{'POS_NAME'} = $pos;
                        $oligo{$seq}{'LEN'} = $plength;
                        $oligo{$seq}{'GC'} = $gc;
                        $oligo{$seq}{'RED'} = $redundancy;
                        $oligo{$seq}{'TM'} = $tm;
                        $oligo{$seq}{'STRAND'} = 'anti';
                    }
                }
            }
        }
    }
}

if (($all{'STRAND'} eq 'sense') or ($all{'STRAND'} eq 'both')) {
    foreach my $primer (keys %oligoAll) {

        my ($pos, $seq) = split(/-/, $primer);
        my $gc = gplusc($seq);

        #CHECK 4 - GC content
        if ($gc >= $all{'GC_CONTENT_MIN'} && $gc <= $all{'GC_CONTENT_MAX'}) {
            my $rev = $seq; #reverse seq
            $rev = reverse($rev);
            $rev =~ tr/CATGYRKMBVDH/GTACRYMKVBHD/;

            my $plength = length($seq);
            my $tm;
            if ($all{'TYPE'} eq 'primer') {
                my $dang = substr($all{'CONSENSUS'}, ($pos-2), 1);
                $tm = `$path_cgi/tm_calculation_pp.pl -primer $rev -type primer -sense R -mg $mg_tot -dang $dang -mon $monovalent -oligo $C -dntp $dNTP_tot`; #Tm for oligo
            } else {
                my $dangS = substr($all{'CONSENSUS'}, ($pos-2), 1);
                my $dangE = substr($all{'CONSENSUS'}, ($pos+$plength-1), 1);
                $tm = `$path_cgi/tm_calculation_pp.pl -primer $rev -type probe -sense R -mg $mg_tot -dang ${dangS},${dangE} -mon $monovalent -oligo $C -dntp $dNTP_tot`; #Tm for oligo
            }
            chomp($tm);
            $tm = sprintf("%.2f", $tm);

            #CHECK 5 - melting temperature
            if (($tm >= $all{'TM_MIN'}) && ($tm <= $all{'TM_MAX'})) {
                my $tail = substr($rev, -5);
                my $tail_gc = gplusc($tail);

                #CHECK 6 - GC clamp
                if (($tail_gc >= $gcclamp_min) && ($tail_gc <= $gcclamp_max)) {
                    my $redundancy = redundancy($primer);
                    my $bad = 0;
                    if ($redundancy > 0) {
                        my $tail = substr($rev, -($all{'END3WILD'}));
                        if ($tail =~ /[RYSWKMBDHVN]/) {
                            $bad = 1;
                        }
                        my $head = substr($rev, 0 ,$all{'END5WILD'});
                        if ($head =~ /[RYSWKMBDHVN]/) {
                            $bad = 1;
                        }
                    }
                    if ($bad == 0) {
                        $presence2++;
                        #oligo information
                        $oligo{$rev}{'POS'} = $pos;
                        $oligo{$rev}{'POS_NAME'} = $pos + $plength -1;
                        $oligo{$rev}{'LEN'} = $plength;
                        $oligo{$rev}{'GC'} = $gc;
                        $oligo{$rev}{'RED'} = $redundancy;
                        $oligo{$rev}{'TM'} = $tm;
                        $oligo{$rev}{'STRAND'} = 'sense';

                    }
                }
            }
        }
    }
}

$diffPres = $presence2/$presence1;

if ($diffPres < 0.5) {
    push @diffMessage, '<li>Try to increase dG and melting temperature</li>';
}

#if no suitable oligos were found during the first screening
if ($presence2 == 0) {
    `mv ${folder}/${nameFile}* ${folder}/inputs/`;
    #send email
    my $to = $all{'EMAIL'};
    my $from = 'gv16363@bristol.ac.uk';
    my $subject = 'PhyloPrimer results - ' . $all{'PROJECT'};

    my $message = "Content-Type: text/html; charset=ISO-8859-1\n\n<html><body>Hi,<br>Phyloprimer did find any suitable forward primers. Try to widen the search criteria.<br><br><br>All the best,<br>Gilda</body>";

    open(MAIL, "|/usr/sbin/sendmail -t");

    # Email Header
    print MAIL "To: $to\n";
    print MAIL "From: $from\n";
    print MAIL "Subject: $subject\n";
    # Email Body
    print MAIL $message;

    close(MAIL);
    exit;
}
my %oligo_pass;

#CHECK 7 - hairpins
#oligo
my $size_oligo = ceil((scalar keys %oligo) / 5); #number of forward oligos
open(my $tmp, ">$folder/tmp/Hairpin_oligo_1.tmp") or die;
my $fileName = "Hairpin_oligo_1.tmp";
my $hair;
foreach my $primer (keys %oligo) {
    $hair++;
    if (($hair % $size_oligo) == 0) {
        print $tmp "$primer\n";
        close($tmp);
        system("perl $path_cgi/Hairpin_checking_web_pp.pl -folder $folder -primer $fileName -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius &");
        $fileName = "Hairpin_oligo_" . $hair . ".tmp";
        open($tmp, ">$folder/tmp/Hairpin_oligo_${hair}.tmp") or die;
    } else {
        print $tmp "$primer\n";
    }
}
system("perl $path_cgi/Hairpin_checking_web_pp.pl -folder $folder -primer $fileName -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius &");

my $done = 0;
while ($done == 0) {
    sleep 30;
    my $list = `ls $folder/tmp/Hairpin_*tmp 2> /dev/null`;
    if ($list eq '') {
        $done = 1;
    }
}
`cat $folder/tmp/SecondaryStructure_Hairpin_*.tmp > $folder/tmp/SecondaryStructure_hairpin.txt`;
`cat $folder/tmp/SecondaryStructure_Hairpin_*.tmp.delta > $folder/tmp/Hairpin_delta_oligo.tmp`;

open(IN, "<$folder/tmp/Hairpin_delta_oligo.tmp") or die;

open(my $tmp, ">$folder/tmp/Self_oligo_1.tmp") or die;
my $fileName = "Self_oligo_1.tmp";

my $self;
while(defined(my $input = <IN>)) { #only for forward
    chomp($input);
    my ($oligo, $dG_hair) = split(/\t/, $input); #what is ) or null?!?
    if (($dG_hair > $all{'OLI_HAIRPIN'}) or ($dG_hair eq '')) {
        $self++;
        if (($self % $size_oligo) == 0) {
            print $tmp "$oligo\n";
            close($tmp);
            system("perl $path_cgi/Self_dimer_checking_web_pp.pl -folder $folder -primer $fileName -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius &");
            $fileName = "Self_oligo_" . $self . ".tmp";
            open($tmp, ">$folder/tmp/Self_oligo_${self}.tmp") or die;
        } else {
            print $tmp "$oligo\n";
        }
    }
}
close(IN);
system("perl $path_cgi/Self_dimer_checking_web_pp.pl -folder $folder -primer $fileName -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius &");

$done = 0;
while ($done == 0) {
    sleep 30;
    my $list = `ls $folder/tmp/Self_*tmp 2> /dev/null`;
    if ($list eq '') {
        $done = 1;
    }
}
`cat $folder/tmp/SecondaryStructure_Self_*.tmp > $folder/tmp/SecondaryStructure_selfDimer.txt`;

#CHECK 8 - self dimers
my %oligo_pass;
my %allPrimer;
my $presence3 = 0;

#forward oligos
`cat $folder/tmp/SecondaryStructure_Self_oligo*.tmp.delta > $folder/tmp/Self_delta_oligo.tmp`;

open(IN, "<$folder/tmp/Self_delta_oligo.tmp") or die;

while(defined(my $input = <IN>)) { #only for forward
    chomp($input);
    my ($oligo, $dG_self) = split(/\t/, $input); #what is ) or null?!?
    if (($dG_self > $all{'OLI_SELF'}) or ($dG_self eq ''))  {
        $oligo{$oligo}{'SELF'} = $dG_self;
        $presence3++;
        #        $oligo_pass{$oligo{$oligo}{'POS'}}{$oligo} = ""; #start forward
        $oligo_pass{$oligo} = ""; #start forward
        $allPrimer{$oligo} = '';
    }
}
close(IN);

$diffPres = $presence3/$presence2;

if ($diffPres < 0.5) {
    push @diffMessage, '<li>amplicon size</li>';
}

if ($presence3 == 0) {
    `mv ${folder}/${nameFile}* ${folder}/inputs/`;
    #send email
    my $to = $all{'EMAIL'};
    my $from = 'gv16363@bristol.ac.uk';
    my $subject = 'PhyloPrimer results - ' . $all{'PROJECT'};
    my $message = "Content-Type: text/html; charset=ISO-8859-1\n\n<html><body>Hi,<br>PhyloPrimer did not find any suitable primers with the selected parameters. It looks like the parameters for the following fields were too stringent:<br><ul>@diffMessage</ul><br><br>All the best,<br>Gilda</body>";

    open(MAIL, "|/usr/sbin/sendmail -t");

    # Email Header
    print MAIL "To: $to\n";
    print MAIL "From: $from\n";
    print MAIL "Subject: $subject\n";
    # Email Body
    print MAIL $message;

    close(MAIL);

    exit;
}

#check cross dimers between same primer with degenerate bases!

my $size = ceil($presence3/5); #divided into max 5 files

my $cross = 0;
my $file_num = 0;
my $combo; #array with all oligo variations
open(my $tmp, ">$folder/tmp/CrossSelf_1.tmp") or die;

foreach my $primer (keys %oligo_pass) {
    if ($primer =~ /[RYSWKMBDHVN]/) { #if deg bases: find all possible oligos
        $cross++;
        my ($pRef) = degenerateAlt($primer);
        my @primerDeg = @{$pRef};
        my $combo = $primer; #original
        my $num = scalar(@primerDeg) - 1; #number of oligos in the pair
        foreach my $c1 (0..($num-1)) { #first oligo
            foreach my $c2 (($c1+1)..$num) { #second
                $combo .= "\t" . $primerDeg[$c1] . "," . $primerDeg[$c2];
            }
        }
        if ($cross == 1) {
            print $tmp "$combo\n";
            $file_num = $cross;
        } elsif (($cross % $size) == 0) {
            print $tmp "$combo\n";
            close($tmp);
            system("perl $path_cgi/Cross_dimer_checking_self_web_pp.pl -folder $folder -primer CrossSelf_${file_num}.tmp -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius &");
            open($tmp, ">$folder/tmp/CrossSelf_${cross}.tmp") or die;
            $file_num = $cross;
        } else {
            print $tmp "$combo\n";
        }
    }
}
close($tmp);

if ($cross == 0) { #check if it empty
    system("perl $path_cgi/Cross_dimer_checking_self_web_pp.pl -folder $folder -primer CrossSelf_1.tmp -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius &");
    `rm ${folder}/tmp/CrossSelf_1.tmp`;
} else {
    system("perl $path_cgi/Cross_dimer_checking_self_web_pp.pl -folder $folder -primer CrossSelf_${file_num}.tmp -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius &");
}

my $done = 0;
while ($done == 0) {
    sleep 30;
    my $list = `ls $folder/tmp/CrossSelf_*tmp 2> /dev/null`;
    if ($list eq '') {
        $done = 1;
    }
}

`cat $folder/tmp/SecondaryStructure_crossSelf_*.tmp > $folder/crossSelf.txt`;

#save oligos in hash

open(IN, "<${folder}/crossSelf.txt") or die; #file with all the necessary parameters

my %allPrimerRemove;

my %crossSelf;
my $disc;
my $oligo;
my $dG_cross;
my $one;
my $a;
my $two;
my $dG;
my $count = 0;

while (defined(my $input = <IN>)) {
    chomp($input);
    if ($input =~ /^>/) {
        ($disc, $oligo) = split(/\t/, $input);
    } elsif ($input =~ /^@/) {
        ($disc, $dG_cross) = split(/\t/, $input);
        if ($dG_cross < $all{'OLI_SELF'}) {
            $allPrimerRemove{$oligo} = ''; #does not take in account these primers
        }
    } elsif ($input =~ /^1:/) {
        ($disc, $one) = split(/\t/, $input);
    } elsif ($input =~ /^A:/) {
        $count++;
        ($disc, $a) = split(/\t/, $input);
    } elsif ($input =~ /^2:/) {
        ($disc, $two) = split(/\t/, $input);
    } elsif ($input =~ /^dG:/) {
        ($disc, $dG) = split(/\t/, $input);
        push @{$crossSelf{$oligo}{$count}}, ($one,$a,$two,$dG);
    }
}
close(IN);

#re-print self file
open(IN, "<$folder/tmp/SecondaryStructure_selfDimer.txt") or die; #I'll have to add /tmp/
open(my $file, ">$folder/selfDimer.txt") or die;
my $oligo;
my $in = 0;
my $lenSelf1;
my $lenSelfA;
my %allCount;
my %alldG;

while(defined(my $input = <IN>)) {
    chomp($input);
    if ($input =~ /^>/) {
        ($disc, $oligo) = split(/\t/, $input);
        if ((defined($oligo_pass{$oligo})) && (!(defined($allPrimerRemove{$oligo})))) {
            $in = 1;
            print $file "\n>\t$oligo\n";
        } else {
            $in = 0;
        }
        $lenSelf1 = 0;
        $lenSelfA = 0;
        undef %allCount;
    } elsif ($input =~ /^1:/) {
        ($disc, $one) = split(/\t/, $input);
        $lenSelf1 = length($one);
    } elsif ($input =~ /^A:/) {
        ($disc, $a) = split(/\t/, $input);
        $lenSelfA = length($a);
    } elsif ($input =~ /^2:/) {
        ($disc, $two) = split(/\t/, $input);
    } elsif ($input =~ /^dG:/) {
        ($disc, $dG) = split(/\t/, $input);
        foreach my $c (sort {$a <=> $b} keys %{$crossSelf{$oligo}}) {
            my $len1 = length($crossSelf{$oligo}{$c}[0]);
            my $len2 = length($crossSelf{$oligo}{$c}[1]);
            if (($lenSelf1 != length($crossSelf{$oligo}{$c}[0])) or ($lenSelfA != length($crossSelf{$oligo}{$c}[1])) or ($dG != $crossSelf{$oligo}{$c}[3])) { #if same
                $allCount{$c} = '';
            }
        }
        if ($in == 1) {
            print $file "\n$one\n$a\n$two\ndG:\t$dG\tkcal/mol\n";
        }
    } elsif ($input =~ /^@/) {
        if ($in == 1) {
            foreach my $c (sort {$a <=> $b} keys %allCount) {
                print $file "\n$crossSelf{$oligo}{$c}[0]\n$crossSelf{$oligo}{$c}[1]\n$crossSelf{$oligo}{$c}[2]\ndG:\t$crossSelf{$oligo}{$c}[3]\tkcal/mol*\n";
                $alldG{$crossSelf{$oligo}{$c}[3]} = '';
            }
            my $dG_cross;
            ($disc, $dG_cross) = split(/\t/, $input);
            $alldG{$dG_cross} = '';
            my $dG_def = 0;
            foreach my $dG (sort {$a <=> $b} keys %alldG) {
                if ($dG < $dG_def) {
                    $dG_def = $dG;
                }
            }
            print $file "\n@\t$dG_def\tkcal/mol\n";
        }
    }
    undef %alldG;
}
close(IN);
close($file);

#re-print hairpin file
open(IN, "<$folder/tmp/SecondaryStructure_hairpin.txt") or die;
open(my $file, ">$folder/hairpin.txt") or die;
my $in = 0;
while(defined(my $input = <IN>)) {
    chomp($input);
    if ($input =~ /^>/) {
        my ($disc, $oligo) = split(/\t/, $input);
        if ((defined($oligo_pass{$oligo})) && (!(defined($allPrimerRemove{$oligo})))) {
            print $file "$input\n";
            $in = 1;
        }
    } else {
        if ($in == 1) {
            print $file "$input\n";
        }
    }
}
close(IN);
close($file);

my %scorePair;

##assign a score to the primer pairs
##taking in consideration:
#have a differed base inside
#have a differed base in the tail
#tm difference betweeen F anf R of less than 1
#haipin dG is higher than -1
#self dimer dG is higher than -1
#cross dimer dG is higher than -1

foreach my $oligo (keys %oligo_pass) {
    if (!(defined($allPrimerRemove{$oligo}))) {
        if (($all{'MAXIMIZE_SEL1'} eq 'yes') or ($all{'MAXIMIZE_SEL2'} eq 'yes')) {
            #start and end of primers
            my $startF = $oligo{$oligo}{'POS'};
            my $endF = $oligo{$oligo}{'POS'} + $oligo{$oligo}{'LEN'};

            #start and end of primer tails

            if ($oligo{$oligo}{'STRAND'}== 'anti') {
                my $startF_tail = $endF - 4; #last 5 bases at 3' end
                my $endF1 = $endF - 1;
                foreach my $d (keys %difference) {
                    if ($all{'MAXIMIZE_SEL2'} eq 'yes') {
                        if (($d >= $startF) && ($d <= $endF)) { #if different bases in F primer
                            $scorePair{$oligo}+=2;
                        }
                    }
                    if ($all{'MAXIMIZE_SEL1'} eq 'yes') {
                        if (($d >= $startF_tail) && ($d <= $endF)) { #if different bases in the F primer tail
                            $scorePair{$oligo} += 2;
                        }
                        #first
                        if ($d == $endF) { #if different bases in the F primer tail
                            $scorePair{$oligo} += 20;
                        }
                        #second
                        if ($d == $endF1) { #if different bases in the F primer tail
                            $scorePair{$oligo} += 10;
                        }
                    }
                }
            } else {
                my $endR_tail = $startF + 4; #last 5 bases at 3' end
                my $startR1 = $startF + 1;
                foreach my $d (keys %difference) {
                    if ($all{'MAXIMIZE_SEL2'} eq 'yes') {
                        if (($d >= $startF) && ($d <= $endR_tail)) { #if different bases in the R primer tail
                            $scorePair{$oligo} += 2;
                        }

                    }
                    if ($all{'MAXIMIZE_SEL1'} eq 'yes') {
                        if (($d >= $endR_tail) && ($d <= $endF)) { #if different bases in the F primer tail
                            $scorePair{$oligo} += 2;
                        }
                        #first
                        if ($d == $startF) { #if different bases in the R primer tail
                            $scorePair{$oligo} += 20;
                        }
                        #second
                        if ($d == $startR1) { #if different bases in the R primer tail
                            $scorePair{$oligo} += 10;
                        }
                    }
                }
            }
        }

        if ($all{'DG_SEL'} eq 'yes') {
            if (($oligo{$oligo}{'HAIR'} >= -1) or ($oligo{$oligo}{'HAIR'} eq '')) {
                $scorePair{$oligo}++;
            }
            if (($oligo{$oligo}{'SELF'} >= -1) or ($oligo{$oligo}{'SELF'} eq '')) {
                $scorePair{$oligo}++;
            }
        }

        if ($all{'DEG_SEL'} eq 'yes') {
            my $deg2 = () = $oligo =~ /R|Y|S|W|K|M/g;
            my $deg3 = () = $oligo =~ /B|D|H|V/g;
            my $deg4 = () = $oligo =~ /N/g;
            $scorePair{$oligo} -= $deg2*2;
            $scorePair{$oligo} -= $deg3*3;
            $scorePair{$oligo} -= $deg4*4;
        }
    }
}

#get the range of score
my %scorePairDef;
my $maxScore = 0;
my $minScore = 1;
foreach my $oligo (keys %scorePair) {
    $scorePairDef{$scorePair{$oligo}}{$oligo} = '';
    if ($scorePair{$oligo} > $maxScore) {
        $maxScore = $scorePair{$oligo};
    }
    if ($scorePair{$oligo} < $minScore) {
        $minScore = $scorePair{$oligo};
    }
}

#create file with all the primers for the blast search
open(my $tmp, ">$folder/tmp/oligo.txt") or die;

my $oligoNum = 500;

my $index = 0;
my %topPair;
foreach (my $score=$maxScore; $score>=$minScore; $score--) {
    foreach my $seq (keys %{$scorePairDef{$score}}) {
        $index++;
        if ($index <= $oligoNum) {
            print $tmp "$index\t$seq\t$oligo{$seq}{'LEN'}\t$oligo{$seq}{'POS'}\toligo\n";
            $topPair{$seq} = '';
        }
    }
}
close($tmp);

#create file containing consensus sequence
open($tmp, ">$folder/tmp/consensus.txt") or die;
print $tmp "$all{'CONSENSUS'}";
close($tmp);

#perform blast + bowtie check
`perl $path_cgi/blast_bowtie_pp.pl -folder $folder -type nt`;

my %accessionMySQL;
my %tableBlast;

my $checkFile = $folder . "/tmp/inSilico_nt.txt";

my $emptyBLAST = 0;
my $tableNt = ">";
my $pieChartTableF;

my %scorePairNew;
my $maxScoreNew = 1;
my $minScoreNew = 1;
my %foundSp;

if (-z $checkFile) { #if file is empty
    $emptyBLAST = 1;
    $tableNt = "none";
    $pieChartTableF = "none";
} else { #if file is not empty

    #retrieve information from BLAST file
    checkBLAST('inSilico_nt.txt');

    #connect to mysql database
    my $dbh;
    my $sth;

    $dbh = DBI->connect ($dsn, $user_name, $password, { RaiseError => 1 });

    #print "start\n";

    my $entry = 0;
    my $ask;
    foreach my $accession (keys %accessionMySQL) {
        if ($entry > 0) {
            $ask .= " OR ";
        }
        $ask .= "(acc='" . $accession . "')";
        $entry++;
    }

    $sth = $dbh->prepare("SELECT * FROM DB2_acc_taxid_pp WHERE ($ask)"); ####nt in mysql

    #execute the prepared statement handle:
    $sth->execute();
    #read results of a query, then clean up
    my %taxid_mysql;
    while (my @ary = $sth->fetchrow_array()) {
        $taxid_mysql{$ary[1]}{$ary[0]} = ''; #taxid - acc
    }
    $sth->finish;
    $entry = 0;
    $ask= '';
    my $insideTaxid = 0;
    foreach my $taxid (keys %taxid_mysql) { ###need to do a unique array - not hash
        $insideTaxid = 1;
        if ($entry > 0) {
            $ask .= " OR ";
        }
        $ask .= "(taxid='" . $taxid . "')";
        $entry++;
    }

    my %taxonomy;
    my %taxonomyAll;

    if ($insideTaxid == 1) { #if at least one corrispondance between accession numbers and taxids
        $sth = $dbh->prepare("SELECT * FROM taxid_taxonomy_pp WHERE ($ask)"); ####create in mysql
        #execute the prepared statement handle:
        $sth->execute();
        #read results of a query, then clean up
        while (my @ary = $sth->fetchrow_array()) {
            foreach my $acc (keys %{$taxid_mysql{$ary[0]}}) { #taxid - acc
                $taxonomy{$acc}{'DOMAIN'} = $ary[1]; #accession - rank = taxonomy
                $taxonomy{$acc}{'PHYLUM'} = $ary[2];
                $taxonomy{$acc}{'CLASS'} = $ary[3];
                $taxonomy{$acc}{'ORDER'} = $ary[4];
                $taxonomy{$acc}{'FAMILY'} = $ary[5];
                $taxonomy{$acc}{'GENUS'} = $ary[6];
                $taxonomy{$acc}{'SPECIES'} = $ary[7];

                $taxonomyAll{'DOMAIN'}{$ary[1]} = '';
                $taxonomyAll{'PHYLUM'}{$ary[2]} = '';
                $taxonomyAll{'CLASS'}{$ary[3]} = '';
                $taxonomyAll{'ORDER'}{$ary[4]} = '';
                $taxonomyAll{'FAMILY'}{$ary[5]} = '';
                $taxonomyAll{'GENUS'}{$ary[6]} = '';
                $taxonomyAll{'SPECIES'}{$ary[7]} = '';
            }
        }
        $sth->finish;
    } else { #if no taxids - all accession numbers are Unclassified
        foreach my $acc (keys %accessionMySQL) {
            $taxonomy{$acc}{'DOMAIN'} = 'Unclassified'; #accession - rank = taxnomy
            $taxonomy{$acc}{'PHYLUM'} = 'Unclassified';
            $taxonomy{$acc}{'CLASS'} = 'Unclassified';
            $taxonomy{$acc}{'ORDER'} = 'Unclassified';
            $taxonomy{$acc}{'FAMILY'} = 'Unclassified';
            $taxonomy{$acc}{'GENUS'} = 'Unclassified';
            $taxonomy{$acc}{'SPECIES'} = 'Unclassified';

            $taxonomyAll{'DOMAIN'}{'Unclassified'} = '';
            $taxonomyAll{'PHYLUM'}{'Unclassified'} = '';
            $taxonomyAll{'CLASS'}{'Unclassified'} = '';
            $taxonomyAll{'ORDER'}{'Unclassified'} = '';
            $taxonomyAll{'FAMILY'}{'Unclassified'} = '';
            $taxonomyAll{'GENUS'}{'Unclassified'} = '';
            $taxonomyAll{'SPECIES'}{'Unclassified'} = '';
        }
    }

    #Javascript colors
    my %colour;
    $colour{1} = 'Blue';
    $colour{2} = 'Wheat';
    $colour{3} = 'BurlyWood';
    $colour{4} = 'CadetBlue';
    $colour{5} = 'Chartreuse';
    $colour{6} = 'Coral';
    $colour{7} = 'CornflowerBlue';
    $colour{8} = 'Crimson';
    $colour{9} = 'DarkBlue';
    $colour{10} = 'DarkCyan';
    $colour{11} = 'DarkGreen';
    $colour{12} = 'DarkMagenta';
    $colour{13} = 'DarkOliveGreen';
    $colour{14} = 'DarkOrange';
    $colour{15} = 'DarkOrchid';
    $colour{16} = 'DarkRed';
    $colour{17} = 'DarkSalmon';
    $colour{18} = 'DarkSeaGreen';
    $colour{19} = 'DarkSlateBlue';
    $colour{20} = 'DarkSlateGray';
    $colour{21} = 'DarkSlateGrey';
    $colour{22} = 'DarkTurquoise';
    $colour{23} = 'DarkViolet';
    $colour{24} = 'DeepSkyBlue';
    $colour{25} = 'DodgerBlue';
    $colour{26} = 'FireBrick';
    $colour{27} = 'ForestGreen';
    $colour{28} = 'Gold';
    $colour{29} = 'GoldenRod';
    $colour{30} = 'Green';
    $colour{31} = 'GreenYellow';
    $colour{32} = 'HotPink';
    $colour{33} = 'IndianRed';
    $colour{34} = 'Indigo';
    $colour{35} = 'Khaki';
    $colour{36} = 'Lavender';
    $colour{37} = 'LightSalmon';
    $colour{38} = 'LawnGreen';
    $colour{39} = 'LemonChiffon';
    $colour{40} = 'LightBlue';
    $colour{41} = 'LightCoral';
    $colour{42} = 'LightCyan';
    $colour{43} = 'LightGreen';
    $colour{44} = 'LightPink';
    $colour{45} = 'Aquamarine';
    $colour{46} = 'LightSeaGreen';
    $colour{47} = 'LightSkyBlue';
    $colour{48} = 'LightSteelBlue';
    $colour{49} = 'Lime';
    $colour{50} = 'LimeGreen';
    $colour{51} = 'Maroon';
    $colour{52} = 'MediumAquaMarine';
    $colour{53} = 'MediumBlue';
    $colour{54} = 'MediumOrchid';
    $colour{55} = 'MediumPurple';
    $colour{56} = 'MediumSeaGreen';
    $colour{57} = 'MediumSlateBlue';
    $colour{58} = 'MediumSpringGreen';
    $colour{59} = 'MediumTurquoise';
    $colour{60} = 'MediumVioletRed';
    $colour{61} = 'MidnightBlue';
    $colour{62} = 'Moccasin';
    $colour{63} = 'Navy';
    $colour{64} = 'BlueViolet';
    $colour{65} = 'OldLace';
    $colour{66} = 'OliveDrab';
    $colour{67} = 'Orange';
    $colour{68} = 'OrangeRed';
    $colour{69} = 'Orchid';
    $colour{70} = 'PaleGoldenRod';
    $colour{71} = 'PaleGreen';
    $colour{72} = 'PaleTurquoise';
    $colour{73} = 'PaleVioletRed';
    $colour{74} = 'Pink';
    $colour{75} = 'Plum';
    $colour{76} = 'PowderBlue';
    $colour{77} = 'Purple';
    $colour{78} = 'RebeccaPurple';
    $colour{79} = 'Red';
    $colour{80} = 'RosyBrown';
    $colour{81} = 'RoyalBlue';
    $colour{82} = 'SaddleBrown';
    $colour{83} = 'Salmon';
    $colour{84} = 'SandyBrown';
    $colour{85} = 'SeaGreen';
    $colour{86} = 'Sienna';
    $colour{87} = 'Silver';
    $colour{88} = 'SkyBlue';
    $colour{89} = 'SlateBlue';
    $colour{90} = 'SlateGrey';
    $colour{91} = 'SpringGreen';
    $colour{92} = 'SteelBlue';
    $colour{93} = 'Tan';
    $colour{94} = 'Teal';
    $colour{95} = 'Thistle';
    $colour{96} = 'Tomato';
    $colour{97} = 'Turquoise';
    $colour{98} = 'Violet';
    $colour{99} = 'Yellow';
    $colour{100} = 'YellowGreen';

    my $count = 0;
    my %taxonomyColour;
    foreach my $rank ('DOMAIN','PHYLUM','CLASS','ORDER','FAMILY','GENUS','SPECIES') {
        foreach my $taxon (keys %{$taxonomyAll{$rank}}) {
            $count++;
            $taxonomyColour{$taxon} = $colour{$count};
            if ($count == 100) { #I have got a list of 100 colours (plus Brown for Unclassified)
                $count = 0;
            }
        }
    }

    foreach (my $score=$maxScore; $score>=1; $score--) {
        foreach my $oligo (keys %{$scorePairDef{$score}}) {

            if (defined($topPair{$oligo})) {
                my $scoreNew = 0;

                if ($all{'SPECIES_SEL'} eq 'yes') {
                    my $selectedRank = 'SPECIES';

                    foreach my $mis (@mismatch) { #print first alignments with less mismatches
                        if (defined($tableBlast{$oligo}{$mis})) {
                            foreach my $acc (keys %{$tableBlast{$oligo}{$mis}}) {
                                if (defined($selectedTaxa{'7'}{$taxonomy{$oligo}{$selectedRank}})) {
                                    if (defined($foundSp{$oligo}{$taxonomy{$oligo}{$selectedRank}})) {
                                        $scoreNew += 5;
                                    } else {
                                        $scoreNew += 20; #incentivate variety
                                        $foundSp{$oligo}{$taxonomy{$acc}{$selectedRank}} = '';
                                    }
                                } else {
                                    $scoreNew -= 40;
                                    $foundSp{$oligo}{$taxonomy{$acc}{$selectedRank}} = '';
                                }
                            }
                        }
                    }
                } else {
                    foreach my $mis (@mismatch) { #print first alignments with less mismatches
                        if (defined($tableBlast{$oligo}{$mis})) {
                            foreach my $acc (keys %{$tableBlast{$oligo}{$mis}}) {
                                $foundSp{$oligo}{$taxonomy{$acc}{'SPECIES'}} = '';
                            }
                        }
                    }
                }

                if ($all{'GENUS_SEL'} eq 'yes') {
                    my $selectedRank = 'GENUS';
                    foreach my $mis (@mismatch) { #print first alignments with less mismatches
                        if (defined($tableBlast{$oligo}{$mis})) {
                            foreach my $acc (keys %{$tableBlast{oligo}{$mis}}) {
                                if (defined($selectedTaxa{'6'}{$taxonomy{$acc}{$selectedRank}})) {
                                    $scoreNew += 20;
                                } else {
                                    $scoreNew -= 40;
                                }
                            }
                        }
                    }
                }

                if ($all{'FAMILY_SEL'} eq 'yes') {
                    my $selectedRank = 'FAMILY';
                    foreach my $mis (@mismatch) { #print first alignments with less mismatches
                        if (defined($tableBlast{$oligo}{$mis})) {
                            foreach my $acc (keys %{$tableBlast{$oligo}{$mis}}) {
                                if (defined($selectedTaxa{'5'}{$taxonomy{$acc}{$selectedRank}})) {
                                    $scoreNew += 20; #incentivate variety
                                } else {
                                    $scoreNew -= 40;
                                }
                            }
                        }
                    }
                }

                if ($all{'ORDER_SEL'} eq 'yes') {
                    my $selectedRank = 'ORDER';
                    foreach my $mis (@mismatch) { #print first alignments with less mismatches
                        if (defined($tableBlast{$oligo}{$mis})) {
                            foreach my $acc (keys %{$tableBlast{$oligo}{$mis}}) {
                                if (defined($selectedTaxa{'4'}{$taxonomy{$acc}{$selectedRank}})) {
                                    $scoreNew += 20; #incentivate variety
                                } else {
                                    $scoreNew -= 40;
                                }
                            }
                        }
                    }
                }

                if ($all{'CLASS_SEL'} eq 'yes') {
                    my $selectedRank = 'CLASS';
                    foreach my $mis (@mismatch) { #print first alignments with less mismatches
                        if (defined($tableBlast{$oligo}{$mis})) {
                            foreach my $acc (keys %{$tableBlast{$oligo}{$mis}}) {
                                if (defined($selectedTaxa{'3'}{$taxonomy{$acc}{$selectedRank}})) {
                                    $scoreNew += 20; #incentivate variety
                                } else {
                                    $scoreNew -= 40;
                                }
                            }
                        }
                    }
                }

                if ($all{'PHYLUM_SEL'} eq 'yes') {
                    my $selectedRank = 'PHYLUM';
                    foreach my $mis (@mismatch) { #print first alignments with less mismatches
                        if (defined($tableBlast{$oligo}{$mis})) {
                            foreach my $acc (keys %{$tableBlast{$oligo}{$mis}}) {
                                if (defined($selectedTaxa{'2'}{$taxonomy{$acc}{$selectedRank}})) {
                                    $scoreNew += 20; #incentivate variety
                                } else {
                                    $scoreNew -= 40;
                                }
                            }
                        }
                    }
                }

                if ($all{'DOMAIN_SEL'} eq 'yes') {
                    my $selectedRank = 'DOMAIN';
                    foreach my $mis (@mismatch) { #print first alignments with less mismatches
                        if (defined($tableBlast{$oligo}{$mis})) {
                            foreach my $acc (keys %{$tableBlast{$oligo}{$mis}}) {
                                if (defined($selectedTaxa{'1'}{$taxonomy{$acc}{$selectedRank}})) {
                                    $scoreNew += 20; #incentivate variety
                                } else {
                                    $scoreNew -= 40;
                                }
                            }
                        }
                    }
                }

                if ($scoreNew == 0) {
                    $scoreNew = -100000;
                }

                $scoreNew += $score;
                $scorePairNew{$scoreNew}{$oligo} = '';
            }
        }
    }

    #get the range of score
    foreach my $score (keys %scorePairNew) {
        if ($score > $maxScoreNew) {
            $maxScoreNew = $score;
        }

        if ($score < $minScoreNew) {
            $minScoreNew = $score;
        }
    }

    #create file with all the primers for the blast search - new for user
    open(my $tmp, ">$folder/tmp/oligo.txt") or die;

    my %pieF;
    my $index = 0;
    my $countAll = 0;
    my $countF = 0;

    foreach (my $score=$maxScoreNew; $score>=$minScoreNew; $score--) {
        foreach my $oligo (keys %{$scorePairNew{$score}}) {
            $index++;

            if ($index <= 100) {

                #print on new oligo.txt
                print $tmp "$index\t$oligo\t$oligo{$oligo}{'LEN'}\t$oligo{$oligo}{'POS'}\toligo\n";

                #create BLAST table & pieCharts
                my $tableF;
                foreach my $mis (@mismatch) { #print first alignments with less mismatches
                    if (defined($tableBlast{$oligo}{$mis})) {
                        foreach my $acc (keys %{$tableBlast{$oligo}{$mis}}) { ###I need to order them somehow

                            $countF++;

                            #data for pieChart
                            $pieF{'DOMAIN'}{$taxonomy{$acc}{'DOMAIN'}}++;
                            $pieF{'PHYLUM'}{$taxonomy{$acc}{'PHYLUM'}}++;
                            $pieF{'CLASS'}{$taxonomy{$acc}{'CLASS'}}++;
                            $pieF{'ORDER'}{$taxonomy{$acc}{'ORDER'}}++;
                            $pieF{'FAMILY'}{$taxonomy{$acc}{'FAMILY'}}++;
                            $pieF{'GENUS'}{$taxonomy{$acc}{'GENUS'}}++;
                            $pieF{'SPECIES'}{$taxonomy{$acc}{'SPECIES'}}++;

                            #accessions present in only oligo primer
                            $tableF .= "<tr><td>" . $acc . "</td>";
                            $tableF .= "<td>" . $oligo . "<br>" .  $tableBlast{$oligo}{$mis}{$acc}{'AL'} . "</td>";
                            $tableF .= "<td>" . $tableBlast{$oligo}{$mis}{$acc}{'START'} . "</td>";
                            $tableF .= "<td>" . $tableBlast{$oligo}{$mis}{$acc}{'END'} . "</td>";
                            #taxonomy
                            if ($taxonomy{$acc}{'DOMAIN'} eq "") {
                                $tableF .= "<td class ='D'>Unclassified</td>";
                            } else {
                                $tableF .= "<td class ='D'>" . $taxonomy{$acc}{'DOMAIN'} . "</td>";
                            }
                            if ($taxonomy{$acc}{'PHYLUM'} eq "") {
                                $tableF .= "<td class ='P'>Unclassified</td>";
                            } else {
                                $tableF .= "<td class ='P'>" . $taxonomy{$acc}{'PHYLUM'} . "</td>";
                            }
                            if ($taxonomy{$acc}{'CLASS'} eq "") {
                                $tableF .= "<td class ='C'>Unclassified</td>";
                            } else {
                                $tableF .= "<td class ='C'>" . $taxonomy{$acc}{'CLASS'} . "</td>";
                            }
                            if ($taxonomy{$acc}{'ORDER'} eq "") {
                                $tableF .= "<td class ='O'>Unclassified</td>";
                            } else {
                                $tableF .= "<td class ='O'>" . $taxonomy{$acc}{'ORDER'} . "</td>";
                            }
                            if ($taxonomy{$acc}{'FAMILY'} eq "") {
                                $tableF .= "<td class ='F'>Unclassified</td>";
                            } else {
                                $tableF .= "<td class ='F'>" . $taxonomy{$acc}{'FAMILY'} . "</td>";
                            }
                            if ($taxonomy{$acc}{'GENUS'} eq "") {
                                $tableF .= "<td class ='G'>Unclassified</td>";
                            } else {
                                $tableF .= "<td class ='G'>" . $taxonomy{$acc}{'GENUS'} . "</td>";
                            }
                            if ($taxonomy{$acc}{'SPECIES'} eq "") {
                                $tableF .= "<td class ='S'>Unclassified</td></tr>";
                            } else {
                                $tableF .= "<td class ='S'>" . $taxonomy{$acc}{'SPECIES'} . "</td>";
                            }
                        }
                    }
                }

                #compose table
                if ($tableF ne '') {
                    $tableNt .= $oligo . "<table>";
                    if ($tableF ne '') {
                        $tableNt .= $tableF;
                    }
                    $tableNt .= "</table>";
                }

                #assemble pieChartTable - F
                $pieChartTableF .= "<" . $oligo . ">";

                foreach my $rank ('DOMAIN','PHYLUM','CLASS','ORDER','FAMILY','GENUS','SPECIES') {
                    foreach my $taxon (keys %{$pieF{$rank}}) {
                        if ($taxon eq "") {
                            $pieChartTableF .= "Unclassified-" . $pieF{$rank}{$taxon} . "-Brown,";
                        } else {
                            $pieChartTableF .= $taxon . "-" . $pieF{$rank}{$taxon} . "-" . $taxonomyColour{$taxon} . ",";
                        }
                    }
                    chop($pieChartTableF);
                    $pieChartTableF .= ":";
                }
                chop($pieChartTableF);
                $pieChartTableF .= "<>";
                undef %pieF;
            }
        }
    }
}
close($tmp);
undef %tableBlast;
undef %accessionMySQL;

#blast if negative user file is present
my $tableUserCheck = ">";

if ($all{'NEGATIVE_FILE'} eq 'yes') {

    #perform blast + bowtie check
    `perl $path_cgi/blast_bowtie_pp.pl -folder $folder -type user`;
    my $checkFile = $folder . "/tmp/inSilico_user.txt";

    if (-z $checkFile) { #if file is empty
        $tableUserCheck = "none";
    } else { #if file is not empty

        #retrieve information from BLAST file
        checkBLAST('inSilico_user.txt');

        my $negative = (substr($folder, 56, 8)) . ".negativefasta";
        my $countAll = `grep -c "^>" $folder/$negative`;
        chomp($countAll);

        foreach (my $score=$maxScoreNew; $score>=$minScoreNew; $score--) {
            foreach my $oligo (keys %{$scorePairNew{$score}}) {

                #create BLAST table & pieCharts
                my $tableF;

                my $countF = 0;

                foreach my $mis (@mismatch) { #print first alignments with less mismatches
                    if (defined($tableBlast{$oligo}{$mis})) {
                        foreach my $acc (keys %{$tableBlast{$oligo}{$mis}}) { ###I need to order them somehow
                            $countF++;
                            #accessions present in only oligo primer
                            $tableF .= "<tr><td>" . $acc . "</td>";
                            $tableF .= "<td>" . $oligo . "<br>" .  $tableBlast{$oligo}{$mis}{$acc}{'AL'} . "</td>";
                            $tableF .= "<td>" . $tableBlast{$oligo}{$mis}{$acc}{'START'} . "</td>";
                            $tableF .= "<td>" . $tableBlast{$oligo}{$mis}{$acc}{'END'} . "</td>";
                        }
                    }
                }

                #compose table
                if ($tableF ne '') {
                    $tableUserCheck .= $oligo . "<" . $countAll . ";" . $countF . "><table>";
                    if ($tableF ne '') {
                        $tableUserCheck .= $tableF;
                    }
                    $tableUserCheck .= "</table>";
                }
            }
        }
    }
} else {
    $tableUserCheck = "no";
}

#print the first 100 primer pairs
open(my $file, ">$folder/info.primer") or die; #file with all the necessary parameters - result page
open(my $file1, ">$folder/results/oligoList.txt") or die; #file with all the necessary parameters - user data

print $file "PROJECT\t$all{'PROJECT'}\n";
print $file "CONSENSUS\t$all{'CONSENSUS'}\n";
print $file "DIFFERENT_POS\t$all{'DIFFERENT_POS'}\n";
#print selection criteria
my @vis_sel = ('DG_SEL', 'DEG_SEL', 'MAXIMIZE_SEL1', 'MAXIMIZE_SEL2', 'SPECIES_SEL', 'GENUS_SEL', 'FAMILY_SEL', 'ORDER_SEL', 'CLASS_SEL', 'PHYLUM_SEL', 'DOMAIN_SEL');
foreach my $v (@vis_sel) {
    if ($all{$v} eq "yes") {
        print $file "$v\t$all{$v}\n";
    }
}
print $file "\n";

my %bestPair;
my $index = 0;
foreach (my $score=$maxScoreNew; $score>=$minScoreNew; $score--) {
    foreach my $oligo (keys %{$scorePairNew{$score}}) {
        $index++;
        $bestPair{$oligo} = '';
        if ($index <= 100) {
            print $file "INDEX\t$index\n";
            print $file "OLIGO\t$oligo\n";
            print $file "SCORE\t$score\n";
            print $file "POSITION\t$oligo{$oligo}{'POS_NAME'}\n";
            print $file "LENGTH\t$oligo{$oligo}{'LEN'}\n";
            print $file "GC\t$oligo{$oligo}{'GC'}\n";
            print $file "TM\t$oligo{$oligo}{'TM'}\n";
            print $file "STRAND\t$oligo{$oligo}{'STRAND'}\n";
            if ($oligo{$oligo}{'SELF'} == 0) {
                print $file "SELF\t>=0\n";
            } else {
                print $file "SELF\t$oligo{$oligo}{'SELF'}\n";
            }
            if ($oligo{$oligo}{'HAIR'} == 0) {
                print $file "HAIR\t>=0\n";
            } else {
                print $file "HAIR\t$oligo{$oligo}{'HAIR'}\n";
            }
            my $allSp;
            foreach my $sp (keys %{$foundSp{$oligo}}) {
                $allSp .= $sp . ";";
            }
            print $file "SPECIES\t$allSp\n\\\\\n";
        }

        print $file "INDEX\t$index\n";
        print $file1 "OLIGO\t$oligo\n";
        print $file1 "SCORE\t$score\n";
        print $file1 "POSITION\t$oligo{$oligo}{'POS_NAME'}\n";
        print $file1 "LENGTH\t$oligo{$oligo}{'LEN'} bases\n";
        print $file1 "GC\t$oligo{$oligo}{'GC'} %\n";
        print $file1 "TM\t$oligo{$oligo}{'TM'} °C\n";
        print $file1 "STRAND\t$oligo{$oligo}{'STRAND'}\n";
        if ($oligo{$oligo}{'SELF'} == 0) {
            print $file1 "SELF\t>=0 kcal/mol\n";
        } else {
            print $file1 "SELF\t$oligo{$oligo}{'SELF'} kcal/mol\n";
        }
        if ($oligo{$oligo}{'HAIR'} == 0) {
            print $file1 "HAIR\t>=0 kcal/mol\n";
        } else {
            print $file1 "HAIR\t$oligo{$oligo}{'HAIR'} kcal/mol\n";
        }
        my $allSp;
        foreach my $sp (keys %{$foundSp{$oligo}}) {
            $allSp .= $sp . ";";
        }
        $allSp =~ s/\[sub1\]/\'/g; #single quote
        $allSp =~ s/\[sub2\]/,/g; #comma
        $allSp =~ s/\[sub3\]/\(/g; #bracket (
        $allSp =~ s/\[sub4\]/\)/g; #bracket )
        $allSp =~ s/\[sub5\]/:/g; #column
        $allSp =~ s/\[sub6\]/;/g; #semi column
        $allSp =~ s/\[sub7\]/\*/g; #semi column
        $allSp =~ s/\[sub8\]/</g; #lower
        $allSp =~ s/\[sub9\]/>/g; #higher
        $allSp =~ s/\[sub10\]/-/g; #minus
        $allSp =~ s/\[sub11\]/\+/g; #plus
        $allSp =~ s/\[sub12\]/\`/g; #hyphen`
        $allSp =~ s/\[sub13\]/\#/g; #
        $allSp =~ s/\[sub14\]/&/g; #&
        $allSp =~ s/\[sub15\]/\^/g; #&
        $allSp =~ s/\[sub16\]/\//g; #/
        $allSp =~ s/\[sub17\]/_/g; #_
        print $file1 "SPECIES\t$allSp\n\\\\\n";
    }
}

foreach (my $score=$maxScore; $score>=$minScore; $score--) {
    foreach my $oligo (keys %{$scorePairDef{$score}}) {
        if (!(defined($bestPair{$oligo}))) {
            print $file "INDEX\t$index\n";
            print $file1 "OLIGO\t$oligo\n";
            print $file1 "SCORE\t$score\n";
            print $file1 "POSITION\t$oligo{$oligo}{'POS_NAME'}\n";
            print $file1 "LENGTH\t$oligo{$oligo}{'LEN'} bases\n";
            print $file1 "GC\t$oligo{$oligo}{'GC'} %\n";
            print $file1 "TM\t$oligo{$oligo}{'TM'} °C\n";
            print $file1 "STRAND\t$oligo{$oligo}{'STRAND'}\n";
            if ($oligo{$oligo}{'SELF'} == 0) {
                print $file1 "SELF\t>=0 kcal/mol\n";
            } else {
                print $file1 "SELF\t$oligo{$oligo}{'SELF'} kcal/mol\n";
            }
            if ($oligo{$oligo}{'HAIR'} == 0) {
                print $file1 "HAIR\t>=0 kcal/mol\n";
            } else {
                print $file1 "HAIR\t$oligo{$oligo}{'HAIR'} kcal/mol\n";
            }
            print $file1 "SPECIES\tnot BLASTed\n";
            print $file1 "\\\\\n";
        }
    }
}

print $file "BLAST_TABLE\t$tableNt\n";
print $file "PIECHART_F\t$pieChartTableF\n";
print $file "USER_TABLE\t$tableUserCheck\n";
close($file);
close($file1);

#send email
my $to = $all{'EMAIL'};
my $from = 'gv16363@bristol.ac.uk';
my $subject = 'PhyloPrimer results - ' . $all{'PROJECT'};

#deconstruct folder name

$folder =~ s/\/var\/www\/cerealgenomics\/phyloprimer\/analysesPhyloprimer\///g;

my $input_kind = chop($folder); #get last letter and understand if a, b or c
my $input_ST = chop($folder); #get last letter and understand if T or S
my $defSet = $folder . $input_kind . $input_ST;
$folder = $folder . $input_ST . $input_kind;

`$path_cgi/Create_README_pp.pl -folder $folder -file $fileInfo`; #create README.txt

`cp $path_html/analysesPhyloprimer/${folder}/selfDimer.txt $path_html/analysesPhyloprimer/${folder}/results/selfDimer.txt`;
`cp $path_html/analysesPhyloprimer/${folder}/hairpin.txt $path_html/analysesPhyloprimer/${folder}/results/hairpin.txt`;

`mv $path_html/analysesPhyloprimer/${folder}/${nameFile}* $path_html/analysesPhyloprimer/${folder}/inputs/`;

#for phylogenetic tree visualisation
`mv $path_html/analysesPhyloprimer/${folder}/inputs/${nameFile}.treeInfo $path_html/analysesPhyloprimer/${folder}/`;
`mv $path_html/analysesPhyloprimer/${folder}/inputs/${nameFile}.allAccession $path_html/analysesPhyloprimer/${folder}/`;

my $folderAll = $path_html . "/analysesPhyloprimer/" . $folder;

#move bowtie files
`mv $path_html/analysesPhyloprimer/${folder}/tmp/oligo_bowtie*.sam $path_html/analysesPhyloprimer/${folder}/results/`;

#remove tmp folder
#`rm -r ${folderAll}/tmp`;

#zip the folder
`zip -r ${folderAll}/PhyloPrimer_${nameFile}.zip ${folder} -x ${folder}/*txt -x ${folder}/info.primer -x ${folder}/${nameFile}.treeInfo -x ${folder}/${nameFile}.allAccession -x "${folder}/tmp/*" -x ${folder}/results/${nameFile}.loadAcc`;

`chown www-data:www-data ${folderAll}/info.primer`;
`chown www-data:www-data ${folderAll}/*txt`;
`chown www-data:www-data ${folderAll}/PhyloPrimer_${nameFile}.zip`;

my $message = "Hi,\n\nPlease find the link to the PhyloPrimer results: https:\/\/www.cerealsdb.uk.net\/cerealgenomics\/cgi-bin\/phyloprimerResultsOligo.cgi?defSet=" . $defSet . "\n\n\nAll the best,\nPhyloPrimer team";

open(MAIL, "|/usr/sbin/sendmail -t");

# Email Header
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: $subject\n\n";
# Email Body
print MAIL $message;

close(MAIL);

#subroutines
sub gplusc {
    my ($seq) = $_[0];
    my $n = () = $seq =~ /G|C/g; #1
    if ($seq =~ /[RYSWKMBDHVN]/) { #if degerenate bases
        my $n1 = () = $seq =~ /S/g; #+1
        my $n2 = () = $seq =~ /R|Y|K|M|N/g; #+0.5
        my $n3 = () = $seq =~ /D|H/g; #+0.33
        my $n4 = () = $seq =~ /B|V/g; #+0.66
        $n +=  $n1 + $n2*0.5 + $n3*0.33 + $n4*0.66;
    }
    my $gc = sprintf("%.1f", (100*($n/length($seq))));
    return($gc);
}

sub redundancy {
    my $primer = shift;
    my $fold = 0;
    my $num = 0;
    #Multiply the fold redundancy by the factor relevent to each ambiguity code found. e.g. Y = C/T = x2 whilst N = CATG = x4.
    while($primer =~ /[YRMWKS]/g){
        $fold *= 2;
        $num++;
    }
    while($primer =~ /[VHDB]/g){
        $fold *= 3;
        $num++;
    }
    while($primer =~ /N/g){
        $fold *= 4;
        $num++;
    }
    return ($fold,$num); #fold and number of wildcards
}

sub degenerateAlt { #if $primer_input has degenerate bases I need to retrieve all the possible alternatives
    my ($primer) = $_[0];

    my %degenerate;
    my %degenerateNew;
    my $count = 0;
    my $inside = 0;
    my @all;
    my $primer0;

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
    foreach my $c (keys %degenerate) {
        push @all, $degenerate{$c};
    }
    my $allRef = \@all;
    return($allRef);
}

#retrieve information from BLAST+BOWTIE results
#populate %accessionMySQL
#populate %tableBlast
sub checkBLAST {

    my ($fileBLAST) = $_[0];

    open(IN, "<$folder/tmp/$fileBLAST") or die; #blast+bowtie file

    my %tableBlastTMP;

    while(defined(my $input = <IN>)) {
        chomp($input); ####not same results
        my ($mis, $pair, $acc, $al, $start, $end) = split(/\t/, $input);
        $tableBlastTMP{$pair}{$mis}{$acc}{'START'} = $start;
        $tableBlastTMP{$pair}{$mis}{$acc}{'END'} = $end;
        $tableBlastTMP{$pair}{$mis}{$acc}{'AL'} = $al;
    }
    close(IN);

    #same a maximumum of 100 entries for each oligo pairs
    foreach my $pair (keys %tableBlastTMP) {
      my $indexTable = 0;
      foreach my $mis (@mismatch) { #print first alignments with less mismatches
        if (defined($tableBlastTMP{$pair}{$mis})) {
          foreach my $acc (sort keys %{$tableBlastTMP{$pair}{$mis}}) { ###I need to order them somehow
              $indexTable++;
              if ($indexTable <= 200) { #show only the first 100 BLAST+bowtie matches
                $accessionMySQL{$acc} = "";
                $tableBlast{$pair}{$mis}{$acc}{'START'} = $tableBlastTMP{$pair}{$mis}{$acc}{'START'};
                $tableBlast{$pair}{$mis}{$acc}{'END'} = $tableBlastTMP{$pair}{$mis}{$acc}{'END'};
                $tableBlast{$pair}{$mis}{$acc}{'AL'} = $tableBlastTMP{$pair}{$mis}{$acc}{'AL'};
              }
            }
          }
        }
      }
}
