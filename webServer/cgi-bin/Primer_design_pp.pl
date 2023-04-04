#!/usr/bin/perl

use strict;
use DBD::mysql;
use POSIX;

#connected to backgroundDesignCheck.pl : primer design pipeline

#set path to folders:
my $path_cgi = 'path_to_cgi_folder';
my $path_html = 'path_to_html_folder';

#set mySQL parameters
my $dsn = "mysql_database";
my $user_name = "mysql_user";
my $password = "mysql_password";

#Usage: ./Primer_design.pl -folder $path_html/analysesPhyloprimer/AkuVricWqWoYjmWDMStmTa -file AkuVricW.info

##NOTES:
#Tm would be the same for forward and reverse if I was not considering dangling ends
#Indicated position of forward and reverse must always be the position at 5' end so F = position and R = position+length

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
  if (($all{'HIGHF'} ne "") && ($all{'HIGHR'} eq "")) { #if only forward area defined
  $area{'F'} = $all{'HIGHF'};
  $all{'HIGHR_START'} = $all{'HIGHF_START'} + $all{'LEN_MIN'} + $all{'AMPL_SIZE_MIN'};
  $all{'HIGHR_END'} = $all{'HIGHF_END'} + $all{'AMPL_SIZE_MAX'} + $all{'LEN_MAX'};
  if ($all{'HIGHR_END'} > $lenConsensus) {
    $all{'HIGHR_END'} = $lenConsensus;
  }
  my $diff = $all{'HIGHR_END'} - $all{'HIGHR_START'};
  $area{'R'} = substr($all{'CONSENSUS'}, ($all{'HIGHR_START'} -1), ($diff+1));
}
if (($all{'HIGHF'} eq "") && ($all{'HIGHR'} ne "")) { #if only reverse is defined
$all{'HIGHF_START'} = $all{'HIGHR_START'} - $all{'AMPL_SIZE_MAX'} - $all{'LEN_MAX'};
$all{'HIGHF_END'} = $all{'HIGHR_END'} - $all{'LEN_MIN'} - $all{'AMPL_SIZE_MIN'};
my $diff = $all{'HIGHF_END'} - $all{'HIGHF_START'};
$area{'F'} = substr($all{'CONSENSUS'}, ($all{'HIGHF_START'} -1), ($diff+1));
$area{'R'} = $all{'HIGHR'};
}
if (($all{'HIGHF'} ne "") && ($all{'HIGHR'} ne "")) { #if both forward and reverse area defined
if (($all{'HIGHR_START'} - $all{'HIGHF_START'} - $all{'LEN_MAX_PRIMER'}) > $all{'AMPL_SIZE_MAX'}) {
  $all{'HIGHF_START'} = $all{'HIGHR_START'} - $all{'AMPL_SIZE_MAX'} - $all{'LEN_MAX'};
  my $diff = $all{'HIGHF_END'} - $all{'HIGHF_START'} + 1;
  $all{'HIGHF'} = substr($all{'HIGHF'}, -($diff));
}
if (($all{'HIGHR_END'} - $all{'HIGHF_END'} - $all{'LEN_MAX_PRIMER'}) > $all{'AMPL_SIZE_MAX'}) {
  $all{'HIGHR_END'} = $all{'HIGHF_END'} + $all{'AMPL_SIZE_MAX'} + $all{'LEN_MAX'};
  my $diff = $all{'HIGHR_END'} - $all{'HIGHR_START'} + 1;
  $all{'HIGHR'} = substr($all{'HIGHR'}, 0, $diff);
}
$area{'F'} = $all{'HIGHF'};
$area{'R'} = $all{'HIGHR'};
}
} elsif ($all{'CONSERVED'} eq "no") {
  $area{'ALL'} = $all{'CONSENSUS'}; #if none of the areas was specified by the user
}

my %selectedTaxa;

my $nameFileSpecies = $nameFile . ".taxonomy";
my $check = $folder . "/" . $nameFileSpecies;
if (-e $check) {
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

my %oligo;
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
    if ($countEX > 10000) { #exit this loop if more than 10000 repetitions, otherwise %exclude takes too much RAM. PhyloPrimer gets many repetitions when many degenerate bases in consensus
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
  if ($countEX > 10000) { #exit this loop if more than 10000 repetitions, otherwise %exclude takes too much RAM. PhyloPrimer gets many repetitions when many degenerate bases in consensus
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
    } elsif ($type eq 'R') {
      $pos = $all{'HIGHR_START'} -1;
      $l = $all{'HIGHR_START'} + length($area{$type});
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
                    $oligo{$type}{$final} = ""; #type - primer = ""
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

my %forwardTMP;
my %reverseTMP;
my %oligoAll;

open(my $fileTM, ">$folder/tmp/tm.tmp") or die;

my $presence2 = 0;
foreach my $type (sort keys %oligo) {
  if ($type eq "ALL") {
    foreach my $primer (keys %{$oligo{$type}}) {
      my ($pos, $seq) = split(/-/, $primer);
      my $gc = gplusc($seq); #check for gc content

      #CHECK 4 - GC content
      if ($gc >= $all{'GC_CONTENT_MIN'} && $gc <= $all{'GC_CONTENT_MAX'}) {
        my $rev = $seq; #reverse seq
        $rev = reverse($rev);
        $rev =~ tr/CATGYRKMBVDH/GTACRYMKVBHD/;

        #CHECK 5 - GC clamp (for)
        my $tail = substr($seq, -5);
        my $tail_gc = gplusc($tail);
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
            #CHECK 6 - melting temperature
            my $plength = length($seq);
            my $dang = substr($all{'CONSENSUS'}, ($pos+$plength-1), 1);
            my $plength = length($seq);
            print $fileTM "$seq\tprimer\tF\t$dang\n";
            $presence2++;
            #forward information
            $forwardTMP{$seq}{'POS'} = $pos;
            $forwardTMP{$seq}{'POS_NAME'} = $pos;
            $forwardTMP{$seq}{'LEN'} = $plength;
            $forwardTMP{$seq}{'GC'} = $gc;
            $forwardTMP{$seq}{'RED'} = $redundancy;
          }
        }

        #CHECK 5 - GC clamp (rev)
        my $tail = substr($rev, -5);
        my $tail_gc = gplusc($tail);
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
            #CHECK 6 - melting temperature
            my $dang = substr($all{'CONSENSUS'}, ($pos-2), 1);
            my $plength = length($rev);
            print $fileTM "$rev\tprimer\tR\t$dang\n";
            $presence2++;
            #reverse information
            $reverseTMP{$rev}{'POS'} = $pos;
            $reverseTMP{$rev}{'POS_NAME'} = $pos + $plength -1;
            $reverseTMP{$rev}{'LEN'} = $plength;
            $reverseTMP{$rev}{'GC'} = $gc;
            $reverseTMP{$rev}{'RED'} = $redundancy;
          }
        }
      }
    }
  } elsif ($type eq "F") {
    foreach my $primer (keys %{$oligo{$type}}) {
      my ($pos, $seq) = split(/-/, $primer);
      my $gc = gplusc($seq);

      #CHECK 4 - GC content
      if ($gc >= $all{'GC_CONTENT_MIN'} && $gc <= $all{'GC_CONTENT_MAX'}) {

        #CHECK 5 - GC clamp
        my $tail = substr($seq, -5);
        my $tail_gc = gplusc($tail);
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
              #CHECK 6 - melting temperature
              my $plength = length($seq);
              my $dang = substr($all{'CONSENSUS'}, ($pos+$plength-1), 1);
              my $plength = length($seq);
              print $fileTM "$seq\tprimer\tF\t$dang\n";
              $presence2++;
              #forward information
              $forwardTMP{$seq}{'POS'} = $pos;
              $forwardTMP{$seq}{'POS_NAME'} = $pos;
              $forwardTMP{$seq}{'LEN'} = $plength;
              $forwardTMP{$seq}{'GC'} = $gc;
              $forwardTMP{$seq}{'RED'} = $redundancy;
          }
        }
      }
    }
  } elsif ($type eq "R") {
    foreach my $primer (keys %{$oligo{$type}}) {
      my ($pos, $seq) = split(/-/, $primer);
      my $gc = gplusc($seq);

      #CHECK 4 - GC content
      if ($gc >= $all{'GC_CONTENT_MIN'} && $gc <= $all{'GC_CONTENT_MAX'}) {
        my $rev = $seq; #reverse seq
        $rev = reverse($rev);
        $rev =~ tr/CATGYRKMBVDH/GTACRYMKVBHD/;

        #CHECK 5 - GC clamp
        my $tail = substr($rev, -5);
        my $tail_gc = gplusc($tail);


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
              #CHECK 6 - melting temperature
              my $dang = substr($all{'CONSENSUS'}, ($pos-2), 1);
              my $plength = length($rev);
              print $fileTM "$rev\tprimer\tR\t$dang\n";
              $presence2++;
              #reverse information
              $reverseTMP{$rev}{'POS'} = $pos;
              $reverseTMP{$rev}{'POS_NAME'} = $pos + $plength -1;
              $reverseTMP{$rev}{'LEN'} = $plength;
              $reverseTMP{$rev}{'GC'} = $gc;
              $reverseTMP{$rev}{'RED'} = $redundancy;
          }
        }
      }
    }
  }
}
close($fileTM);

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

#start tm calculation
my $size = ceil($presence2/10); #divided into max 5 files

my $cross = 0;
my $in = 0;

open(IN, "<$folder/tmp/tm.tmp") or die; #information on tm
open(my $tmp, ">$folder/tmp/tm_1.tmp") or die; #first subset file
my $fileName = "tm_1.tmp";

#read file
while (defined(my $input = <IN>)) {
  chomp($input);
  $cross++;
  if (($cross % $size) == 0) {
    print $tmp "$input\n";
    close($tmp);
    system("perl $path_cgi/tm_calculation_file_pp.pl -file $fileName -folder $folder -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot &");
    $fileName = "tm_" . ${cross} . ".tmp";
    open($tmp, ">$folder/tmp/$fileName") or die;
    $in = 0;
  } else {
    print $tmp "$input\n";
    $in = 1;
  }
}
close(IN);
close($tmp);

#remove original tm file
`rm $folder/tmp/tm.tmp`;

if ($in == 1) {
    system("perl $path_cgi/tm_calculation_file_pp.pl -file $fileName -folder $folder -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot &");
} else {
    #remove last empty file
    `rm $folder/tmp/$fileName`;
}

#stall till all tm parallel calculations are not finished
my $done = 0;
while ($done == 0) {
    sleep 30;
    my $list = `ls $folder/tmp/tm_*.tmp 2> /dev/null`;
    if ($list eq '') {
        $done = 1;
    }
}

#merge aal tm files
`cat $folder/tmp/result_tm_*.tmp > $folder/tmp/result_tm.tmp`;

my $forCheck = keys %forwardTMP;
my $revCheck = keys %reverseTMP;
my $dataCheck = localtime();
print "how many $forCheck and $revCheck and $dataCheck\n";

#open tm file
my %forward;
my %reverse;
my %forward_pass;
my %reverse_pass;
open(IN, "<$folder/tmp/result_tm.tmp") or die;
while(defined(my $input = <IN>)) {
    chomp($input);
    my @info = split(/\t/, $input);
    if ($info[2] eq 'F') {
        my $tm = sprintf("%.2f", $info[4]);
        if (($tm >= $all{'TM_MIN'}) && ($tm <= $all{'TM_MAX'})) {
            $forward{$info[0]}{'POS'} = $forwardTMP{$info[0]}{'POS'};
            $forward{$info[0]}{'POS_NAME'} = $forwardTMP{$info[0]}{'POS_NAME'};
            $forward{$info[0]}{'LEN'} = $forwardTMP{$info[0]}{'LEN'};
            $forward{$info[0]}{'GC'} = $forwardTMP{$info[0]}{'GC'};
            $forward{$info[0]}{'RED'} = $forwardTMP{$info[0]}{'RED'};
            $forward{$info[0]}{'TM'} = $tm;
            $forward_pass{$forwardTMP{$info[0]}{'POS'}}{$info[0]} = '';
        }
    } elsif ($info[2] eq 'R') {
        my $tm = sprintf("%.2f", $info[4]);
        if (($tm >= $all{'TM_MIN'}) && ($tm <= $all{'TM_MAX'})) {
            $reverse{$info[0]}{'POS'} = $reverseTMP{$info[0]}{'POS'};
            $reverse{$info[0]}{'POS_NAME'} = $reverseTMP{$info[0]}{'POS_NAME'};
            $reverse{$info[0]}{'LEN'} = $reverseTMP{$info[0]}{'LEN'};
            $reverse{$info[0]}{'GC'} = $reverseTMP{$info[0]}{'GC'};
            $reverse{$info[0]}{'RED'} = $reverseTMP{$info[0]}{'RED'};
            $reverse{$info[0]}{'TM'} = $tm;
            $reverse_pass{$reverseTMP{$info[0]}{'POS'}}{$info[0]} = '';
        }
    }
}
close(IN);

my $forCheck1 = keys %forward;
my $revCheck1 = keys %reverse;
my $dataCheck1 = localtime();
print "how many1 $forCheck1 and $revCheck1 and $dataCheck1\n";

#Create pair $pairs (check Ta and diff Tm)
my $pair =0;
my $Ta = 0;
my $tm_diff;
my $presence3 = 0;
my $presence4 = 0;

my %for_to_rev;
my %rev_to_for;
my %for_to_count;
my %rev_to_count;

foreach my $pos_F (sort {$a<=>$b} keys %forward_pass) {
    foreach my $pos_R (sort {$a<=>$b} keys %reverse_pass) {
        $presence3++;
        my $diff = $pos_R - $pos_F;
        #CHECK 7 - amplicon size
        if (($diff > $all{'AMPL_SIZE_MIN'}) && ($diff < $all{'AMPL_SIZE_MAX'})) {
            foreach my $primer_F (sort {$a<=>$b} keys %{$forward_pass{$pos_F}}) {
                foreach my $primer_R (sort {$a<=>$b} keys %{$reverse_pass{$pos_R}}) {
                    if ($reverse{$primer_R}{'TM'} >= $forward{$primer_F}{'TM'}) {
                        $Ta = $forward{$primer_F}{'TM'} - 5;
                        $tm_diff = $reverse{$primer_R}{'TM'} - $forward{$primer_F}{'TM'};
                    } else {
                        $Ta = $reverse{$primer_R}{'TM'} - 5;
                        $tm_diff = $reverse{$primer_R}{'TM'} - $forward{$primer_F}{'TM'};
                    }
                    #CHECK 8 - annealing temperature and Tm difference between primers
                    if (($Ta >= $all{'TA_MIN'}) && ($Ta <= $all{'TA_MAX'}) && ($tm_diff <= $all{'TM_DIFF'})) {
                        $presence4++;
                        $for_to_rev{$primer_F}{$primer_R} = ""; #forward - reverse
                        $rev_to_for{$primer_R}{$primer_F} = ""; #reverse - forward
                        $for_to_count{$primer_F}++; #reverse = count
                        $rev_to_count{$primer_R}++; #reverse = count
                    }
                }
            }
        }
    }
}

if ($presence3 == 0) {
    $diffPres = 0;
} else {
    $diffPres = $presence4/$presence3;
}

if ($diffPres < 0.5) {
    push @diffMessage, '<li>amplicon size</li>';
}

if (($presence3 == 0) or ($presence4 == 0)) {
    `mv ${folder}/${nameFile}* ${folder}/inputs/`;
    #send email
    my $to = $all{'EMAIL'};
    my $from = 'gv16363@bristol.ac.uk';
    my $subject = 'PhyloPrimer results - ' . $all{'PROJECT'};
    my $message = "Content-Type: text/html; charset=ISO-8859-1\n\n<html><body>Hi,<br>PhyloPrimer did not find any suitable primer pairs with the selected parameters. It looks like the parameters for the following fields were too stringent:<br><ul>@diffMessage</ul><br><br>All the best,<br>Gilda</body>";

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

my $forCheck2 = keys %for_to_count;
my $revCheck2 = keys %rev_to_count;
my $dataCheck2 = localtime();
print "how many2 $forCheck2 and $revCheck2 and $dataCheck2\n";

#CHECK 9 - hairpins for forward oligo
my $size_for = ceil((scalar keys %for_to_count) / 10); #number of forward oligos

my $hair = 0;
my $in = 0;

open(my $tmp, ">$folder/tmp/Hairpin_for_1.tmp") or die;
my $fileName = "Hairpin_for_1.tmp";

foreach my $primer (keys %for_to_count) {
    $hair++;
    if (($hair % $size_for) == 0) {
        print $tmp "$primer\n";
        close($tmp);
        system("perl $path_cgi/Hairpin_checking_web_pp.pl -folder $folder -primer $fileName -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius &");
        $fileName = "Hairpin_for_" . $hair . ".tmp";
        open($tmp, ">$folder/tmp/Hairpin_for_${hair}.tmp") or die;
        $in = 0;
    } else {
        print $tmp "$primer\n";
        $in = 1;
    }
}
close($tmp);
if ($in == 1) {
    system("perl $path_cgi/Hairpin_checking_web_pp.pl -folder $folder -primer $fileName -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius &");
} else {
    #remove last empty file
    `rm $folder/tmp/$fileName`;
}

#stall till has not finished
my $done = 0;
while ($done == 0) {
    sleep 30;
    my $list = `ls $folder/tmp/Hairpin_for_*tmp 2> /dev/null`;
    if ($list eq '') {
        $done = 1;
    }
}
`cat $folder/tmp/SecondaryStructure_Hairpin_*.tmp > $folder/tmp/SecondaryStructure_hairpin.txt`;

#read forward oligos & retrieve only reverse oligos
`cat $folder/tmp/SecondaryStructure_Hairpin_for*.tmp.delta > $folder/tmp/Hairpin_delta_for.tmp`;
open(IN, "<$folder/tmp/Hairpin_delta_for.tmp") or die;
while(defined(my $input = <IN>)) { #only for forward
    chomp($input);
    my ($primer_F, $dG_hair) = split(/\t/, $input); #what is ) or null?!?
    if ($dG_hair < $all{'OLI_HAIRPIN'}) {
        #retrieve only reverse primers that have valid hairpin values for correspondent forward oligo
        foreach my $primer_R (keys %{$for_to_rev{$primer_F}}) {
            $for_to_count{$primer_F} = 0;
            $rev_to_count{$primer_R}--;
        }
    }
}
close(IN);
my %oligo_retrieve;
foreach my $primer (keys %rev_to_count) {
    if ($rev_to_count{$primer} > 0) {
        $oligo_retrieve{$primer} = '';
    }
}
`rm $folder/tmp/SecondaryStructure_Hairpin*`;

my $revCheck3 = keys %oligo_retrieve;
my $dataCheck3 = localtime();
print "how many3 $revCheck3 and $dataCheck3\n";

#CHECK 10 - hairpins for reverse oligos
my $size_rev = ceil((scalar keys %oligo_retrieve) / 10); #number of reverse oligos

my $hair = 0;
my $in = 0;

open(my $tmp, ">$folder/tmp/Hairpin_rev_1.tmp") or die;
my $fileName = "Hairpin_rev_1.tmp";

foreach my $primer (keys %oligo_retrieve) {
    $hair++;
    if (($hair % $size_rev) == 0) {
        print $tmp "$primer\n";
        close($tmp);
        system("perl $path_cgi/Hairpin_checking_web_pp.pl -folder $folder -primer $fileName -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius &");
        $fileName = "Hairpin_rev_" . $hair . ".tmp";
        open($tmp, ">$folder/tmp/Hairpin_rev_${hair}.tmp") or die;
        $in = 0;
    } else {
        print $tmp "$primer\n";
        $in = 1;
    }
}
close($tmp);
if ($in == 1) {
    system("perl $path_cgi/Hairpin_checking_web_pp.pl -folder $folder -primer $fileName -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius &");
} else {
    #remove last empty file
    `rm $folder/tmp/$fileName`;
}

#stall till has not finished
my $done = 0;
while ($done == 0) {
    sleep 30;
    my $list = `ls $folder/tmp/Hairpin_rev_*tmp 2> /dev/null`;
    if ($list eq '') {
        $done = 1;
    }
}
`cat $folder/tmp/SecondaryStructure_Hairpin_*.tmp >> $folder/tmp/SecondaryStructure_hairpin.txt`;

#read reverse oligos & retrieve only forward oligos
`cat $folder/tmp/SecondaryStructure_Hairpin_rev*.tmp.delta > $folder/tmp/Hairpin_delta_rev.tmp`;
open(IN, "<$folder/tmp/Hairpin_delta_rev.tmp") or die;
while(defined(my $input = <IN>)) { #only for reverse
    chomp($input);
    my ($primer_R, $dG_hair) = split(/\t/, $input); #what is ) or null?!?
    if ($dG_hair < $all{'OLI_HAIRPIN'}) {
        #retrieve only reverse primers that have valid hairpin values for correspondent forward oligo
        foreach my $primer_F (keys %{$rev_to_for{$primer_R}}) {
            $for_to_count{$primer_F}--;
            $rev_to_count{$primer_R} = 0;
        }
    }
}
close(IN);
undef %oligo_retrieve;
foreach my $primer (keys %for_to_count) {
    if ($for_to_count{$primer} > 0) {
        $oligo_retrieve{$primer} = '';
    }
}
`rm $folder/tmp/SecondaryStructure_Hairpin*`;

my $revCheck4 = keys %oligo_retrieve;
my $dataCheck4 = localtime();
print "how many4 $revCheck4 and $dataCheck4\n";

#CHECK 11 - self dimers for forward oligos
my $size_for;
if ($all{'NUMBERWILD_MAX'} < 3) {
    $size_for = ceil((scalar keys %oligo_retrieve) / 7); #number of forward oligos
} elsif ($all{'NUMBERWILD_MAX'} < 5) {
    $size_for = ceil((scalar keys %oligo_retrieve) / 10); #number of forward oligos
} else {
    $size_for = ceil((scalar keys %oligo_retrieve) / 20); #number of forward oligos
}

my $self = 0;
my $in = 0;

open(my $tmp, ">$folder/tmp/Self_for_1.tmp") or die;
my $fileName = "Self_for_1.tmp";

foreach my $primer (keys %oligo_retrieve) { #only for forward
    $self++;
    if (($self % $size_for) == 0) {
            print $tmp "$primer\n";
            close($tmp);
            system("perl $path_cgi/Self_dimer_checking_web_pp.pl -folder $folder -primer $fileName -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius &");
            $fileName = "Self_for_" . $self . ".tmp";
            open($tmp, ">$folder/tmp/Self_for_${self}.tmp") or die;
            $in = 0;
        } else {
            print $tmp "$primer\n";
            $in = 1;
        }
}
close(IN);
if ($in == 1) {
    system("perl $path_cgi/Self_dimer_checking_web_pp.pl -folder $folder -primer $fileName -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius &");
} else {
    #remove last empty file
    `rm $folder/tmp/$fileName`;
}

#stall till has not finished
my $done = 0;
while ($done == 0) {
    sleep 30;
    my $list = `ls $folder/tmp/Self_for_*tmp 2> /dev/null`;
    if ($list eq '') {
        $done = 1;
    }
}
`cat $folder/tmp/SecondaryStructure_Self_*.tmp > $folder/tmp/SecondaryStructure_selfDimer.txt`;

#read forward oligos & retrieve only reverse oligos
`cat $folder/tmp/SecondaryStructure_Self_for*.tmp.delta > $folder/tmp/Self_delta_for.tmp`;
open(IN, "<$folder/tmp/Self_delta_for.tmp") or die;
while(defined(my $input = <IN>)) { #only for forward
    chomp($input);
    my ($primer_F, $dG_self) = split(/\t/, $input); #what is ) or null?!?
    if ($dG_self < $all{'OLI_SELF'}) {
        #retrieve only reverse primers that have valid hairpin values for correspondent forward oligo
        foreach my $primer_R (keys %{$for_to_rev{$primer_F}}) {
            $for_to_count{$primer_F} = 0;
            $rev_to_count{$primer_R}--;
        }
    }
}
close(IN);
undef %oligo_retrieve;
foreach my $primer (keys %rev_to_count) {
    if ($rev_to_count{$primer} > 0) {
        $oligo_retrieve{$primer} = '';
    }
}
`rm $folder/tmp/SecondaryStructure_Self*`;

my $revCheck5 = keys %oligo_retrieve;
my $dataCheck5 = localtime();
print "how many5 $revCheck5 and $dataCheck5\n";

#CHECK 12 - self dimers for reverse oligos
my $size_rev;
if ($all{'NUMBERWILD_MAX'} < 3) {
    $size_rev = ceil((scalar keys %oligo_retrieve) / 7); #number of reverse oligos
} elsif ($all{'NUMBERWILD_MAX'} < 5) {
    $size_rev = ceil((scalar keys %oligo_retrieve) / 10); #number of reverse oligos
} else {
    $size_rev = ceil((scalar keys %oligo_retrieve) / 20); #number of reverse oligos
}

my $self = 0;
my $in = 0;

open(my $tmp, ">$folder/tmp/Self_rev_1.tmp") or die;
my $fileName = "Self_rev_1.tmp";

foreach my $primer (keys %oligo_retrieve) { #only for reverse
    $self++;
    if (($self % $size_rev) == 0) {
            print $tmp "$primer\n";
            close($tmp);
            system("perl $path_cgi/Self_dimer_checking_web_pp.pl -folder $folder -primer $fileName -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius &");
            $fileName = "Self_rev_" . $self . ".tmp";
            open($tmp, ">$folder/tmp/Self_rev_${self}.tmp") or die;
            $in = 0;
        } else {
            print $tmp "$primer\n";
            $in = 1;
        }
}
close(IN);
if ($in == 1) {
    system("perl $path_cgi/Self_dimer_checking_web_pp.pl -folder $folder -primer $fileName -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius &");
} else {
    #remove last empty file
    `rm $folder/tmp/$fileName`;
}

#stall till has not finished
my $done = 0;
while ($done == 0) {
    sleep 30;
    my $list = `ls $folder/tmp/Self_rev_*tmp 2> /dev/null`;
    if ($list eq '') {
        $done = 1;
    }
}
`cat $folder/tmp/SecondaryStructure_Self_*.tmp >> $folder/tmp/SecondaryStructure_selfDimer.txt`;

#read reverse oligos & retrieve only forward oligos
`cat $folder/tmp/SecondaryStructure_Self_rev*.tmp.delta > $folder/tmp/Self_delta_rev.tmp`;
open(IN, "<$folder/tmp/Self_delta_rev.tmp") or die;
while(defined(my $input = <IN>)) { #only for reverse
    chomp($input);
    my ($primer_R, $dG_hair) = split(/\t/, $input); #what is ) or null?!?
    if ($dG_hair < $all{'OLI_SELF'}) {
        #retrieve only reverse primers that have valid hairpin values for correspondent forward oligo
        foreach my $primer_F (keys %{$rev_to_for{$primer_R}}) {
            $for_to_count{$primer_F}--;
            $rev_to_count{$primer_R} = 0;
        }
    }
}
close(IN);
undef %oligo_retrieve;
foreach my $primer (keys %for_to_count) { #not necessary
    if ($for_to_count{$primer} > 0) {
        $oligo_retrieve{$primer} = '';
    }
}
`rm $folder/tmp/SecondaryStructure_Self*`;

my $revCheck6 = keys %oligo_retrieve;
my $dataCheck6 = localtime();
print "how many6 $revCheck6 and $dataCheck6\n";

# #CHECK 8 - self dimers
# my %forward_pass;
# my %reverse_pass;
# my %allPrimer;
# my $presence3 = 0;
#
# #forward oligos
# `cat $folder/tmp/SecondaryStructure_Self_for*.tmp.delta > $folder/tmp/Self_delta_for.tmp`;
#
# open(IN, "<$folder/tmp/Self_delta_for.tmp") or die;
#
# while(defined(my $input = <IN>)) { #only for forward
#     chomp($input);
#     my ($oligo, $dG_self) = split(/\t/, $input); #what is ) or null?!?
#     if (($dG_self > $all{'OLI_SELF'}) or ($dG_self eq ''))  {
#         $forward{$oligo}{'SELF'} = $dG_self;
#         $presence3++;
#         $forward_pass{$forward{$oligo}{'POS'}}{$oligo} = ""; #start forward
#         $allPrimer{$oligo} = '';
#     }
# }
# close(IN);
#
# #reverse oligos
# `cat $folder/tmp/SecondaryStructure_Self_rev*.tmp.delta > $folder/tmp/Self_delta_rev.tmp`;
#
# open(IN, "<$folder/tmp/Self_delta_rev.tmp") or die;
#
# while(defined(my $input = <IN>)) { #only for forward
#     chomp($input);
#     my ($oligo, $dG_self) = split(/\t/, $input); #what is ) or null?!?
#     if (($dG_self > $all{'OLI_SELF'}) or ($dG_self eq ''))  {
#         $reverse{$oligo}{'SELF'} = $dG_self;
#         $presence3++;
#         $reverse_pass{$reverse{$oligo}{'POS'}}{$oligo} = ""; #start reverse
#         $allPrimer{$oligo} = '';
#     }
# }
# close(IN);
#


#run
my %allPrimer;
my %allPrimerRemove;
my $presence3 = 0;
foreach my $primer (keys %for_to_count) {
    if ($for_to_count{$primer} > 0) {
        $allPrimer{$primer} = '';
        $presence3++;
    } else {
        $allPrimerRemove{$primer} = '';
    }
}
foreach my $primer (keys %rev_to_count) {
    if ($rev_to_count{$primer} > 0) {
        $allPrimer{$primer} = '';
        $presence3++;
    } else {
        $allPrimerRemove{$primer} = '';
    }
}

$diffPres = $presence3/$presence2;

if ($diffPres < 0.5) {
    push @diffMessage, '<li>dG value for hairpin and self dimers</li>';
}

#if no suitable oligos were found during the first screening
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
my $size_rev;
if ($all{'NUMBERWILD_MAX'} < 3) {
    $size = ceil($presence3/7);
} elsif ($all{'NUMBERWILD_MAX'} < 5) {
    $size = ceil($presence3/10);
} else {
    $size = ceil($presence3/10);
}

my $cross = 0;
my $file_num = 0;
my $combo; #array with all oligo variations
open(my $tmp, ">$folder/tmp/CrossSelf_1.tmp") or die;

foreach my $primer (keys %allPrimer) {
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

$done = 0;
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

my $pair =0;
my $Ta = 0;

my %pairPass;
my %pairPassCross;
my %pairPassTa;

my %printCross;

my %pairPass;
my $pair=0;

my %forwardPrint;
my %reversePrint;

my $tm_diff;
$presence3 = 0;
$presence4 = 0;
foreach my $pos_F (sort {$a<=>$b} keys %forward_pass) {
    my $num;
    foreach my $pos_R (sort {$a<=>$b} keys %reverse_pass) {
        $presence3++;
        my $diff = $pos_R - $pos_F;

        #CHECK 9 - amplicon size
        if (($diff > $all{'AMPL_SIZE_MIN'}) && ($diff < $all{'AMPL_SIZE_MAX'})) {
            foreach my $primer_F (sort {$a<=>$b} keys %{$forward_pass{$pos_F}}) {
                if (!(defined($allPrimerRemove{$primer_F}))) {
                    foreach my $primer_R (sort {$a<=>$b} keys %{$reverse_pass{$pos_R}}) {
                        if (!(defined($allPrimerRemove{$primer_R}))) {
                            if ($reverse{$primer_R}{'TM'} >= $forward{$primer_F}{'TM'}) {
                                $Ta = $forward{$primer_F}{'TM'} - 5;
                                $tm_diff = $reverse{$primer_R}{'TM'} - $forward{$primer_F}{'TM'};
                            } else {
                                $Ta = $reverse{$primer_R}{'TM'} - 5;
                                $tm_diff = $reverse{$primer_R}{'TM'} - $forward{$primer_F}{'TM'};
                            }
                            #CHECK 10 - annealing temperature and Tm difference between primers
                            if (($Ta >= $all{'TA_MIN'}) && ($Ta <= $all{'TA_MAX'}) && ($tm_diff <= $all{'TM_DIFF'})) {
                                $presence4++;
                                $pair++;
                                $pairPass{$pair}{'F'} = $primer_F; #pairnumber - F = primer ####new
                                $pairPass{$pair}{'R'} = $primer_R; #pairnumber - R = primer ####new
                                $forwardPrint{$primer_F} = '';
                                $reversePrint{$primer_R} = '';
                            }
                        }
                    }
                }
            }
        }
    }
}

if ($presence3 == 0) {
    $diffPres = 0;
} else {
    $diffPres = $presence4/$presence3;
}

if ($diffPres < 0.5) {
    push @diffMessage, '<li>amplicon size</li>';
}

if (($presence3 == 0) or ($presence4 == 0)) {
    `mv ${folder}/${nameFile}* ${folder}/inputs/`;
    #send email
    my $to = $all{'EMAIL'};
    my $from = 'gv16363@bristol.ac.uk';
    my $subject = 'PhyloPrimer results - ' . $all{'PROJECT'};
    my $message = "Content-Type: text/html; charset=ISO-8859-1\n\n<html><body>Hi,<br>PhyloPrimer did not find any suitable primer pairs with the selected parameters. It looks like the parameters for the following fields were too stringent:<br><ul>@diffMessage</ul><br><br>All the best,<br>Gilda</body>";

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
        if (defined($forwardPrint{$oligo}) or defined($reversePrint{$oligo})) {
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
my $oligo;
my $in = 0;
while(defined(my $input = <IN>)) {
    chomp($input);
    if ($input =~ /^>/) {
        ($disc, $oligo) = split(/\t/, $input);
        if (defined($forwardPrint{$oligo}) or defined($reversePrint{$oligo})) {
            print $file "$input\n";
            $in = 1;
        } else {
            $in = 0;
        }
    } else {
        if ($in == 1) {
            print $file "$input\n";
        }
    }
}
close(IN);
close($file);

#print forward file
open(my $file, ">$folder/results/forwardList.txt") or die;
foreach my $for (keys %forwardPrint) {
    print $file "FORWARD\t$for\n";
    print $file "POSITION\t$forward{$for}{'POS_NAME'}\n";
    print $file "LENGTH\t$forward{$for}{'LEN'} bases\n";
    print $file "GC\t$forward{$for}{'GC'} %\n";
    print $file "TM\t$forward{$for}{'TM'} °C\n";
    if ($forward{$for}{'SELF'} == 0) {
        print $file "SELF DIMER\t>=0 kcal/mol\n";
    } else {
        print $file "SELF DIMER\t$forward{$for}{'SELF'} kcal/mol\n";
    }
    if ($forward{$for}{'HAIR'} == 0) {
        print $file "HAIRPIN\t>=0 kcal/mol\n";
    } else {
        print $file "HAIRPIN\t$forward{$for}{'HAIR'} kcal/mol\n";
    }
    print $file "\\\\\n";
}
close($file);

#print reverse file
open(my $file, ">$folder/results/reverseList.txt") or die;
foreach my $rev (keys %reversePrint) {
    print $file "REVERSE\t$rev\n";
    print $file "POSITION\t$reverse{$rev}{'POS_NAME'}\n";
    print $file "LENGTH\t$reverse{$rev}{'LEN'} bases\n";
    print $file "GC\t$reverse{$rev}{'GC'} %\n";
    print $file "TM\t$reverse{$rev}{'TM'} °C\n";
    if ($reverse{$rev}{'SELF'} == 0) {
        print $file "SELF DIMER\t>=0 kcal/mol\n";
    } else {
        print $file "SELF DIMER\t$reverse{$rev}{'SELF'} kcal/mol\n";
    }
    if ($reverse{$rev}{'HAIR'} == 0) {
        print $file "HAIRPIN\t>=0 kcal/mol\n";
    } else {
        print $file "HAIRPIN\t$reverse{$rev}{'HAIR'} kcal/mol\n";
    }
    print $file "\\\\\n";
}
close($file);

my %scorePair;
my %scorePairEx;

##assign a score to the primer pairs
##taking in consideration:
#have a differed base inside
#have a differed base in the tail
#tm difference betweeen F anf R of less than 1
#haipin dG is higher than -1
#self dimer dG is higher than -1
#cross dimer dG is higher than -1
foreach my $pair (keys %pairPass) {

    $scorePair{$pair} = 1; #so if not any selection, all the pairs will equal to 1 and will be randomly selected

    if (($all{'MAXIMIZE_SEL1'} eq 'yes') or ($all{'MAXIMIZE_SEL2'} eq 'yes')) {
        #start and end of primers
        my $startF = $forward{$pairPass{$pair}{'F'}}{'POS'};
        my $startR = $reverse{$pairPass{$pair}{'R'}}{'POS'};
        my $endF = $forward{$pairPass{$pair}{'F'}}{'POS'} + $forward{$pairPass{$pair}{'F'}}{'LEN'} - 1;
        my $endR = $reverse{$pairPass{$pair}{'R'}}{'POS'} + $reverse{$pairPass{$pair}{'R'}}{'LEN'} - 1;

        #start and end of primer tails
        my $startF_tail = $endF - 4; #last 5 bases at 3' end
        my $endR_tail = $startR + 4; #last 5 bases at 3' end

        my $endF1 = $endF - 1;
        my $startR1 = $startR + 1;

        foreach my $d (keys %difference) {

            if ($all{'MAXIMIZE_SEL2'} eq 'yes') {
                if (($d >= $startF) && ($d <= $endF)) { #if different bases in F primer
                    $scorePair{$pair}+=2;
                }
                if (($d >= $startR) && ($d <= $endR)) { #if different bases in R primer
                    $scorePair{$pair}+=2;
                }

            }
            if ($all{'MAXIMIZE_SEL1'} eq 'yes') {

                if (($d >= $startF_tail) && ($d <= $endF)) { #if different bases in the F primer tail
                    $scorePair{$pair} += 2;
                }
                if (($d >= $startR) && ($d <= $endR_tail)) { #if different bases in the R primer tail
                    $scorePair{$pair} += 2;
                }

                #first
                if ($d == $endF) { #if different bases in the F primer tail
                    $scorePair{$pair} += 20;
                }
                if ($d == $startR) { #if different bases in the R primer tail
                    $scorePair{$pair} += 20;
                }

                #second
                if ($d == $endF1) { #if different bases in the F primer tail
                    $scorePair{$pair} += 10;
                }
                if ($d == $startR1) { #if different bases in the R primer tail
                    $scorePair{$pair} += 10;
                }

            }
        }
    }

    if ($all{'TM_SEL'} eq 'yes') {
        my $tm_diff = $forward{$pairPass{$pair}{'F'}}{'TM'} - $forward{$pairPass{$pair}{'R'}}{'TM'};
        if (($tm_diff >= -1) && ($tm_diff <= 1)) {
            $scorePair{$pair}++;
        }
    }

    if ($all{'DG_SEL'} eq 'yes') {
        if (($forward{$pairPass{$pair}{'F'}}{'HAIR'} >= -1) or ($forward{$pairPass{$pair}{'F'}}{'HAIR'} eq '')) {
            $scorePair{$pair}++;
        }
        if (($reverse{$pairPass{$pair}{'R'}}{'HAIR'} >= -1) or ($reverse{$pairPass{$pair}{'R'}}{'HAIR'} eq '')) {
            $scorePair{$pair}++;
        }
        if (($forward{$pairPass{$pair}{'F'}}{'SELF'} >= -1) or ($forward{$pairPass{$pair}{'F'}}{'SELF'} eq '')) {
            $scorePair{$pair}++;
        }
        if (($reverse{$pairPass{$pair}{'R'}}{'SELF'} >= -1) or ($reverse{$pairPass{$pair}{'R'}}{'SELF'} eq '')) {
            $scorePair{$pair}++;
        }
    }

    if ($all{'DEG_SEL'} eq 'yes') {
        my $deg2 = () = $pairPass{$pair}{'F'} =~ /R|Y|S|W|K|M/g;
        my $deg3 = () = $pairPass{$pair}{'F'} =~ /B|D|H|V/g;
        my $deg4 = () = $pairPass{$pair}{'F'} =~ /N/g;
        $scorePair{$pair} -= $deg2*2;
        $scorePair{$pair} -= $deg3*3;
        $scorePair{$pair} -= $deg4*4;
        my $deg2 = () = $pairPass{$pair}{'R'} =~ /R|Y|S|W|K|M/g;
        my $deg3 = () = $pairPass{$pair}{'R'} =~ /B|D|H|V/g;
        my $deg4 = () = $pairPass{$pair}{'R'} =~ /N/g;
        $scorePair{$pair} -= $deg2*2;
        $scorePair{$pair} -= $deg3*3;
        $scorePair{$pair} -= $deg4*4;
    }
}

#get the range of score
my %scorePairDef;
my $maxScore = 0;
my $minScore = 1;
foreach my $pair (keys %scorePair) {
    $scorePairDef{$scorePair{$pair}}{$pair} = ''; #score - pair
    if ($scorePair{$pair} > $maxScore) {
        $maxScore = $scorePair{$pair};
    }
    if ($scorePair{$pair} < $minScore) {
        $minScore = $scorePair{$pair};
    }
}

my %top;
my %topPair;

my %doneF;
my %doneR;
my %cluster;


my $oligoNum = 500; #500 only if differences between positive and negative consensus. If no differences, the design of the oligos will be spread on all the consensus and the oligos will be all different between each other and the mySQL search takes too long.

if (($all{'DIFFERENT_POS'} eq "no") or (($all{'MAXIMIZE_SEL1'} ne 'yes') && ($all{'MAXIMIZE_SEL2'} ne 'yes'))) {
    $oligoNum = 250;
}

###Cross_dimer calculated on 1000 pairs.
###even if not 1000 calculate on all
open(my $tmp1, ">$folder/tmp/Cross_1.tmp") or die;
open(my $tmp2, ">$folder/tmp/Cross_2.tmp") or die;
open(my $tmp3, ">$folder/tmp/Cross_3.tmp") or die;
open(my $tmp4, ">$folder/tmp/Cross_4.tmp") or die;
open(my $tmp5, ">$folder/tmp/Cross_5.tmp") or die;

my $index = 0;
my @pairIndex;
foreach (my $score=$maxScore; $score>=$minScore; $score--) {
    foreach my $pair (keys %{$scorePairDef{$score}}) {
        $index++;
        push @pairIndex, $pair;
        if ($index <= 1000) { #check 1000 primer pairs for cross dimers
            my $combo = $pairPass{$pair}{'F'} . "\t" . $pairPass{$pair}{'R'};
            if ($index > 800) {
                print $tmp1 "$combo\n";
            } elsif ($index > 600) {
                print $tmp2 "$combo\n";
            } elsif ($index > 400) {
                print $tmp3 "$combo\n";
            } elsif ($index > 200) {
                print $tmp4 "$combo\n";
            } else {
                print $tmp5 "$combo\n";
            }
        }
    }
}
close($tmp1);
close($tmp2);
close($tmp3);
close($tmp4);
close($tmp5);

system("perl $path_cgi/Cross_dimer_checking_web_pp.pl -folder $folder -primer Cross_1.tmp -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius &"); #Tm is the same for forward and reverse.
system("perl $path_cgi/Cross_dimer_checking_web_pp.pl -folder $folder -primer Cross_2.tmp -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius &"); #Tm is the same for forward and reverse.
system("perl $path_cgi/Cross_dimer_checking_web_pp.pl -folder $folder -primer Cross_3.tmp -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius &"); #Tm is the same for forward and reverse.
system("perl $path_cgi/Cross_dimer_checking_web_pp.pl -folder $folder -primer Cross_4.tmp -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius &"); #Tm is the same for forward and reverse.
system("perl $path_cgi/Cross_dimer_checking_web_pp.pl -folder $folder -primer Cross_5.tmp -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius &"); #Tm is the same for forward and reverse.

my $done = 0;
while ($done == 0) {
    sleep 10;
    my $list = `ls $folder/tmp/Cross_*tmp 2> /dev/null`;
    if ($list eq '') {
        $done = 1;
    }
}

`cat $folder/tmp/SecondaryStructure_crossDimer_*.tmp > $folder/crossDimer.txt`;

my $disc;
my $for;
my $rev;
my $dG_cross;
my $presence5 = 0;

my %forwardPrint;
my %reversePrint;

my $pair;
my $index = 0;
my %pairPassCross;
my %pairPassTa;

open(IN, "<$folder/crossDimer.txt") or die;
while(defined(my $input = <IN>)) {
    chomp($input);
    if ($input =~ /^>/) {
        ($disc, $for, $rev) = split(/\t/, $input);
        $pair = $pairIndex[$index];
        $index++;
    } elsif ($input =~ /^@/) {
        ($disc, $dG_cross) = split(/\t/, $input);
        $dG_cross =~ s/ kcal\/mol//g;

        #CHECK 11 - cross dimers
        if ($dG_cross > $all{'OLI_CROSS'}) {
            $presence5++; ####check that OK

            $pairPassCross{$pair} = $dG_cross;

            #re-calculate Ta to be saved
            if ($reverse{$rev}{'TM'} >= $forward{$for}{'TM'}) {
                $Ta = $forward{$for}{'TM'} - 5;
            } else {
                $Ta = $reverse{$rev}{'TM'} - 5;
            }

            $pairPassTa{$pair} = $Ta;

        }
    }
}
close(IN);

if ($presence5 == 0) {
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

#create file with all the primers for the blast search
open(my $tmp, ">$folder/tmp/oligo.txt") or die;
my $index = 0;
foreach (my $score=$maxScore; $score>=$minScore; $score--) {
    foreach my $pair (keys %{$scorePairDef{$score}}) {
        if (defined($pairPassCross{$pair})) { ###only the oligo passed cross dimer check
            $index++;
            if ($index <= $oligoNum) {
                print $tmp "$pair\t$pairPass{$pair}{'F'}\t$forward{$pairPass{$pair}{'F'}}{'LEN'}\t$forward{$pairPass{$pair}{'F'}}{'POS'}\tfor\n";
                print $tmp "$pair\t$pairPass{$pair}{'R'}\t$reverse{$pairPass{$pair}{'R'}}{'LEN'}\t$reverse{$pairPass{$pair}{'R'}}{'POS'}\trev\n";
            }
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

my %scorePairNew;
my $maxScoreNew = 1;
my $minScoreNew = 1;
my %foundSp;

my $emptyBLAST = 0;
my $tableNt = ">";
my $pieChartTableFR;

if (-z $checkFile) { #if file is empty
    $emptyBLAST = 1;
    $tableNt = "none";
    $pieChartTableFR = "none";

} else { #if file is not empty

    #retrieve information from BLAST file
    checkBLAST('inSilico_nt.txt');

    #connect to mysql database
    my $dbh;
    my $sth;

    $dbh = DBI->connect ($dsn, $user_name, $password, { RaiseError => 1 });

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

    foreach (my $score=$maxScore; $score>=$minScore; $score--) {
        foreach my $pair (keys %{$scorePairDef{$score}}) {

            if (defined($pairPass{$pair})) {
                my $scoreNew = 0;
                my $for = $pairPass{$pair}{'F'};
                my $rev = $pairPass{$pair}{'R'};

                if ($all{'SPECIES_SEL'} eq 'yes') {
                    my $selectedRank = 'SPECIES';

                    foreach my $mis (@mismatch) { #print first alignments with less mismatches
                        if (defined($tableBlast{$pair}{$mis})) {
                            foreach my $acc (keys %{$tableBlast{$pair}{$mis}}) {
                                if (defined($selectedTaxa{'7'}{$taxonomy{$acc}{$selectedRank}})) {
                                    if (defined($foundSp{$pair}{$taxonomy{$acc}{$selectedRank}})) {
                                        $scoreNew += 5;
                                    } else {
                                        $scoreNew += 20; #incentivate variety
                                        $foundSp{$pair}{$taxonomy{$acc}{$selectedRank}} = '';
                                    }
                                } else {
                                    $scoreNew -= 40;
                                    $foundSp{$pair}{$taxonomy{$acc}{$selectedRank}} = '';
                                }
                            }
                        }
                    }
                } else {
                    foreach my $mis (@mismatch) { #print first alignments with less mismatches
                        if (defined($tableBlast{$pair}{$mis})) {
                            foreach my $acc (keys %{$tableBlast{$pair}{$mis}}) {
                                $foundSp{$pair}{$taxonomy{$acc}{'SPECIES'}} = '';
                            }
                        }
                    }
                }

                if ($all{'GENUS_SEL'} eq 'yes') {
                    my $selectedRank = 'GENUS';
                    foreach my $mis (@mismatch) { #print first alignments with less mismatches
                        if (defined($tableBlast{$pair}{$mis})) {
                            foreach my $acc (keys %{$tableBlast{$pair}{$mis}}) {
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
                        if (defined($tableBlast{$pair}{$mis})) {
                            foreach my $acc (keys %{$tableBlast{$pair}{$mis}}) {
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
                        if (defined($tableBlast{$pair}{$mis})) {
                            foreach my $acc (keys %{$tableBlast{$pair}{$mis}}) {
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
                        if (defined($tableBlast{$pair}{$mis})) {
                            foreach my $acc (keys %{$tableBlast{$pair}{$mis}}) {
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
                        if (defined($tableBlast{$pair}{$mis})) {
                            foreach my $acc (keys %{$tableBlast{$pair}{$mis}}) {
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
                        if (defined($tableBlast{$pair}{$mis})) {
                            foreach my $acc (keys %{$tableBlast{$pair}{$mis}}) {
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
                $scorePairNew{$scoreNew}{$pair} = '';
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

    #create table
    my %pieFR;
    my $index = 0;
    foreach (my $score=$maxScoreNew; $score>=$minScoreNew; $score--) {
        foreach my $pair (keys %{$scorePairNew{$score}}) {
            $index++;
            if ($index <= 100) {

                #print on new oligo.txt
                print $tmp "$pair\t$pairPass{$pair}{'F'}\t$forward{$pairPass{$pair}{'F'}}{'LEN'}\t$forward{$pairPass{$pair}{'F'}}{'POS'}\tfor\n";
                print $tmp "$pair\t$pairPass{$pair}{'R'}\t$reverse{$pairPass{$pair}{'R'}}{'LEN'}\t$reverse{$pairPass{$pair}{'R'}}{'POS'}\trev\n";

                #create BLAST table & pieCharts
                my $for = $pairPass{$pair}{'F'};
                my $rev = $pairPass{$pair}{'R'};
                my $combined = $pairPass{$pair}{'F'} . "-" . $pairPass{$pair}{'R'};
                my $tableFR;

                foreach my $mis (@mismatch) { #print first alignments with less mismatches
                    if (defined($tableBlast{$pair}{$mis})) {
                        foreach my $acc (keys %{$tableBlast{$pair}{$mis}}) { ###I need to order them somehow

                            #data for pieChart
                            $pieFR{'DOMAIN'}{$taxonomy{$acc}{'DOMAIN'}}++;
                            $pieFR{'PHYLUM'}{$taxonomy{$acc}{'PHYLUM'}}++;
                            $pieFR{'CLASS'}{$taxonomy{$acc}{'CLASS'}}++;
                            $pieFR{'ORDER'}{$taxonomy{$acc}{'ORDER'}}++;
                            $pieFR{'FAMILY'}{$taxonomy{$acc}{'FAMILY'}}++;
                            $pieFR{'GENUS'}{$taxonomy{$acc}{'GENUS'}}++;
                            $pieFR{'SPECIES'}{$taxonomy{$acc}{'SPECIES'}}++;

                            #accessions present in both forward and reverse primers
                            $tableFR .= "<tr><td rowspan='2'>" . $acc . "</td>";
                            $tableFR .= "<td>F</td>";
                            $tableFR .= "<td>" . $for . "<br>" . $tableBlast{$pair}{$mis}{$acc}{'AL_F'} . "</td>";
                            $tableFR .= "<td>" . $tableBlast{$pair}{$mis}{$acc}{'START_F'} . "</td>";
                            $tableFR .= "<td>" . $tableBlast{$pair}{$mis}{$acc}{'END_F'} . "</td>";
                            $tableFR .= "<td rowspan='2'>" . $tableBlast{$pair}{$mis}{$acc}{'LEN'} . "</td>";

                            #taxonomy
                            if ($taxonomy{$acc}{'DOMAIN'} eq "") {
                                $tableFR .= "<td rowspan='2' class ='D'>Unclassified</td>";
                            } else {
                                $tableFR .= "<td rowspan='2' class ='D'>" . $taxonomy{$mis}{'DOMAIN'} . "</td>";
                            }
                            if ($taxonomy{$acc}{'PHYLUM'} eq "") {
                                $tableFR .= "<td rowspan='2' class ='P'>Unclassified</td>";
                            } else {
                                $tableFR .= "<td rowspan='2' class ='P'>" . $taxonomy{$acc}{'PHYLUM'} . "</td>";
                            }
                            if ($taxonomy{$acc}{'CLASS'} eq "") {
                                $tableFR .= "<td rowspan='2' class ='C'>Unclassified</td>";
                            } else {
                                $tableFR .= "<td rowspan='2' class ='C'>" . $taxonomy{$acc}{'CLASS'} . "</td>";
                            }
                            if ($taxonomy{$acc}{'ORDER'} eq "") {
                                $tableFR .= "<td rowspan='2' class ='O'>Unclassified</td>";
                            } else {
                                $tableFR .= "<td rowspan='2' class ='O'>" . $taxonomy{$acc}{'ORDER'} . "</td>";
                            }
                            if ($taxonomy{$acc}{'FAMILY'} eq "") {
                                $tableFR .= "<td rowspan='2' class ='F'>Unclassified</td>";
                            } else {
                                $tableFR .= "<td rowspan='2' class ='F'>" . $taxonomy{$acc}{'FAMILY'} . "</td>";
                            }
                            if ($taxonomy{$acc}{'GENUS'} eq "") {
                                $tableFR .= "<td rowspan='2' class ='G'>Unclassified</td>";
                            } else {
                                $tableFR .= "<td rowspan='2' class ='G'>" . $taxonomy{$acc}{'GENUS'} . "</td>";
                            }
                            if ($taxonomy{$acc}{'SPECIES'} eq "") {
                                $tableFR .= "<td rowspan='2' class ='S'>Unclassified</td></tr>";
                            } else {
                                $tableFR .= "<td rowspan='2' class ='S'>" . $taxonomy{$acc}{'SPECIES'} . "</td></tr>";
                            }

                            $tableFR .= "<tr><td>R</td>";
                            $tableFR .= "<td>" . $rev . "<br>" . $tableBlast{$pair}{$mis}{$acc}{'AL_R'} . "</td>";
                            $tableFR .= "<td>" . $tableBlast{$pair}{$mis}{$acc}{'START_R'} . "</td>";
                            $tableFR .= "<td>" . $tableBlast{$pair}{$mis}{$acc}{'END_R'} . "</td></tr>";
                        }
                    }
                }

                #compose table
                if ($tableFR ne '') {
                    $tableNt .= $combined . "<table>";
                    $tableNt .= $tableFR;
                    $tableNt .= "<tr><td colspan='8' id='dividingCell'></td></tr>";
                    $tableNt =~ s/<tr><td colspan='8' id='dividingCell'><\/td><\/tr>$//g;
                    $tableNt .= "</table>";
                }

                #assemble pieChartTable - FR
                $pieChartTableFR .= "<" . $for . "-" . $rev . ">";

                foreach my $rank ('DOMAIN','PHYLUM','CLASS','ORDER','FAMILY','GENUS','SPECIES') {
                    foreach my $taxon (keys %{$pieFR{$rank}}) {
                        if ($taxon eq "") {
                            $pieChartTableFR .= "Unclassified-" . $pieFR{$rank}{$taxon} . "-Brown,";
                        } else {
                            $pieChartTableFR .= $taxon . "-" . $pieFR{$rank}{$taxon} . "-" . $taxonomyColour{$taxon} . ",";
                        }
                    }
                    chop($pieChartTableFR);
                    $pieChartTableFR .= ":";
                }
                chop($pieChartTableFR);
                $pieChartTableFR .= "<>";

                undef %pieFR;
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
            foreach my $pair (keys %{$scorePairNew{$score}}) {

                #create BLAST table for additional user BLAST search
                my $for = $pairPass{$pair}{'F'};
                my $rev = $pairPass{$pair}{'R'};
                my $combined = $pairPass{$pair}{'F'} . "-" . $pairPass{$pair}{'R'};

                my $tableFR;

                my $countFR = 0;

                foreach my $mis (@mismatch) { #print first alignments with less mismatches
                    if (defined($tableBlast{$pair}{$mis})) {
                        foreach my $acc (keys %{$tableBlast{$pair}{$mis}}) { ###I need to order them somehow
                            $countFR++;
                            #accessions present in both forward and reverse primers
                            $tableFR .= "<tr><td rowspan='2'>" . $acc . "</td>";
                            $tableFR .= "<td>F</td>";
                            $tableFR .= "<td>" . $for . "<br>" . $tableBlast{$pair}{$mis}{$acc}{'AL_F'} . "</td>";
                            $tableFR .= "<td>" . $tableBlast{$pair}{$mis}{$acc}{'START_F'} . "</td>";
                            $tableFR .= "<td>" . $tableBlast{$pair}{$mis}{$acc}{'END_F'} . "</td>";
                            $tableFR .= "<td rowspan='2'>" . $tableBlast{$pair}{$mis}{$acc}{'LEN'} . "</td>";

                            $tableFR .= "<tr><td>R</td>";
                            $tableFR .= "<td>" . $rev . "<br>" . $tableBlast{$pair}{$mis}{$acc}{'AL_R'} . "</td>";
                            $tableFR .= "<td>" . $tableBlast{$pair}{$mis}{$acc}{'START_R'} . "</td>";
                            $tableFR .= "<td>" . $tableBlast{$pair}{$mis}{$acc}{'END_R'} . "</td></tr>";
                        }
                    }
                }

                #compose table
                if ($tableFR ne '') {
                    $tableUserCheck .= $combined . "<" . $countAll . ";" . $countFR . "><table>";
                    $tableUserCheck .= $tableFR;
                    $tableUserCheck .= "<tr><td colspan='8' id='dividingCell'></td></tr>";
                    $tableUserCheck =~ s/<tr><td colspan='8' id='dividingCell'><\/td><\/tr>$//g;
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
open(my $file1, ">$folder/results/assayList.txt") or die; #file with all the necessary parameters - user data

print $file "PROJECT\t$all{'PROJECT'}\n";
print $file "CONSENSUS\t$all{'CONSENSUS'}\n";
print $file "DIFFERENT_POS\t$all{'DIFFERENT_POS'}\n";

#print selection criteria
my @vis_sel = ('TM_SEL', 'DG_SEL', 'DEG_SEL', 'MAXIMIZE_SEL1', 'MAXIMIZE_SEL2', 'SPECIES_SEL', 'GENUS_SEL', 'FAMILY_SEL', 'ORDER_SEL', 'CLASS_SEL', 'PHYLUM_SEL', 'DOMAIN_SEL');
foreach my $v (@vis_sel) {
    if ($all{$v} eq "yes") {
        print $file "$v\t$all{$v}\n";
    }
}

print $file "\n";

my %bestPair;
my $index = 0;
foreach (my $score=$maxScoreNew; $score>=$minScoreNew; $score--) {
    foreach my $pair (keys %{$scorePairNew{$score}}) {
        $index++;
        $bestPair{$pair} = '';
        if ($index <= 100) {

            print $file "INDEX\t$index\n";
            print $file "SCORE\t$score\n";
            print $file "FORWARD\t$pairPass{$pair}{'F'}\n";
            print $file "POSITION\t$forward{$pairPass{$pair}{'F'}}{'POS_NAME'}\n";
            print $file "LENGTH\t$forward{$pairPass{$pair}{'F'}}{'LEN'}\n";
            print $file "GC\t$forward{$pairPass{$pair}{'F'}}{'GC'}\n";
            print $file "TM\t$forward{$pairPass{$pair}{'F'}}{'TM'}\n";
            if ($forward{$pairPass{$pair}{'F'}}{'SELF'} == 0) {
                print $file "SELF\t>=0\n";
            } else {
                print $file "SELF\t$forward{$pairPass{$pair}{'F'}}{'SELF'}\n";
            }
            if ($forward{$pairPass{$pair}{'F'}}{'HAIR'} == 0) {
                print $file "HAIR\t>=0\n";
            } else {
                print $file "HAIR\t$forward{$pairPass{$pair}{'F'}}{'HAIR'}\n";
            }
            print $file "REVERSE\t$pairPass{$pair}{'R'}\n";
            print $file "POSITION\t$reverse{$pairPass{$pair}{'R'}}{'POS_NAME'}\n";
            print $file "LENGTH\t$reverse{$pairPass{$pair}{'R'}}{'LEN'}\n";
            print $file "GC\t$reverse{$pairPass{$pair}{'R'}}{'GC'}\n";
            print $file "TM\t$reverse{$pairPass{$pair}{'R'}}{'TM'}\n";
            if ($reverse{$pairPass{$pair}{'R'}}{'SELF'} == 0) {
                print $file "SELF\t>=0\n";
            } else {
                print $file "SELF\t$reverse{$pairPass{$pair}{'R'}}{'SELF'}\n";
            }
            if ($reverse{$pairPass{$pair}{'R'}}{'HAIR'} == 0) {
                print $file "HAIR\t>=0\n";
            } else {
                print $file "HAIR\t$reverse{$pairPass{$pair}{'R'}}{'HAIR'}\n";
            }
            print $file "TA\t$pairPassTa{$pair}\n";
            if (defined($pairPassCross{$pair})) {
                if ($pairPassCross{$pair} == 0) {
                    print $file "CROSS\t>=0\n";
                } else {
                    print $file "CROSS\t$pairPassCross{$pair}\n";
                }
            } else {
                print $file "CROSS\tnot calculated\n";
            }
            my $allSp;
            foreach my $sp (keys %{$foundSp{$pair}}) {
                $allSp .= $sp . ";";
            }
            print $file "SPECIES\t$allSp\n\\\\\n";
        }

        print $file1 "INDEX\t$index\n";
        print $file1 "SCORE\t$score\n";
        print $file1 "FORWARD\t$pairPass{$pair}{'F'}\n";
        print $file1 "POSITION\t$forward{$pairPass{$pair}{'F'}}{'POS_NAME'}\n";
        print $file1 "LENGTH\t$forward{$pairPass{$pair}{'F'}}{'LEN'} bases\n";
        print $file1 "GC\t$forward{$pairPass{$pair}{'F'}}{'GC'} %\n";
        print $file1 "TM\t$forward{$pairPass{$pair}{'F'}}{'TM'} °C\n";
        if ($forward{$pairPass{$pair}{'F'}}{'SELF'} == 0) {
            print $file1 "SELF DIMER\t>=0 kcal/mol\n";
        } else {
            print $file1 "SELF DIMER\t$forward{$pairPass{$pair}{'F'}}{'SELF'} kcal/mol\n";
        }
        if ($forward{$pairPass{$pair}{'F'}}{'HAIR'} == 0) {
            print $file1 "HAIRPIN\t>=0 kcal/mol\n";
        } else {
            print $file1 "HAIRPIN\t$forward{$pairPass{$pair}{'F'}}{'HAIR'} kcal/mol\n";
        }
        print $file1 "REVERSE\t$pairPass{$pair}{'R'}\n";
        print $file1 "POSITION\t$reverse{$pairPass{$pair}{'R'}}{'POS_NAME'}\n";
        print $file1 "LENGTH\t$reverse{$pairPass{$pair}{'R'}}{'LEN'} bases\n";
        print $file1 "GC\t$reverse{$pairPass{$pair}{'R'}}{'GC'} %\n";
        print $file1 "TM\t$reverse{$pairPass{$pair}{'R'}}{'TM'} °C\n";
        if ($reverse{$pairPass{$pair}{'R'}}{'SELF'} == 0) {
            print $file1 "SELF DIMER\t>=0\n";
        } else {
            print $file1 "SELF DIMER\t$reverse{$pairPass{$pair}{'R'}}{'SELF'}\n";
        }
        if ($reverse{$pairPass{$pair}{'R'}}{'HAIR'} == 0) {
            print $file1 "HAIRPIN\t>=0 kcal/mol\n";
        } else {
            print $file1 "HAIRPIN\t$reverse{$pairPass{$pair}{'R'}}{'HAIR'} kcal/mol\n";
        }
        print $file1 "TA\t$pairPassTa{$pair} °C\n";
        if (defined($pairPassCross{$pair})) {
            if ($pairPassCross{$pair} == 0) {
                print $file1 "CROSS DIMER\t>=0 kcal/mol\n";
            } else {
                print $file1 "CROSS DIMER\t$pairPassCross{$pair} kcal/mol\n";
            }
        } else {
            print $file1 "CROSS DIMER\tnot calculated\n";
        }
        my $allSp;
        foreach my $sp (keys %{$foundSp{$pair}}) {
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
    foreach my $pair (keys %{$scorePairDef{$score}}) {
        if (!(defined($bestPair{$pair}))) {
            print $file1 "INDEX\t$index\n";
            print $file1 "SCORE\t$score\n";
            print $file1 "FORWARD\t$pairPass{$pair}{'F'}\n";
            print $file1 "POSITION\t$forward{$pairPass{$pair}{'F'}}{'POS_NAME'}\n";
            print $file1 "LENGTH\t$forward{$pairPass{$pair}{'F'}}{'LEN'} bases\n";
            print $file1 "GC\t$forward{$pairPass{$pair}{'F'}}{'GC'} %\n";
            print $file1 "TM\t$forward{$pairPass{$pair}{'F'}}{'TM'} °C\n";
            if ($forward{$pairPass{$pair}{'F'}}{'SELF'} == 0) {
                print $file1 "SELF DIMER\t>=0 kcal/mol\n";
            } else {
                print $file1 "SELF DIMER\t$forward{$pairPass{$pair}{'F'}}{'SELF'} kcal/mol\n";
            }
            if ($forward{$pairPass{$pair}{'F'}}{'HAIR'} == 0) {
                print $file1 "HAIRPIN\t>=0 kcal/mol\n";
            } else {
                print $file1 "HAIRPIN\t$forward{$pairPass{$pair}{'F'}}{'HAIR'} kcal/mol\n";
            }
            print $file1 "REVERSE\t$pairPass{$pair}{'R'}\n";
            print $file1 "POSITION\t$reverse{$pairPass{$pair}{'R'}}{'POS_NAME'}\n";
            print $file1 "LENGTH\t$reverse{$pairPass{$pair}{'R'}}{'LEN'} bases\n";
            print $file1 "GC\t$reverse{$pairPass{$pair}{'R'}}{'GC'} %\n";
            print $file1 "TM\t$reverse{$pairPass{$pair}{'R'}}{'TM'} °C\n";
            if ($reverse{$pairPass{$pair}{'R'}}{'SELF'} == 0) {
                print $file1 "SELF DIMER\t>=0\n";
            } else {
                print $file1 "SELF DIMER\t$reverse{$pairPass{$pair}{'R'}}{'SELF'}\n";
            }
            if ($reverse{$pairPass{$pair}{'R'}}{'HAIR'} == 0) {
                print $file1 "HAIRPIN\t>=0 kcal/mol\n";
            } else {
                print $file1 "HAIRPIN\t$reverse{$pairPass{$pair}{'R'}}{'HAIR'} kcal/mol\n";
            }
            print $file1 "TA\t$pairPassTa{$pair} °C\n";
            if (defined($pairPassCross{$pair})) {
                if ($pairPassCross{$pair} == 0) {
                    print $file1 "CROSS DIMER\t>=0 kcal/mol\n";
                } else {
                    print $file1 "CROSS DIMER\t$pairPassCross{$pair} kcal/mol\n";
                }
            } else {
                print $file1 "CROSS DIMER\tnot calculated\n";
            }
            print $file1 "SPECIES\tnot BLASTed\n";
            print $file1 "\\\\\n";
        }
    }
}

print $file "BLAST_TABLE\t$tableNt\n";
print $file "PIECHART_FR\t$pieChartTableFR\n";
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
`cp $path_html/analysesPhyloprimer/${folder}/crossDimer.txt $path_html/analysesPhyloprimer/${folder}/results/crossDimer.txt`;
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

#send email
my $message = "Hi,\n\nPlease find the link to the PhyloPrimer results: https:\/\/www.cerealsdb.uk.net\/cerealgenomics\/cgi-bin\/phyloprimerResultsPrimer.cgi?defSet=" . $defSet . "\n\n\nAll the best,\nPhyloPrimer team";

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
  if ($seq =~ /[RYSWKMBDHVN]/) { #if degenerate bases
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
        my ($mis, $pair, $acc, $al_f, $al_r, $start_f, $end_f, $start_r, $end_r, $len) = split(/\t/, $input);
        $tableBlastTMP{$pair}{$mis}{$acc}{'LEN'} = $len;
        $tableBlastTMP{$pair}{$mis}{$acc}{'START_F'} = $start_f;
        $tableBlastTMP{$pair}{$mis}{$acc}{'END_F'} = $end_f;
        $tableBlastTMP{$pair}{$mis}{$acc}{'START_R'} = $start_r;
        $tableBlastTMP{$pair}{$mis}{$acc}{'END_R'} = $end_r;
        $tableBlastTMP{$pair}{$mis}{$acc}{'AL_F'} = $al_f;
        $tableBlastTMP{$pair}{$mis}{$acc}{'AL_R'} = $al_r;
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
                $tableBlast{$pair}{$mis}{$acc}{'LEN'} = $tableBlastTMP{$pair}{$mis}{$acc}{'LEN'};
                $tableBlast{$pair}{$mis}{$acc}{'START_F'} = $tableBlastTMP{$pair}{$mis}{$acc}{'START_F'};
                $tableBlast{$pair}{$mis}{$acc}{'END_F'} = $tableBlastTMP{$pair}{$mis}{$acc}{'END_F'};
                $tableBlast{$pair}{$mis}{$acc}{'START_R'} = $tableBlastTMP{$pair}{$mis}{$acc}{'START_R'};
                $tableBlast{$pair}{$mis}{$acc}{'END_R'} = $tableBlastTMP{$pair}{$mis}{$acc}{'END_R'};
                $tableBlast{$pair}{$mis}{$acc}{'AL_F'} = $tableBlastTMP{$pair}{$mis}{$acc}{'AL_F'};
                $tableBlast{$pair}{$mis}{$acc}{'AL_R'} = $tableBlastTMP{$pair}{$mis}{$acc}{'AL_R'};
              }
            }
          }
        }
      }
}
