#!/usr/bin/perl

use strict;
use DBD::mysql;
use POSIX;

#connected to backgroundDesignCheck.pl : probe check pipeline

#set path to folders:
my $path_cgi = 'path_to_cgi_folder';
my $path_html = 'path_to_html_folder';

#set mySQL parameters
my $dsn = "mysql_database";
my $user_name = "mysql_user";
my $password = "mysql_password";

#Usage: ./Probe_check_pp.pl -folder $path_html/analysesPhyloprimer/BZpTYlKPtufEzMrMRTylCa -file BZpTYlKP.info5

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
my $fileFasta = $nameFile . ".fasta";

#save input into array
my %all;
open(IN, "<$folder/$fileInfo") or die; #file with all the necessary parameters
while (defined(my $input = <IN>)) {
    chomp($input);
    my @data = split(/\t/, $input);
    $all{$data[0]} = $data[1];
}
close(IN);

my @mismatch = (0,1,2,3,4,5,6); #number of mismatches - max allowed is 3 per oligo (6 per oligo pair

##parameters for dG and dT calculation
my $monovalent =  $all{'MON_DG'}/1000;
my $mg_tot = $all{'MG_DG'}/1000;
my $C = $all{'OLIGO_DG'}/1000000;
my $dNTP_tot = $all{'DNTP_DG'}/1000;
my $temperature_celsius = $all{'T_DG'};

#make temporary directory
`mkdir $folder/tmp 2> /dev/null`;

#make directory that will be compressed at the end for the user
`mkdir $folder/results 2> /dev/null`;
`mkdir $folder/inputs 2> /dev/null`;

#open oligo file
my %check;
my %checkPaired;

open(IN, "<$folder/$fileFasta") or die; #file with all the necessary parameters
my $pair = 0;
my $lenMin = 300;
my $len;
while (defined(my $input = <IN>)) {
    $pair++;
    chomp($input);
    my @data = split(/\t/, $input);
    if ($data[0] =~ /-/) {
        my @oligo = split(/-/, $data[0]);
        $check{$oligo[1]} = "";
        $checkPaired{$pair}{'F'}{'TITLE'} = $data[0];
        $checkPaired{$pair}{'F'}{'SEQ'} = $oligo[1];
        $len = length($oligo[1]);
    } else {
        $check{$data[0]} = "";
        $checkPaired{$pair}{'F'}{'TITLE'} = $data[0];
        $checkPaired{$pair}{'F'}{'SEQ'} = $data[0];
        $len = length($data[0]);
    }
    if ($data[1] =~ /-/) {
        my @oligo = split(/-/, $data[1]);
        $check{$oligo[1]} = "";
        $checkPaired{$pair}{'R'}{'TITLE'} = $data[1];
        $checkPaired{$pair}{'R'}{'SEQ'} = $oligo[1];
        $len = length($oligo[1]);
    } else {
        $check{$data[1]} = "";
        $checkPaired{$pair}{'R'}{'TITLE'} = $data[1];
        $checkPaired{$pair}{'R'}{'SEQ'} = $data[1];
        $len = length($data[1]);
    }
    if ($data[2] =~ /-/) {
        my @oligo = split(/-/, $data[2]);
        $check{$oligo[1]} = "";
        $checkPaired{$pair}{'P'}{'TITLE'} = $data[2];
        $checkPaired{$pair}{'P'}{'SEQ'} = $oligo[1];
        $len = length($oligo[1]);
    } else {
        $check{$data[2]} = "";
        $checkPaired{$pair}{'P'}{'TITLE'} = $data[2];
        $checkPaired{$pair}{'P'}{'SEQ'} = $data[2];
        $len = length($data[2]);
    }
    if ($len < $lenMin) {
        $lenMin = $len;
    }
}
close(IN);

my %checkInfo;

#tm on all -- how to do with dangl???
#secondary hairpin
#secondary self
open(my $tmp1, ">$folder/tmp/Hairpin_1.tmp") or die;
open(my $tmp2, ">$folder/tmp/Self_1.tmp") or die;

foreach my $oligo (keys %check) {

    #create files for secondary structure check
    print $tmp1 "$oligo\n";
    print $tmp2 "$oligo\n";

    #GC%
    $checkInfo{$oligo}{'GC'} = gplusc($oligo); #check for gc content

    #len
    $checkInfo{$oligo}{'LEN'} = length($oligo);

    #Tm
    my $tm = `$path_cgi/tm_calculation_pp.pl -primer $oligo -type primer -sense F -mg $mg_tot -dang X -mon $monovalent -oligo $C -dntp $dNTP_tot`; #Tm is the same for forward and reverse
    chomp($tm);
    $tm = sprintf("%.2f", $tm);
    $checkInfo{$oligo}{'TM'} = $tm;

}
close($tmp1);
close($tmp2);

#Hairpin formation
`$path_cgi/Hairpin_checking_web_pp.pl -folder $folder -primer Hairpin_1.tmp -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius`; #Tm is the same for forward and reverse. tm is not ok

#Self dimer formation
`$path_cgi/Self_dimer_checking_web_pp.pl -folder $folder -primer Self_1.tmp -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius`; #Tm is the same for forward and reverse. tm is not ok

#Hairpin data
open(IN, "<$folder/tmp/SecondaryStructure_Hairpin_1.tmp.delta") or die;
while(defined(my $input = <IN>)) { #only for forward
    chomp($input);
    my ($oligo, $dG_hair) = split(/\t/, $input);
    if (($dG_hair eq '') or ($dG_hair eq '0.00')) {
        $dG_hair = 0;
    }
    $checkInfo{$oligo}{'HAIR'} = $dG_hair;
}

#Self dimer data
open(IN, "<$folder/tmp/SecondaryStructure_Self_1.tmp.delta") or die;
while(defined(my $input = <IN>)) { #only for forward
    chomp($input);
    my ($oligo, $dG_self) = split(/\t/, $input);
    if (($dG_self eq '') or ($dG_self eq '0.00')) {
        $dG_self = 0;
    }
    $checkInfo{$oligo}{'SELF'} = $dG_self;
}

#calculate cross self dimers
open(my $tmp, ">$folder/tmp/CrossSelf_1.tmp") or die;
foreach my $oligo (keys %check) {
    if ($oligo =~ /[RYSWKMBDHVN]/) { #if deg bases: find all possible oligos
        my ($pRef) = degenerateAlt($oligo);
        my @primerDeg = @{$pRef};
        print $tmp "$oligo"; #original
        my $num = scalar(@primerDeg) - 1; #number of oligos in the pair
        foreach my $c1 (0..($num-1)) { #first oligo
            foreach my $c2 (($c1+1)..$num) { #second
                print $tmp "\t$primerDeg[$c1],$primerDeg[$c2]";
            }
        }
    }
    print $tmp "\n";
}
close($tmp);
`$path_cgi/Cross_dimer_checking_self_web_pp.pl -folder $folder -primer CrossSelf_1.tmp -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius`;
`cat $folder/tmp/SecondaryStructure_crossSelf_*.tmp > $folder/crossSelf.txt`;

#read crossSelf file
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

#print file for cross self dimers
open(IN, "<$folder/tmp/SecondaryStructure_Self_1.tmp") or die; #I'll have to add /tmp/
open(my $file, ">$folder/selfDimer.txt") or die;
my $oligo;
my $lenSelf1;
my $lenSelfA;
my %allCount;
my %alldG;

while(defined(my $input = <IN>)) {
    chomp($input);
    if ($input =~ /^>/) {
        ($disc, $oligo) = split(/\t/, $input);
        print $file "\n>\t$oligo\n";
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
        print $file "\n$one\n$a\n$two\ndG:\t$dG\tkcal/mol\n";
    } elsif ($input =~ /^@/) {
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
    undef %alldG;
}
close(IN);
close($file);

#cross dimers calculation
open(my $tmp, ">$folder/tmp/Cross_1.tmp") or die;
foreach my $pair (sort {$a<=>$b} keys %checkPaired) { # check that in order - need to modyfy Cross script
    print $tmp "$checkPaired{$pair}{'F'}{'SEQ'}\t$checkPaired{$pair}{'R'}{'SEQ'},";
    print $tmp "$checkPaired{$pair}{'R'}{'SEQ'}\t$checkPaired{$pair}{'P'}{'SEQ'},";
    print $tmp "$checkPaired{$pair}{'F'}{'SEQ'}\t$checkPaired{$pair}{'P'}{'SEQ'}\n";
}
close($tmp);
`$path_cgi/Cross_dimer_checking_probe_web_pp.pl -folder $folder -primer Cross_1.tmp -mg $mg_tot -mon $monovalent -oligo $C -dntp $dNTP_tot -t $temperature_celsius`;
`cat $folder/tmp/SecondaryStructure_crossDimer_*.tmp > $folder/crossDimer.txt`;
my %checkCross;
my %checkTa;
my $pair = 0;
open(IN, "<$folder/crossDimer.txt") or die;
while(defined(my $input = <IN>)) {
    chomp($input);
    if ($input =~ /^>/) {
        $pair++;
    } elsif ($input =~ /^@/) {
        my $Ta;
        my ($disc, $dG_cross, $unit) = split(/\t/, $input); # unit of measure
        $checkCross{$pair} = $dG_cross;
        if ($checkInfo{$checkPaired{$pair}{'R'}{'SEQ'}}{'TM'} >= $checkInfo{$checkPaired{$pair}{'F'}{'SEQ'}}{'TM'}) { #need to check how is for probes
            $Ta = $checkInfo{$checkPaired{$pair}{'F'}{'SEQ'}}{'TM'} - 5;
        } else {
            $Ta = $checkInfo{$checkPaired{$pair}{'R'}{'SEQ'}}{'TM'} - 5;
        }
        $checkTa{$pair} = $Ta;
    }
}
close(IN);

#create file with all the primers for the blast search
open(my $tmp, ">$folder/tmp/oligo.txt") or die;
foreach my $pair (sort keys %checkPaired) {
    my $len_f = length($checkPaired{$pair}{'F'}{'SEQ'});
    my $len_r = length($checkPaired{$pair}{'R'}{'SEQ'});
    my $len_p = length($checkPaired{$pair}{'P'}{'SEQ'});
    print $tmp "$pair\t$checkPaired{$pair}{'F'}{'SEQ'}\t$len_f\tfor\n";
    print $tmp "$pair\t$checkPaired{$pair}{'R'}{'SEQ'}\t$len_r\trev\n";
    print $tmp "$pair\t$checkPaired{$pair}{'P'}{'SEQ'}\t$len_p\tprobe\n";
}
close($tmp);

#perform blast + bowtie check
`perl $path_cgi/blast_bowtie_check_pp.pl -folder $folder -type nt`;

my %accessionMySQL;
my %tableBlast;

my $checkFile = $folder . "/tmp/inSilico_nt.txt";

my %foundSp;

my $emptyBLAST = 0;
my $tableNt = ">";
my $pieChartTableFRP;

if (-z $checkFile) { #if file is empty
    $emptyBLAST = 1;
    $tableNt = "none";
    $pieChartTableFRP = "none";
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
        $sth = $dbh->prepare("SELECT * FROM taxid_taxonomy_pp WHERE ($ask)");
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

    foreach my $pair (sort {$a <=> $b} keys %checkPaired) {

        my %pieFRP;

        #create BLAST table for additional user BLAST search
        my $for = $checkPaired{$pair}{'F'}{'SEQ'};
        my $rev = $checkPaired{$pair}{'R'}{'SEQ'};
        my $probe = $checkPaired{$pair}{'P'}{'SEQ'};
        my $combined = $checkPaired{$pair}{'F'}{'TITLE'} . "-" . $checkPaired{$pair}{'R'}{'TITLE'} . "-" . $checkPaired{$pair}{'P'}{'TITLE'};

        my $tableFRP;
        my $tableFR;

        my $countFRP=0;
        my $countFR=0;
        my $countFP=0;
        my $countRP=0;
        my $countF=0;
        my $countR=0;
        my $countP=0;

        foreach my $mis (@mismatch) { #print first alignments with less mismatches
            foreach my $acc (sort {$tableBlast{$pair}{$mis}{$a}{'MATCH'} <=> $tableBlast{$pair}{$mis}{$b}{'MATCH'}} keys %{$tableBlast{$pair}{$mis}}) {
                $countFRP++;

                ##save all species
                $foundSp{$pair}{$taxonomy{$acc}{'SPECIES'}} = '';

                #data for pieChart
                $pieFRP{'DOMAIN'}{$taxonomy{$acc}{'DOMAIN'}}++;
                $pieFRP{'PHYLUM'}{$taxonomy{$acc}{'PHYLUM'}}++;
                $pieFRP{'CLASS'}{$taxonomy{$acc}{'CLASS'}}++;
                $pieFRP{'ORDER'}{$taxonomy{$acc}{'ORDER'}}++;
                $pieFRP{'FAMILY'}{$taxonomy{$acc}{'FAMILY'}}++;
                $pieFRP{'GENUS'}{$taxonomy{$acc}{'GENUS'}}++;
                $pieFRP{'SPECIES'}{$taxonomy{$acc}{'SPECIES'}}++;

                if ($tableBlast{$pair}{$mis}{$acc}{'START_P'} ne 'no') { #match between forward + reverse + probe

                    #accessions present in both forward and reverse primers
                    $tableFRP .= "<tr><td rowspan='3'>" . $acc . "</td>";
                    $tableFRP .= "<td>F</td>";
                    $tableFRP .= "<td>" . $for . "<br>" . $tableBlast{$pair}{$mis}{$acc}{'AL_F'} . "</td>";
                    $tableFRP .= "<td>" . $tableBlast{$pair}{$mis}{$acc}{'START_F'} . "</td>";
                    $tableFRP .= "<td>" . $tableBlast{$pair}{$mis}{$acc}{'END_F'} . "</td>";
                    $tableFRP .= "<td rowspan='3'>" . $tableBlast{$pair}{$mis}{$acc}{'LEN'} . "</td>";

                    #taxonomy
                    if ($taxonomy{$acc}{'DOMAIN'} eq "") {
                        $tableFRP .= "<td rowspan='3' class ='D'>Unclassified</td>";
                    } else {
                        $tableFRP .= "<td rowspan='3' class ='D'>" . $taxonomy{$mis}{'DOMAIN'} . "</td>";
                    }
                    if ($taxonomy{$acc}{'PHYLUM'} eq "") {
                        $tableFRP .= "<td rowspan='3' class ='P'>Unclassified</td>";
                    } else {
                        $tableFRP .= "<td rowspan='3' class ='P'>" . $taxonomy{$acc}{'PHYLUM'} . "</td>";
                    }
                    if ($taxonomy{$acc}{'CLASS'} eq "") {
                        $tableFRP .= "<td rowspan='3' class ='C'>Unclassified</td>";
                    } else {
                        $tableFRP .= "<td rowspan='3' class ='C'>" . $taxonomy{$acc}{'CLASS'} . "</td>";
                    }
                    if ($taxonomy{$acc}{'ORDER'} eq "") {
                        $tableFRP .= "<td rowspan='3' class ='O'>Unclassified</td>";
                    } else {
                        $tableFRP .= "<td rowspan='3' class ='O'>" . $taxonomy{$acc}{'ORDER'} . "</td>";
                    }
                    if ($taxonomy{$acc}{'FAMILY'} eq "") {
                        $tableFRP .= "<td rowspan='3' class ='F'>Unclassified</td>";
                    } else {
                        $tableFRP .= "<td rowspan='3' class ='F'>" . $taxonomy{$acc}{'FAMILY'} . "</td>";
                    }
                    if ($taxonomy{$acc}{'GENUS'} eq "") {
                        $tableFRP .= "<td rowspan='3' class ='G'>Unclassified</td>";
                    } else {
                        $tableFRP .= "<td rowspan='3' class ='G'>" . $taxonomy{$acc}{'GENUS'} . "</td>";
                    }
                    if ($taxonomy{$acc}{'SPECIES'} eq "") {
                        $tableFRP .= "<td rowspan='3' class ='S'>Unclassified</td></tr>";
                    } else {
                        $tableFRP .= "<td rowspan='3' class ='S'>" . $taxonomy{$acc}{'SPECIES'} . "</td></tr>";
                    }

                    $tableFRP .= "<tr><td>R</td>";
                    $tableFRP .= "<td>" . $rev . "<br>" . $tableBlast{$pair}{$mis}{$acc}{'AL_R'} . "</td>";
                    $tableFRP .= "<td>" . $tableBlast{$pair}{$mis}{$acc}{'START_R'} . "</td>";
                    $tableFRP .= "<td>" . $tableBlast{$pair}{$mis}{$acc}{'END_R'} . "</td></tr>";

                    $tableFRP .= "<tr><td>P</td>";
                    $tableFRP .= "<td>" . $probe . "<br>" . $tableBlast{$pair}{$mis}{$acc}{'AL_P'} . "</td>";
                    $tableFRP .= "<td>" . $tableBlast{$pair}{$mis}{$acc}{'START_P'} . "</td>";
                    $tableFRP .= "<td>" . $tableBlast{$pair}{$mis}{$acc}{'END_P'} . "</td></tr>";

                } else { #match between forward + reverse

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
        if (($tableFRP ne '') or ($tableFR ne '')) { # if there is at least one that is defined
            $tableNt .= $combined . "<table>";

            if ($tableFRP ne '') {
                $tableNt .= $tableFRP;
                $tableNt .= "<tr><td colspan='8' id='dividingCell'></td></tr>";
            }
            if ($tableFR ne '') {
                $tableNt .= $tableFR;
                $tableNt .= "<tr><td colspan='8' id='dividingCell'></td></tr>";
            }
            $tableNt =~ s/<tr><td colspan='8' id='dividingCell'><\/td><\/tr>$//g;
            $tableNt .= "</table>";
        }

        #assemble pieChartTable - FRP
        if (defined($pieFRP{'DOMAIN'})) {
            $pieChartTableFRP .= "<" . $combined . ">";

            foreach my $rank ('DOMAIN','PHYLUM','CLASS','ORDER','FAMILY','GENUS','SPECIES') {
                foreach my $taxon (keys %{$pieFRP{$rank}}) {
                    if ($taxon eq "") {
                        $pieChartTableFRP .= "Unclassified-" . $pieFRP{$rank}{$taxon} . "-Brown,";
                    } else {
                        $pieChartTableFRP .= $taxon . "-" . $pieFRP{$rank}{$taxon} . "-" . $taxonomyColour{$taxon} . ",";
                    }
                }
                chop($pieChartTableFRP);
                $pieChartTableFRP .= ":";
            }
            chop($pieChartTableFRP);
            $pieChartTableFRP .= "<>";
        }
    }
}
undef %tableBlast;
undef %accessionMySQL;

#blast if negative user file is present
my $tableUserCheck = ">";

if ($all{'NEGATIVE_FILE'} eq 'yes') {

    #perform blast + bowtie check
    `perl $path_cgi/blast_bowtie_check_pp.pl -folder $folder -type user`;
    my $checkFile = $folder . "/tmp/inSilico_user.txt";

    if (-z $checkFile) { #if file is empty
        $tableUserCheck = "none";
    } else { #if file is not empty

        #retrieve information from BLAST file
        checkBLAST('inSilico_user.txt');

        my $negative = (substr($folder, 56, 8)) . ".negativefasta";
        my $countAll = `grep -c "^>" $folder/$negative`;
        chomp($countAll);

        foreach my $pair (sort {$a <=> $b} keys %checkPaired) {

            #create BLAST table & pieCharts
            my $for = $checkPaired{$pair}{'F'}{'SEQ'};
            my $rev = $checkPaired{$pair}{'R'}{'SEQ'};
            my $probe = $checkPaired{$pair}{'P'}{'SEQ'};
            my $combined = $checkPaired{$pair}{'F'}{'TITLE'} . "-" . $checkPaired{$pair}{'R'}{'TITLE'} . "-" . $checkPaired{$pair}{'P'}{'TITLE'};

            my $tableFRP;
            my $tableFR;

            my $countFRP = 0;
            my $countFR = 0;

            foreach my $mis (@mismatch) { #print first alignments with less mismatches
                if (defined($tableBlast{$pair}{$mis})) {
                    foreach my $acc (keys %{$tableBlast{$pair}{$mis}}) { ###I need to order them somehow
                        if ($tableBlast{$pair}{$mis}{$acc}{'START_P'} ne 'no') { #match between forward + reverse + probe

                            $countFRP++;

                            #accessions present in both forward and reverse primers
                            $tableFRP .= "<tr><td rowspan='3'>" . $acc . "</td>";
                            $tableFRP .= "<td>F</td>";
                            $tableFRP .= "<td>" . $for . "<br>" . $tableBlast{$pair}{$mis}{$acc}{'AL_F'} . "</td>";
                            $tableFRP .= "<td>" . $tableBlast{$pair}{$mis}{$acc}{'START_F'} . "</td>";
                            $tableFRP .= "<td>" . $tableBlast{$pair}{$mis}{$acc}{'END_F'} . "</td>";
                            $tableFRP .= "<td rowspan='3'>" . $tableBlast{$pair}{$mis}{$acc}{'LEN'} . "</td>";

                            $tableFRP .= "<tr><td>R</td>";
                            $tableFRP .= "<td>" . $rev . "<br>" . $tableBlast{$pair}{$mis}{$acc}{'AL_R'} . "</td>";
                            $tableFRP .= "<td>" . $tableBlast{$pair}{$mis}{$acc}{'START_R'} . "</td>";
                            $tableFRP .= "<td>" . $tableBlast{$pair}{$mis}{$acc}{'END_R'} . "</td></tr>";

                            $tableFRP .= "<tr><td>P</td>";
                            $tableFRP .= "<td>" . $probe . "<br>" . $tableBlast{$pair}{$mis}{$acc}{'AL_P'} . "</td>";
                            $tableFRP .= "<td>" . $tableBlast{$pair}{$mis}{$acc}{'START_P'} . "</td>";
                            $tableFRP .= "<td>" . $tableBlast{$pair}{$mis}{$acc}{'END_P'} . "</td></tr>";

                        } else { #match between forward + reverse

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
            }
            #compose table
            if (($tableFRP ne '') or ($tableFR ne '')) { # if there is at least one that is defined
                $tableUserCheck .= $combined . "<" . $countAll . ";" . $countFRP . ";" . $countFR . "><table>";

                if ($tableFRP ne '') {
                    $tableUserCheck .= $tableFRP;
                    $tableUserCheck .= "<tr><td colspan='8' id='dividingCell'></td></tr>";
                }
                if ($tableFR ne '') {
                    $tableUserCheck .= $tableFR;
                    $tableUserCheck .= "<tr><td colspan='8' id='dividingCell'></td></tr>";
                }
                $tableUserCheck =~ s/<tr><td colspan='8' id='dividingCell'><\/td><\/tr>$//g;
                $tableUserCheck .= "</table>";
            }
        }
    }
} else {
    $tableUserCheck = "no";
}

#print the data
open(my $file, ">$folder/info.primer") or die; #file with all the necessary parameters - result page
open(my $file1, ">$folder/results/oligoList.txt") or die; #file with all the necessary parameters - user data

print $file "PROJECT\t$all{'PROJECT'}\n";

foreach my $pair (sort {$a <=> $b} keys %checkPaired) {
    print $file "INDEX\t$pair\n";
    print $file1 "INDEX\t$pair\n";
    print $file "FORWARD\t$checkPaired{$pair}{'F'}{'TITLE'}\n";
    print $file1 "FORWARD\t$checkPaired{$pair}{'F'}{'TITLE'}\n";
    print $file "LENGTH\t$checkInfo{$checkPaired{$pair}{'F'}{'SEQ'}}{'LEN'}\n";
    print $file1 "LENGTH\t$checkInfo{$checkPaired{$pair}{'F'}{'SEQ'}}{'LEN'} bases\n";
    print $file "GC\t$checkInfo{$checkPaired{$pair}{'F'}{'SEQ'}}{'GC'}\n";
    print $file1 "GC\t$checkInfo{$checkPaired{$pair}{'F'}{'SEQ'}}{'GC'} %\n";
    print $file "TM\t$checkInfo{$checkPaired{$pair}{'F'}{'SEQ'}}{'TM'}\n";
    print $file1 "TM\t$checkInfo{$checkPaired{$pair}{'F'}{'SEQ'}}{'TM'} °C\n";
    if ($checkInfo{$checkPaired{$pair}{'F'}{'SEQ'}}{'SELF'} == 0) {
        print $file "SELF\t>=0\n";
        print $file1 "SELF\t>=0 kcal/mol\n";
    } else {
        print $file "SELF\t$checkInfo{$checkPaired{$pair}{'F'}{'SEQ'}}{'SELF'}\n";
        print $file1 "SELF\t$checkInfo{$checkPaired{$pair}{'F'}{'SEQ'}}{'SELF'} kcal/mol\n";
    }
    if ($checkInfo{$checkPaired{$pair}{'F'}{'SEQ'}}{'HAIR'} == 0) {
        print $file "HAIR\t>=0\n";
        print $file1 "HAIR\t>=0 kcal/mol\n";
    } else {
        print $file "HAIR\t$checkInfo{$checkPaired{$pair}{'F'}{'SEQ'}}{'HAIR'}\n";
        print $file1 "HAIR\t$checkInfo{$checkPaired{$pair}{'F'}{'SEQ'}}{'HAIR'} kcal/mol\n";
    }
    print $file "REVERSE\t$checkPaired{$pair}{'R'}{'TITLE'}\n";
    print $file1 "REVERSE\t$checkPaired{$pair}{'R'}{'TITLE'}\n";
    print $file "LENGTH\t$checkInfo{$checkPaired{$pair}{'R'}{'SEQ'}}{'LEN'}\n";
    print $file1 "LENGTH\t$checkInfo{$checkPaired{$pair}{'R'}{'SEQ'}}{'LEN'} bases\n";
    print $file "GC\t$checkInfo{$checkPaired{$pair}{'R'}{'SEQ'}}{'GC'}\n";
    print $file1 "GC\t$checkInfo{$checkPaired{$pair}{'R'}{'SEQ'}}{'GC'} %\n";
    print $file "TM\t$checkInfo{$checkPaired{$pair}{'R'}{'SEQ'}}{'TM'}\n";
    print $file1 "TM\t$checkInfo{$checkPaired{$pair}{'R'}{'SEQ'}}{'TM'} °C\n";

    if ($checkInfo{$checkPaired{$pair}{'R'}{'SEQ'}}{'SELF'} == 0) {
        print $file "SELF\t>=0\n";
        print $file1 "SELF\t>=0 kcal/mol\n";
    } else {
        print $file "SELF\t$checkInfo{$checkPaired{$pair}{'R'}{'SEQ'}}{'SELF'}\n";
        print $file1 "SELF\t$checkInfo{$checkPaired{$pair}{'R'}{'SEQ'}}{'SELF'} kcal/mol\n";
    }
    if ($checkInfo{$checkPaired{$pair}{'R'}{'SEQ'}}{'HAIR'} == 0) {
        print $file "HAIR\t>=0\n";
        print $file1 "HAIR\t>=0 kcal/mol\n";
    } else {
        print $file "HAIR\t$checkInfo{$checkPaired{$pair}{'R'}{'SEQ'}}{'HAIR'}\n";
        print $file1 "HAIR\t$checkInfo{$checkPaired{$pair}{'R'}{'SEQ'}}{'HAIR'} kcal/mol\n";
    }
    print $file "PROBE\t$checkPaired{$pair}{'P'}{'TITLE'}\n";
    print $file1 "PROBE\t$checkPaired{$pair}{'P'}{'TITLE'}\n";
    print $file "LENGTH\t$checkInfo{$checkPaired{$pair}{'P'}{'SEQ'}}{'LEN'}\n";
    print $file1 "LENGTH\t$checkInfo{$checkPaired{$pair}{'P'}{'SEQ'}}{'LEN'} bases\n";
    print $file "GC\t$checkInfo{$checkPaired{$pair}{'P'}{'SEQ'}}{'GC'}\n";
    print $file1 "GC\t$checkInfo{$checkPaired{$pair}{'P'}{'SEQ'}}{'GC'} %\n";
    print $file "TM\t$checkInfo{$checkPaired{$pair}{'P'}{'SEQ'}}{'TM'}\n";
    print $file1 "TM\t$checkInfo{$checkPaired{$pair}{'P'}{'SEQ'}}{'TM'} °C\n";

    if ($checkInfo{$checkPaired{$pair}{'P'}{'SEQ'}}{'SELF'} == 0) {
        print $file "SELF\t>=0\n";
        print $file1 "SELF\t>=0 kcal/mol\n";
    } else {
        print $file "SELF\t$checkInfo{$checkPaired{$pair}{'P'}{'SEQ'}}{'SELF'}\n";
        print $file1 "SELF\t$checkInfo{$checkPaired{$pair}{'P'}{'SEQ'}}{'SELF'} kcal/mol\n";
    }
    if ($checkInfo{$checkPaired{$pair}{'P'}{'SEQ'}}{'HAIR'} == 0) {
        print $file "HAIR\t>=0\n";
        print $file1 "HAIR\t>=0 kcal/mol\n";
    } else {
        print $file "HAIR\t$checkInfo{$checkPaired{$pair}{'P'}{'SEQ'}}{'HAIR'}\n";
        print $file1 "HAIR\t$checkInfo{$checkPaired{$pair}{'P'}{'SEQ'}}{'HAIR'} kcal/mol\n";
    }
    print $file "TA\t$checkTa{$pair}\n";
    print $file1 "TA\t$checkTa{$pair} °C\n";
    if ($checkCross{$pair} == 0) {
        print $file "CROSS\t>=0\n";
        print $file1 "CROSS\t>=0 kcal/mol\n";
    } else {
        print $file "CROSS\t$checkCross{$pair}\n";
        print $file1 "CROSS\t$checkCross{$pair} kcal/mol\n";
    }
    my $allSp;
    foreach my $sp (keys %{$foundSp{$pair}}) {
        $allSp .= $sp . ";";
    }
    print $file "SPECIES\t$allSp\n\\\\\n";
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

print $file "BLAST_TABLE\t$tableNt\n";
print $file "PIECHART_FRP\t$pieChartTableFRP\n";
print $file "USER_TABLE\t$tableUserCheck\n";
close($file);

#send email
my $to = $all{'EMAIL'};
my $from = 'gv16363@bristol.ac.uk';
my $subject = 'Phyloprimer results - ' . $all{'PROJECT'};

#deconstruct folder name
$folder =~ s/\/var\/www\/cerealgenomics\/phyloprimer\/analysesPhyloprimer\///g;

my $input_kind = chop($folder); #get last letter and understand if a, b or c
my $input_ST = chop($folder); #get last letter and understand if T or S
my $defSet = $folder . $input_kind . $input_ST;
$folder = $folder . $input_ST . $input_kind;

`$path_cgi/Create_README_pp.pl -folder $folder -file $fileInfo`; #create README.txt

`cp $path_html/analysesPhyloprimer/${folder}/selfDimer.txt $path_html/analysesPhyloprimer/${folder}/results/selfDimer.txt`;

`mv $path_html/analysesPhyloprimer/${folder}/tmp/SecondaryStructure_Hairpin_1.tmp $path_html/analysesPhyloprimer/${folder}/hairpin.txt`;
`cp $path_html/analysesPhyloprimer/${folder}/hairpin.txt $path_html/analysesPhyloprimer/${folder}/results/hairpin.txt`;

`cp $path_html/analysesPhyloprimer/${folder}/crossDimer.txt $path_html/analysesPhyloprimer/${folder}/results/crossDimer.txt`;

`mv $path_html/analysesPhyloprimer/${folder}/${nameFile}* $path_html/analysesPhyloprimer/${folder}/inputs/`;

my $folderAll = $path_html . "/analysesPhyloprimer/" . $folder;

#move bowtie files
`mv $path_html/analysesPhyloprimer/${folder}/tmp/oligo_bowtie*.sam $path_html/analysesPhyloprimer/${folder}/results/`;

#remove tmp folder
#`rm -r ${folderAll}/tmp`;

#zip the folder
`zip -r ${folderAll}/PhyloPrimer_${nameFile}.zip ${folder} -x ${folder}/*txt -x ${folder}/info.primer -x "${folder}/tmp/*"`;

`chown www-data:www-data ${folderAll}/info.primer`;
`chown www-data:www-data ${folderAll}/*txt`;
`chown www-data:www-data ${folderAll}/PhyloPrimer_${nameFile}.zip`;

#Send email
my $message = "Hi,\n\nPlease find the link to the PhyloPrimer results: https:\/\/www.cerealsdb.uk.net\/cerealgenomics\/cgi-bin\/phyloprimerResultsCheckProbe.cgi?defSet=" . $defSet . "\n\n\nAll the best,\nPhyloPrimer team";

open(MAIL, "|/usr/sbin/sendmail -t");

## Email Header
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
        chomp($input);
        my ($mis, $pair, $acc, $al_f, $al_r, $start_f, $end_f, $start_r, $end_r, $len, $al_p, $start_p, $end_p) = split(/\t/, $input);
        $tableBlastTMP{$pair}{$mis}{$acc}{'LEN'} = $len;
        $tableBlastTMP{$pair}{$mis}{$acc}{'START_F'} = $start_f;
        $tableBlastTMP{$pair}{$mis}{$acc}{'END_F'} = $end_f;
        $tableBlastTMP{$pair}{$mis}{$acc}{'START_R'} = $start_r;
        $tableBlastTMP{$pair}{$mis}{$acc}{'END_R'} = $end_r;
        $tableBlastTMP{$pair}{$mis}{$acc}{'AL_F'} = $al_f;
        $tableBlastTMP{$pair}{$mis}{$acc}{'AL_R'} = $al_r;
        $tableBlastTMP{$pair}{$mis}{$acc}{'START_P'} = $start_p;
        $tableBlastTMP{$pair}{$mis}{$acc}{'END_P'} = $end_p;
        $tableBlastTMP{$pair}{$mis}{$acc}{'AL_P'} = $al_p;
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
                $tableBlast{$pair}{$mis}{$acc}{'START_P'} = $tableBlastTMP{$pair}{$mis}{$acc}{'START_P'};
                $tableBlast{$pair}{$mis}{$acc}{'END_P'} = $tableBlastTMP{$pair}{$mis}{$acc}{'END_P'};
                $tableBlast{$pair}{$mis}{$acc}{'AL_P'} = $tableBlastTMP{$pair}{$mis}{$acc}{'AL_P'};
              }
            }
          }
        }
      }
}
