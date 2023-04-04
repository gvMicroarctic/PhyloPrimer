#! /usr/bin/perl -w

use strict;
use POSIX;

#connected to Primer_check_pp.pl, Probe_check_pp.pl and Oligo_check_pp.pl: perform blast + bowtie search on Phyloprimer checked oligos

#usage: ./blast_bowtie_pp.pl -folder folder -type [nt-user]

#set path to folders:
my $path_db =  'path_to_DB2';

my %args = @ARGV;

my $folder = $args{-folder};
my $type = $args{-type};

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

#open file with all oligo infomation
open(FILEIN, "<$folder/tmp/oligo.txt") or die "Couldn't open the file: $!";

my %oligoFile;
my %oligoPair;
my $pair = 0;

my $typeAll = 'primer';

#define oligo pairs
while (defined(my $input = <FILEIN>)) {
    chomp($input);
    my ($pair, $seq, $len, $type) = split(/\t/, $input);
    $oligoFile{$seq} = $len; #oligo = length sequence
    if ($type eq 'for') {
        $oligoPair{$pair}{'for'} = $seq; #pair - type = seq
    } elsif ($type eq 'rev') {
        $oligoPair{$pair}{'rev'} = $seq; #pair - type = seq
    } elsif ($type eq 'probe') {
        $oligoPair{$pair}{'probe'} = $seq; #pair - type = seq
        $typeAll = 'probe';
    } elsif ($type eq 'oligo') {
        $oligoPair{$pair}{'oligo'} = $seq; #pair - type = seq - so then it treated as for
        $typeAll = 'oligo';
    }
}
close(FILEIN);

#create file with all oligo sequences without degenerate bases
open(my $file, ">$folder/tmp/bowtieOligo.fasta") or die "Couldn't open the file: $!";
foreach my $oligo (keys %oligoFile) {
    if ($oligo =~ /[RYSWKMBDHVN]/) { #if oligo sequence has degenearte bases
        my ($pRef) = degenerateAlt($oligo);
        my @primerDeg = @{$pRef};
        foreach my $pD (@primerDeg) {
            print $file ">$oligo\n$pD\n"; #title is the same primer and seq is with substitutive bases from degenerate bases
        }
    } else {
        print $file ">$oligo\n$oligo\n";
    }
}
close($file);

#open/close files in case blast does not recover any match and this script exits
#I create these empty files so then the main script can find them even if empty
open(my $fileTMP, ">$folder/tmp/inSilico_nt.txt") or die "Couldn't open the file: $!";
close($fileTMP);
open(my $fileTMP1, ">$folder/tmp/inSilico_user.txt") or die "Couldn't open the file: $!";
close($fileTMP1);

if ($typeAll ne 'oligo') { #probe or primer

    my %add;
    #defined which taxonomy check to perform
    if ($type eq 'nt') {

        open(my $fileF, ">$folder/tmp/forward.fasta") or die "Couldn't open the file: $!"; #forward oligos
        open(my $fileR, ">$folder/tmp/reverse.fasta") or die "Couldn't open the file: $!"; #reverse oligos
        foreach my $pair (sort keys %oligoPair) {
            print $fileF ">$pair\n";
            print $fileF "$oligoPair{$pair}{'for'}\n";
            print $fileR ">$pair\n";
            print $fileR "$oligoPair{$pair}{'rev'}\n";
        }
        close($fileF);
        close($fileR);

        #BLAST of forward oligos
        `blastn -db $path_db/DB2 -query $folder/tmp/forward.fasta -word_size 7 -perc_identity 60 -qcov_hsp_perc 50 -evalue 10000 -outfmt 6 -num_threads 20 -max_target_seqs 40000 -ungapped > $folder/tmp/genbankBLAST_for.m8`;

        #retrieve only accessions that BLASTed against forward oligos
        open(FILEIN, "<$folder/tmp/genbankBLAST_for.m8") or die "Couldn't open the file: $!"; #input file with aligned sequences

        my %m8;
        my %m8Prov;

        while (defined(my $input = <FILEIN>)) {
            chomp($input);
            my @info = split /\t/, $input;
            my $pos = ceil(($info[9] + $info[8]) / 2);
            if ($info[8] < $info[9]) {
                $m8{$info[0]}{'for'}{$info[1]}{$pos} = '+'; #pair - for - accession - position = sense
                $m8Prov{$info[1]}{$pos} = '+'; #pair - position = sense
            } else {
                $m8{$info[0]}{'for'}{$info[1]}{$pos} = '-'; #pair - for - accession - position = sense
                $m8Prov{$info[1]}{$pos} = '+'; #accession - position = sense
            }
        }
        close(FILEIN);

        open($file, ">$folder/tmp/recover.txt") or die "Couldn't open the file: $!"; #input file with aligned sequences
        foreach my $acc (keys %m8Prov) {
            print $file "$acc\n";
        }
        close($file);

        undef(%m8Prov);

        #BLAST of reverse oligos
        `blastn -db $path_db/DB2 -query $folder/tmp/reverse.fasta -word_size 7 -perc_identity 60 -qcov_hsp_perc 40 -evalue 10000 -outfmt 6 -num_threads 20 -max_target_seqs 50000 -ungapped -seqidlist $folder/tmp/recover.txt > $folder/tmp/genbankBLAST_rev.m8`;

        open(FILEIN, "<$folder/tmp/genbankBLAST_rev.m8") or die "Couldn't open the file: $!"; #input file with aligned sequences
        while (defined(my $input = <FILEIN>)) {
            chomp($input);
            my @info = split /\t/, $input;
            my $pos = ceil(($info[9] + $info[8]) / 2);
            if ($info[8] < $info[9]) {
                $m8{$info[0]}{'rev'}{$info[1]}{$pos} = '+'; #pair - rev - accession - position = sense
            } else {
                $m8{$info[0]}{'rev'}{$info[1]}{$pos} = '-'; #pair - rev - accession - position = sense
            }
        }
        close(FILEIN);

        open(FILEIN, "<$folder/tmp/oligo.txt") or die "Couldn't open the file: $!";
        my $for;
        my $rev;
        my %finalAcc;
        while (defined(my $input = <FILEIN>)) {
            chomp($input);
            my ($pair, $seq, $len, $type) = split(/\t/, $input);
            #correspondance between $seq and cluster - from cluster I have info from m8For and m8Rev
            if ($type eq 'for') {
                $for = $seq;
            } elsif ($type eq 'rev') {
                $rev = $seq;
                foreach my $acc (keys %{$m8{$pair}{'for'}}) {
                    if (defined($m8{$pair}{'rev'}{$acc})) {
                        foreach my $posFor (keys %{$m8{$pair}{'for'}{$acc}}) {
                            foreach my $posRev (keys %{$m8{$pair}{'rev'}{$acc}}) {
                                my $diff = $posFor - $posRev;
                                #                                if ($m8{$pair}{'for'}{$acc}{$posFor} ne $m8{$pair}{'rev'}{$acc}{$posRev}) { #opposite directions - sense and antisense
                                if (($diff < 3000) && ($diff > -3000)) {
                                    $finalAcc{$acc}{$posFor} = ''; #accession - position = ''
                                    $finalAcc{$acc}{$posRev} = ''; #accession - position = ''
                                }
                                #                                }
                            }
                        }
                    }
                }
            }
        }
        close(FILEIN);

        open($file, ">$folder/tmp/fastaN.txt") or die "Couldn't open the file: $!"; #input file with aligned sequences
        foreach my $acc (keys %finalAcc) {
            my $min = 100000000000000000000000000000;
            my $max = 0;

            foreach my $pos (keys %{$finalAcc{$acc}}) {
                if ($pos < $min) {
                    $min = $pos;
                }
                if ($pos > $max) {
                    $max = $pos;
                }
            }
            $min = $min - 50;
            $max = $max + 50;

            if ($min < 0) {
                $min = 1;
            }
            $add{$acc} = $min; #accession = min
            print $file "$acc\t${min}-${max}\n"; #retrieve only sequence portion with valid blast matches
        }

        #retrieve sequences from nt database to use for blast of reverse oligo clusters
        `blastdbcmd -db $path_db/DB2 -entry_batch $folder/tmp/fastaN.txt > $folder/tmp/blast.fasta`;

        my $check = $folder . "/tmp/blast.fasta";

        #check size of sequence file for bowtie database construction
        my $size = -s $check;
        if ($size > 2500000000) { #2.5 Gb
            `bowtie-build $folder/tmp/blast.fasta $folder/tmp/blastBow --large-index --threads 15`;
        } else {
            `bowtie-build $folder/tmp/blast.fasta $folder/tmp/blastBow --threads 15`;
        }

        #run bowtie
        `bowtie $folder/tmp/blastBow -f $folder/tmp/bowtieOligo.fasta -l 7 -k 100000 -v 3 -I 50 -X 3000 --sam-nosq -S $folder/tmp/oligo_bowtie_nt.sam --sam-nohead --sam-nosq`;
        #N.B. bowtie can search matches with a maximum of 2 mismatches but it is accurate and performs global alignment. BBmap, bowtie2 and bwa do not perform as well as bowtie.

        open(FILEIN, "<$folder/tmp/oligo_bowtie_nt.sam") or die "Couldn't open the file: $!"; #input file with aligned sequences

    } else {
        my $negative = (substr($folder, 56, 8)) . ".negativefasta";

        my $check = $folder . "/" . $negative;

        #check size of sequence file for bowtie database construction
        my $size = -s $check;
        if ($size > 2500000000) { #2.5 Gb
            `bowtie-build $folder/$negative $folder/tmp/blastBow --large-index --threads 10`;
        } else {
            `bowtie-build $folder/$negative $folder/tmp/blastBow --threads 10`;
        }

        #run bowtie
        `bowtie $folder/tmp/blastBow -f $folder/tmp/bowtieOligo.fasta -l 7 -k 100000 -v 3 -I 50 -X 3000 --sam-nosq -S $folder/tmp/oligo_bowtie_user.sam --sam-nohead --sam-nosq`;
        #N.B. bowtie can search matches with a maximum of 2 mismatches but it is accurate and performs global alignment. BBmap, bowtie2 and bwa do not perform as well as bowtie.

        open(FILEIN, "<$folder/tmp/oligo_bowtie_user.sam") or die "Couldn't open the file: $!"; #input file with aligned sequences
    }

    my %bowtie;

    my $count = 0;

    while (defined(my $input = <FILEIN>)) {
        chomp($input);
        $count++;
        my @info = split /\t/, $input;
        my ($nm) = $info[13] =~ /NM:i:(.*)/;
        my ($md) = $info[12] =~ /MD:Z:(.*)/;
        $bowtie{$info[0]}{$info[2]}{$count}{'nm'} = $nm; #sequence - accession - acc_count - 'nm' = nm
        $bowtie{$info[0]}{$info[2]}{$count}{'md'} = $md; #sequence - accession - acc_count - 'md' = md
        $bowtie{$info[0]}{$info[2]}{$count}{'al'} = $info[3]; #sequence - accession - acc_count - 'al' = al
        $bowtie{$info[0]}{$info[2]}{$count}{'flag'} = $info[1]; #sequence - accession - acc_count - 'flag' = bitwise flag
    }
    close(FILEIN);

    my %score;

    if ($typeAll eq 'primer') {
        foreach my $pair (keys %oligoPair) {
            my $for = $oligoPair{$pair}{'for'};
            my $rev = $oligoPair{$pair}{'rev'};
            foreach my $acc (keys %{$bowtie{$for}}) {

                if (defined($bowtie{$rev}{$acc})) {
                    foreach my $countFor (keys %{$bowtie{$for}{$acc}}) {
                        foreach my $countRev (keys %{$bowtie{$rev}{$acc}}) {
                            my $diff;
                            if ($bowtie{$for}{$acc}{$countFor}{'al'} < $bowtie{$rev}{$acc}{$countRev}{'al'}) {
                                $diff = $bowtie{$rev}{$acc}{$countRev}{'al'} - $bowtie{$for}{$acc}{$countFor}{'al'} + $oligoFile{$rev}; # + 1
                            } else {
                                $diff = $bowtie{$for}{$acc}{$countFor}{'al'} - $bowtie{$rev}{$acc}{$countRev}{'al'} + $oligoFile{$for}; # + 1
                            }
                            if ($diff < 3000) { #primers are close enough
                                if ($bowtie{$for}{$acc}{$countFor}{'flag'} != $bowtie{$rev}{$acc}{$countRev}{'flag'}) { #check one is forward and one reverse
                                    if ($type eq 'nt') {
                                        if (defined($add{$acc})) {
                                            $score{$pair}{$acc}{'startFor'} = $add{$acc} + $bowtie{$for}{$acc}{$countFor}{'al'};
                                            $score{$pair}{$acc}{'endFor'} = $add{$acc} + $bowtie{$for}{$acc}{$countFor}{'al'} + $oligoFile{$for} - 1;
                                            $score{$pair}{$acc}{'startRev'} = $add{$acc} + $bowtie{$rev}{$acc}{$countRev}{'al'};
                                            $score{$pair}{$acc}{'endRev'} = $add{$acc} + $bowtie{$rev}{$acc}{$countRev}{'al'} + $oligoFile{$rev} - 1;
                                            $score{$pair}{$acc}{'diff'} = $diff;
                                            $score{$pair}{$acc}{'nmFor'} = $bowtie{$for}{$acc}{$countFor}{'md'};
                                            $score{$pair}{$acc}{'nmRev'} = $bowtie{$rev}{$acc}{$countRev}{'md'};
                                        }
                                    } else {
                                        $score{$pair}{$acc}{'startFor'} = $bowtie{$for}{$acc}{$countFor}{'al'};
                                        $score{$pair}{$acc}{'endFor'} = $bowtie{$for}{$acc}{$countFor}{'al'} + $oligoFile{$for} - 1;
                                        $score{$pair}{$acc}{'startRev'} = $bowtie{$rev}{$acc}{$countRev}{'al'};
                                        $score{$pair}{$acc}{'endRev'} = $bowtie{$rev}{$acc}{$countRev}{'al'} + $oligoFile{$rev} - 1;
                                        $score{$pair}{$acc}{'diff'} = $diff;
                                        $score{$pair}{$acc}{'nmFor'} = $bowtie{$for}{$acc}{$countFor}{'md'};
                                        $score{$pair}{$acc}{'nmRev'} = $bowtie{$rev}{$acc}{$countRev}{'md'};
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
        foreach my $pair (keys %oligoPair) {
            my $for = $oligoPair{$pair}{'for'};
            my $rev = $oligoPair{$pair}{'rev'};
            my $probe = $oligoPair{$pair}{'probe'};
            foreach my $acc (keys %{$bowtie{$for}}) {
                if (defined($bowtie{$rev}{$acc})) {
                    foreach my $countFor (keys %{$bowtie{$for}{$acc}}) {
                        foreach my $countRev (keys %{$bowtie{$rev}{$acc}}) {
                            my $diff;
                            my $min;
                            my $max;
                            if ($bowtie{$for}{$acc}{$countFor}{'al'} < $bowtie{$rev}{$acc}{$countRev}{'al'}) {
                                $diff = $bowtie{$rev}{$acc}{$countRev}{'al'} - $bowtie{$for}{$acc}{$countFor}{'al'} + $oligoFile{$rev}; # + 1
                                $min = $bowtie{$for}{$acc}{$countFor}{'al'};
                                $max = $bowtie{$rev}{$acc}{$countRev}{'al'};
                            } else {
                                $diff = $bowtie{$for}{$acc}{$countFor}{'al'} - $bowtie{$rev}{$acc}{$countRev}{'al'} + $oligoFile{$for}; # + 1
                                $min = $bowtie{$rev}{$acc}{$countRev}{'al'};
                                $max = $bowtie{$for}{$acc}{$countFor}{'al'};
                            }
                            if ($diff < 3000) { #primers are close enough
                                if ($bowtie{$for}{$acc}{$countFor}{'flag'} != $bowtie{$rev}{$acc}{$countRev}{'flag'}) { #check one is forward and one reverse
                                    if ($type eq 'nt') {
                                        if (defined($add{$acc})) {
                                            $score{$pair}{$acc}{'startFor'} = $add{$acc} + $bowtie{$for}{$acc}{$countFor}{'al'};
                                            $score{$pair}{$acc}{'endFor'} = $add{$acc} + $bowtie{$for}{$acc}{$countFor}{'al'} + $oligoFile{$for} - 1;
                                            $score{$pair}{$acc}{'startRev'} = $add{$acc} + $bowtie{$rev}{$acc}{$countRev}{'al'};
                                            $score{$pair}{$acc}{'endRev'} = $add{$acc} + $bowtie{$rev}{$acc}{$countRev}{'al'} + $oligoFile{$rev} - 1;
                                            $score{$pair}{$acc}{'diff'} = $diff;
                                            $score{$pair}{$acc}{'nmFor'} = $bowtie{$for}{$acc}{$countFor}{'md'};
                                            $score{$pair}{$acc}{'nmRev'} = $bowtie{$rev}{$acc}{$countRev}{'md'};
                                            my $inside = 0;
                                            foreach my $countProbe (keys %{$bowtie{$probe}{$acc}}) { #check probe falls inside forward and reverse
                                                if (($bowtie{$probe}{$acc}{$countProbe}{'al'} > $min) && ($bowtie{$probe}{$acc}{$countProbe}{'al'} < $max)) {
                                                    $score{$pair}{$acc}{'startProbe'} = $add{$acc} + $bowtie{$probe}{$acc}{$countProbe}{'al'};
                                                    $score{$pair}{$acc}{'endProbe'} = $add{$acc} + $bowtie{$probe}{$acc}{$countProbe}{'al'} + $oligoFile{$probe} - 1;
                                                    $score{$pair}{$acc}{'nmProbe'} = $bowtie{$probe}{$acc}{$countProbe}{'md'};
                                                    $inside = 1;
                                                }
                                            }
                                            if ($inside == 0) {
                                                $score{$pair}{$acc}{'startProbe'} = 'no';
                                                $score{$pair}{$acc}{'endProbe'} = 'no';
                                                $score{$pair}{$acc}{'nmProbe'} = 'no';
                                            }
                                        }
                                    } else {
                                        $score{$pair}{$acc}{'startFor'} = $bowtie{$for}{$acc}{$countFor}{'al'};
                                        $score{$pair}{$acc}{'endFor'} = $bowtie{$for}{$acc}{$countFor}{'al'} + $oligoFile{$for} - 1;
                                        $score{$pair}{$acc}{'startRev'} = $bowtie{$rev}{$acc}{$countRev}{'al'};
                                        $score{$pair}{$acc}{'endRev'} = $bowtie{$rev}{$acc}{$countRev}{'al'} + $oligoFile{$rev} - 1;
                                        $score{$pair}{$acc}{'diff'} = $diff;
                                        $score{$pair}{$acc}{'nmFor'} = $bowtie{$for}{$acc}{$countFor}{'md'};
                                        $score{$pair}{$acc}{'nmRev'} = $bowtie{$rev}{$acc}{$countRev}{'md'};
                                        my $inside = 0;
                                        foreach my $countProbe (keys %{$bowtie{$probe}{$acc}}) { #check probe falls inside forward and reverse
                                            if (($bowtie{$probe}{$acc}{$countProbe}{'al'} > $min) && ($bowtie{$probe}{$acc}{$countProbe}{'al'} < $max)) {
                                                $score{$pair}{$acc}{'startProbe'} = $add{$acc} + $bowtie{$probe}{$acc}{$countProbe}{'al'};
                                                $score{$pair}{$acc}{'endProbe'} = $add{$acc} + $bowtie{$probe}{$acc}{$countProbe}{'al'} + $oligoFile{$probe} - 1;
                                                $score{$pair}{$acc}{'nmProbe'} = $bowtie{$probe}{$acc}{$countProbe}{'md'};
                                                $inside = 1;
                                            }
                                        }
                                        if ($inside == 0) {
                                            $score{$pair}{$acc}{'startProbe'} = 'no';
                                            $score{$pair}{$acc}{'endProbe'} = 'no';
                                            $score{$pair}{$acc}{'nmProbe'} = 'no';
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    #create a file with start alignment, end alignment, size amplicon and alignment
    if ($type eq 'nt') {
        open($file, ">$folder/tmp/inSilico_nt.txt") or die "Couldn't open the file: $!"; #input file with aligned sequences
    } else {
        open($file, ">$folder/tmp/inSilico_user.txt") or die "Couldn't open the file: $!"; #input file with aligned sequences
    }

    foreach my $pair (sort {$a <=> $b} keys %oligoPair) {
        $count = 0;
        foreach my $acc (keys %{$score{$pair}}) {
            $count++;
            my @seq = split(/(A|C|G|T)/, $score{$pair}{$acc}{'nmFor'});
            my $alignment;
            my $mis = 0;
            foreach my $base (@seq) {
                if ($base =~ /[A|C|G|T]+/) {
                    $alignment .= $base;
                    $mis++;
                } elsif ($base > 0) {
                    $alignment .= "*" x $base;
                }
            }
            @seq = split(/(A|C|G|T)/, $score{$pair}{$acc}{'nmRev'});
            my $alignment_f = $alignment;
            $alignment = '';
            foreach my $base (@seq) {
                if ($base =~ /[A|C|G|T]+/) {
                    $alignment .= $base;
                    $mis++;
                } elsif ($base > 0) {
                    $alignment .= "*" x $base;
                }
            }
            my $alignment_r = $alignment;
            if ($typeAll eq 'primer') {
                print $file "$mis\t$pair\t$acc\t$alignment_f\t$alignment_r\t$score{$pair}{$acc}{'startFor'}\t$score{$pair}{$acc}{'endFor'}\t$score{$pair}{$acc}{'startRev'}\t$score{$pair}{$acc}{'endRev'}\t$score{$pair}{$acc}{'diff'}\n";
            } else {
                if ($score{$pair}{$acc}{'nmProbe'} ne 'no') { #if also match with probe
                    @seq = split(/(A|C|G|T)/, $score{$pair}{$acc}{'nmProbe'});
                    $alignment = '';
                    foreach my $base (@seq) {
                        if ($base =~ /[A|C|G|T]+/) {
                            $alignment .= $base;
                            $mis++;
                        } elsif ($base > 0) {
                            $alignment .= "*" x $base;
                        }
                    }
                    my $alignment_p = $alignment;
                    print $file "$mis\t$pair\t$acc\t$alignment_f\t$alignment_r\t$score{$pair}{$acc}{'startFor'}\t$score{$pair}{$acc}{'endFor'}\t$score{$pair}{$acc}{'startRev'}\t$score{$pair}{$acc}{'endRev'}\t$score{$pair}{$acc}{'diff'}\t$alignment_p\t$score{$pair}{$acc}{'startProbe'}\t$score{$pair}{$acc}{'endProbe'}\n";
                } else {
                    print $file "$mis\t$pair\t$acc\t$alignment_f\t$alignment_r\t$score{$pair}{$acc}{'startFor'}\t$score{$pair}{$acc}{'endFor'}\t$score{$pair}{$acc}{'startRev'}\t$score{$pair}{$acc}{'endRev'}\t$score{$pair}{$acc}{'diff'}\tno\t$score{$pair}{$acc}{'startProbe'}\t$score{$pair}{$acc}{'endProbe'}\n";
                }
            }
        }
    }
    close($file);

} else { #oligo
    my %add;

    #defined which taxonomy check to perform
    if ($type eq 'nt') {

        open(my $fileO, ">$folder/tmp/oligo.fasta") or die "Couldn't open the file: $!"; #forward oligos
        foreach my $pair (sort keys %oligoPair) {
            print $fileO ">$pair\n";
            print $fileO "$oligoPair{$pair}{'oligo'}\n";
        }
        close($fileO);

        #BLAST of forward oligo clusters
        `blastn -db $path_db/DB2 -query $folder/tmp/oligo.fasta -word_size 7 -perc_identity 60 -qcov_hsp_perc 50 -evalue 10000 -outfmt 6 -num_threads 20 -max_target_seqs 500 -ungapped > $folder/tmp/genbankBLAST.m8`;

        #retrieve only accessions that BLASTed against forward oligo cluster
        open(FILEIN, "<$folder/tmp/genbankBLAST.m8") or die "Couldn't open the file: $!"; #input file with aligned sequences

        my %finalAcc;

        while (defined(my $input = <FILEIN>)) {
            chomp($input);
            my @info = split /\t/, $input;
            my $pos = ceil(($info[9] + $info[8]) / 2);

            if ($info[8] < $info[9]) {
                $finalAcc{$info[1]}{$pos} = '+'; #accession - position = sense
            } else {
                $finalAcc{$info[1]}{$pos} = '-'; #accession - position = sense
            }
        }
        close(FILEIN);

        open($file, ">$folder/tmp/fastaN.txt") or die "Couldn't open the file: $!"; #input file with aligned sequences
        foreach my $acc (keys %finalAcc) {
            my $min = 100000000000000000000000000000;
            my $max = 0;

            foreach my $pos (keys %{$finalAcc{$acc}}) {
                if ($pos < $min) {
                    $min = $pos;
                }
                if ($pos > $max) {
                    $max = $pos;
                }
            }

            $min = $min - 30;
            $max = $max + 30;

            if ($min < 0) {
                $min = 1;
            }
            my $diff = $max - $min;

            if ($diff < 100000) { #ony sequences long less than 10 Kb
                $add{$acc} = $min; #accession = min
                print $file "$acc\t${min}-${max}\n"; #retrieve only sequence portion with valid blast matches
            }
        }

        #retrieve sequences from nt database to use for blast of reverse oligo clusters
        `blastdbcmd -db $path_db/DB2 -entry_batch $folder/tmp/fastaN.txt > $folder/tmp/blast.fasta`;

        my $check = $folder . "/tmp/blast.fasta";

        #check size of sequence file for bowtie database construction
        my $size = -s $check;
        if ($size > 2500000000) { #2.5 Gb
            `bowtie-build $folder/tmp/blast.fasta $folder/tmp/blastBow --large-index --threads 15`;
        } else {
            `bowtie-build $folder/tmp/blast.fasta $folder/tmp/blastBow --threads 15`;
        }

        #run bowtie
        `bowtie $folder/tmp/blastBow -f $folder/tmp/bowtieOligo.fasta -l 7 -k 100000 -v 3 -I 50 -X 3000 --sam-nosq -S $folder/tmp/oligo_bowtie_nt.sam --sam-nohead --sam-nosq`;
        #N.B. bowtie can search matches with a maximum of 2 mismatches but it is accurate and performs global alignment. BBmap, bowtie2 and bwa do not perform as well as bowtie.

        open(FILEIN, "<$folder/tmp/oligo_bowtie_nt.sam") or die "Couldn't open the file: $!"; #input file with aligned sequences

    } else {
        my $negative = (substr($folder, 56, 8)) . ".negativefasta";

        my $check = $folder . "/" . $negative;

        #check size of sequence file for bowtie database construction
        my $size = -s $check;
        if ($size > 2500000000) { #2.5 Gb
            `bowtie-build $folder/$negative $folder/tmp/blastBow --large-index --threads 10`;
        } else {
            `bowtie-build $folder/$negative $folder/tmp/blastBow --threads 10`;
        }

        #run bowtie
        `bowtie $folder/tmp/blastBow -f $folder/tmp/bowtieOligo.fasta -l 7 -k 100000 -v 3 -I 50 -X 3000 --sam-nosq -S $folder/tmp/oligo_bowtie_user.sam --sam-nohead --sam-nosq`;
        #N.B. bowtie can search matches with a maximum of 2 mismatches but it is accurate and performs global alignment. BBmap, bowtie2 and bwa do not perform as well as bowtie.

        open(FILEIN, "<$folder/tmp/oligo_bowtie_user.sam") or die "Couldn't open the file: $!"; #input file with aligned sequences
    }

    my %bowtie;

    my $count = 0;

    while (defined(my $input = <FILEIN>)) {
        chomp($input);
        $count++;
        my @info = split /\t/, $input;
        my ($nm) = $info[13] =~ /NM:i:(.*)/;
        my ($md) = $info[12] =~ /MD:Z:(.*)/;
        $bowtie{$info[0]}{$info[2]}{$count}{'nm'} = $nm; #sequence - accession - acc_count - 'nm' = nm
        $bowtie{$info[0]}{$info[2]}{$count}{'md'} = $md; #sequence - accession - acc_count - 'md' = md
        $bowtie{$info[0]}{$info[2]}{$count}{'al'} = $info[3]; #sequence - accession - acc_count - 'al' = al
        $bowtie{$info[0]}{$info[2]}{$count}{'flag'} = $info[1]; #sequence - accession - acc_count - 'al' = flag
    }
    close(FILEIN);

    my %score;

    foreach my $pair (keys %oligoPair) {
        my $oligo = $oligoPair{$pair}{'oligo'};
        foreach my $acc (keys %{$bowtie{$oligo}}) {
            foreach my $countOligo (keys %{$bowtie{$oligo}{$acc}}) {
                if ($bowtie{$oligo}{$acc}{$countOligo}{'flag'} != 4) {
                if ($type eq 'nt') {
                    if (defined($add{$acc})) {
                        $score{$pair}{$acc}{'startOligo'} = $add{$acc} + $bowtie{$oligo}{$acc}{$countOligo}{'al'};
                        $score{$pair}{$acc}{'endOligo'} = $add{$acc} + $bowtie{$oligo}{$acc}{$countOligo}{'al'} + $oligoFile{$oligo} - 1;
                        $score{$pair}{$acc}{'nmOligo'} = $bowtie{$oligo}{$acc}{$countOligo}{'md'};
                    }
                } else {
                    $score{$pair}{$acc}{'startOligo'} = $bowtie{$oligo}{$acc}{$countOligo}{'al'};
                    $score{$pair}{$acc}{'endOligo'} = $bowtie{$oligo}{$acc}{$countOligo}{'al'} + $oligoFile{$oligo} - 1;
                    $score{$pair}{$acc}{'nmOligo'} = $bowtie{$oligo}{$acc}{$countOligo}{'md'};
                }
                }
            }
        }
    }

    #create a file with start alignment, end alignment, size amplicon and alignment
    if ($type eq 'nt') {
        open($file, ">$folder/tmp/inSilico_nt.txt") or die "Couldn't open the file: $!"; #input file with aligned sequences
    } else {
        open($file, ">$folder/tmp/inSilico_user.txt") or die "Couldn't open the file: $!"; #input file with aligned sequences
    }

    foreach my $pair (sort {$a <=> $b} keys %oligoPair) {
        $count = 0;
        foreach my $acc (keys %{$score{$pair}}) {
            $count++;
            my @seq = split(/(A|C|G|T)/, $score{$pair}{$acc}{'nmOligo'});
            my $alignment;
            my $mis = 0;
            foreach my $base (@seq) {
                if ($base =~ /[A|C|G|T]+/) {
                    $alignment .= $base;
                    $mis++;
                } elsif ($base > 0) {
                    $alignment .= "*" x $base;
                }
            }
            print $file "$mis\t$pair\t$acc\t$alignment\t$score{$pair}{$acc}{'startOligo'}\t$score{$pair}{$acc}{'endOligo'}\n";
        }
    }
    close($file);
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
