#!/usr/bin/perl

#Primer design

use strict;
use POSIX;

#Usage: ./Create_README_pp.pl -folder $path_html/analysesPhyloprimer/BZpTYlKPtufEzMrMRTylCa -file BZpTYlKP.info5

#connected to Primer_design_pp.pl, Probe_design_pp.pl, Oligo_design_pp.pl, Primer_check_pp.pl, Probe_check_pp.pl and Oligo_check_pp.pl : create README file

#set path to folders:
my $path_html = 'path_to_html_folder';

my %args = @ARGV;

my $fileInfo = $args{-file};
my $index = substr($fileInfo, -1);
my $user = substr($fileInfo, 0, 8);

my $folder = $args{-folder};
my $index1 = substr($folder, -2);

my $path = "$path_html/analysesPhyloprimer/" . $folder;

open(my $file, ">${path}/README") or die; #file with all the necessary parameters

print $file "#########################################\n###############PhyloPrimer###############\n#########################################\n";
print $file "There are two different folders:\n";

my $checkFile = $folder . "/" . $user . ".fasta.out.alignment.negative";
my $checkFile1 = $folder . "/" . $user . ".tax.txt";
my $checkFile2 = $folder . "/" . $user . ".pro.txt";

print $file "\n\tinputs\n";

if ($index1 eq 'Ta') {
    
    print $file "\t\t$user.fasta : DNA sequences\n";
    print $file "\t\t$user.fasta.out : BLAST search output\n";
    print $file "\t\t$user.fasta.out.tree : guide tree\n";
    print $file "\t\t$user.fasta.out.alignment : DNA alignment\n";
    print $file "\t\t$user.fasta.out.alignment.positive : selected DNA sequences from dynamic tree\n";
    if (-e $checkFile) {
        print $file "\t\t$user.fasta.out.alignment.negative : DNA sequences not selected in the dynamic tree\n";
    }
    print $file "\t\t$user.consensus.positive : consensus calculated from the selected sequences\n";
    if (-e $checkFile) {
        print $file "\t\t$user.consensus.negative : consensus calculated from sequences that were not selected in the dynamic tree\n";
    }
    print $file "\t\t$user.info : oligo design parameter file\n";
    my $out_file_species = $path . "/" . (substr($folder, 0, 8)) . ".taxonomy";
    
    if (-e $out_file_species) {
        print $file "\t\t$user.taxonomy : taxonomy of the selected tree sequences where 1 corresponds to Domain, 2 to Phylum, 3 to Class, 4 to Order, 5 to Family, 6 to Genus and 7 to Species\n";
    }
    
} elsif ($index1 eq 'Tb') {
    print $file "\t\t$user.fasta : DNA sequences\n";
    if (-e $checkFile1) {
        print $file "\t\t$user.tax.txt : tab delimited taxid file\n";
    }
    if (-e $checkFile2) {
        print $file "\t\t$user.pro.txt : tab delimited gene file\n";
    }
    print $file "\t\t$user.alignment : DNA alignment\n";
    print $file "\t\t$user.fasta.tree  : guide tree\n";
    print $file "\t\t$user.alignment.positive : selected DNA sequences from dynamic tree\n";
    if (-e $checkFile) {
        print $file "\t\t$user.fasta.out.alignment.negative : DNA sequences not selected in the dynamic tree\n";
    }
    print $file "\t\t$user.consensus.positive : consensus calculated from the selected sequences\n";
    if (-e $checkFile) {
        print $file "\t\t$user.consensus.negative : consensus calculated from sequences that were not selected in the dynamic tree\n";
    }
    print $file "\t\t$user.info : oligo design parameter file\n";
    
    my $out_file_species = $path . "/" . (substr($folder, 0, 8)) . ".taxonomy";
    
    if (-e $out_file_species) {
         print $file "\t\t$user.taxonomy : taxonomy of the selected tree sequences where 1 corresponds to Domain, 2 to Phylum, 3 to Class, 4 to Order, 5 to Family, 6 to Genus and 7 to Species\n";
    }
    
    
    
} elsif ($index1 eq 'Tc') {

    print $file "\t\t$user.tree : newick tree\n";
    print $file "\t\t$user.alignment : DNA alignment\n";
    if (-e $checkFile1) {
        print $file "\t\t$user.tax.txt : tab delimited taxid file\n";
    }
    if (-e $checkFile2) {
        print $file "\t\t$user.pro.txt : tab delimited gene file\n";
    }
    
    print $file "\t\t$user.alignment.positive : selected DNA sequences from dynamic tree\n";
    if (-e $checkFile) {
        print $file "\t\t$user.fasta.out.alignment.negative : DNA sequences not selected in the dynamic tree\n";
    }
    print $file "\t\t$user.consensus.positive : consensus calculated from the selected sequences\n";
    if (-e $checkFile) {
        print $file "\t\t$user.consensus.negative : consensus calculated from sequences that were not selected in the dynamic tree\n";
    }
    print $file "\t\t$user.info : oligo design parameter file\n";
    my $out_file_species = $path . "/" . (substr($folder, 0, 8)) . ".taxonomy";
    
    if (-e $out_file_species) {
        print $file "\t\t$user.taxonomy : taxonomy of the selected tree sequences where 1 corresponds to Domain, 2 to Phylum, 3 to Class, 4 to Order, 5 to Family, 6 to Genus and 7 to Species\n";
    }
    
} elsif ($index1 eq 'Sa') {
    print $file "\t\t$user.fasta : DNA sequences\n";
    print $file "\t\t$user.alignment : DNA alignment\n";
    print $file "\t\t$user.consensus.positive : consensus sequence\n";
    print $file "\t\t$user.info${index} : oligo design parameter file\n";
} elsif ($index1 eq 'Sb') {
    print $file "\t\t$user.alignment : DNA alignment\n";
    print $file "\t\t$user.consensus.positive : consensus sequence\n";
    print $file "\t\t$user.info${index} : oligo design parameter file\n";
} elsif ($index1 eq 'Sc') {
    print $file "\t\t$user.consensus : consensus sequence\n";
    print $file "\t\t$user.consensus.positive : consensus sequence\n";
    print $file "\t\t$user.info${index} : oligo design parameter file\n";
} elsif ($index1 eq 'Ca') {
    print $file "\t\t$user.fasta : oligo DNA sequence\n";
    print $file "\t\t$user.info${index} : oligo design parameter file\n";
}


print $file "\n\tresults\n\n";

my $checkFile = $folder . "/tmp/oligo_bowtie_user.sam";

if ($index eq '1') {
    
    print $file "\t\tforwardList.txt : list of forward primers\n";
    print $file "\t\treverseList.txt : list of reverse primers\n";
    print $file "\t\tassayList.txt : list of primer pairs\n";
    print $file "\t\tselfDimer.txt : self dimer values\n";
    print $file "\t\tcrossDimer.txt : cross dimer values\n";
    print $file "\t\thairpin.txt : hairpin values\n";
    if (-e $checkFile) {
        print $file "\t\toligo_bowtie_user.sam : oligo mapping to your sequences\n";
    }
    print $file "\t\toligo_bowtie_nt.sam : oligo mapping to the DB2 database\n";
    
} elsif ($index eq '2') {
    
    print $file "\t\tforwardList.txt : list of forward primers\n";
    print $file "\t\treverseList.txt : list of reverse primers\n";
    print $file "\t\tprobeList.txt : list of probes\n";
    print $file "\t\tassayList.txt : list of primer pairs\n";
    print $file "\t\tselfDimer.txt : self dimer values\n";
    print $file "\t\tcrossDimer.txt : cross dimer values\n";
    print $file "\t\thairpin.txt : hairpin values\n";
    if (-e $checkFile) {
        print $file "\t\toligo_bowtie_user.sam : oligo mapping to your sequences\n";
    }
    print $file "\t\toligo_bowtie_nt.sam : oligo mapping to the DB2 database\n";
    
} elsif ($index eq '3') {
    
    print $file "\t\toligoList.txt : list of oligos\n";
    print $file "\t\tselfDimer.txt : self dimer values\n";
    print $file "\t\thairpin.txt : hairpin values\n";
    if (-e $checkFile) {
        print $file "\t\toligo_bowtie_user.sam : oligo mapping to your sequences\n";
    }
    print $file "\t\toligo_bowtie_nt.sam : oligo mapping to the DB2 database\n";
    
} elsif ($index eq '4') {
    
    print $file "\t\toligoList.txt : list of oligos\n";
    print $file "\t\tselfDimer.txt : self dimer values\n";
    print $file "\t\tcrossDimer.txt : cross dimer values\n";
    print $file "\t\thairpin.txt : hairpin values\n";
    if (-e $checkFile) {
        print $file "\t\toligo_bowtie_user.sam : oligo mapping to your sequences\n";
    }
    print $file "\t\toligo_bowtie_nt.sam : oligo mapping to the DB2 database\n";
    
} elsif ($index eq '5') {
    
    print $file "\t\toligoList.txt : list of oligos\n";
    print $file "\t\tselfDimer.txt : self dimer values\n";
    print $file "\t\tcrossDimer.txt : cross dimer values\n";
    print $file "\t\thairpin.txt : hairpin values\n";
    if (-e $checkFile) {
        print $file "\t\toligo_bowtie_user.sam : oligo mapping to your sequences\n";
    }
    print $file "\t\toligo_bowtie_nt.sam : oligo mapping to the DB2 database\n";
    
} elsif ($index eq '6') {
    
    print $file "\t\toligoList.txt : list of oligos\n";
    print $file "\t\tselfDimer.txt : self dimer values\n";
    print $file "\t\thairpin.txt : hairpin values\n";
    if (-e $checkFile) {
        print $file "\t\toligo_bowtie_user.sam : oligo mapping to your sequences\n";
    }
    print $file "\t\toligo_bowtie_nt.sam : oligo mapping to the DB2 database\n";
}
