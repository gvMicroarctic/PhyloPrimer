#!/usr/bin/perl

use JSON; #if not already installed, just run "cpan JSON"
use strict;
use warnings;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser );
use File::Basename;

#connected to tree_input.html : check data uploaded in Dynamic selection page

#set path to folders:
my $path_html = 'path_to_html_folder';

#Set maximum size of uploaded file
$CGI::POST_MAX = 1024 * 50000; #50 MegaBytes

#header needs to start with > and all the white spaces are changed with _
#folder name is passed only if the file is correct!

my $newhtml;
my $folderOld;
my $pass;
my %fasta;
my $fastaLen = 0;
my %replace;

#create folder with 20 random characters (N.B. it needs to be sent back to html)
my @chars = ("A".."Z", "a".."z");
$folderOld .= $chars[rand @chars] for 1..20;

#open cgi connection
my $cgi = CGI->new;
print $cgi->header('application/json;charset=UTF-8');

#import variables from html
my $input_tree = $cgi->param("input_tree"); #radio button. it can be DNA, alignment or newickh

my $in_sequence = $cgi->param("in_sequence"); #fasta sequence
my $in_fasta_phylo = $cgi->param("in_fasta_phylo"); #fasta file

my $in_al_phylo = $cgi->param("in_al_phylo"); #DNA sequence file
my $in_al_phylo_tax = $cgi->param("in_al_phylo_tax"); #tab delimited taxonomy
my $in_al_phylo_pro = $cgi->param("in_al_phylo_pro"); #tab delimited protein data

my $in_newick_phylo = $cgi->param("in_newick_phylo"); #newick file
my $in_newick_phylo_seq = $cgi->param("in_newick_phylo_seq"); #alignment file for newick
my $in_newick_phylo_tax = $cgi->param("in_newick_phylo_tax"); #tab delimited taxonomy
my $in_newick_phylo_pro = $cgi->param("in_newick_phylo_pro"); #tab delimited protein data

#convert  data to JSON
my $op = JSON -> new -> utf8 -> pretty(1);

if ($input_tree eq "gene") { #upload of few genes/sequences of interest that will be used as input for BLASt search
    my $folder = $folderOld . "Ta";
    $pass = $folderOld . "a";
    my $upload_dir = $path_html . "/uploadsPhyloprimer/" . $folder . "/"; #path to the folder
    `mkdir $upload_dir`; #create new directory
    my $filenameNew = (substr($folder, 0, 8)) . ".fasta";
    my $inside = 0;
    #fasta sequence text area
    if ($in_sequence ne "") {
        $inside++;
        my ($fastaSeq, $minL, $maxL) = sequenceload($in_sequence, $filenameNew, $upload_dir);
        if ($fastaSeq eq 'CONTENT') {
            $newhtml = "There was an error uploading your DNA sequences. All the DNA sequences must be in fasta format. Please check the correct format <a href='../phyloprimer/format_page.html#DNA_seq'>here</a>.<br>";
        } elsif ($fastaSeq eq 'UNIQUE') {
            $newhtml = "There was an error uploading your DNA sequences. Some of the sequence headings are not unique.<br>";
            
        } else {
            if ($minL == $maxL) {
                $fastaLen = $minL;
            } else {
                $fastaLen = $minL . " - " . $maxL;
            }
            if ($maxL <= 6000) {
                if ($fastaSeq == 1) {
                    $newhtml = "You uploaded $fastaSeq DNA sequence. The length of the sequence is $fastaLen bp.<br>";
                } elsif ($fastaSeq <= 10) {
                    $newhtml = "You uploaded $fastaSeq DNA sequences. The length of the sequences is $fastaLen bp.<br>";
                    
                    #mafft alignment to check gene similarity
                    `mafft --localpair --maxiterate 1000 --thread 30 --quiet ${upload_dir}/${filenameNew}  > ${upload_dir}/${filenameNew}.check`;
                    #consensus check
                    my $geneCheckPerc = `geneCheck_pp.pl -file ${upload_dir}/${filenameNew}.check`;
                    `rm ${upload_dir}/${filenameNew}.check`; #remove check file
                    
                    #if gap value too high: give warning
                    my @geneCheck = split(/-/, $geneCheckPerc);
                    if ($geneCheck[0] > 50) {
                        $newhtml .= "<b style='color:blue'\;>WARNING:</b> We checked your genes and looks like one or more sequences do not represent the same gene portion as the others.<br>";
                    }
                    
                } else {
                    $newhtml = "There was an error uploading your DNA sequences. The maximum number of allowed sequences is 10.<br>";
                }
            } else {
                $newhtml = "There was an error uploading your DNA sequences. The maximum allowed length is 6000 bp.<br>";
            }
        }
    }
    #fasta file
    if ($in_fasta_phylo ne "") {
        $inside++;
        my $arg = "in_fasta_phylo";
        my ($fastafile, $headRef, $minL, $maxL) = fastaload($in_fasta_phylo, $filenameNew, $upload_dir, $arg);
        if ($fastafile eq 'SIZE') {
            $newhtml = "There was an error uploading your DNA sequences. The file must be smaller than 50 Mb.<br>";
        } elsif ($fastafile eq 'NAME') {
            $newhtml = "There was an error uploading your DNA sequences. The file name contains forbidden characters. The only allowed characters are a-z, A-Z,  ., - and  _.<br>";
        } elsif ($fastafile eq 'CONTENT') {
            $newhtml = "There was an error uploading your DNA sequences. All the DNA sequences must be in fasta format. Please check the correct format <a href='../phyloprimer/format_page.html#DNA_seq'>here</a>.<br>";
            #`rm $upload_dir/$filename`;
        } elsif ($fastafile eq 'UNIQUE') {
            $newhtml = "There was an error uploading your DNA sequences. Some of the sequence headings are not unique.<br>";
        } else {
            if ($minL == $maxL) {
                $fastaLen = $minL;
            } else {
                $fastaLen = $minL . " - " . $maxL;
            }
            if ($maxL <= 6000) {
                if ($fastafile == 1) {
                    $newhtml = "You uploaded $fastafile DNA sequence. The length of the sequence is $fastaLen bp.<br>";
                } elsif ($fastafile <= 10) {
                    $newhtml = "You uploaded $fastafile DNA sequences. The length of the sequences is $fastaLen bp<br>";
                    
                    #mafft alignment to check gene similarity
                    `mafft --localpair --maxiterate 1000 --thread 30 --quiet ${upload_dir}/${filenameNew}  > ${upload_dir}/${filenameNew}.check`;
                    #consensus check
                    my $geneCheckPerc = `geneCheck_pp.pl -file ${upload_dir}/${filenameNew}.check`;
                    `rm ${upload_dir}/${filenameNew}.check`; #remove check file
                    
                    #if gap value too high: give warning
                    my @geneCheck = split(/-/, $geneCheckPerc);
                    if ($geneCheck[0] > 50) {
                        $newhtml .= "<b style='color:blue'\;>WARNING:</b> We checked your genes and looks like one or more sequences do not represent the same gene portion as the others.<br>";
                    }
                    
                } else {
                    $newhtml = "There was an error uploading your DNA sequences. The maximum number of allowed sequences is 10.<br>";
                }
            } else {
                $newhtml = "There was an error uploading your DNA sequences. The maximum allowed length is 6000 bp.<br>";
            }
        }
    }
    if ($inside > 1) { #if the user submitted both sequences and file
        $newhtml = "There was an error uploading your DNA sequences. Please input your DNA sequences either pasting them in the text area input field or uploading a fasta file.<br>";
    } elsif ($inside == 0) {
        $newhtml = "Please input your DNA sequences.<br>";
    }
    
} elsif ($input_tree eq "sequence") { #upload of DNA sequences
    my $folder = $folderOld . "Tb";
    $pass = $folderOld . "b";
    my $upload_dir = $path_html . "/uploadsPhyloprimer/" . $folder . "/"; #path to the folder
    `mkdir $upload_dir`; #create new directory
    my @headFasta = ();
    #fasta file
    if ($in_al_phylo ne "") {
        my $arg = "in_al_phylo";
        my $filenameNew = (substr($folder, 0, 8)) . ".fasta";
        my ($fastafile, $headRef, $minL, $maxL) =  fastaload($in_al_phylo, $filenameNew, $upload_dir,$arg);
        @headFasta = @{$headRef};
        
        foreach my $h (@headFasta) {
            $fasta{$h} = "";
        }
        if ($fastafile eq 'SIZE') {
            $newhtml = "There was an error uploading your DNA sequences. The file must be smaller than 50 Mb.<br>";
        } elsif ($fastafile eq 'NAME') {
            $newhtml = "There was an error uploading your DNA sequences. The file name contains forbiden characters. The only allowed characters are a-z, A-Z,  ., - and  _.<br>";
        } elsif ($fastafile eq 'CONTENT') {
            $newhtml = "There was an error uploading your DNA sequences. All the DNA sequences must be in fasta format. Please check the correct format <a href='../phyloprimer/format_page.html#DNA_seq'>here</a>.<br>";
        } elsif ($fastafile eq 'UNIQUE') {
            $newhtml = "There was an error uploading your DNA sequences. Some of the sequence headings are not unique.<br>";
        } else {
            if ($minL == $maxL) {
                $fastaLen = $minL;
            } else {
                $fastaLen = $minL . " - " . $maxL;
            }
            if ($fastafile == 1) { #at least 4 sequences
                $newhtml =  "There was an error uploading your DNA sequences. You uploaded only $fastafile sequence and the minimum requirement is 4 sequences.<br>";
            } elsif ($fastafile < 3) { #at least 4 sequences
                $newhtml =  "There was an error uploading you fasta file. You uploaded only $fastafile sequence and the minimum requirement is 4 sequences.<br>";
            } elsif ($fastafile > 500) { #at least 4 sequences
                $newhtml =  "There was an error uploading your DNA sequences. The maximum number of allowed sequences is 500.<br>";
            } else {
                $newhtml =  "You uploaded $fastafile sequences. The DNA sequence length is $fastaLen bp.<br>";
            }
        }
        #taxonomy file - check only if sequence file!!!
        if ($in_al_phylo_tax ne "") {
            my $arg = "in_al_phylo_tax";
            my $filenameNew = (substr($folder, 0, 8)) . ".tax.txt";
            my $fastaRef = \@headFasta;
            my ($datafile, $headRef) = dataload($in_al_phylo, $filenameNew, $upload_dir, $arg, $fastaRef); ##nedd to pass also array.
            if ($datafile eq 'SIZE') {
                $newhtml = $newhtml . "There was an error uploading your taxonomy reference file. The file must be smaller than 50 Mb.<br>";
            } elsif ($datafile eq 'NAME') {
                $newhtml = $newhtml . " There was an error uploading your taxonomy reference file. The file name contains forbiden characters. The only allowed characters are a-z, A-Z,  ., - and  _.<br>";
            } elsif ($datafile eq 'CONTENT1') {
                $newhtml = $newhtml . "There was an error uploading your taxonomy reference file. Some of the accession numbers present in the taxonomic reference file are not present in fasta file.<br>";
            } elsif ($datafile eq 'CONTENT2') {
                $newhtml = $newhtml . "There was an error uploading your taxonomy reference file. The taxonomy information must be reported with a taxid number. Please check the correct format <a href='../phyloprimer/format_page.html#tab_file'>here</a>.<br>";
            } else {
                my @headData = @{$headRef};
                my %tax;
                foreach my $h (@headData) {
                    $tax{$h} = '';
                }
                #Check if all the accessions in the fasta file are in taxonomical file
                my @excluded = ();
                foreach my $f (keys %fasta) {
                    if (!(defined($tax{$f}))) { #some of the accession in the taxonomy file are not present in the fasta file
                        push @excluded, $f;
                    }
                }
                if (defined($excluded[1])) {
                    my $excl = join ', ', @excluded;
                    $newhtml = $newhtml . "Your taxonomy reference file was successfully submitted. The accessions $excl do not have a taxonomical correspondence.<br>";
                } elsif (defined($excluded[0])) {
                    $newhtml = $newhtml . "Your taxonomy reference file was successfully submitted. The accession $excluded[0] does not have a taxonomical correspondence.<br>";
                } else {
                    $newhtml = $newhtml . "Your taxonomy reference file was successfully submitted. All the accessions in the fasta file are present in the taxonomy reference file.<br>";
                }
            }
        } else {
            $newhtml = $newhtml . "No reference taxonomy file was uploaded so there will not be any taxonomy information in the guide tree.<br>";
        }
        #protein file - check only if sequence file!!!
        if ($in_al_phylo_pro ne "") {
            my $arg = "in_al_phylo_pro";
            my $filenameNew = (substr($folder, 0, 8)) . ".pro.txt";
            my $fastaRef = \@headFasta;
            my ($datafile, $headRef) = dataload($in_al_phylo, $filenameNew, $upload_dir, $arg, $fastaRef); #need to pass also array.
            if ($datafile eq 'SIZE') {
                $newhtml = $newhtml . "There was an error uploading your protein-gene reference file. The file must be smaller than 50 Mb.<br>";
            } elsif ($datafile eq 'NAME') {
                $newhtml = $newhtml . "There was an error uploading your protein-gene reference file. The file name contains forbiden characters. The only allowed characters are a-z, A-Z,  ., - and  _.<br>";
            } elsif ($datafile eq 'CONTENT1') {
                $newhtml = $newhtml . "There was an error uploading your protein-gene reference file. Some of the accession numbers present in the taxonomic reference file are not present in fasta file.<br>";
            } elsif ($datafile eq 'CONTENT2') {
                $newhtml = $newhtml . "There was an error uploading your protein-gene reference file. Please check the correct format <a href='../phyloprimer/format_page.html#tab_file'>here</a>.";
            } else {
                my @headData = @{$headRef};
                my %pro;
                foreach my $h (@headData) {
                    $pro{$h} = '';
                }
                #Check if all the accessions in the fasta file are in taxonomical file
                my @excluded = ();
                foreach my $f (keys %fasta) {
                    if (!(defined($pro{$f}))) { #some of the accession in the taxonomy file are not present in the fasta file
                        push @excluded, $f;
                    }
                }
                if (defined($excluded[1])) {
                    my $excl = join ', ', @excluded;
                    $newhtml = $newhtml . "Your protein-gene reference file was successfully submitted. The accessions $excl do not have a protein-gene correspondence.<br>";
                } elsif (defined($excluded[0])) {
                    $newhtml = $newhtml . "Your protein-gene reference file was successfully submitted. The accession $excluded[0] does not have a protein-gene correspondence.<br>";
                } else {
                    $newhtml = $newhtml . "Your protein-gene reference file was successfully submitted. All the accessions in the fasta file are present in the protein-gene reference file.<br>";
                }
            }
        } else {
            $newhtml = $newhtml . "No reference protein-gene reference file was uploaded so there will not be any protein-gene information in the guide tree.<br>";
        }
    } else {
        $newhtml = "Please input your DNA sequences.<br>";
    }
} elsif ($input_tree eq "newick") { #upload of newick tree
    my $folder = $folderOld . "Tc";
    $pass = $folderOld . "c";
    my $upload_dir = $path_html . "/uploadsPhyloprimer/" . $folder . "/"; #path to the folder
    `mkdir $upload_dir`; #create new directory
    #fasta file
    if ($in_newick_phylo_seq ne "") {
        my $arg = "in_newick_phylo_seq";
        my $filenameNew = (substr($folder, 0, 8)) . ".alignment";
        my ($fastafile, $headRef, $minL, $maxL) =  fastaload($in_newick_phylo_seq, $filenameNew, $upload_dir,$arg);
        my @headFasta = @{$headRef};
        foreach my $h (@headFasta) {
            $fasta{$h} = "";
        }
        if ($fastafile eq 'SIZE') {
            $newhtml = "There was an error uploading your alignment sequences. The file must be smaller than 50 Mb.<br>";
        } elsif ($fastafile eq 'NAME') {
            $newhtml = "There was an error uploading your alignment sequences. The file name contains forbiden characters. The only allowed characters are a-z, A-Z,  ., - and  _.<br>";
        } elsif ($fastafile eq 'LEN') {
            $newhtml = "There was an error uploading your alignment sequences. All the alignment sequences should be the same length.<br>";
        } elsif ($fastafile eq 'CONTENT') {
            $newhtml = "There was an error uploading your alignment sequences. All the alignment sequences must be in fasta format. Please check the correct format <a href='../phyloprimer/format_page.html#alignment'>here</a>.<br>";
        } elsif ($fastafile eq 'UNIQUE') {
            $newhtml = "There was an error uploading your alignment sequences. Some of the sequence headings were not unique.<br>";
        } else {
            if ($minL == $maxL) {
                $fastaLen = $minL;
            } else {
                $fastaLen = $minL . " - " . $maxL;
            }
            $newhtml = "You uploaded $fastafile DNA sequences. The length range is $fastaLen bp.<br>";
        }
    } else {
        $newhtml = "Please input your alignment sequences.<br>";
    }
    #newick tree
    if ($in_newick_phylo ne "") {
        my $arg = "in_newick_phylo";
        my $filenameNew = (substr($folder, 0, 8)) . ".tree";
        my ($newickfile, $headRef)=  newickload($in_newick_phylo, $filenameNew, $upload_dir,$arg);
        my @headNew = @{$headRef};
        if ($newickfile eq 'SIZE') {
            $newhtml = "There was an error uploading your newick file. The file must be smaller than 50 Mb.<br>" . $newhtml;
        } elsif ($newickfile eq 'NAME') { #I may remove it
            $newhtml = "There was an error uploading your newick file. The file name contains forbiden characters. The only allowed characters are a-z, A-Z,  ., - and  _.<br>" . $newhtml;
        } elsif ($newickfile eq 'CONTENT') {
            $newhtml = "There was an error uploading your newick file. The tree must be in newick format. Please check the correct format <a href='../phyloprimer/format_page.html#newick'>here</a>.<br>" . $newhtml;
        } elsif ($newickfile eq 'UNIQUE') {
            $newhtml = "There was an error uploading your newick file. Some of the sequence headings were not unique.<br>";
        } else {
            #Check there is correspondence between the two files
            my @excluded = ();
            foreach my $h (@headNew) {
                if (!(defined($fasta{$h}))) { #some of the accession in the alignment are not present in the fasta file
                    push @excluded, $h;
                }
            }
            if (defined($excluded[1])) {
                my $excl = join ', ', @excluded;
                $newhtml = "There was an error uploading your newick file. The entries $excl found in the newick file are missing in the fasta file.<br>" . $newhtml;
            } elsif (defined($excluded[0])) {
                $newhtml = "There was an error uploading your newick file. The entry $excluded[0] found in the newick file is missing in the fasta file.<br>" . $newhtml;
            } else {
                if ($newickfile > 500) {
                    $newhtml = "There was an error uploading your newick file. The maximum number of allowed entries is 500.<br>" . $newhtml;
                } else {
                    $newhtml = "Your newick tree was successfully uploaded. The newick tree has $newickfile entries.<br>" . $newhtml;
                }
            }
        }
        
        #taxonomy file - check only if sequence file!!!
        if ($in_newick_phylo_tax ne "") {
            my $arg = "in_newick_phylo_tax";
            my $filenameNew = (substr($folder, 0, 8)) . ".tax.txt";
            my $fastaRef = \@headNew;
            my ($datafile, $headRef) = dataload($in_newick_phylo_tax, $filenameNew, $upload_dir, $arg, $fastaRef); ##nedd to pass also array.
            if ($datafile eq 'SIZE') {
                $newhtml = $newhtml . "There was an error uploading your taxonomy reference file. The file must be smaller than 50 Mb.<br>";
            } elsif ($datafile eq 'NAME') {
                $newhtml = $newhtml . "There was an error uploading your taxonomy reference file. The file name contains forbiden characters. The only allowed characters are a-z, A-Z,  ., - and  _.<br>";
            } elsif ($datafile eq 'CONTENT1') {
                $newhtml = $newhtml . "There was an error uploading your taxonomy reference file. Some of the accession numbers present in the taxonomic reference file are not present in fasta file.<br>";
            } elsif ($datafile eq 'CONTENT2') {
                $newhtml = $newhtml . "There was an error uploading your taxonomy reference file. The taxonomy information must be reported with a taxid number. Please check the correct format <a href='../phyloprimer/format_page.html#tab_file'>here</a>.<br>";
            } else {
                my @headData = @{$headRef};
                my %tax;
                foreach my $h (@headData) {
                    $tax{$h} = '';
                }
                #Check if all the accessions in the fasta file are in taxonomical file
                my @excluded = ();
                foreach my $f (keys %fasta) {
                    if (!(defined($tax{$f}))) { #some of the accession in the taxonomy file are not present in the fasta file
                        push @excluded, $f;
                    }
                }
                if (defined($excluded[1])) {
                    my $excl = join ', ', @excluded;
                    $newhtml = $newhtml . "Your taxonomy reference file was successfully submitted. The accessions $excl do not have a taxonomical correspondence.<br>";
                } elsif (defined($excluded[0])) {
                    $newhtml = $newhtml . "Your taxonomy reference file was successfully submitted. The accession $excluded[0] does not have a taxonomical correspondence.<br>";
                } else {
                    $newhtml = $newhtml . "Your taxonomy reference file was successfully submitted. All the accessions in the fasta file are present in the taxonomy reference file.<br>";
                }
            }
        } else {
            $newhtml = $newhtml . "No reference taxonomy file was uploaded so there will not be any taxonomy information in the guide tree.<br>";
        }
        #protein file - check only if sequence file!!!
        if ($in_newick_phylo_pro ne "") {
            my $arg = "in_newick_phylo_pro";
            my $filenameNew = (substr($folder, 0, 8)) . ".pro.txt";
            my $fastaRef = \@headNew;
            my ($datafile, $headRef) = dataload($in_newick_phylo_pro, $filenameNew, $upload_dir, $arg, $fastaRef); ##need to pass also array.
            if ($datafile eq 'SIZE') {
                $newhtml = $newhtml . "There was an error uploading your protein-gene reference file. The file must be smaller than 50 Mb.<br>";
            } elsif ($datafile eq 'NAME') {
                $newhtml = $newhtml . "There was an error uploading your protein-gene reference file. The file name contains forbiden characters. The only allowed characters are a-z, A-Z,  ., - and  _.<br>";
            } elsif ($datafile eq 'CONTENT1') {
                $newhtml = $newhtml . "There was an error uploading your protein-gene reference file. Some of the accession numbers present in the protein-gene reference file are not present in fasta file.<br>";
            } elsif ($datafile eq 'CONTENT2') {
                $newhtml = $newhtml . "There was an error uploading your protein-gene reference file. Please check the correct format <a href='../phyloprimer/format_page.html#tab_file'>here</a>.";
            } else {
                my @headData = @{$headRef};
                my %pro;
                foreach my $h (@headData) {
                    $pro{$h} = '';
                }
                #Check if all the accessions in the fasta file are in taxonomical file
                my @excluded = ();
                foreach my $f (keys %fasta) {
                    if (!(defined($pro{$f}))) { #some of the accession in the taxonomy file are not present in the fasta file
                        push @excluded, $f;
                    }
                }
                if (defined($excluded[1])) {
                    my $excl = join ', ', @excluded;
                    $newhtml = $newhtml . "Your protein-gene reference file was successfully submitted. The accessions $excl do not have a protein-gene correspondence.<br>";
                } elsif (defined($excluded[0])) {
                    $newhtml = $newhtml . "Your protein-gene reference file was successfully submitted. The accession $excluded[0] does not have a protein-gene correspondence.<br>";
                } else {
                    $newhtml = $newhtml . "Your protein-gene reference file was successfully submitted. All the accessions in the fasta file are present in the protein-gene reference file.<br>";
                }
            }
        } else {
            $newhtml = $newhtml . "No reference protein-gene reference file was uploaded so there will not be any protein-gene information in the guide tree.<br>";
        }
        
    } else {
        $newhtml = "Please input your newick tree.<br>" . $newhtml;
    }
} else {
    $newhtml = "Something went wrong. Please reload the page.<br>";  #if the input_tree value was none of the expected
}

#check and load sequence
sub sequenceload {
    my ($content) = $_[0];
    my ($filenameNew) = $_[1];
    my ($upload_dir) = $_[2];
    my $fastaSeq;
    my @all = split(/\n/, $content);
    my $header = 0;
    my $count = 0;
    my $sequence = 0;
    my $badEntry = 0;
    my $last;
    my %print;
    my $seq;
    my $accession;
    my %all_len;
    my $in_count = 0;
    my $min;
    my $max;
    my %unique;
    foreach my $input (@all) {
        chomp($input);
        $input =~ s/^\s+|\s+$//g;
        if ($input =~ /^>/) {
            $header++; #how many headers
            $count = 0;
            $last = "header";
            ($accession) = $input =~ /^>(.*)/;
            #substitute accession special characters
            $accession =~ s/\'/\[sub1\]/g; #single quote
            $accession =~ s/,/\[sub2\]/g; #comma
            $accession =~ s/\(/\[sub3\]/g; #bracket (
            $accession =~ s/\)/\[sub4\]/g; #bracket )
            $accession =~ s/:/\[sub5\]/g; #column
            $accession =~ s/;/\[sub6\]/g; #semi column
            $accession =~ s/\*/\[sub7\]/g; #semi column
            $accession =~ s/</\[sub8\]/g; #lower
            $accession =~ s/>/\[sub9\]/g; #higher
            $accession =~ s/-/\[sub10\]/g; #minus
            $accession =~ s/\+/\[sub11\]/g; #plus
            $accession =~ s/\`/\[sub12\]/g; #hyphen`
            $accession =~ s/\#/\[sub13\]/g; #
            $accession =~ s/&/\[sub14\]/g; #&
            $accession =~ s/\^/\[sub15\]/g; #&
            $accession =~ s/\//\[sub16\]/g; #/
            $accession =~ s/_/\[sub17\]/g; #_
            $accession =~ s/\s/_/g;#all spaces to _
            if (defined($unique{$accession})) {
                $badEntry = 2;
            } else {
                $unique{$accession} = '';
            }
            $accession = "USER_INPUT" . $header . " -userHeader- " . $accession;
            if ($header > 1) {
                if (defined($print{$accession})) {
                    $accession = $accession . $header;
                    $print{$header}{$accession} = $seq; #accession = sequence
                } else {
                    $print{$header}{$accession} = $seq; #accession = sequence
                }
            }
            $seq = "";
        } else {
            if ($input ne "") {
                $last = "sequence";
                $count++;
                if ($input !~ /\A[acgtnryswkmbdhvn]+\z/i) { #safe characters for sequence
                    $badEntry = 1;
                    last;
                } else {
                    $seq .= $input;
                }
                if ($count == 1) { #how many sequences
                    $sequence++;
                }
            }
        }
        $print{$header}{$accession} = $seq; #accession = sequence
    }
    
    if (($last eq "sequence") && ($header == 0)) { #onlyone sequence and no header
        $header = 1;
        undef %print;
        $print{'1'}{'USER_INPUT1 -userHeader- '} = $seq;
    }
    
    if (($badEntry == 0) && ($last eq "sequence") && (($header == $sequence) or (($header == 0) && ($sequence == 1)))) {
        $fastaSeq = $header;
        open (UPLOADFILE, ">$upload_dir/$filenameNew" ) or die "$!"; #print only if entry is okay
        foreach my $c (sort {$a <=> $b} keys %print) {
            foreach my $head (keys %{$print{$c}}) {
                print UPLOADFILE ">$head\n$print{$c}{$head}\n";
                my $len = length($print{$c}{$head});
                $all_len{$len}++;
            }
        }
        close(UPLOADFILE);
        foreach my $l (sort {$a <=> $b} keys %all_len) {
            $in_count++;
            if ($in_count == 1) {
                $min = $l;
                $max = $l;
            } else {
                if ($l > $max) {
                    $max = $l;
                } elsif ($l < $min) {
                    $min = $l;
                }
            }
        }
    } elsif ($badEntry == 2) {
        $fastaSeq = "UNIQUE";
    } else {
        $fastaSeq = "CONTENT";
    }
    return ($fastaSeq, $min, $max);
}

#check and load file
sub fastaload {
    my ($filename) = $_[0];
    my ($filenameNew) = $_[1];
    my ($upload_dir) = $_[2];
    my ($arg) = $_[3];
    
    my $fastafile;
    my @head;
    
    my %all_len;
    my $in_count = 0;
    my $min;
    my $max;
    my %unique;
    
    if (!(defined($filename))) { #if the file is too big
        $fastafile = "SIZE";
    } else { #if the file size is lower than the threshold
        
        my ( $name, $path, $extension ) = fileparse ( $filename, '..*' );
        $filename = $name . $extension;
        
        if ( $filename =~ /^([a-zA-Z0-9_\.\-]+)$/ ) { #safe characters for file name
            my $upload_filehandle = $cgi->upload($arg);
            my $header = 0;
            my $count = 0;
            my $sequence = 0;
            my $badEntry = 0;
            my %all_seq;
            my $last;
            my %print;
            my $seq;
            my $accession;
            
            while (defined(my $input = <$upload_filehandle>)) { #read the file
                chomp($input);
                $input =~ s/\r//g;
                if ($input =~ /^>/) {
                    $header++; #how many headers
                    $count = 0;
                    $last = "header";
                    ($accession) = $input =~ /^>(.*)/;
                    #substitute accession special characters
                    $accession =~ s/^\s+|\s+$//g; #taxon name
                    $accession =~ s/\'/\[sub1\]/g; #single quote
                    $accession =~ s/,/\[sub2\]/g; #comma
                    $accession =~ s/\(/\[sub3\]/g; #bracket (
                    $accession =~ s/\)/\[sub4\]/g; #bracket )
                    $accession =~ s/:/\[sub5\]/g; #column
                    $accession =~ s/;/\[sub6\]/g; #semi column
                    $accession =~ s/\*/\[sub7\]/g; #semi column
                    $accession =~ s/</\[sub8\]/g; #lower
                    $accession =~ s/>/\[sub9\]/g; #higher
                    $accession =~ s/-/\[sub10\]/g; #minus
                    $accession =~ s/\+/\[sub11\]/g; #plus
                    $accession =~ s/\`/\[sub12\]/g; #hyphen`
                    $accession =~ s/\#/\[sub13\]/g; #
                    $accession =~ s/&/\[sub14\]/g; #&
                    $accession =~ s/\^/\[sub15\]/g; #&
                    $accession =~ s/\//\[sub16\]/g; #/
                    $accession =~ s/_/\[sub17\]/g; #_
                    $accession =~ s/\s/_/g;#all spaces to _
                    if (defined($unique{$accession})) {
                        $badEntry = 2;
                    } else {
                        $unique{$accession} = '';
                    }
                    if ($arg eq 'in_fasta_phylo') { #if gene - modify user sequences in any case (not only if > 20 characters)
                        push @head, $accession; #save all accession
                        $accession = "USER_INPUT" . $header . " -userHeader- " . $accession;
                    } else { #if others
                        if (length($accession) > 20) {
                            push @head, $accession; #save all accession
                            my $accession0 = $accession;
                            $accession = "USER_INPUT" . $header . " -userHeader- " . $accession;
                            ###########REPLACE
                            $replace{$accession0} = $accession;
                        } else {
                            push @head, $accession; #save all accession
                        }
                    }
                    $seq = "";
                } else {
                    if ($input ne "") {
                        $last = "sequence";
                        $count++;
                        if ((($input !~ /\A[acgtnryswkmbdhvn]+\z/i) && ($arg ne "in_newick_phylo_seq")) or (($input !~ /\A[acgtnryswkmbdhvn\.\-]+\z/i) && ($arg eq "in_newick_phylo_seq"))) { #safe characters for sequence / alignments
                            $badEntry = 1;
                            last;
                        } else {
                            $seq .= $input;
                        }
                        if ($count == 1) { #how many sequences
                            $sequence++;
                        }
                    }
                }
                
                $print{$header}{$accession} = $seq; #accession = sequence
            }
            
            if (($last eq "sequence") && ($header == 0)) { #onlyone sequence and no header
                $header = 1;
                undef %print;
                $print{'1'}{'USER_INPUT1 -userHeader- '} = $seq;
            }
            
            #check if all sequence alignment are long the same
            if ($arg eq 'in_newick_phylo_seq') {
                my $c = 0;
                my $lenAlOld;
                foreach my $head (keys %print) {
                    $c++;
                    if ($c == 1) {
                        $lenAlOld = length($print{$head});
                    } else {
                        my $lenAl = length($print{$head});
                        if ($lenAl != $lenAlOld) {
                            $badEntry = 3;
                            last;
                        }
                    }
                }
            }
            
            #tolerate if no header so only one sequence
            #BUT if more than one header, same number between headers and sequences
            if (($badEntry == 0) && ($last eq "sequence") && (($header == $sequence) or (($header == 0) && ($sequence == 1)))) {
                $fastafile = $header;
                open (UPLOADFILE, ">$upload_dir/$filenameNew" ) or die "$!"; #print only if entry is okay
                binmode UPLOADFILE;
                foreach my $c (sort {$a <=> $b} keys %print) {
                    foreach my $head (keys %{$print{$c}}) {
                        print UPLOADFILE ">$head\n$print{$c}{$head}\n";
                        my $len = length($print{$c}{$head});
                        $all_len{$len}++;
                    }
                }
                close(UPLOADFILE);
                foreach my $l (sort {$a <=> $b} keys %all_len) {
                    $in_count++;
                    if ($in_count == 1) {
                        $min = $l;
                        $max = $l;
                    } else {
                        if ($l > $max) {
                            $max = $l;
                        } elsif ($l < $min) {
                            $min = $l;
                        }
                    }
                }
                
            } elsif ($badEntry == 2) {
                $fastafile = "UNIQUE";
            } elsif ($badEntry == 3) {
                $fastafile = "LEN"; #not all alignment have same length
            } else {
                $fastafile = "CONTENT";
            }
        } else { #uncorrect file name
            $fastafile = "NAME";
        }
    }
    my $headRef = \@head;
    return ($fastafile,$headRef,$min,$max);
}

#newick load
sub newickload {
    #subrountine inputs
    my ($filename) = $_[0];
    my ($filenameNew) = $_[1];
    my ($upload_dir) = $_[2];
    my ($arg) = $_[3];
    
    my $newickfile;
    my $count = 0;
    my @head; #
    my $newick; #
    my $badEntry = 0; #
    my $header = 0;
    my $acc1;
    
    my %unique;
    
    
    if (!(defined($filename))) { #if the file is too big
        $newickfile = "SIZE";
    } else { #if the file size is lower than the threshold
        my ($name, $path, $extension) = fileparse ( $filename, '..*' );
        $filename = $name . $extension;
        
        if ($filename =~ /^([a-zA-Z0-9_\.\-]+)$/) { #safe characters for file name
            my $upload_filehandle = $cgi->upload($arg);
            while (defined(my $input = <$upload_filehandle>)) { #read the file
                chomp($input);
                $input =~ s/^\s+|\s+$//g;
                $newick .= $input;
                
            }
            my @all = split(/:/,$newick); #divide newick tree where :
            my $len = scalar @all; #number of lines
            ##check newick tree
            foreach my $all (@all) {
                $count++;
                #quick check that file starts with right characters
                if ($count == 1) {
                    if ($all =~ /^\(/) {
                        my ($acc) = $all =~ /[\(]+(.*)$/;
                        #substitute accession special characters
                        $acc =~ s/^\s+|\s+$//g; #taxon name
                        $acc =~ s/\'/\[sub1\]/g; #single quote
                        $acc =~ s/,/\[sub2\]/g; #comma
                        $acc =~ s/\(/\[sub3\]/g; #bracket (
                        $acc =~ s/\)/\[sub4\]/g; #bracket )
                        $acc =~ s/:/\[sub5\]/g; #column
                        $acc =~ s/;/\[sub6\]/g; #semi column
                        $acc =~ s/\*/\[sub7\]/g; #semi column
                        $acc =~ s/</\[sub8\]/g; #lower
                        $acc =~ s/>/\[sub9\]/g; #higher
                        $acc =~ s/-/\[sub10\]/g; #minus
                        $acc =~ s/\+/\[sub11\]/g; #plus
                        $acc =~ s/\`/\[sub12\]/g; #hyphen`
                        $acc =~ s/\#/\[sub13\]/g; #
                        $acc =~ s/&/\[sub14\]/g; #&
                        $acc =~ s/\^/\[sub15\]/g; #&
                        $acc =~ s/\//\[sub16\]/g; #/
                        $acc =~ s/_/\[sub17\]/g; #_
                        $acc =~ s/\s/_/g;#all spaces to _
                        if (defined($unique{$acc})) {
                            $badEntry = 2;
                        } else {
                            $unique{$acc} = '';
                        }
                        push @head, $acc; #save all accession - accession are the originals
                        if (length($acc) > 20) {
                            $header++;
                            $acc1 = "USER_INPUT" . $header . " -userHeader- " . $acc;
                            ###########REPLACE
                            $replace{$acc} = $acc1;
                        }
                        $newickfile++;
                    } else {
                        $badEntry = 1; #if tree does not start with ( -> something wrong
                    }
                } elsif ($count == $len) { #end tree
                    
                    
                    if (($all !~ /^([0-9\.:]+)\)$/) && ($all !~ /^([0-9\.:]+)\);$/)) {
                        $badEntry = 1; #if tree does not finish with ); -> something wrong
                    } else {
                        if ($all !~ /;$/) {
                            $newick .= ";";
                        }
                    }
                } else { #middle of the tree
                    if ($all =~ /^([0-9\.]+),/) { #when the branch length is folled by , -> accession number
                        if ($all =~ /\(/) { #,((
                            my ($acc) = $all =~ /[\(]+(.*)$/;
                            #substitute accession special characters
                            $acc =~ s/^\s+|\s+$//g; #taxon name
                            $acc =~ s/\'/\[sub1\]/g; #single quote
                            $acc =~ s/,/\[sub2\]/g; #comma
                            $acc =~ s/\(/\[sub3\]/g; #bracket (
                            $acc =~ s/\)/\[sub4\]/g; #bracket )
                            $acc =~ s/:/\[sub5\]/g; #column
                            $acc =~ s/;/\[sub6\]/g; #semi column
                            $acc =~ s/\*/\[sub7\]/g; #semi column
                            $acc =~ s/</\[sub8\]/g; #lower
                            $acc =~ s/>/\[sub9\]/g; #higher
                            $acc =~ s/-/\[sub10\]/g; #minus
                            $acc =~ s/\+/\[sub11\]/g; #plus
                            $acc =~ s/\`/\[sub12\]/g; #hyphen`
                            $acc =~ s/\#/\[sub13\]/g; #
                            $acc =~ s/&/\[sub14\]/g; #&
                            $acc =~ s/\^/\[sub15\]/g; #&
                            $acc =~ s/\//\[sub16\]/g; #/
                            $acc =~ s/_/\[sub17\]/g; #_
                            $acc =~ s/\s/_/g;#all spaces to _
                            if (defined($unique{$acc})) {
                                $badEntry = 2;
                            } else {
                                $unique{$acc} = '';
                            }
                            push @head, $acc; #save all accession
                            if (length($acc) > 20) {
                                $header++;
                                $acc1 = "USER_INPUT" . $header . " -userHeader- " . $acc;
                                ###########REPLACE
                                $replace{$acc} = $acc1;
                            }
                            $newickfile++;
                        } else { #,
                            my ($acc) = $all =~ /,(.*)/;
                            #substitute accession special characters
                            $acc =~ s/^\s+|\s+$//g; #taxon name
                            $acc =~ s/\'/\[sub1\]/g; #single quote
                            $acc =~ s/,/\[sub2\]/g; #comma
                            $acc =~ s/\(/\[sub3\]/g; #bracket (
                            $acc =~ s/\)/\[sub4\]/g; #bracket )
                            $acc =~ s/:/\[sub5\]/g; #column
                            $acc =~ s/;/\[sub6\]/g; #semi column
                            $acc =~ s/\*/\[sub7\]/g; #semi column
                            $acc =~ s/</\[sub8\]/g; #lower
                            $acc =~ s/>/\[sub9\]/g; #higher
                            $acc =~ s/-/\[sub10\]/g; #minus
                            $acc =~ s/\+/\[sub11\]/g; #plus
                            $acc =~ s/\`/\[sub12\]/g; #hyphen`
                            $acc =~ s/\#/\[sub13\]/g; #
                            $acc =~ s/&/\[sub14\]/g; #&
                            $acc =~ s/\^/\[sub15\]/g; #&
                            $acc =~ s/\//\[sub16\]/g; #/
                            $acc =~ s/_/\[sub17\]/g; #_
                            $acc =~ s/\s/_/g;#all spaces to _
                            if (defined($unique{$acc})) {
                                $badEntry = 2;
                            } else {
                                $unique{$acc} = '';
                            }
                            push @head, $acc; #save all accession
                            if (length($acc) > 20) {
                                $header++;
                                $acc1 = "USER_INPUT" . $header . " -userHeader- " . $acc;
                                ###########REPLACE
                                $replace{$acc} = $acc1;
                            }
                            $newickfile++;
                        }
                    } elsif ($all =~ /^([0-9\.]+)\)/) { #when branch length is followed by ) -> bootstrap or end of the line
                        if (($all !~ /\)$/) && ($all !~ /\)[0-9]+$/)) {
                            $badEntry = 1; #if after ) it is not the end of the line or bootstrap
                        }
                    } else {
                        $badEntry = 1; #if the branch length is followed by other characters other then , and )
                    }
                }
            }
            if ($badEntry == 0) {
                foreach my $r (keys %replace) {
                    $newick =~ s/$r/$replace{$r}/g;
                }
                open (UPLOADFILE, ">$upload_dir/$filenameNew" ) or die "$!";
                binmode UPLOADFILE;
                print UPLOADFILE "$newick";
                close(UPLOADFILE);
            } elsif ($badEntry == 2) {
                $newickfile = "UNIQUE";
            } else {
                $newickfile = "CONTENT";
            }
        } else { #uncorrect file name
            $newickfile = "NAME";
        }
    }
    my $headRef = \@head;
    return ($newickfile,$headRef);
}



#check and load data file (taxonomy and protein referenrce file)
sub dataload {
    my ($filename) = $_[0];
    my ($filenameNew) = $_[1];
    my ($upload_dir) = $_[2];
    my ($arg) = $_[3];
    my @fastaAcc = @{$_[4]}; #array with fasta accession numbers
    my %fastaAcc = ();
    
    foreach my $a (@fastaAcc) {
        $fastaAcc{$a} = '';
    }
    
    my $fastafile;
    my @head;
    
    my %all_len;
    my $in_count = 0;
    my $min;
    my $max;
    
    if (!(defined($filename))) { #if the file is too big
        $fastafile = "SIZE";
    } else { #if the file size is lower than the threshold
        
        my ( $name, $path, $extension ) = fileparse ( $filename, '..*' );
        $filename = $name . $extension;
        
        if ( $filename =~ /^([a-zA-Z0-9_\.\-]+)$/) { #safe characters for file name
            my $upload_filehandle = $cgi->upload($arg);
            my $badEntry = 0;
            my %print;
            while (defined(my $input = <$upload_filehandle>)) { #read the file
                chomp($input);
                my @info = split(/\t/, $input);
                if ($info[1] ne '') {
                    if ((($arg =~ /_tax/) && ($info[1] =~ /^([0-9]+)$/)) or ($arg =~ /_pro/)) {
                        
                        $info[0] =~ tr/A-Za-z0-9.-_/_/sc; #substitute all the characters not included in A-Z a-z 0-9 . - _ ( ) with _
                        $info[1] =~ tr/A-Za-z0-9.-_/_/sc; #substitute all the characters not included in A-Z a-z 0-9 . - _ ( ) with _
                        if (defined($fastaAcc{$info[0]})) {
                            push @head, $info[0]; #save all accession
                            $print{$info[0]} = $info[1];
                        } else {
                            $badEntry = 1;
                            $fastafile = "CONTENT1";
                        }
                    } else {
                        $badEntry = 1;
                        $fastafile = "CONTENT2";
                    }
                } else {
                    $badEntry = 1;
                    $fastafile = "CONTENT2";
                }
            }
            if ($badEntry == 0) {
                open (UPLOADFILE, ">$upload_dir/$filenameNew" ) or die "$!"; #print only if the file does not have errors
                binmode UPLOADFILE;
                foreach my $acc (keys %print) {
                    if (defined($replace{$acc})) {
                        my ($acc1) = $replace{$acc} =~ /(.*) -userHeader-/;
                        print UPLOADFILE "$acc1\t$print{$acc}\n";
                    } else {
                        print UPLOADFILE "$acc\t$print{$acc}\n";
                    }
                }
                close(UPLOADFILE);
            }
        } else { #uncorrect file name
            $fastafile = "NAME";
        }
    }
    my $headRef = \@head;
    return ($fastafile,$headRef);
}

#send data back to html

if (($newhtml =~ /error/) or ($newhtml =~/Please input/)) {
    $newhtml = "<h3>" . $newhtml . "<br><b>There was an error in your uploads. Please check and resubmit again.</b></h3>";
} else {
    if ($input_tree eq "gene") {
        $newhtml = "<h3>" . $newhtml . "<br><b>The uploaded data is correct. If you are happy with the upload, please go ahead and select the BLAST parameters you want Phyloprimer to use during the BLAST search against the modified GenBank database. We suggest not to set the identity and the coverage percentage too low because the BLAST search may retrieve sequences too distant related from your query and the output tree may be too populated for having an easy visualization. In case you want to modify your inputs, please restart from the Inputs section and then press Check Uploads button again. Once that you are ready, press the Create Phylogenetic Tree button. The tree construction may take few minutes.</b></h3>";
    } else {
        $newhtml = "<h3>" . $newhtml . "<br><b>The uploaded data is correct. If you are happy with the upload, go ahead and press the Create Phylogenetic Tree button. In case you want to modify your inputs, please restart from the Inputs section and then press Check Uploads button again. Once that you are ready, press the Create Phylogenetic Tree button. The tree construction may take few minutes.</b></h3>";
    }
    
    
}

my @result;
push @result, "result";
push @result, $newhtml;

push @result, "set";
push @result, $pass;

my $json = $op -> encode({
    @result
});
print $json;
