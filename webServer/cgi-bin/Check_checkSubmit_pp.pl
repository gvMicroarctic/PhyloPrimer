#!/usr/bin/perl

use JSON; #if not already installed, just run "cpan JSON"
use strict;
use warnings;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser );
use File::Basename;

#connected to check_input.html : check parameters set in Oligo check page

#set path to folders:
my $path_html = 'path_to_html_folder';

#Set maximum size of uploaded file
$CGI::POST_MAX = 1024 * 50000; #50 MegaBytes

#header needs to start with > and all the white spaces are changed with _
#folder name is passed only if the file is correct!

my $newhtml;
my $folderOld;
my $folder;
my $upload_dir;
my $infoFile;
my $pass;
my %fasta;
my $fastaLen = 0;
my %replace;

#open cgi connection
my $cgi = CGI->new;
print $cgi->header('application/json;charset=UTF-8');

#Section 0
my $project = $cgi->param('project'); #word
my $email = $cgi->param('email'); #word - needs to have @
my $mailing_list = $cgi->param('mailing_list'); #user subscribed to mailing list

#Section 1
my $use = $cgi->param('use'); #how to use the script
my $defSet = $cgi->param('defSet'); #how to use the script

#Section 2
my $in_sequence = $cgi->param("in_sequence"); #oligo sequence
my $in_oligo = $cgi->param("in_oligo"); #oligo file

#Section 3
my $in_negative = $cgi->param("in_negative"); #negative file
my $spec_mis = $cgi->param('spec_mis'); #positive integer
my $spec3_mis = $cgi->param('spec3_mis'); #positive integer
my $spec3_length = $cgi->param('spec3_length'); #positive integer

#Section 4
my $mon_dG = $cgi->param('mon_dG'); #n
my $mg_dG = $cgi->param('mg_dG'); #n
my $oligo_dG = $cgi->param('oligo_dG'); #n
my $dNTP_dG = $cgi->param('dNTP_dG'); #n
my $t_dG = $cgi->param('t_dG'); #n

#convert  data to JSON
my $op = JSON -> new -> utf8 -> pretty(1);

if ($defSet eq "new") {

    ##Section 0
    $newhtml = $newhtml . "<h4><b>USER DATA</b></h4>";

    #project
    if ($project eq '') {
       $newhtml = $newhtml . "<li><b style='color:red'\;>Job ID:</b> please check. Insert a Job ID.</li><br>";
    } else {
        if ($project  =~ /\A[0-9A-Za-z_\-\. ]+\z/) {
            $newhtml = $newhtml . "<li>Job ID: $project.</li><br>";
            $project =~ tr/ /_/;
        } else {
            $newhtml = $newhtml . "<li><b style='color:red'\;>Job ID:</b> please check. No special characters allowed.</li><br>";
        }
    }

    #email
    if ($email =~ /\@/) {
        $newhtml = $newhtml . "<li>Email: $email.</li>";
    } else {
        $newhtml = $newhtml . "<li><b style='color:red'\;>Email:</b> please check. We need a valid email address to send you the results.</li>";
    }

    $mailing_list =~ tr/ /_/;
    if ($mailing_list eq '') {
        $mailing_list = "no";
    } else {
        $mailing_list = "yes";
    }

    $newhtml = $newhtml . "<br>";

    ##Section 2
    $newhtml .= "<h4><b>INPUTS</b></h4>";
    my @chars = ("A".."Z", "a".."z"); #create folder with 20 random characters (N.B. it needs to be sent back to html)
    $folderOld .= $chars[rand @chars] for 1..20;
    $folder = $folderOld . "Ca";
    $pass = $folderOld . "a";
    $upload_dir = $path_html . "/uploadsPhyloprimer/" . $folder . "/"; #path to the folder
    `mkdir $upload_dir`; #create new directory
    my $filenameNew = (substr($folder, 0, 8)) . ".fasta";
    my $inside = 0;
    #oligo sequence text area
    if ($in_sequence ne "") {
        $inside++;
        my ($fastaSeq, $minL, $maxL) = sequenceload($in_sequence, $filenameNew, $upload_dir);
        if ($fastaSeq eq 'CONTENT') {
            $newhtml .= "<li><b style='color:red'\;>Oligo upload:</b> There was an error uploading your DNA sequences. All the DNA sequences must be in fasta format. Please check the correct format <a href='https://www.cerealsdb.uk.net/cerealgenomics/phyloprimer/format_page.html#oligo'>here</a>.</li>";
        } elsif ($fastaSeq eq 'LEN') {
            $newhtml .= "<li><b style='color:red'\;>Oligo upload:</b> There was an error uploading your DNA sequences. The maximum allowed length is 200 bp.</li>";
        } else {
            if ($minL == $maxL) {
                $fastaLen = $minL;
            } else {
                $fastaLen = $minL . " - " . $maxL;
            }
            my @fasta = split(/,/, $fastaSeq); #details
            if (($fasta[1] > 0) && ($fasta[2] == 0) && ($fasta[3] == 0)) { #single oligo
                if ($fasta[0] == 1) {
                    $newhtml .= "<li>Oligo upload: You uploaded $fasta[0] DNA oligo. The length of the sequence is $fastaLen bp.</li>";
                } elsif ($fasta[0] <= 10) {
                    $newhtml .= "<li>Oligo upload: You uploaded $fasta[0] DNA oligos. The length of the sequences is $fastaLen bp.</li>";
                } else {
                    $newhtml .= "<li><b style='color:red'\;>Oligo upload:</b> There was an error uploading your oligos. The maximum number of allowed sequences is 10.</li>";
                }
            } elsif (($fasta[1] == 0) && ($fasta[2] > 0) && ($fasta[3] == 0)) { #forward and reverse oligo
                if ($fasta[0] == 1) {
                    $newhtml .= "<li>Oligo upload: You uploaded $fasta[0] oligo pair. The length of the sequences is $fastaLen bp.</li>";
                } elsif ($fasta[0] <= 10) {
                    $newhtml = "<li>Oligo upload: You uploaded $fasta[0] oligo pairs. The length of the sequences is $fastaLen bp.</li>";
                } else {
                    $newhtml .= "<li><b style='color:red'\;>Oligo upload:</b> There was an error uploading your oligo pairs. The maximum number of allowed sequences is 10.</li>";
                }
            } elsif (($fasta[1] == 0) && ($fasta[2] == 0) && ($fasta[3] > 0)) { #forward and reverse and probe
                if ($fasta[0] == 1) {
                    $newhtml .= "<li>Oligo upload: You uploaded $fasta[0] oligo and probe pair. The length of the sequences is $fastaLen bp.</li>";
                } elsif ($fasta[0] <= 10) {
                    $newhtml .= "<li>Oligo upload: You uploaded $fasta[0] oligo and probe pairs. The length of the sequences is $fastaLen bp.</li>";
                } else {
                    $newhtml .= "<li><b style='color:red'\;>Oligo upload:</b> There was an error uploading your oligo and probe pairs. The maximum number of allowed sequences is 10.</li>";
                }
            } else {
                if ($fasta[0] == 1) {
                    $newhtml .= "<li>Oligo upload: You uploaded $fasta[0] combo. The length of the sequences is $fastaLen bp.</li>";
                } elsif ($fasta[0] <= 10) {
                    $newhtml .= "<li>Oligo upload: You uploaded $fasta[0] combos. The length of the sequences is $fastaLen bp.</li>";
                } else {
                    $newhtml .= "<li><b style='color:red'\;>Oligo upload:</b> There was an error uploading your combos. The maximum number of allowed sequences is 10.</li>";
                }
            }
        }
    }
    #fasta file
    if ($in_oligo ne "") {
        $inside++;
        my $arg = "in_oligo";
	my ($fastafile,$minL, $maxL,$type) = oligoload($in_oligo, $filenameNew, $upload_dir, $arg);
        if ($fastafile eq 'SIZE') {
            $newhtml .= "<li><b style='color:red'\;>Oligo upload:</b> There was an error uploading your DNA oligos. The file must be smaller than 50 Mb.</li>";
        } elsif ($fastafile eq 'NAME') { #I may remove this
            $newhtml .= "<li><b style='color:red'\;>Oligo upload:</b> There was an error uploading your DNA oligos. The file name contains forbidden characters. The only allowed characters are a-z, A-Z,  ., - and  _.</li>";
        } elsif ($fastafile eq 'CONSISTENCY') {
            $newhtml .= "<li><b style='color:red'\;>Oligo upload:</b> There was an error uploading your DNA oligos. All the oligos need to be consist so you can upload or all three, all two, all one.</li>";
        } elsif ($fastafile eq 'CONTENT') {
            $newhtml .= "<li><b style='color:red'\;>Oligo upload:</b> There was an error uploading your DNA oligos. All the DNA sequences must be in fasta format. Please check the correct format <a href='https://www.cerealsdb.uk.net/cerealgenomics/phyloprimer/format_page.html#oligo'>here</a>.</li>";
        } elsif ($fastafile eq 'LEN') {
            $newhtml .= "<li><b style='color:red'\;>Oligo upload:</b> There was an error uploading your DNA sequences. The maximum allowed length is 200 bp.</li>";
        } else {
            if ($minL == $maxL) {
                $fastaLen = $minL;
            } else {
                $fastaLen = $minL . " - " . $maxL;
            }
            my @fasta = split(/,/, $fastafile); #details
            if ($type == 1) { #single oligo
                if ($fasta[0] == 1) {
                    $newhtml .= "<li>Oligo upload: You uploaded $fasta[0] DNA oligo. The length of the sequence is $fastaLen bp.</li>";
                } elsif ($fasta[0] <= 10) {
                    $newhtml .= "<li>Oligo upload: You uploaded $fasta[0] DNA oligos. The length of the sequences is $fastaLen bp.</li>";
                } else {
                    $newhtml .= "<li><b style='color:red'\;>Oligo upload:</b> there was an error uploading your oligos. The maximum number of allowed sequences is 10.</li>";
                }
            } elsif ($type == 2) { #forward and reverse oligo
                if ($fasta[0] == 1) {
                    $newhtml .= "<li>Oligo upload: you uploaded $fasta[0] oligo pair. The length of the sequences is $fastaLen bp.</li>";
                } elsif ($fasta[0] <= 10) {
                    $newhtml .= "<li>Oligo upload: you uploaded $fasta[0] oligo pairs. The length of the sequences is $fastaLen bp.</li>";
                } else {
                    $newhtml .= "<li><b style='color:red'\;>Oligo upload:</b> there was an error uploading your oligo pairs. The maximum number of allowed sequences is 10.</li>";
                }
            } elsif ($type == 3) { #forward and reverse and probe
                if ($fasta[0] == 1) {
                    $newhtml .= "<li>Oligo upload: You uploaded $fasta[0] oligo and probe pair. The length of the sequences is $fastaLen bp.</li>";
                } elsif ($fasta[0] <= 10) {
                    $newhtml .= "<li>Oligo upload: You uploaded $fasta[0] oligo and probe pairs. The length of the sequences is $fastaLen bp.</li>";
                } else {
                    $newhtml .= "<li><b style='color:red'\;>Oligo upload:</b> there was an error uploading your oligo and probe pairs. The maximum number of allowed sequences is 10.</li>";
                }
            }
        }
    }
    if ($inside > 1) { #if the user submitted both sequences and file
        $newhtml .= "<li><b style='color:red'\;>Oligo upload:</b> there was an error uploading your DNA sequences. Please input your DNA sequences either pasting them in the text area input field or uploading a fasta file.</li>";
    } elsif ($inside == 0) {
        $newhtml .= "<li><b style='color:red'\;>Oligo upload:</b> please input your oligos.</li>";
    }
    $newhtml .= "<br>";

    ##Section 3
    $newhtml .= "<h4><b>BLAST CHECK</b></h4>";
    #negative file
    if ($in_negative ne "") {
        my $arg = "in_negative";
        my $filenameNew = (substr($folder, 0, 8)) . ".negativefasta";
        my ($fastafile, $minL, $maxL) =  fastaload($in_negative, $filenameNew, $upload_dir, $arg);
        if ($fastafile eq 'SIZE') {
            $newhtml .= "<li><b style='color:red'\;>Secondary BLAST check:</b> there was an error uploading your DNA sequences. The file must be smaller than 50 Mb.</li>";
        } elsif ($fastafile eq 'NAME') {
            $newhtml .= "<li><b style='color:red'\;>Secondary BLAST check:</b> there was an error uploading your DNA sequences. The file name contains forbiden characters. The only allowed characters are a-z, A-Z,  ., - and  _.</li>";
        } elsif ($fastafile eq 'CONTENT') {
            $newhtml .= "<li><b style='color:red'\;>Secondary BLAST check:</b> there was an error uploading your DNA sequences. All the DNA sequences must be in fasta format. Please check the correct format <a href='https://www.cerealsdb.uk.net/cerealgenomics/phyloprimer/format_page.html#DNA_seq'>here</a>.</li>";
        } elsif ($fastafile eq 'UNIQUE') {
            $newhtml .= "<li><b style='color:red'\;>Secondary BLAST check:</b> there was an error uploading your DNA sequences. Some of the sequence headings are not unique.</li>";
        } else {
            if ($minL == $maxL) {
                $fastaLen = $minL;
            } else {
                $fastaLen = $minL . " - " . $maxL;
            }
            if ($fastafile == 1) { #at least 4 sequences
                $newhtml .=  "<li>Secondary BLAST check: you uploaded $fastafile DNA sequence. The length of the sequence is $fastaLen bp.</li>";
            } elsif ($fastafile > 1500) { #at least 4 sequences
                $newhtml .=  "<li><b style='color:red'\;>Secondary BLAST check:</b> there was an error uploading your DNA sequences. The maximum number of allowed sequences is 1500.</li>";
            } else {
                $newhtml .=  "<li>Secondary BLAST check: you uploaded $fastafile sequences. The DNA sequence length is $fastaLen bp.</li>";
            }
        }
    } else {
        $newhtml .= "<li>Secondary BLAST check: no negative fasta was uploaded.</li>";
    }
    $newhtml .= "<br>";

    my $badEntry = 0;

    #BLAST mismatch settings
    if ($spec_mis !~ /\A[0-9]+\z/) {
        $newhtml = $newhtml . "<li><b style='color:red'\;>Maximum number of mismatches in the entire oligo sequence:</b> please check. The value must be a positive integer.</li>";
        $badEntry = 1;
    } else {
        $newhtml = $newhtml . "<li>Maximum number of mismatches in the entire oligo sequence: $spec_mis.</li>";
    }

    $newhtml = $newhtml . "<br>";

    if ($spec3_mis !~ /\A[0-9]+\z/) {
        $newhtml = $newhtml . "<li><b style='color:red'\;>Maximum number of mismatches at the 3' oligo end:</b> please check. The value must be a positive integer.</li>";
        $badEntry = 1;
    } else {
        $newhtml = $newhtml . "<li>Maximum number of mismatches at the 3' oligo end: $spec3_mis.</li>";
    }

    $newhtml = $newhtml . "<br>";

    if ($spec3_length !~ /\A[0-9]+\z/) {
        $newhtml = $newhtml . "<li><b style='color:red'\;>Number of oligo bases intended as the 3' oligo end:</b> please check. The value must be a positive integer.</li>";
        $badEntry = 1;
    } else {
        $newhtml = $newhtml . "<li>Number of oligo bases intended as the 3' oligo end: $spec3_length.</li>";
    }

    $newhtml = $newhtml . "<br>";

    ##Section 4
    $newhtml .= "<h4><b>PCR CONDITIONS</b></h4>";

    #na content
    if ($mon_dG !~ /\A[0-9.]+\z/i) {
        $newhtml .= "<li><b style='color:red'\;>Monovalent ion content:</b> please check. This value must be a positive number.</li>";
        $badEntry = 1;
    } else {
        $newhtml .= "<li>Monovalent ion content: $mon_dG mM.</li>";
    }

    $newhtml .= "<br>";

    #mg content
    if ($mg_dG !~ /\A[0-9.]+\z/i) {
        $newhtml .= "<li><b style='color:red'\;>Mg<sup>2+</sup> content:</b> please check. This value must be a positive number.</li>";
    } else {
        $newhtml .= "<li>Mg<sup>2+</sup> content: $mg_dG mM.</li>";
    }

    $newhtml .= "<br>";

    #oligo content
    if ($oligo_dG !~ /\A[0-9.]+\z/i) {
        $newhtml .= "<li><b style='color:red'\;>Oligo content:</b> please check. This value must be a positive number.</li>";
        $badEntry = 1;
    } else {
        $newhtml .= "<li>Oligo content: $oligo_dG uM.</li>";
    }

    $newhtml .= "<br>";

    #dNTP content
    if ($dNTP_dG !~ /\A[0-9.]+\z/i) {
        $newhtml .= "<li><b style='color:red'\;>dNTP content:</b> please check. This value must me a positive number.</li>";
        $badEntry = 1;
    } else {
        $newhtml .= "<li>dNTP content: $dNTP_dG uM.</li>";
    }

    $newhtml .= "<br>";

    #temperature
    if ($t_dG !~ /\A[0-9.]+\z/i) {
        $newhtml .= "<li><b style='color:red'\;>Temperature for secondary structure calculation:</b> please check. This value must be a positive number.</li>";
        $badEntry = 1;
    } else {
        $newhtml .= "<li>Temperature for secondary structure calculation: $t_dG &deg;C.</li>";
    }

    $newhtml .= "<br>";

} else {
    my $input_kind = chop($defSet); #get last letter and understand if a, b or c - here it a
    $folder = $defSet . "C" . $input_kind; #re-construct folder name
    $upload_dir = $path_html . "/uploadsPhyloprimer/" . $folder; #path to the folder

    my $in_file = (substr($folder, 0, 8)) ;

    open(IN, "<${upload_dir}/${in_file}.fasta") or die;
    my $type;
    my $count=0;
    while(defined(my $input=<IN>)) {
	chomp($input);
	$count++;
	if ($count == 1) {
		my @data = split(/\t/, $input);
		$type = scalar(@data);
	}
    }
    close(IN);

    if ($type == 2) { #primer
        $in_file .= ".info4"; #open folder from previous pages
    } elsif ($type == 3) { #probe
        $in_file .= ".info5"; #open folder from previous pages
    } elsif ($type == 1) { #oligo
        $in_file .= ".info6"; #open folder from previous pages
    }

    open(my $file, ">${upload_dir}/${in_file}");

    #Section 0
    print $file "PROJECT\t$project\n";
    print $file "EMAIL\t$email\n";
    print $file "MAILING_LIST\t$mailing_list\n";

    #Section 3
    my $checkNegative = ${upload_dir} . "/" . (substr($folder, 0, 8)) . ".negativefasta";
    if (-e $checkNegative) {
        print $file "NEGATIVE_FILE\tyes\n";
    } else {
        print $file "NEGATIVE_FILE\tno\n";
    }
    #print $file "SPEC_MIS\t$spec_mis\n";
    #print $file "SPEC3_MIS\t$spec3_mis\n";
    #print $file "SPEC3_LENGTH\t$spec3_length\n";

    #Section 4
    print $file "MON_DG\t$mon_dG\n";
    print $file "MG_DG\t$mg_dG\n";
    print $file "OLIGO_DG\t$oligo_dG\n";
    print $file "DNTP_DG\t$dNTP_dG\n";
    print $file "T_DG\t$t_dG\n";

    close($file);


    #open folder from previous pages - user
    my $in_file = (substr($folder, 0, 8)) . ".info";

    open(my $file, ">${upload_dir}/${in_file}");
    ##Section 1
    print $file "\nPhyloPrimer mode: Oligo check\n";

    #Section 0
    print $file "\nUSER DATA\n";
    print $file "Job ID: $project\n";
    print $file "Email: $email\n";
    print $file "Mailing list: $mailing_list\n";

    #Section 3

    print $file "\nBLAST CHECK\n";
    if (-e $checkNegative) {
        print $file "Secondary BLAST check: yes\n";
    } else {

        print $file "Secondary BLAST check: no\n";
    }
    #print $file "Maximum number of mismatches in the entire oligo sequence: $spec_mis\n";
    #print $file "Maximum number of mismatches at the 3' oligo end: $spec3_mis\n";
    #print $file "Number of oligo bases intended as the 3' oligo end: $spec3_length\n";

    #Section 4
    print $file "\nPCR CONDITIONS\n";
    print $file "Monovalent ion concentration: $mon_dG mM\n";
    print $file "Mg2+ concentration: $mg_dG mM\n";
    print $file "Oligo concentration: $oligo_dG uM\n";
    print $file "dNTP concentration: $dNTP_dG uM\n";
    print $file "Temperature for secondary structure calculation: $t_dG °C\n";

    close($file);

}


#oligo loaded in the text area
sub sequenceload {
    my ($content) = $_[0];
    my ($filenameNew) = $_[1];
    my ($upload_dir) = $_[2];
    my $fastafile;
    my @all = split(/\n/, $content);
    my $badEntry = 0;
    my %all_len;
    my $min;
    my $max;
    my $fastafile1 = 0;
    my $fastafile2 = 0;
    my $fastafile3 = 0;
    my %print;

    foreach my $input (@all) {
        chomp($input);
        $input =~ s/^\s+|\s+$//g;
        $fastafile++;
        my @data = split(/\s/, $input);
        my $len = scalar(@data);
        my $entry ='';
        if ($len == 0) {
            $badEntry = 1;
        } elsif ($len == 1) {
            $fastafile1++;
            if ($data[0] =~ /^[acgtnryswkmbdhvn]+$/i) {
                $entry .= $data[0]; #save
                $all_len{length($data[0])} = ''; #seq length
            } else {
                my @oligo = split(/-/, $data[0]);
                my $len = scalar(@oligo);
                if (($oligo[1] =~ /\A[acgtnryswkmbdhvn]+\z/i) && ($len == 2)) {
                    $oligo[0] =~ s/\s/_/g;
                    $entry .= $oligo[0] . "-" . $oligo[1]; #save
                    $all_len{length($oligo[1])} = ''; #seq length
                } else {
                    $badEntry = 1;
                }
            }
            $print{$fastafile}= $entry;
        } elsif ($len == 2) {
            $fastafile2++;
            foreach $a (0,1) {
                if ($data[$a] =~ /\A[acgtnryswkmbdhvn]+\z/i) {
                    $entry .= $data[$a]; #save
                    $all_len{length($data[$a])} = ''; #seq length
                } else {
                    my @oligo = split(/-/, $data[$a]);
                    my $len = scalar(@oligo);
                    if (($oligo[1] =~ /\A[acgtnryswkmbdhvn]+\z/i) && ($len == 2)) {
                        $oligo[0] =~ s/\s/_/g;
                        $entry .= $oligo[0] . "-" . $oligo[1]; #save
                        $all_len{length($oligo[1])} = ''; #seq length
                    } else {
                        $badEntry = 1;
                    }
                }
                if ($a != 1) {
                    $entry .= "\t"; #save
                }
            }
            $print{$fastafile}= $entry;
        } elsif ($len == 3) {
            $fastafile3++;
            foreach $a (0,1,2) {
                if ($data[$a] =~ /\A[acgtnryswkmbdhvn]+\z/i) {
                    $entry .= $data[$a]; #save
                    $all_len{length($data[$a])} = ''; #seq length
                } else {
                    my @oligo = split(/-/, $data[$a]);
                    my $len = scalar(@oligo);
                    if (($oligo[1] =~ /\A[acgtnryswkmbdhvn]+\z/i) && ($len == 2)) {
                        $oligo[0] =~ s/\s/_/g;
                        $entry .= $oligo[0] . "-" . $oligo[1]; #save
                        $all_len{length($oligo[1])} = ''; #seq length
                    } else {
                        $badEntry = 1;
                    }
                }
                if ($a != 2) {
                    $entry .= "\t"; #save
                }
            }
            $print{$fastafile}= $entry;
        } else {
            $badEntry = 1;
        }
    }
    if ($badEntry == 0) {
        my $in_count = 0;
        foreach my $l (sort {$a <=> $b} keys %all_len) { #select min and max length
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
        if ($max > 200) {
            $badEntry = 2;
            $fastafile = "LEN";
        } else {
            open (UPLOADFILE, ">$upload_dir/$filenameNew" ) or die "$!"; #print only if entry is okay
            binmode UPLOADFILE;
            foreach my $oligo (sort keys %print) {
                print UPLOADFILE "$print{$oligo}\n";
            }
            close(UPLOADFILE);
            $fastafile .= "," . $fastafile1 . "," . $fastafile2 . "," . $fastafile3;
        }
    } elsif ($badEntry == 1) {
        $fastafile = "CONTENT";
    }
    return ($fastafile, $min, $max);
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
                        $accession = "SEQUENCE_pp" . $header . " -userHeader- " . $accession;
                    } else { #if others
                        if (length($accession) > 20) {
                            push @head, $accession; #save all accession
                            my $accession0 = $accession;
                            $accession = "SEQUENCE_pp" . $header . " -userHeader- " . $accession;
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
                $print{$accession} = $seq; #accession = sequence
            }

            #tolerate if no header so only one sequence
            #BUT if more than one header, same number between headers and sequences
            if (($badEntry == 0) && ($last eq "sequence") && (($header == $sequence) or (($header == 0) && ($sequence == 1)))) {
                $fastafile = $header;
                open (UPLOADFILE, ">$upload_dir/$filenameNew" ) or die "$!"; #print only if entry is okay
                binmode UPLOADFILE;
                foreach my $head (keys %print) {
                    print UPLOADFILE ">$head\n$print{$head}\n";
                    my $len = length($print{$head});
                    $all_len{$len}++;
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
            } else {
                $fastafile = "CONTENT";
            }
        } else { #uncorrect file name
            $fastafile = "NAME";
        }
    }
    return ($fastafile,$min,$max);
}


#check and load oligo file
sub oligoload {
    my ($filename) = $_[0];
    my ($filenameNew) = $_[1];
    my ($upload_dir) = $_[2];
    my ($arg) = $_[3];

    my $fastafile = 0;
    my $min;
    my $max;
    my $type;

    if (!(defined($filename))) { #if the file is too big
        $fastafile = "SIZE";
    } else { #if the file size is lower than the threshold

        my ( $name, $path, $extension ) = fileparse ( $filename, '..*' );
        $filename = $name . $extension;

        if ( $filename =~ /^([a-zA-Z0-9_\.\-]+)$/ ) { #safe characters for file name
            my $badEntry = 0;
            my $fastafile1 = 0;
            my $fastafile2 = 0;
            my $fastafile3 = 0;
            my %print;
            my %all_len;
            my $upload_filehandle = $cgi->upload($arg);

            while (defined(my $input = <$upload_filehandle>)) { #read the file
                chomp($input);
                $input =~ s/\r//g;
                $fastafile++;
                my @data = split(/\s/, $input); #both space and tabs are valid
                my $len = scalar(@data);
                my $entry ='';
                if ($len == 0) {
                    $badEntry = 1;
                } elsif ($len == 1) {
                    $fastafile1++;
                    if ($data[0] =~ /^[acgtnryswkmbdhvn]+$/i) {
                        $entry .= $data[0]; #save
                        $all_len{length($data[0])} = ''; #seq length
                    } else {
                        my @oligo = split(/-/, $data[0]);
                        my $len = scalar(@oligo);
                        if (($oligo[1] =~ /\A[acgtnryswkmbdhvn]+\z/i) && ($len == 2)) {
                            $oligo[0] =~ s/\s/_/g;
                            $entry .= $oligo[0] . "-" . $oligo[1]; #save
                            $all_len{length($oligo[1])} = ''; #seq length
                        } else {
                            $badEntry = 1;
                        }
                    }
                    $print{$fastafile}= $entry;
                } elsif ($len == 2) {
                    $fastafile2++;
                    foreach $a (0,1) {
                        if ($data[$a] =~ /\A[acgtnryswkmbdhvn]+\z/i) {
                            $entry .= $data[$a]; #save
                            $all_len{length($data[$a])} = ''; #seq length
                        } else {
                            my @oligo = split(/-/, $data[$a]);
                            my $len = scalar(@oligo);
                            if (($oligo[1] =~ /\A[acgtnryswkmbdhvn]+\z/i) && ($len == 2)) {
                                $oligo[0] =~ s/\s/_/g;
                                $entry .= $oligo[0] . "-" . $oligo[1]; #save
                                $all_len{length($oligo[1])} = ''; #seq length
                            } else {
                                $badEntry = 1;
                            }
                        }
                        if ($a != 1) {
                            $entry .= "\t"; #save
                        }
                    }
                    $print{$fastafile}= $entry;
                } elsif ($len == 3) {
                    $fastafile3++;
                    foreach $a (0,1,2) {
                        if ($data[$a] =~ /\A[acgtnryswkmbdhvn]+\z/i) {
                            $entry .= $data[$a]; #save
                            $all_len{length($data[$a])} = ''; #seq length
                        } else {
                            my @oligo = split(/-/, $data[$a]);
                            my $len = scalar(@oligo);
                            if (($oligo[1] =~ /\A[acgtnryswkmbdhvn]+\z/i) && ($len == 2)) {
                                $oligo[0] =~ s/\s/_/g;
                                $entry .= $oligo[0] . "-" . $oligo[1]; #save
                                $all_len{length($oligo[1])} = ''; #seq length
                            } else {
                                $badEntry = 1;
                            }
                        }
                        if ($a != 2) {
                            $entry .= "\t"; #save
                        }
                    }
                    $print{$fastafile}= $entry;
                } else {
                    $badEntry = 1;
                }
            }
            if (($fastafile1 > 0) && ($fastafile2 == 0) && ($fastafile3 == 0)) {
                $type = 1;
            } elsif (($fastafile1 == 0) && ($fastafile2 > 0) && ($fastafile3 == 0)) {
                $type = 2;
            } elsif (($fastafile1 == 0) && ($fastafile2 == 0) && ($fastafile3 > 0)) {
                $type = 3;
            } else {
                $badEntry = 2;
            }
            if ($badEntry == 0) {
                my $in_count = 0;
                foreach my $l (sort {$a <=> $b} keys %all_len) { #select min and max length
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
                if ($max > 200) {
                    $badEntry = 2;
                    $fastafile = "LEN";
                } else {
                    open (UPLOADFILE, ">$upload_dir/$filenameNew" ) or die "$!"; #print only if entry is okay
                    binmode UPLOADFILE;
                    foreach my $oligo (sort keys %print) {
                        print UPLOADFILE "$print{$oligo}\n";
                    }
                    close(UPLOADFILE);
                    $fastafile .= "," . $fastafile1 . "," . $fastafile2 . "," . $fastafile3;
                }
            } elsif ($badEntry == 1) {
                $fastafile = "CONTENT";
            } elsif ($badEntry == 2) {
                $fastafile = "CONSISTENCY";
            }
        } else { #uncorrect file name
            $fastafile = "NAME";
        }
    }
    return ($fastafile,$min,$max,$type);
}

#send data back to html

if (($newhtml =~ /error/) or ($newhtml =~/Please input/) or ($newhtml =~/please check/)) {
    $newhtml = "<h3>" . $newhtml . "<br><b>There was an error in your uploads. Please check and resubmit again.</b></h3><br><br>";
} else {
    $newhtml = "<h3>" . $newhtml . "<br><b>The uploaded data is correct. If you are happy with the upload, please go ahead click on the Check Oligos button so that PhyloPrimer will start with the analysis of your oligos.</b></h3><br><br>";
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
