#!/usr/bin/perl

use JSON; #if not already installed, just run "cpan JSON"
use strict;
use warnings;
use CGI qw(:standard);
use CGI::Carp qw ( fatalsToBrowser );
use File::Basename;

#connected to sequence_input.html : check data uploaded in Premade selection page

#set path to folders:
my $path_html = 'path_to_html_folder';

#Set maximum size of uploaded file
$CGI::POST_MAX = 1024 * 50000; #50 MegaBytes

#header needs to start with > and all the white spaces are changed with _
my $newhtml;
my $folderOld;

#create folder with 20 random characters (N.B. it needs to be sent back to html)
my @chars = ("A".."Z", "a".."z");
$folderOld .= $chars[rand @chars] for 1..20;

#open cgi connection
my $cgi = CGI->new;
print $cgi->header('application/json;charset=UTF-8');

#import variables from html
my $input_kind = $cgi->param('input_kind'); #radio button. it can be DNA, alignment or consensus

my $in_sequence = $cgi->param("in_sequence"); #fasta sequence
my $in_fasta = $cgi->param("in_fasta"); #fasta file

my $in_alignment = $cgi->param("in_alignment"); #alignment file

my $in_consensus = $cgi->param('in_consensus'); #consensus file

my $pass;
my $fastaLen = 0;

my $buttonValue;

#folder name is passed only if the file is correct!

#convert  data to JSON
my $op = JSON -> new -> utf8 -> pretty(1);

if ($input_kind eq "DNA") { #DNA upload
    my $folder = $folderOld . "Sa";
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
            $newhtml = "There was an error uploading your DNA sequences. All the DNA sequences must be in fasta format. Please check the correct format <a href='https://www.cerealsdb.uk.net/cerealgenomics/phyloprimer/format_page.html#DNA_seq'>here</a>.<br>";
        } elsif ($fastaSeq eq 'UNIQUE') {
            $newhtml = "There was an error uploading your DNA sequences. Some of the sequence headings are not unique.<br>";
        } else {
            if ($minL == $maxL) {
                $fastaLen = $minL;
            } else {
                $fastaLen = $minL . " - " . $maxL;
            }
                if ($fastaSeq == 1) {
                    $buttonValue = "Check Consensus";
                    $newhtml = "You uploaded $fastaSeq DNA sequence and it will directly be used as consensus sequence. The length of the sequence is $fastaLen bp.<br>";
                } elsif ($fastaSeq <= 1500) {
                    $buttonValue = "Create Consensus";
                    $newhtml = "You uploaded $fastaSeq DNA sequences. The length of the sequences is $fastaLen bp.<br>";
                } else {
                    $newhtml = "There was an error uploading your DNA sequences. The maximum number of allowed sequences is 1500.<br>";
                }

        }
    }
    #fasta file
    if ($in_fasta ne "") {
        $inside++;
        my $arg = "in_fasta";
        my ($fastafile, $minL, $maxL) = fastaload($in_fasta, $filenameNew, $upload_dir, $arg);
        if ($fastafile eq 'SIZE') {
            $newhtml = "There was an error uploading your DNA sequences. The file must be smaller than 50 Mb.<br>";
        } elsif ($fastafile eq 'NAME') {
            $newhtml = "There was an error uploading your DNA sequences. The file name contains forbidden characters. The only allowed characters are a-z, A-Z,  ., - and  _.<br>";
        } elsif ($fastafile eq 'CONTENT') {
            $newhtml = "There was an error uploading your DNA sequences. All the DNA sequences must be in fasta format. Please check the correct format <a href='https://www.cerealsdb.uk.net/cerealgenomics/phyloprimer/format_page.html#DNA_seq'>here</a>.<br>";
        } elsif ($fastafile eq 'UNIQUE') {
            $newhtml = "There was an error uploading your DNA sequences. Some of the sequence headings are not unique.<br>";
        } else {
            if ($minL == $maxL) {
                $fastaLen = $minL;
            } else {
                $fastaLen = $minL . " - " . $maxL;
            }
                if ($fastafile == 1) {
                    $newhtml = "You uploaded $fastafile DNA sequence and it will directly be used as consensus sequence. The length of the sequence is $fastaLen bp.<br>";
                    $buttonValue = "Check Consensus";
                } elsif ($fastafile <= 1500) {
                    $buttonValue = "Create Consensus";
                    $newhtml = "You uploaded $fastafile DNA sequences. The length of the sequences is $fastaLen bp.<br>";
                } else {
                    $newhtml = "There was an error uploading your DNA sequences. The maximum number of allowed sequences is 1500.<br>";
                }
        }
    }
    
    
    if ($inside > 1) { #if the user submitted both sequences and file
        $newhtml = "There was an error uploading your DNA sequences. Please input your DNA sequences either pasting them in the text area input field or uploading a fasta file.<br>";
    } elsif ($inside == 0) {
        $newhtml = "Please input your DNA sequences.<br>";
    }
    
} elsif ($input_kind eq "alignment") { #alignment upload
    my $folder = $folderOld . "Sb";
    $pass = $folderOld . "b";
    my $upload_dir = $path_html . "/uploadsPhyloprimer/" . $folder . "/"; #path to the folder
    `mkdir $upload_dir`; #create new directory
    
    #alignment file
    if ($in_alignment ne "") {
        my $arg = "in_alignment";
        my $filenameNew = (substr($folder, 0, 8)) . ".alignment";
        my ($fastafile, $minL, $maxL) = fastaload($in_alignment, $filenameNew, $upload_dir, $arg);
        if ($fastafile eq 'SIZE') {
            $newhtml = "There was an error uploading your alignment sequences. The file must be smaller than 50 Mb.<br>";
        } elsif ($fastafile eq 'NAME') {
            $newhtml = "There was an error uploading your alignment sequences. The file name contains forbiden characters. The only allowed characters are a-z, A-Z,  ., - and  _.<br>";
        } elsif ($fastafile eq 'LEN') {
            $newhtml = "There was an error uploading your alignment sequences. All the alignment sequences should be the same length.<br>";
        } elsif ($fastafile eq 'CONTENT') {
            $newhtml = "There was an error uploading your alignment sequences. All the DNA sequences must be in fasta format. Please check the correct format <a href='https://www.cerealsdb.uk.net/cerealgenomics/phyloprimer/format_page.html#alignment'>here</a>.<br>";
        } elsif ($fastafile eq 'UNIQUE') {
            $newhtml = "There was an error uploading your alignment sequences. Some of the sequence headings are not unique.<br>";
        } else {
            if ($minL == $maxL) {
                $fastaLen = $minL;
            } else {
                $fastaLen = $minL . " - " . $maxL;
            }
            if ($fastafile == 1) { #if one seq
                $newhtml =  "You uploaded $fastafile alignment sequence and it will directly be used as consensus sequence. The length of the sequence is $fastaLen bp.<br>";
                $buttonValue = "Check Consensus";
            } elsif ($fastafile <= 1500) { #maximum 1500 sequence
                $buttonValue = "Create Consensus";
                $newhtml = "You uploaded $fastafile alignment sequences. The length of the sequences is $fastaLen bp.<br>";
            } else {
                $newhtml = "There was an error uploading your alignment sequences. The maximum number of allowed sequences is 1500.<br>";
            }
        }
    } else {
        $newhtml = "Please input your alignment sequences.<br>";
    }
    
} elsif ($input_kind eq "consensus") { #consensus upload
    my $folder = $folderOld . "Sc";
    $pass = $folderOld . "c";
    my $upload_dir = $path_html . "/uploadsPhyloprimer/" . $folder . "/"; #path to the folder
    `mkdir $upload_dir`; #create new directory
    #fasta file
    if ($in_consensus ne "") {
        my $arg = "in_consensus";
        my $filenameNew = (substr($folder, 0, 8)) . ".consensus";
        my ($fastafile, $minL, $maxL) = fastaload($in_consensus, $filenameNew, $upload_dir, $arg);
        if ($fastafile eq 'SIZE') {
            $newhtml = "There was an error uploading your DNA sequences. The file must be smaller than 50 Mb.<br>";
        } elsif ($fastafile eq 'NAME') {
            $newhtml = "There was an error uploading your DNA sequences. The file name contains forbiden characters. The only allowed characters are a-z, A-Z,  ., - and  _.<br>";
        } elsif ($fastafile eq 'CONTENT') {
            $newhtml = "There was an error uploading your DNA sequences. All the DNA sequences must be in fasta format. Please check the correct format <a href='https://www.cerealsdb.uk.net/cerealgenomics/phyloprimer/format_page.html#consensus'>here</a>.<br>";
        } else {
            if ($minL == $maxL) {
                $fastaLen = $minL;
            } else {
                $fastaLen = $minL . " - " . $maxL;
            }
            if ($fastafile == 1) { #if one seq
                $buttonValue = "Check Consensus";
                $newhtml =  "You uploaded $fastafile DNA sequence and it will directly be used as consensus sequence. The length of the sequence is $fastaLen bp.<br>";
            } elsif ($fastafile > 1) { #maximum 1500 sequence
                $newhtml = "There was an error uploading your DNA sequences. The maximum number of allowed sequences is 1.<br>";
            }
        }
    } else {
        $newhtml = "Please input your consensus sequence.<br>";
    }
} else {
    $newhtml = "Something went wrong";
}


#send data back to html
if (($newhtml =~ /error/) or ($newhtml =~/Please input/)) {
    $newhtml = "<h3>" . $newhtml . "<br><b>There was an error in your uploads. Please check and resubmit again.</b></h3>";
} else {
    $newhtml = "<h3>" . $newhtml . "<br><b>The uploaded data is correct. If you are happy with the upload, go ahead and press the " . $buttonValue . " button.</b></h3>";
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
        $print{$accession} = $seq; #accession = sequence
    }
    if (($badEntry == 0) && ($last eq "sequence") && (($header == $sequence) or (($header == 0) && ($sequence == 1)))) {
        $fastaSeq = $header;
        open (UPLOADFILE, ">$upload_dir/$filenameNew" ) or die "$!"; #print only if entry is okay
	if ($sequence > 1) {
        foreach my $head (keys %print) {
            print UPLOADFILE ">$head\n$print{$head}\n";
            my $len = length($print{$head});
            $all_len{$len}++;
        }
	} else {
	         foreach my $head (keys %print) {
            print UPLOADFILE ">$head\n$print{$head}";
            my $len = length($print{$head});
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
                    $seq = "";
                } else {
                    if ($input ne "") {
                        $last = "sequence";
                        $count++;
                        if ((($input =~ /\A[acgtnryswkmbdhvn]+\z/i) && ($arg eq "in_fasta")) or (($input =~ /\A[acgtnryswkmbdhvn\.\-]+\z/i) && (($arg eq "in_alignment") or ($arg eq "in_consensus")))) { #safe characters for sequence / alignments
                            $seq .= $input;
                        } else {
                            $badEntry = 1;
                            last;
                        }
                        if ($count == 1) { #how many sequences
                            $sequence++;
                        }
                    }
                }
                $print{$accession} = $seq; #accession = sequence
            }
            
            
            #chack if all sequence alignment are long the same
            if ($arg eq 'in_alignment') {
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
		if ($sequence > 1) {
                foreach my $head (keys %print) {
                    print UPLOADFILE ">$head\n$print{$head}\n";
                    my $len = length($print{$head});
                    $all_len{$len}++;
                }
		} else {
		                foreach my $head (keys %print) {
                    print UPLOADFILE ">$head\n$print{$head}";
                    my $len = length($print{$head});
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
    return ($fastafile,$min,$max);
}

#send data back to html
my @result;
push @result, "result";
push @result, $newhtml;

push @result, "button";
push @result, $buttonValue;

push @result, "set";
push @result, $pass;

my $json = $op -> encode({
    @result
});
print $json;




