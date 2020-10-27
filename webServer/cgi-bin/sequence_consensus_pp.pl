#!/usr/bin/perl

use strict;
use warnings;
use JSON; #if not already installed, just run "cpan JSON"
use CGI;
use File::Basename;

#connected to sequence_input.html : calculate consensus sequence in Premade Selection page

#set path to folders:
my $path_cgi = 'path_to_cgi_folder';
my $path_html = 'path_to_html_folder';

my $cgi = CGI->new;

print $cgi->header('application/json;charset=UTF-8');

##input section
my $defSet = $cgi->param('defSet'); #folder

my $newhtml;
my $consensus;
my $origin;

#convert  data to JSON
my $op = JSON -> new -> utf8 -> pretty(1);

my $input_kind = chop($defSet); #get last letter and understand kind of folder
my $folder = $defSet . "S" . $input_kind; #"S" because sequence input
my $upload_dir = $path_html . "/uploadsPhyloprimer/" . $folder; #path to the folder

if ($input_kind eq "a") { #if DNA
    my $in_file = (substr($folder, 0, 8)) . ".fasta";
    my $out_file = (substr($folder, 0, 8)) . ".alignment";
    #check file (again)
    my $seqCheck;
    ($seqCheck,$consensus) =  fastaload($in_file,$upload_dir); #pass file and directory and get back number of sequence/error
    if ($seqCheck =~ /ERROR/) {
        $newhtml = $newhtml . $seqCheck;
    } else {
        if ($seqCheck > 1) {
            #run mafft
            `mafft --localpair --maxiterate 1000 --quiet --thread 35 ${upload_dir}/${in_file} > ${upload_dir}/${out_file}`; #adjust threads
            #consensus
            $consensus = `$path_cgi/consensus_pp.pl ${upload_dir}/${out_file}`;
            my $out_file = (substr($folder, 0, 8)) . ".consensus.positive";
            open(my $file, ">${upload_dir}/${out_file}");
            print $file "$consensus";
            close($file);
            $newhtml = "SUCCESS";
            $origin = "calculated";
        } else {
            my $out_file = (substr($folder, 0, 8)) . ".consensus.positive";
            `cp ${upload_dir}/${in_file} ${upload_dir}/${out_file}`;
            $newhtml = "SUCCESS";
            $origin = "checked";
        }
    }
} elsif ($input_kind eq "b") { #if alignment
    my $in_file = (substr($folder, 0, 8)) . ".alignment";
    #check file (again)
    my $seqCheck;
    ($seqCheck,$consensus) =  fastaload($in_file,$upload_dir); #pass file and directory and get back number of sequence/error
    
    if ($seqCheck ne 'ERROR') {
        if ($seqCheck > 1) {
            #consensus
            $consensus = `$path_cgi/consensus_pp.pl ${upload_dir}/${in_file}`;
            my $out_file = (substr($folder, 0, 8)) . ".consensus.positive";
            open(my $file, ">${upload_dir}/${out_file}");
            print $file "$consensus";
            close($file);
            $newhtml = "SUCCESS";
            $origin = "calculated";
        } else {
            my $out_file = (substr($folder, 0, 8)) . ".consensus.positive";
            `cp ${upload_dir}/${in_file} ${upload_dir}/${out_file}`;
            $newhtml = "SUCCESS";
            $origin = "checked";
        }
    } else {
        $newhtml = "ERROR";
    }
    
} elsif ($input_kind eq "c") { #if consensus
    my $in_file = (substr($folder, 0, 8)) . ".consensus";
    #check file (again)
    my $seqCheck;
    ($seqCheck,$consensus) =  fastaload($in_file,$upload_dir); #pass file and directory and get back number of sequence/error
    if ($seqCheck ne 'ERROR') {
        if ($seqCheck == 1) {
            $newhtml = "SUCCESS";
            $origin = "checked";
        } else {
            $newhtml = "ERROR";
        }
    } else {
        $newhtml = "ERROR";
    }
} else {
    $newhtml = "ERROR54 " . $input_kind . " and " . $defSet;
}


$consensus =~ tr/-//d;
$consensus =~ tr/.//d;

my $all = () = $consensus =~ /R|Y|S|W|K|M|B|D|H|V|N|A|T|G|C/g; #all bases
my $deg = () = $consensus =~ /R|Y|S|W|K|M|B|D|H|V|N/g; #degenerate bases

my $perc = ($deg/$all)*100;
$perc = sprintf("%.0f", $perc);

if ($perc > 0) {
    $newhtml = "ERROR_PERC";
}


#result array
my @result;

push @result, "result";
push @result, $newhtml;

push @result, "origin";
push @result, $origin;

push @result, "pos";
push @result, $consensus;

push @result, "perc";
push @result, $perc;


my $json = $op -> encode({
    @result
});
print $json;

sub fastaload {
    my ($in_file) = $_[0]; #filename
    my ($upload_dir) = $_[1]; #directory
    my $fastafile;
    my $consensusF;
    if ($in_file =~ /^([a-zA-Z0-9_\.\-]+)$/ ) { #safe characters for file name
        my $header = 0;
        my $count = 0;
        my $sequence = 0;
        my $badEntry = 0;
        my $last;
        
        open(IN_CHECK, "<${upload_dir}/${in_file}");
        while (defined(my $input = <IN_CHECK>)) { #read the file
            chomp($input);
            if ($input =~ /^>/) {
                $header++; #how many headers
                $count = 0;
                $last = "header";
            } else {
                if ($input ne "") {
                    $last = "sequence";
                    $count++;
                    $consensusF .= $input;
                    if ($input !~ /\A[acgtnryswkmbdhvn\-\.]+\z/i) { #safe characters for sequence
                        $badEntry = 1;
                    }
                    if ($count == 1) { #how many sequences
                        $sequence++;
                    }
                }
            }
        }
        close(IN_CHECK);
        #tolerate if no header so only one sequence
        #BUT if more than one header, same number between headers and sequences
        if (($badEntry == 0) && ($last eq "sequence") && (($header == $sequence) or (($header == 0) && ($sequence == 1)))) {
            $fastafile = $sequence;
            if ($fastafile == 1) {
                my $out_file = (substr($folder, 0, 8)) . ".consensus.positive";
                open(my $file, ">${upload_dir}/${out_file}");
                print $file "$consensusF";
                close($file);
            }
        } else {
            $fastafile = "ERROR";
        }
    } else { #uncorrect file name
        $fastafile = "ERROR";
    }
    return ($fastafile,$consensusF);
}
