#!/usr/bin/perl

use strict;

############################
####Different modalities####
############################

#1) primer design - Primer_design_pp.pl
#2) probe design - Probe_design_pp.pl
#3) oligo design - Oligo_design_pp.pl

#4) primer check - Primer_check_pp.pl
#5) probe check - Probe_check_pp.pl
#6) oligo check - Oligo_check_pp.pl

#set path to folders:
my $path_cgi = 'path_to_cgi_folder';
my $path_html = 'path_to_html_folder';

my $done = 0;

while ($done == 0) {
    my $listAnalyses1 = `ls $path_html/uploadsPhyloprimer/*/*info1 2> /dev/null`; #completed updates ready to be analysed
    my $listAnalyses2 = `ls $path_html/uploadsPhyloprimer/*/*info2 2> /dev/null`; #completed updates ready to be analysed
    my $listAnalyses3 = `ls $path_html/uploadsPhyloprimer/*/*info3 2> /dev/null`; #completed updates ready to be analysed
    my $listAnalyses4 = `ls $path_html/uploadsPhyloprimer/*/*info4 2> /dev/null`; #completed updates ready to be analysed
    my $listAnalyses5 = `ls $path_html/uploadsPhyloprimer/*/*info5 2> /dev/null`; #completed updates ready to be analysed
    my $listAnalyses6 = `ls $path_html/uploadsPhyloprimer/*/*info6 2> /dev/null`; #completed updates ready to be analysed
    
    my $listProcess = `ls $path_html/analysesPhyloprimer/*/*.info 2> /dev/null`; #how many processes
    my @listProcess = split(/\n/, $listProcess);
    
    my $active = 0;
    foreach my $l (@listProcess) {
        $active++;
    }
    my $allowed = 4 - $active; #maximum of processes I can launched
    
    print "Free processes: $allowed\n";
    
    if ($listAnalyses1 ne '') {
        my @listAnalyses = split(/\n/, $listAnalyses1);
        
        foreach my $l (@listAnalyses) {
            if ($allowed > 0) {
                my ($folder) = $l =~ /uploadsPhyloprimer\/(.*?)\//;
                my $folder0 = $path_html . "/uploadsPhyloprimer/" . $folder;
                my ($file) = $l =~ /${folder}\/(.*)/;
                my $datestring = localtime();
                print "$datestring : moved $folder0\n";
                `mv $folder0 $path_html/analysesPhyloprimer`;
                sleep 5;
                system("perl $path_cgi/Primer_design_pp.pl -file $file -folder $path_html/analysesPhyloprimer/$folder &");
                $allowed--;
            }
        }
        
    }
    
    if ($listAnalyses2 ne '') {
        my @listAnalyses = split(/\n/, $listAnalyses2);
        foreach my $l (@listAnalyses) {
            if ($allowed > 0) {
                my ($folder) = $l =~ /uploadsPhyloprimer\/(.*?)\//;
                my $folder0 = $path_html . "/uploadsPhyloprimer/" . $folder;
                my ($file) = $l =~ /${folder}\/(.*)/;
                my $datestring = localtime();
                print "$datestring : moved $folder0\n";
                `mv $folder0 $path_html/analysesPhyloprimer`;
                sleep 5;
                system("perl $path_cgi/Probe_design_pp.pl -file $file -folder $path_html/analysesPhyloprimer/$folder &");
                $allowed--;
            }
        }
    }
    
    if ($listAnalyses3 ne '') {
        my @listAnalyses = split(/\n/, $listAnalyses3);
        foreach my $l (@listAnalyses) {
            if ($allowed > 0) {
                my ($folder) = $l =~ /uploadsPhyloprimer\/(.*?)\//;
                my $folder0 = $path_html . "/uploadsPhyloprimer/" . $folder;
                my ($file) = $l =~ /${folder}\/(.*)/;
                my $datestring = localtime();
                print "$datestring : moved $folder0\n";
                `mv $folder0 $path_html/analysesPhyloprimer`;
                sleep 5;
                system("perl $path_cgi/Oligo_design_pp.pl -file $file -folder $path_html/analysesPhyloprimer/$folder &");
                $allowed--;
            }
        }
    }
    
    if ($listAnalyses4 ne '') {
        my @listAnalyses = split(/\n/, $listAnalyses4);
        foreach my $l (@listAnalyses) {
            if ($allowed > 0) {
                my ($folder) = $l =~ /uploadsPhyloprimer\/(.*?)\//;
                my $folder0 = $path_html . "/uploadsPhyloprimer/" . $folder;
                my ($file) = $l =~ /${folder}\/(.*)/;
                my $datestring = localtime();
                print "$datestring : moved $folder0\n";
                `mv $folder0 $path_html/analysesPhyloprimer`;
                sleep 5;
                system("perl $path_cgi/Primer_check_pp.pl -file $file -folder $path_html/analysesPhyloprimer/$folder &");
                $allowed--;
            }
        }
    }
    
    if ($listAnalyses5 ne '') {
        my @listAnalyses = split(/\n/, $listAnalyses5);
        foreach my $l (@listAnalyses) {
            if ($allowed > 0) {
                my ($folder) = $l =~ /uploadsPhyloprimer\/(.*?)\//;
                my $folder0 = $path_html . "/uploadsPhyloprimer/" . $folder;
                my ($file) = $l =~ /${folder}\/(.*)/;
                my $datestring = localtime();
                print "$datestring : moved $folder0\n";
                `mv $folder0 $path_html/analysesPhyloprimer`;
                sleep 5;
                system("perl $path_cgi/Probe_check_pp.pl -file $file -folder $path_html/analysesPhyloprimer/$folder &");
                $allowed--;
            }
        }
    }
    
    if ($listAnalyses6 ne '') {
        my @listAnalyses = split(/\n/, $listAnalyses6);
        foreach my $l (@listAnalyses) {
            if ($allowed > 0) {
                my ($folder) = $l =~ /uploadsPhyloprimer\/(.*?)\//;
                my $folder0 = $path_html . "/uploadsPhyloprimer/" . $folder;
                my ($file) = $l =~ /${folder}\/(.*)/;
                my $datestring = localtime();
                print "$datestring : moved $folder0\n";
                `mv $folder0 $path_html/analysesPhyloprimer`;
                sleep 5;
                system("perl $path_cgi/Oligo_check_pp.pl -file $file -folder $path_html/analysesPhyloprimer/$folder &");
                $allowed--;
            }
        }
    }
    
    sleep 60;
}
