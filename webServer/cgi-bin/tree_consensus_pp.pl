#!/usr/bin/perl

use strict;
use warnings;
use JSON; #if not already installed, just run "cpan JSON"
use CGI;
use DBD::mysql;
use File::Basename;

#connected to tree_input.html and tree_input.cgi : calculate consensus sequence in Dynamic Selection page

#set path to folders:
my $path_cgi = 'path_to_cgi_folder';
my $path_html = 'path_to_html_folder';

#set mySQL parameters
my $dsn = "mysql_database";
my $user_name = "mysql_user";
my $password = "mysql_password";

my $cgi = CGI->new;

print $cgi->header('application/json;charset=UTF-8');

##input section
my $defSet = $cgi->param('defSet');

my $accession = $cgi->param('accession'); #phylocanvas.selectedNodes

#convert  data to JSON
my $op = JSON -> new -> utf8 -> pretty(1);

my @accession = split(/,/,$accession);

my $folder;

my $badEntry = 0;
my $newhtml;
my $consensus_pos = 'X';
my $consensus_neg = 'X';
my $perc = 0; #percentage of degenerate bases


if ($accession ne '') {
    
    #check defSet
    if ((length($defSet) < 21) or (length($defSet) > 23) or ($defSet !~ /\A[A-Za-z0-9_]+\z/)) {
        $badEntry = 1;
    }
    
    #check accession numbers
    my @accesionNew = @accession;
    foreach my $a (@accession) {
        if ($a =~ /<</) { #user entry
            my ($new) = $a =~ /<<(.*)/;
            $a = $new;
        }
        if ((length($a) > 20) or ($a !~ /\A[A-Za-z0-9_\-\.]+\z/)) {
            $badEntry = 1;
        }
    }
    
    $newhtml = "ERROR_UNK";
    
    if ($badEntry == 0) {
        #read blast file
        my $upload_dir; #path to the folder
        my $in_file;
        my $out_file_positive;
        my $out_file_negative;
        
        my $entry = chop($defSet); #a if DNA, b if alignment and c if newick
        if ($entry eq "a") {
            $folder = $defSet . "Ta";
            #read alignment file
            $upload_dir = $path_html . "/uploadsPhyloprimer/" . $folder;

            $in_file = (substr($folder, 0, 8)) . ".fasta.out.alignment"; #Ta
            
            $out_file_positive = $in_file . ".positive";
            $out_file_negative = $in_file . ".negative";
            
            #remove if previously created
            my $out_file_species = $upload_dir . "/" . (substr($folder, 0, 8)) . ".taxonomy"; #Tb
            
            if (-e $out_file_species) {
                `rm $out_file_species`;
            }
            
            #taxonomy information for visualization criteria
            #connect to mysql database
            my $dbh;
            my $sth;
            
            $dbh = DBI->connect ($dsn, $user_name, $password, { RaiseError => 1 });
            
            my $entry = 0;
            my $ask;
            foreach my $accession (@accession) {
                if ($entry > 0) {
                    $ask .= " OR ";
                }
                $ask .= "(acc='" . $accession . "')";
                $entry++;
            }
            
            $sth = $dbh->prepare("SELECT * FROM DB1_acc_taxid_pp WHERE ($ask)");
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
            
            my %taxaPhylo;
            if ($insideTaxid == 1) { #if at least one corrispondance between accession numbers and taxids
                $sth = $dbh->prepare("SELECT * FROM taxid_taxonomy_pp WHERE ($ask)");
                #execute the prepared statement handle:
                $sth->execute();
                #read results of a query, then clean up
                while (my @ary = $sth->fetchrow_array()) {
                    foreach my $accession (keys %{$taxid_mysql{$ary[0]}}) { #taxid - acc
                        foreach my $i (1..7) {
                            $ary[$i] =~ s/\[sub1\]/\'/g;
                            $ary[$i] =~ s/\[sub2\]/,/g;
                            $ary[$i] =~ s/\[sub3\]/\(/g;
                            $ary[$i] =~ s/\[sub4\]/\)/g;
                            $ary[$i] =~ s/\[sub5\]/:/g;
                            $ary[$i] =~ s/\[sub7\]/\*/g;
                            $ary[$i] =~ s/\[sub8\]/</g;
                            $ary[$i] =~ s/\[sub9\]/>/g;
                            $ary[$i] =~ s/\[sub10\]/-/g;
                            $ary[$i] =~ s/\[sub11\]/\+/g;
                            $ary[$i] =~ s/\[sub12\]/\'/g;
                            $ary[$i] =~ s/\[sub13\]/\#/g;
                            $ary[$i] =~ s/\[sub14\]/&/g;
                            $ary[$i] =~ s/\[sub15\]/\^/g;
                            $ary[$i] =~ s/\[sub16\]/\//g;
                            $taxaPhylo{$i}{$ary[$i]} = ''; #rank - taxon = ''
                        }
                    }
                }
                $sth->finish;
            }
            open(my $file, ">$out_file_species");
            
            foreach my $i (1..7) {
                foreach my $t (sort keys %{$taxaPhylo{$i}}) {
                    print $file "$i\t$t\n";
                }
            }
            close($file);
            
        } elsif ($entry eq "b") {
            $folder = $defSet . "Tb";
            #read alignment file
            $upload_dir = $path_html . "/uploadsPhyloprimer/" . $folder;
            $in_file = (substr($folder, 0, 8)) . ".alignment"; #Tb
            $out_file_positive = $in_file . ".positive";
            $out_file_negative = $in_file . ".negative";
            
            #remove if previously created
            my $out_file_species = $upload_dir . "/" . (substr($folder, 0, 8)) . ".taxonomy"; #Tb
            
            if (-e $out_file_species) {
                `rm $out_file_species`;
            }
            
            #taxonomy information for visualization criteria
            my %userTax;
            
            my $taxUser = (substr($folder, 0, 8)) . ".tax.txt"; #user taxonomy information
            open(TAX, "<${upload_dir}/${taxUser}");
            while (defined(my $input= <TAX> )) { #open taxonomy info
                chomp($input);
                my ($acc,$taxid) = split(/\t/, $input);
                $userTax{$acc} = $taxid;
            }
            close(TAX);
            
            my $entry = 0;
            my $ask;
            my $insideTaxid = 0;
            
            foreach my $acc (@accession) {
                if (defined($userTax{$acc})) {
                    $insideTaxid = 1;
                    if ($entry > 0) {
                        $ask .= " OR ";
                    }
                    $ask .= "(taxid='" . $userTax{$acc} . "')";
                    $entry++;
                }
            }
            
            if ($ask ne '') { #if some correspondances
                
                #connect to mysql database
                my $dbh;
                my $sth;
                
                $dbh = DBI->connect ($dsn, $user_name, $password, { RaiseError => 1 });
                
                my %taxaPhylo;
                if ($insideTaxid == 1) { #if at least one corrispondance between accession numbers and taxids
                    $sth = $dbh->prepare("SELECT * FROM taxid_taxonomy_pp WHERE ($ask)");
                    #execute the prepared statement handle:
                    $sth->execute();
                    #read results of a query, then clean up
                    while (my @ary = $sth->fetchrow_array()) {
                        foreach my $accession (keys %userTax) { #taxid - acc
                            foreach my $i (1..7) {
                                $ary[$i] =~ s/\[sub1\]/\'/g;
                                $ary[$i] =~ s/\[sub2\]/,/g;
                                $ary[$i] =~ s/\[sub3\]/\(/g;
                                $ary[$i] =~ s/\[sub4\]/\)/g;
                                $ary[$i] =~ s/\[sub5\]/:/g;
                                $ary[$i] =~ s/\[sub7\]/\*/g;
                                $ary[$i] =~ s/\[sub8\]/</g;
                                $ary[$i] =~ s/\[sub9\]/>/g;
                                $ary[$i] =~ s/\[sub10\]/-/g;
                                $ary[$i] =~ s/\[sub11\]/\+/g;
                                $ary[$i] =~ s/\[sub12\]/\'/g;
                                $ary[$i] =~ s/\[sub13\]/\#/g;
                                $ary[$i] =~ s/\[sub14\]/&/g;
                                $ary[$i] =~ s/\[sub15\]/\^/g;
                                $ary[$i] =~ s/\[sub16\]/\//g;
                                $taxaPhylo{$i}{$ary[$i]} = ''; #rank - taxon = ''
                            }
                        }
                    }
                    $sth->finish;
                }
                
                open(my $file, ">$out_file_species");
                
                foreach my $i (1..7) {
                    foreach my $t (sort keys %{$taxaPhylo{$i}}) {
                        print $file "$i\t$t\n";
                    }
                }
                close($file);
            }
        } elsif ($entry eq "c") {
            $folder = $defSet . "Tc";
            #read alignment file
            $upload_dir = $path_html . "/uploadsPhyloprimer/" . $folder;
            $in_file = (substr($folder, 0, 8)) . ".alignment"; #Tc
            $out_file_positive = $in_file . ".positive";
            $out_file_negative = $in_file . ".negative";
            
            #remove if previously created
            my $out_file_species = $upload_dir . "/" . (substr($folder, 0, 8)) . ".taxonomy"; #Tb
            
            if (-e $out_file_species) {
                `rm $out_file_species`;
            }
            
            #taxonomy information for visualization criteria
            my %userTax;
            
            my $taxUser = (substr($folder, 0, 8)) . ".tax.txt"; #user taxonomy information
            open(TAX, "<${upload_dir}/${taxUser}");
            while (defined(my $input= <TAX> )) { #open taxonomy info
                chomp($input);
                my ($acc,$taxid) = split(/\t/, $input);
                $userTax{$acc} = $taxid;
            }
            close(TAX);
            
            my $entry = 0;
            my $ask;
            my $insideTaxid = 0;
            
            foreach my $acc (@accession) {
                if (defined($userTax{$acc})) {
                    $insideTaxid = 1;
                    if ($entry > 0) {
                        $ask .= " OR ";
                    }
                    $ask .= "(taxid='" . $userTax{$acc} . "')";
                    $entry++;
                }
            }
            
            if ($ask ne '') { #if some correspondances
                
                #connect to mysql database
                my $dbh;
                my $sth;
                
                $dbh = DBI->connect ($dsn, $user_name, $password, { RaiseError => 1 });
                
                my %taxaPhylo;
                if ($insideTaxid == 1) { #if at least one corrispondance between accession numbers and taxids
                    $sth = $dbh->prepare("SELECT * FROM taxid_taxonomy_pp WHERE ($ask)");
                    #execute the prepared statement handle:
                    $sth->execute();
                    #read results of a query, then clean up
                    while (my @ary = $sth->fetchrow_array()) {
                        foreach my $accession (keys %userTax) { #taxid - acc
                            foreach my $i (1..7) {
                                $ary[$i] =~ s/\[sub1\]/\'/g;
                                $ary[$i] =~ s/\[sub2\]/,/g;
                                $ary[$i] =~ s/\[sub3\]/\(/g;
                                $ary[$i] =~ s/\[sub4\]/\)/g;
                                $ary[$i] =~ s/\[sub5\]/:/g;
                                $ary[$i] =~ s/\[sub7\]/\*/g;
                                $ary[$i] =~ s/\[sub8\]/</g;
                                $ary[$i] =~ s/\[sub9\]/>/g;
                                $ary[$i] =~ s/\[sub10\]/-/g;
                                $ary[$i] =~ s/\[sub11\]/\+/g;
                                $ary[$i] =~ s/\[sub12\]/\'/g;
                                $ary[$i] =~ s/\[sub13\]/\#/g;
                                $ary[$i] =~ s/\[sub14\]/&/g;
                                $ary[$i] =~ s/\[sub15\]/\^/g;
                                $ary[$i] =~ s/\[sub16\]/\//g;
                                $taxaPhylo{$i}{$ary[$i]} = ''; #rank - taxon = ''
                            }
                        }
                    }
                    $sth->finish;
                }
                
                open(my $file, ">$out_file_species");
                
                foreach my $i (1..7) {
                    foreach my $t (sort keys %{$taxaPhylo{$i}}) {
                        print $file "$i\t$t\n";
                    }
                }
                close($file);
            }
        } else {
            $newhtml = "ERROR_UNK";
            last;
        }
        
        
        my $negCheckFile = ${upload_dir} . "/" . ${out_file_negative};
        if (-e $negCheckFile) { #check that negative file does not exist already
            `rm $negCheckFile`;
            my $out_file = (substr($folder, 0, 8)) . ".consensus.negative";
            $negCheckFile = ${upload_dir} . "/" . ${out_file};
            `rm $negCheckFile`;
        }
        
        
        
        my %prev_fasta;
        foreach my $acc (@accession) { #selected accession numbers
            $prev_fasta{$acc} = "";
        }
        
        #check fasta -- later
        #retrieve sequences
        $accession = "";
        my $inside = 0;
        open(IN_FASTA, "<${upload_dir}/${in_file}");
        my $seq;
        my %positive;
        my %negative;
        my $pos_check = 0;
        my $neg_check = 0;
        while (defined(my $input = <IN_FASTA>)) {
            chomp($input);
            if ($input =~ /^>/) {
                ($accession) = $input =~ /^>(.*)/;
                if ($accession =~ /-userHeader-/) {
                    my @info = split(/ -userHeader- /, $accession);
                    $accession = $info[0];
                }
                if (defined($prev_fasta{$accession})) {
                    $inside = 1;
                    $seq = '';
                } else {
                    $inside = 0;
                    $seq = '';
                }
            } else {
                $seq .= $input;
                if ($inside == 1) {
                    $pos_check = 1;
                    $positive{$accession} = $seq;
                } else {
                    $neg_check = 1;
                    $negative{$accession} = $seq;
                }
            }
        }
        close(IN_FASTA);
        
        
        if ($pos_check == 1) {
            open(my $file_positive, ">${upload_dir}/${out_file_positive}");
            foreach my $p (keys %positive) {
                print $file_positive ">$p\n$positive{$p}\n";
            }
            close($file_positive);
            #run consensus - positive
            $consensus_pos = `$path_cgi/consensus_pp.pl ${upload_dir}/${out_file_positive}`;
            my $out_file = (substr($folder, 0, 8)) . ".consensus.positive";
            open(my $file, ">${upload_dir}/${out_file}");
            print $file "$consensus_pos";
            close($file);
            
            my $all = () = $consensus_pos =~ /R|Y|S|W|K|M|B|D|H|V|N|A|T|G|C/g; #all bases
            my $deg = () = $consensus_pos =~ /R|Y|S|W|K|M|B|D|H|V|N/g; #degenerate bases
            
            $perc = ($deg/$all)*100;
            
        } else {
            $newhtml = "ERROR_UNK";
        }
        
        if ($neg_check == 1) {
            open(my $file_negative, ">${upload_dir}/${out_file_negative}");
            foreach my $n (keys %negative) {
                print $file_negative ">$n\n$negative{$n}\n";
            }
            close($file_negative);
            #run consensus - negative
            $consensus_neg = `$path_cgi/consensus_pp.pl ${upload_dir}/${out_file_negative}`;
            my $out_file = (substr($folder, 0, 8)) . ".consensus.negative";
            open(my $file, ">${upload_dir}/${out_file}");
            print $file "$consensus_neg";
            close($file);
        }
                
        #html is then SUCCESS or ERROR. if SUCCESS other script is launched
        if ($consensus_pos ne "X") { #what about if no negatives??
            $newhtml = "SUCCESS";
            if ($consensus_neg eq 'X') { #no negative consensus
                $consensus_pos =~ tr/-//d;
                
            } else { #negative consensus is present
                
                my @pos = ($consensus_pos =~ m/./g); #create hashes with position and base
                my @neg = ($consensus_neg =~ m/./g); #create hashes with position and base
                
                $consensus_pos = '';
                $consensus_neg = '';
                
                my %pos_count;
                my %neg_count;
                
                my $count = 0;
                foreach my $p (@pos) {
                    $count++;
                    $pos_count{$count} = $p;
                }
                
                my $count = 0;
                foreach my $n (@neg) {
                    $count++;
                    $neg_count{$count} = $n;
                }
                
                foreach my $pos (sort {$a <=> $b} keys %pos_count) { #remove gaps when present in both consensuses
                    if (($pos_count{$pos} ne '-') or ($neg_count{$pos} ne '-')) {
                        $consensus_pos .= $pos_count{$pos};
                        $consensus_neg .= $neg_count{$pos};
                    }
                }
                
            }
            
        }
        
        
        $perc = sprintf("%.0f", $perc);
        
        if ($perc > 0) {
            $newhtml = "ERROR_PERC";
        }
        
        #print in file all accession numbers
        my $out_file = (substr($folder, 0, 8)) . ".allAccession";
        open(my $fileAcc, ">${upload_dir}/${out_file}") or die;
        print $fileAcc join("\n", @accesionNew);
        print $fileAcc "\n";
        close($fileAcc);
        
    }
    
    
    
} else {
    $newhtml = "ERROR_POS"; #if no selection
}
#result array
my @result;

push @result, "result";
push @result, $newhtml;

push @result, "pos";
push @result, $consensus_pos;

push @result, "neg";
push @result, $consensus_neg;

push @result, "perc";
push @result, $perc;

my $json = $op -> encode({
    @result
});
print $json;




