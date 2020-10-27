#!/usr/bin/perl

use strict;
use warnings;
use JSON; #if not already installed, just run "cpan JSON"
use CGI qw(:standard);
use DBI;
use DBD::mysql;
use File::Basename;
use POSIX;

#connected to tree_input.html : create dynamic tree in Dynamic selection page

#set path to folders:
my $path_html = 'path_to_html_folder';
my $path_db =  'path_to_DB1';

#set mySQL parameters
my $dsn = "mysql_database";
my $user_name = "mysql_user";
my $password = "mysql_password";

my $cgi = CGI->new;
print $cgi->header('application/json;charset=UTF-8');

#radio button that indicates type on input
#the user can input:
#1. DNA sequences - run BLAST, alignment and iq-tree
#2. alignments - run iq-tree
#3. newick - just give newick
my $input_tree = $cgi->param('input_tree');

#folder
my $defSet = $cgi->param('defSet');
#blast parameters
my $identity_blast = $cgi->param('identity_blast');
my $coverage_blast = $cgi->param('coverage_blast');
my $evalue_blast = $cgi->param('evalue_blast');

my $noBLAST = 0;

my $newhtmlSen;
my $newhtml;
my $count_acc;
my $count_acc_tot;
my $taxonomy;

#convert data to JSON
my $op = JSON -> new -> utf8 -> pretty(1);

my @acc_mysql;
my @acc_mysql_tot;

my $newick;

my $match = 0;

my $badEntry = 0;

my $sub = 0;

my $query_count = 0;

#create tables
my $userTable = "<h3>PhyloPrimer substitutes the header of the inputted sequences with USER_INPUT followed by a progressive number. Please find the details in the following table.</h3><table style='width:100%'><tr><th>PhyloPrimer name</th><th>User original name</th></tr>";
my $clusterTable = "<h3>If there are identical sequences represented by more than four entries, PhyloPrimer pools the entries together and renames them as CLUSTER followed by a progressive number. Please find which entries were grouped within each cluster in the following table.</h3><table style='width:100%'><tr><th style='width: 10%;'>Cluster</th><th style='width: 10%;'>Number of entries</th><th style='width: 60%;'>Entries</th><th style='width: 20%;'>Species-level taxonomy</th></tr>";
my %user;
my %clusterAll;
my %cluster;
my %clusterTable_first;
my %clusterTable_count;
my %prev_fasta;
my %data_fasta;


my $multiple = 0; #number of clusters
my $unique = 0; #how many unique sequences retrieved by BLAST
#my %prev_fasta_acc;
#print genes
my %gene;

my $double = 0;

#check BLAST parameters
if ($evalue_blast !~ /\A[0-9\.]+\z/i) {
    $badEntry = 1;
}
if (($identity_blast < 60) or ($identity_blast > 100) or ($identity_blast !~ /\A[0-9]+\z/)) {
    $badEntry = 1;
}
    
if (($coverage_blast < 60) or ($coverage_blast > 100) or ($coverage_blast !~ /\A[0-9]+\z/)) {
    $badEntry = 1;
}


#check that folder name is correct
if ((length($defSet) != 21) or ($defSet !~ /\A[A-Za-z]+\z/)) {
    $badEntry = 2;
}


my $query;
my $bitScore;
my %winner;
my $count_representative = 0;

my %userTable_order;

my $folder;
my $upload_dir;

my %taxAssigned;

if ($badEntry == 0) {
    #selected DNA
    if ($input_tree eq 'gene')  {
        my @blast;
        chop($defSet);
        $folder = $defSet . "Ta";
        $upload_dir = $path_html . "/uploadsPhyloprimer/" . $folder; #path to the folder
        my $in_file = (substr($folder, 0, 8)) . ".fasta";
        my $accession;
        
        open(IN_FASTA, "<${upload_dir}/${in_file}");
        while (defined(my $input = <IN_FASTA>)) {
            chomp($input);
            if ($input =~ /^>/) {
                if ($input =~ / -userHeader- /) { #get new and old name
                    my @new = split(/ -userHeader- /, $input);
                    ($accession) = $new[0] =~ /^>(.*)/;
                    $user{$accession} = $new[1]; #new name = user name
                    #substitute accession special characters
                    $new[1] =~ s/_/ /g;#all spaces to _
                    $new[1] =~ s/^\s+|\s+$//g; #taxon name
                    $new[1] =~ s/\[sub1\]/\'/g; #single quote
                    $new[1] =~ s/\[sub2\]/,/g; #comma
                    $new[1] =~ s/\[sub3\]/\(/g; #bracket (
                    $new[1] =~ s/\[sub4\]/\)/g; #bracket )
                    $new[1] =~ s/\[sub5\]/:/g; #column
                    $new[1] =~ s/\[sub6\]/;/g; #semi column
                    $new[1] =~ s/\[sub7\]/\*/g; #semi column
                    $new[1] =~ s/\[sub8\]/</g; #lower
                    $new[1] =~ s/\[sub9\]/>/g; #higher
                    $new[1] =~ s/\[sub10\]/-/g; #minus
                    $new[1] =~ s/\[sub11\]/\+/g; #plus
                    $new[1] =~ s/\[sub12\]/\`/g; #hyphen`
                    $new[1] =~ s/\[sub13\]/\#/g; #
                    $new[1] =~ s/\[sub14\]/&/g; #&
                    $new[1] =~ s/\[sub15\]/\^/g; #&
                    $new[1] =~ s/\[sub16\]/\//g; #/
                    $new[1] =~ s/\[sub17\]/_/g; #_
                    $userTable_order{"<tr><td>" . $accession . "</td><td>" . $new[1] . "</td></tr>"} = '';
                }
                $query_count++;
            } else {
                if ($input !~ /\A[acgtnryswkmbdhvn]+\z/i) { #safe characters for sequence - bases and wildcards
                } else {
                    push @{$prev_fasta{$input}{'ALL'}}, $accession; #all the accession numbers
                    $prev_fasta{$input}{'COUNT'}++;
                }
            }
        }
        close(IN_FASTA);
        
        my $countBLAST = ceil(400/$query_count); #how many unique sequences must be retrieved per gene
        
        foreach my $u (sort keys %userTable_order) {
            $userTable .= $u;
        }
        
        #blast
        @blast = `blastn -db $path_db/DB1 -query ${upload_dir}/${in_file} -task megablast -max_target_seqs 2500 -perc_identity ${identity_blast} -qcov_hsp_perc ${coverage_blast} -evalue ${evalue_blast} -outfmt 5 -num_threads 30`;
    
        #clean blast result
        my $inside = 0;
        my $seq;
        my %all;
        my @info = (); #accession row
        my $id1;
        
        #BLAST fasta file
        my $out_file = $in_file . ".out";
        open(my $file, ">${upload_dir}/${out_file}");
        
        my $countUnique = 0;
        my %uniqueSeq;
        my $over = 0;
        
        my $pro = ''; #protein
        foreach my $input (@blast) { #get information from BLAST output
            chomp($input);
            if ($input =~ /<Iteration_query-def>/) { #user sequence name
                my ($query0) = $input =~ /<Iteration_query-def>(.*)<\/Iteration_query-def>/;
                my @new = split(/ -userHeader- /, $query0);
                $query = $new[0];
                $countUnique = 0;
                $over = 0;
                undef %uniqueSeq;
            } elsif ($input =~ /<Hit_def>/) { #GenBank sequence accession number
                ($accession) = $input =~ /<Hit_def>(.*?) /;
                #my ($ref) = $input =~ /<Hit_def>${accession}\s+(.*?) /;
                #($pro) = $input =~ /${accession}\s+${ref}\s+(.*)<\/Hit_def>/;
                ($pro) = $input =~ /${accession}(.*)<\/Hit_def>/;
                $pro =~ s/^\s+|\s+$//g;
                $inside = 0;
            } elsif ($input =~ /<Hsp_hseq>/) { #sequence
                $inside++;
                ($seq) = $input =~ /<Hsp_hseq>(.*)<\/Hsp_hseq>/;
                $seq =~ s/-//g;
                if ($inside == 1) {
                    $count_acc_tot++;
                    $double = 0;
                    if ((!(defined($uniqueSeq{$seq})))) {
                        $countUnique++;
                        $uniqueSeq{$seq} = '';
                    }
                    if ($countUnique <= $countBLAST) { #do not enter here after the threshold
                        $gene{$accession} = $pro;
                        $match++;
                        push @{$data_fasta{$accession}{'QUERY'}}, $query; ##I could create an hash
                        #checking if accession is unique because it may be not unique when working with more input sequences
                        if (defined($all{$accession})) {
                            if (length($seq) > length($all{$accession})) {
                                $all{$accession} = $seq;
                            }
                        } else {
                            $all{$accession} = $seq;
                            push @acc_mysql_tot, $accession;
                            push @{$prev_fasta{$seq}{'ALL'}}, $accession; #all the accession numbers
                            $prev_fasta{$seq}{'COUNT'}++;
                        }
                    } else {
                        $over = 1;
                    }
                }
            } elsif ($input =~ /<Hsp_bit-score>/) { #bit-score
                if ($over == 0) {
                    ($bitScore) = $input =~ /<Hsp_bit-score>(.*)<\/Hsp_bit-score>/;
                    if ((defined($data_fasta{$accession}{'BIT'})) && ($bitScore > $data_fasta{$accession}{'BIT'})) {
                        $data_fasta{$accession}{'BIT'} = $bitScore;
                        $double = 1;
                        $winner{$accession} = $query; #reporting BLAST results for this query
                    }
                }
            } elsif ($input =~ /<Hsp_evalue>/) { #e-value
                if ($over == 0) {
                    if ((!(defined($data_fasta{$accession}{'EV'}))) or ($double == 1)) {
                        $data_fasta{$accession}{'BIT'} = $bitScore;
                        my ($evalue) = $input =~ /<Hsp_evalue>(.*)<\/Hsp_evalue>/;
                        $data_fasta{$accession}{'EV'} = $evalue;
                        $winner{$accession} = $query; #reporting BLAST results for this query
                    }
                }
            } elsif ($input =~ /<Hsp_query-from>/) { #user start
                if ($over == 0) {
                    if ((!(defined($data_fasta{$accession}{'Q1'}))) or ($double == 1)) {
                        my ($query1) = $input =~ /<Hsp_query-from>(.*)<\/Hsp_query-from>/;
                        $data_fasta{$accession}{'Q1'} = $query1;
                    }
                }
            } elsif ($input =~ /<Hsp_query-to>/) { #user start
                if ($over == 0) {
                    if ((!(defined($data_fasta{$accession}{'Q2'}))) or ($double == 1)) {
                        my ($query2) = $input =~ /<Hsp_query-to>(.*)<\/Hsp_query-to>/;
                        $data_fasta{$accession}{'Q2'} = $query2;
                    }
                }
            } elsif ($input =~ /<Hsp_hit-from>/) { #hit start
                if ($over == 0) {
                    if ((!(defined($data_fasta{$accession}{'H1'}))) or ($double == 1)) {
                        my ($hit1) = $input =~ /<Hsp_hit-from>(.*)<\/Hsp_hit-from>/;
                        $data_fasta{$accession}{'H1'} = $hit1;
                    }
                }
            } elsif ($input =~ /<Hsp_hit-to>/) { #hit end
                if ($over == 0) {
                    if ((!(defined($data_fasta{$accession}{'H2'}))) or ($double == 1)) {
                        my ($hit2) = $input =~ /<Hsp_hit-to>(.*)<\/Hsp_hit-to>/;
                        $data_fasta{$accession}{'H2'} = $hit2;
                    }
                }
            } elsif ($input =~ /<Hsp_identity>/) { #id
                if ($over == 0) {
                    if ((!(defined($data_fasta{$accession}{'ID'}))) or ($double == 1)) {
                        ($id1) = $input =~ /<Hsp_identity>(.*)<\/Hsp_identity>/;
                    }
                }
            } elsif ($input =~ /<Hsp_align-len>/) { #len
                if ($over == 0) {
                    if ((!(defined($data_fasta{$accession}{'ID'}))) or ($double == 1)) {
                        my ($id2) = $input =~ /<Hsp_align-len>(.*)<\/Hsp_align-len>/;
                        my $perc = sprintf("%.2f", (($id1/$id2) * 100));
                        $data_fasta{$accession}{'ID'} = $perc;
                    }
                }
            }
        }
        
        
        foreach my $seq (keys %prev_fasta) { #all the sequences
            $unique++; #count how many unique sequences
            my $cluster;
            if ($prev_fasta{$seq}{'COUNT'} < 4) {
                
                $count_representative += $prev_fasta{$seq}{'COUNT'};
                foreach my $acc (@{$prev_fasta{$seq}{'ALL'}}) {
                    print $file ">$acc\n$seq\n";
                    push @acc_mysql, $acc;
                    
                }
            } else { #if a sequence is the same for more than 4 entries
                $multiple++; #how many clusters
                $cluster = "CLUSTER" . $multiple;
                print $file ">$cluster\n$seq\n";
                foreach my $acc (@{$prev_fasta{$seq}{'ALL'}}) {
                    $clusterAll{$acc} = $cluster; #acc - cluster
                    if (defined($user{$acc})) {
                        my $accessionHigh = "<b>" . $acc . "</b>";
                        $acc = $accessionHigh;
                        $user{$cluster} = '';
                    }
                }
                $clusterTable_first{$cluster} = "<tr><td>" . $cluster . "</td><td>" . $prev_fasta{$seq}{'COUNT'} . "</td><td>" . join(', ', @{$prev_fasta{$seq}{'ALL'}}) . "</td><td>"; #cluster table
                push @acc_mysql, $cluster;
                $clusterTable_count{$cluster} = $prev_fasta{$seq}{'COUNT'};
                $count_representative += $prev_fasta{$seq}{'COUNT'};
            }
        }
        close($file);
        
        if ($match < 4) { #if the BLAST search retrieves less than 4 sequences
            $noBLAST = 1;
        } elsif (($unique < 4) && ($count_acc_tot == 2500)) { #if BLAST retrieved more than 4 sequences but not more then 4 unique sequences and maximum retrieved
            $noBLAST = 2;
        } elsif ($unique < 4) { #if BLAST retrieved more than 4 sequences but not more then 4 unique sequences
            $noBLAST = 3;
        } else {
            #alignment
            $in_file = $out_file;
            $out_file = $in_file . ".alignment";
            my $out_fileTree = $in_file . ".tree";
            
            if ($match <= 200) { #if sequences are less than 200
                `mafft --localpair --maxiterate 1000 --thread 30 --quiet ${upload_dir}/${in_file}  > ${upload_dir}/${out_file}`;
                `FastTreeMP -nt -gtr -boot 1000 < ${upload_dir}/${out_file} > ${upload_dir}/${out_fileTree}`;
            } else { #if sequences are more than 200
                `mafft --maxiterate 1000 --thread 30 --quiet ${upload_dir}/${in_file}  > ${upload_dir}/${out_file}`;
                `FastTreeMP -nt -gtr -boot 1000 < ${upload_dir}/${out_file} > ${upload_dir}/${out_fileTree}`;
            }
            
            my %al;
            #open alignment from mafft -
            open(AL, "<${upload_dir}/${out_file}");
            while (defined(my $input = <AL>)) {
                chomp($input);
                if ($input =~ /^>/) {
                    my ($old) = $input =~ /^>(.*)/;
                    my $new = $old;
                    $new =~ s/\./_/g; #convert accession as mafft does in guide tree
                    $al{$new} = $old;
                }
            }
            close(AL);
            #open guide tree from mafft
            my $in_tree = $in_file . ".tree";
            open(TREE, "<${upload_dir}/${in_tree}");
            while (defined(my $input = <TREE>)) {
                chomp($input);
                if ($input =~ /^[0-9]+/) {
                    my ($input1) = $input =~ /^[0-9]+_(.*)/; #remove number added by mafft
                    #$newick .= $input1; #save newick removing all new line characters
                    if (defined($al{$input1})) {
                        $newick .= $al{$input1}; #save newick removing all new line characters
                    } else {
                        $newick .= $input1; #save newick removing all new line characters
                    }
                } else {
                    $newick .= $input; #save newick removing all new line characters
                }
            }
            $newick .= ";";
            close(TREE);
        }
        #######################################
        ###########retrieve taxonomy###########
        #######################################
        if ($noBLAST == 0) {
            
            #connect to mysql database
            my $dbh;
            my $sth;
            
            $dbh = DBI->connect ($dsn, $user_name, $password, { RaiseError => 1 });
            
            my $entry = 0;
            my $ask;

            foreach my $accession (@acc_mysql_tot) {
                if ($entry > 0) {
                    $ask .= " OR ";
                }
                $newhtml .= $accession;
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
            foreach my $taxid (keys %taxid_mysql) {
                if ($entry > 0) {
                    $ask .= " OR ";
                }
                $ask .= "(taxid='" . $taxid . "')";
                $entry++;
            }
            my $newick_domain = $newick;
            my $newick_phylum = $newick;
            my $newick_class = $newick;
            my $newick_order = $newick;
            my $newick_family = $newick;
            my $newick_genus = $newick;
            my $newick_species = $newick;
            
            $sth = $dbh->prepare("SELECT * FROM taxid_taxonomy_pp WHERE ($ask)");
            #execute the prepared statement handle:
            $sth->execute();
            #read results of a query, then clean up
            
            my $composed;
            while (my @ary = $sth->fetchrow_array()) {
                foreach my $acc (keys %{$taxid_mysql{$ary[0]}}) { #taxid - acc
                    
                    $taxAssigned{$acc} = ''; #accessions with associated description
                    
                    if ($newick_domain =~ /$acc:/) {
                        
                        if ($query_count == 1) {
                            #domain
                            $composed = $acc . ";" . $ary[1] . ";" . $gene{$acc} . ";;" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_domain =~ s/$acc:/$composed/g;
                            #phylum
                            $composed = $acc . ";" . $ary[2] . ";" . $gene{$acc} . ";;" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_phylum =~ s/$acc:/$composed/g;
                            #class
                            $composed = $acc . ";" . $ary[3] . ";" . $gene{$acc} . ";;" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_class =~ s/$acc:/$composed/g;
                            #order
                            $composed = $acc . ";" . $ary[4] . ";" . $gene{$acc} . ";;" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_order =~ s/$acc:/$composed/g;
                            #family
                            $composed = $acc . ";" . $ary[5] . ";" . $gene{$acc} . ";;" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_family =~ s/$acc:/$composed/g;
                            #genus
                            $composed = $acc . ";" . $ary[6] . ";" . $gene{$acc} . ";;" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_genus =~ s/$acc:/$composed/g;
                            #species
                            $composed = $acc . ";" . $ary[7] . ";" . $gene{$acc} . ";;" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_species =~ s/$acc:/$composed/g;
                        } else {
                            
                            my $match = scalar(@{$data_fasta{$acc}{'QUERY'}});
                            my $query_info = '';
                            if ($match == 1) {
                                $query_info = "This entry got hit only by " . $data_fasta{$acc}{'QUERY'}[0] . ".";
                                
                            } else {
                                $query_info = "This entry got hit by " . join(' - ', @{$data_fasta{$acc}{'QUERY'}}) . ". The reported results are for " . $winner{$acc} . ".";
                            }
                            #domain
                            $composed = $acc . ";" . $ary[1] . ";" . $gene{$acc} . ";" . $query_info . ";" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_domain =~ s/$acc:/$composed/g;
                            #phylum
                            $composed = $acc . ";" . $ary[2] . ";" . $gene{$acc} . ";" . $query_info . ";" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_phylum =~ s/$acc:/$composed/g;
                            #class
                            $composed = $acc . ";" . $ary[3] . ";" . $gene{$acc} . ";" . $query_info . ";" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_class =~ s/$acc:/$composed/g;
                            #order
                            $composed = $acc . ";" . $ary[4] . ";" . $gene{$acc} . ";" . $query_info . ";" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_order =~ s/$acc:/$composed/g;
                            #family
                            $composed = $acc . ";" . $ary[5] . ";" . $gene{$acc} . ";" . $query_info . ";" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_family =~ s/$acc:/$composed/g;
                            #genus
                            $composed = $acc . ";" . $ary[6] . ";" . $gene{$acc} . ";" . $query_info . ";" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_genus =~ s/$acc:/$composed/g;
                            #species
                            $composed = $acc . ";" . $ary[7] . ";" . $gene{$acc} . ";" . $query_info . ";" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_species =~ s/$acc:/$composed/g;
                        }
                    } else {
                        $cluster{$clusterAll{$acc}}{'domain'}{$ary[1]}++; #cluster - rank - taxon
                        $cluster{$clusterAll{$acc}}{'phylum'}{$ary[2]}++; #cluster - rank - taxon
                        $cluster{$clusterAll{$acc}}{'class'}{$ary[3]}++; #cluster - rank - taxon
                        $cluster{$clusterAll{$acc}}{'order'}{$ary[4]}++; #cluster - rank - taxon
                        $cluster{$clusterAll{$acc}}{'family'}{$ary[5]}++; #cluster - rank - taxon
                        $cluster{$clusterAll{$acc}}{'genus'}{$ary[6]}++; #cluster - rank - taxon
                        $cluster{$clusterAll{$acc}}{'species'}{$ary[7]}++; #cluster - rank - taxon
                    }
                    
                }
            }
            $sth->finish;
            
            #add description to taxa without taxonomy
            foreach my $acc (@acc_mysql_tot) {
                if (!(defined($taxAssigned{$acc}))) {
                    if ($newick_domain =~ /$acc:/) {
                        if ($query_count == 1) {
                            #domain
                            $composed = $acc . ";;" . $gene{$acc} . ";;" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_domain =~ s/$acc:/$composed/g;
                            #phylum
                            $composed = $acc . ";;" . $gene{$acc} . ";;" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_phylum =~ s/$acc:/$composed/g;
                            #class
                            $composed = $acc . ";;" . $gene{$acc} . ";;" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_class =~ s/$acc:/$composed/g;
                            #order
                            $composed = $acc . ";;" . $gene{$acc} . ";;" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_order =~ s/$acc:/$composed/g;
                            #family
                            $composed = $acc . ";;" . $gene{$acc} . ";;" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_family =~ s/$acc:/$composed/g;
                            #genus
                            $composed = $acc . ";;" . $gene{$acc} . ";;" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_genus =~ s/$acc:/$composed/g;
                            #species
                            $composed = $acc . ";;" . $gene{$acc} . ";;" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_species =~ s/$acc:/$composed/g;
                        } else {
                            
                            my $match = scalar(@{$data_fasta{$acc}{'QUERY'}});
                            my $query_info = '';
                            if ($match == 1) {
                                $query_info = "This entry got hit only by " . $data_fasta{$acc}{'QUERY'}[0] . ".";
                                
                            } else {
                                $query_info = "This entry got hit by " . join(' - ', @{$data_fasta{$acc}{'QUERY'}}) . ". The reported results are for " . $winner{$acc} . ".";
                            }
                            #domain
                            $composed = $acc . ";;" . $gene{$acc} . ";" . $query_info . ";" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_domain =~ s/$acc:/$composed/g;
                            #phylum
                            $composed = $acc . ";;" . $gene{$acc} . ";" . $query_info . ";" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_phylum =~ s/$acc:/$composed/g;
                            #class
                            $composed = $acc . ";;" . $gene{$acc} . ";" . $query_info . ";" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_class =~ s/$acc:/$composed/g;
                            #order
                            $composed = $acc . ";;" . $gene{$acc} . ";" . $query_info . ";" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_order =~ s/$acc:/$composed/g;
                            #family
                            $composed = $acc . ";;" . $gene{$acc} . ";" . $query_info . ";" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_family =~ s/$acc:/$composed/g;
                            #genus
                            $composed = $acc . ";;" . $gene{$acc} . ";" . $query_info . ";" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_genus =~ s/$acc:/$composed/g;
                            #species
                            $composed = $acc . ";;" . $gene{$acc} . ";" . $query_info . ";" . $data_fasta{$acc}{'BIT'} . ";" . $data_fasta{$acc}{'EV'} . ";" . $data_fasta{$acc}{'Q1'} . " - " . $data_fasta{$acc}{'Q2'} . ";" . $data_fasta{$acc}{'H1'} . " - " . $data_fasta{$acc}{'H2'} . ";" . $data_fasta{$acc}{'ID'} . ":";
                            $newick_species =~ s/$acc:/$composed/g;
                        }
                    }
                }
            }
            
            #substitute clusters
            
            foreach my $clu (sort keys %cluster) {
                $entry = 0;
                foreach my $taxon (sort {$a <=> $b} keys %{$cluster{$clu}{'domain'}}) {
                    if ($entry == 0) {
                        $composed = $clu . ";" . $taxon . ":";
                        $newick_domain =~ s/$clu:/$composed/g;
                    }
                    $entry++;
                }
                $entry = 0;
                foreach my $taxon (sort {$cluster{$clu}{'phylum'}{$b} <=> $cluster{$clu}{'phylum'}{$a}} keys %{$cluster{$clu}{'phylum'}}) {
                    if ($entry == 0) {
                        $composed = $clu . ";" . $taxon . ":";
                        $newick_phylum =~ s/$clu:/$composed/g;
                    }
                    $entry++;
                }
                $entry = 0;
                foreach my $taxon (sort {$cluster{$clu}{'class'}{$b} <=> $cluster{$clu}{'class'}{$a}} keys %{$cluster{$clu}{'class'}}) {
                    if ($entry == 0) {
                        $composed = $clu . ";" . $taxon . ":";
                        $newick_class =~ s/$clu:/$composed/g;
                    }
                    $entry++;
                }
                $entry = 0;
                foreach my $taxon (sort {$cluster{$clu}{'order'}{$b} <=> $cluster{$clu}{'order'}{$a}} keys %{$cluster{$clu}{'order'}}) {
                    if ($entry == 0) {
                        $composed = $clu . ";" . $taxon . ":";
                        $newick_order =~ s/$clu:/$composed/g;
                    }
                    $entry++;
                }
                $entry = 0;
                foreach my $taxon (sort {$cluster{$clu}{'family'}{$b} <=> $cluster{$clu}{'family'}{$a}} keys %{$cluster{$clu}{'family'}}) {
                    if ($entry == 0) {
                        $composed = $clu . ";" . $taxon . ":";
                        $newick_family =~ s/$clu:/$composed/g;
                    }
                    $entry++;
                }
                $entry = 0;
                foreach my $taxon (sort {$cluster{$clu}{'genus'}{$b} <=> $cluster{$clu}{'genus'}{$a}} keys %{$cluster{$clu}{'genus'}}) {
                    if ($entry == 0) {
                        $composed = $clu . ";" . $taxon . ":";
                        $newick_genus =~ s/$clu:/$composed/g;
                    }
                    $entry++;
                }
                
                my @tableSpecies;
                $entry = 0;
                foreach my $taxon (sort {$cluster{$clu}{'species'}{$b} <=> $cluster{$clu}{'species'}{$a}} keys %{$cluster{$clu}{'species'}}) {
                    if ($entry == 0) {
                        $composed = $clu . ";" . $taxon . ":";
                        $newick_species =~ s/$clu:/$composed/g;
                    }
                    my $perc = ($cluster{$clu}{'species'}{$taxon} / $clusterTable_count{$clu}) * 100;
                    $perc = sprintf("%.1f", $perc);
                    $taxon =~ s/\[sub1\]/\'/g; #single quote
                    $taxon =~ s/\[sub2\]/,/g; #comma
                    $taxon =~ s/\[sub3\]/\(/g; #bracket (
                    $taxon =~ s/\[sub4\]/\)/g; #bracket )
                    $taxon =~ s/\[sub5\]/:/g; #column
                    $taxon =~ s/\[sub6\]/;/g; #semi column
                    $taxon =~ s/\[sub7\]/\*/g; #semi column
                    $taxon =~ s/\[sub8\]/</g; #lower
                    $taxon =~ s/\[sub9\]/>/g; #higher
                    $taxon =~ s/\[sub10\]/-/g; #minus
                    $taxon =~ s/\[sub11\]/\+/g; #plus
                    $taxon =~ s/\[sub12\]/\`/g; #hyphen`
                    $taxon =~ s/\[sub13\]/\#/g; #
                    $taxon =~ s/\[sub14\]/&/g; #&
                    $taxon =~ s/\[sub15\]/\^/g; #&
                    $taxon =~ s/\[sub16\]/\//g; #/
                    my $comp = $taxon . " (" . $perc . "%)";
                    push @tableSpecies, $comp;
                    $entry++;
                }
                $clusterTable_first{$clu} = $clusterTable_first{$clu} . join(', ', @tableSpecies) . "</td></tr>";
            }
            foreach my $num (1..$multiple) {
                my $clu = "CLUSTER" . $num;
                $clusterTable = $clusterTable . $clusterTable_first{$clu};
            }
            
            #pass all phylogetic trees to html
            $newhtml = $newick_domain . "&" . $newick_phylum . "&" . $newick_class . "&" . $newick_order . "&" . $newick_family . "&" . $newick_genus . "&" . $newick_species;
            
            foreach my $bold (keys %user) {
                if ($newhtml =~ /$bold;/) {
                    $composed = "<<" . $bold . ";";
                    $newhtml =~ s/$bold;/$composed/g;
                } elsif ($newhtml =~ /$bold:/) {
                    $composed = "<<" . $bold . ":";
                    $newhtml =~ s/$bold:/$composed/g;
                }
            }
            
            #pass number of accession numbers
            $count_acc = scalar @acc_mysql;
            
        } elsif ($noBLAST == 1) { #if no BLAST matches were found
            $newhtml = "ERROR";
            $newhtmlSen = "<h3>The BLAST search found less than 4 matches. Try to relax the BLAST parameters.</h3>";
            $count_acc = 0;
        } elsif ($noBLAST == 2) { #if BLAST retrieved more than 4 sequences but not more then 4 unique sequences and maximum retrieved
            $newhtml = "ERROR";
            $newhtmlSen = "<h3>The BLAST search retrieved $count_acc_tot DNA sequences with an identity percentage higher than $identity_blast\%, a coverage higher than $coverage_blast\% and an e-value lower than $evalue_blast. However, less than 4 sequences were unique and the maximum number of sequence has been retrieved by the BLAST search so we suggest you to create your own sequence dataset and then submit it directly.</h3>";
            $count_acc = 0;
        } elsif ($noBLAST == 3) { #if BLAST retrieved more than 4 sequences but not more then 4 unique sequences
            $newhtml = "ERROR";
            $newhtmlSen = "<h3>The BLAST search retrieved $count_acc_tot DNA sequences with an identity percentage higher than $identity_blast\%, a coverage higher than $coverage_blast\% and an e-value lower than $evalue_blast. However, less than 4 sequences were unique so we could not build the dynamic tree. Try to relax the BLAST parameters.</h3>";
            $count_acc = 0;
        }
        
    } elsif ($input_tree eq 'sequence')  {
        chop($defSet);
        $folder = $defSet . "Tb";
        $upload_dir = $path_html . "/uploadsPhyloprimer/" . $folder; #path to the folder
        my $in_file = (substr($folder, 0, 8)) . ".fasta";
        
        my $accession;
        
        my $seqPP = 0;
        
        open(IN_SEQUENCE, "<${upload_dir}/${in_file}");
        while (defined(my $input = <IN_SEQUENCE>)) {
            chomp($input);
            if ($input =~ /^>/) {
                if ($input =~ / -userHeader- /) { #get new and old name
                    my @new = split(/ -userHeader- /, $input);
                    ($accession) = $new[0] =~ /^>(.*)/;
                    #substitute accession special characters
                    $new[1] =~ s/_/ /g;#all spaces to _
                    $new[1] =~ s/^\s+|\s+$//g; #taxon name
                    $new[1] =~ s/\[sub1\]/\'/g; #single quote
                    $new[1] =~ s/\[sub2\]/,/g; #comma
                    $new[1] =~ s/\[sub3\]/\(/g; #bracket (
                    $new[1] =~ s/\[sub4\]/\)/g; #bracket )
                    $new[1] =~ s/\[sub5\]/:/g; #column
                    $new[1] =~ s/\[sub6\]/;/g; #semi column
                    $new[1] =~ s/\[sub7\]/\*/g; #semi column
                    $new[1] =~ s/\[sub8\]/</g; #lower
                    $new[1] =~ s/\[sub9\]/>/g; #higher
                    $new[1] =~ s/\[sub10\]/-/g; #minus
                    $new[1] =~ s/\[sub11\]/\+/g; #plus
                    $new[1] =~ s/\[sub12\]/\`/g; #hyphen`
                    $new[1] =~ s/\[sub13\]/\#/g; #
                    $new[1] =~ s/\[sub14\]/&/g; #&
                    $new[1] =~ s/\[sub15\]/\^/g; #&
                    $new[1] =~ s/\[sub16\]/\//g; #/
                    $new[1] =~ s/\[sub17\]/_/g; #_
                    $seqPP++;
                    $userTable_order{"<tr><td>" . $accession . "</td><td>" . $new[1] . "</td></tr>"} = '';
                } else {
                    my ($accession) = $input =~ /^>(.*)/;
                }
                push @acc_mysql, $accession;
            }
        }
        close(IN_SEQUENCE);
        
        foreach my $num (1..$seqPP) {
            my $u = "USER_INPUT" . $num;
            $userTable .= $u;
        }
        
        
        #perform alignment
        my $out_file = (substr($folder, 0, 8)) . ".alignment";
        
        if (scalar(@acc_mysql) <= 200) { #if sequences are less than 200
            `mafft --localpair --treeout --maxiterate 1000 --thread 30 --quiet ${upload_dir}/${in_file}  > ${upload_dir}/${out_file}`;
            
        } else { #if sequences are more than 200
            `mafft --treeout --maxiterate 1000 --thread 30 --quiet ${upload_dir}/${in_file}  > ${upload_dir}/${out_file}`;
        }
        
        my %al;
        #open alignment from mafft -
        open(AL, "<${upload_dir}/${out_file}");
        while (defined(my $input = <AL>)) {
            chomp($input);
            if ($input =~ /^>/) {
                my ($old) = $input =~ /^>(.*)/;
                my $new = $old;
                $new =~ s/\./_/g; #convert accession as mafft does in guide tree
                $al{$new} = $old;
            }
        }
        close(AL);
        
        #open guide tree from mafft
        my $in_tree = (substr($folder, 0, 8)) . ".fasta.tree";
        open(TREE, "<${upload_dir}/${in_tree}");
        while (defined(my $input = <TREE>)) {
            chomp($input);
            if ($input =~ /^[0-9]+/) {
                if ($input =~ /-userHeader-/) { #remove double name
                    my ($input1) = $input =~ /^[0-9]+_(.*)_-userHeader-/; #remove number added by mafft
                    $input1 =~ s/^\s+|\s+$//g;
                    if (defined($al{$input1})) {
                        $newick .= $al{$input1}; #save newick removing all new line characters
                    } else {
                        $newick .= $input1; #save newick removing all new line characters
                    }
                } else {
                    my ($input1) = $input =~ /^[0-9]+_(.*)/; #remove number added by mafft
                    $input1 =~ s/^\s+|\s+$//g;
                    if (defined($al{$input1})) {
                        $newick .= $al{$input1}; #save newick removing all new line characters
                    } else {
                        $newick .= $input1; #save newick removing all new line characters
                    }
                }
            } else {
                $input =~ s/^\s+|\s+$//g;
                $newick .= $input; #save newick removing all new line characters
            }
        }
        $newick .= ";";
        close(TREE);
        
        #format newick tree
        my %gene;
        my %taxid_mysql;
        my $proIN = 0;
        
        my $filename = $upload_dir . "/" . (substr($folder, 0, 8)) . ".pro.txt";
        if (-e $filename) { #if file is in the folder
            $sub += 1;
            open(IN_PRO, "<$filename");
            while (defined(my $input = <IN_PRO>)) {
                chomp($input);
                my @info = split(/\t/, $input);
                $gene{$info[0]} = $info[1];
                $proIN = 1;
            }
            close(IN_PRO);
        }
        
        $filename = $upload_dir . "/" . (substr($folder, 0, 8)) . ".tax.txt";
        if (-e $filename) { #if file is in the folder
            $sub += 2;
            open(IN_TAX, "<$filename");
            while (defined(my $input = <IN_TAX>)) {
                chomp($input);
                my @info = split(/\t/, $input);
                $taxid_mysql{$info[1]}{$info[0]} = '';
            }
            close(IN_TAX);
        }
        
        my $entry = 0;
        my $ask= '';
        foreach my $taxid (keys %taxid_mysql) {
            if ($entry > 0) {
                $ask .= " OR ";
            }
            $ask .= "(taxid='" . $taxid . "')";
            $entry++;
        }
        
        
        my $newick_domain = $newick;
        my $newick_phylum = $newick;
        my $newick_class = $newick;
        my $newick_order = $newick;
        my $newick_family = $newick;
        my $newick_genus = $newick;
        my $newick_species = $newick;
        
        
        if ($ask ne '') {
            
            #format tree
            #connect to mysql database
            my $dbh;
            my $sth;
            
            $dbh = DBI->connect ($dsn, $user_name, $password, { RaiseError => 1 });
            
            $sth = $dbh->prepare("SELECT * FROM taxid_taxonomy_pp WHERE ($ask)");
            
            $sth->execute(); #execute the prepared statement handle:
            
            my $composed;
            while (my @ary = $sth->fetchrow_array()) { #read results of a query, then clean up
                foreach my $acc (keys %{$taxid_mysql{$ary[0]}}) { #taxid - acc
                    if ($newick_domain =~ /$acc:/) {
                        #domain
                        $composed = $acc . ";" . $ary[1] . ";" . $gene{$acc} . ":";
                        $newick_domain =~ s/$acc:/$composed/g;
                        #phylum
                        $composed = $acc . ";" . $ary[2] . ";" . $gene{$acc} . ":";
                        $newick_phylum =~ s/$acc:/$composed/g;
                        #class
                        $composed = $acc . ";" . $ary[3] . ";" . $gene{$acc} . ":";
                        $newick_class =~ s/$acc:/$composed/g;
                        #order
                        $composed = $acc . ";" . $ary[4] . ";" . $gene{$acc} . ":";
                        $newick_order =~ s/$acc:/$composed/g;
                        #family
                        $composed = $acc . ";" . $ary[5] . ";" . $gene{$acc} . ":";
                        $newick_family =~ s/$acc:/$composed/g;
                        #genus
                        $composed = $acc . ";" . $ary[6] . ";" . $gene{$acc} . ":";
                        $newick_genus =~ s/$acc:/$composed/g;
                        #species
                        $composed = $acc . ";" . $ary[7] . ";" . $gene{$acc} . ":";
                        $newick_species =~ s/$acc:/$composed/g;
                    }
                    
                }
            }
            $sth->finish;
            
        } elsif ($proIN == 1) {
            my $composed;
            foreach my $acc (keys %gene) { #taxid - acc
                #domain
                $composed = $acc . ";;" . $gene{$acc} . ":";
                $newick_domain =~ s/$acc:/$composed/g;
                #phylum
                $newick_phylum =~ s/$acc:/$composed/g;
                #class
                $newick_class =~ s/$acc:/$composed/g;
                #order
                $newick_order =~ s/$acc:/$composed/g;
                #family
                $newick_family =~ s/$acc:/$composed/g;
                #genus
                $newick_genus =~ s/$acc:/$composed/g;
                #species
                $newick_species =~ s/$acc:/$composed/g;
                
            }
        }
        
        #pass number of accession numbers
        $count_acc = scalar @acc_mysql;
        
        #pass all phylogetic trees to html
        $newhtml = $newick_domain . "&" . $newick_phylum . "&" . $newick_class . "&" . $newick_order . "&" . $newick_family . "&" . $newick_genus . "&" . $newick_species;
        
    } elsif ($input_tree eq 'newick')  {
        chop($defSet);
        $folder = $defSet . "Tc";
        $upload_dir = $path_html . "/uploadsPhyloprimer/" . $folder; #path to the folder
        my $in_file = (substr($folder, 0, 8)) . ".tree";
        
        open(IN_NEWICK, "<${upload_dir}/${in_file}");
        my @all;
        while (defined(my $input = <IN_NEWICK>)) {
            chomp($input);
            @all = split(/:/,$input); #divide newick tree where :
            $input =~ s/^\s+|\s+$//g;
            $newick = $input;
        }
        close(IN_NEWICK);
        
        #count how many acc numbers - and user table construction
        my $count = 0;
        foreach my $all (@all) {
            $count++;
            if ($count == 1) { #start tree
                if ($all =~ /^\(/) {
                    my ($acc) = $all =~ /[\(]+(.*)$/;
                    push @acc_mysql, $acc; #save all accession
                    if ($acc =~ / -userHeader- /) { #get new and old name for user table
                        my @new = split(/ -userHeader- /, $acc);
                        #substitute accession special characters
                        $new[1] =~ s/_/ /g;#all spaces to _
                        $new[1] =~ s/^\s+|\s+$//g; #taxon name
                        $new[1] =~ s/\[sub1\]/\'/g; #single quote
                        $new[1] =~ s/\[sub2\]/,/g; #comma
                        $new[1] =~ s/\[sub3\]/\(/g; #bracket (
                        $new[1] =~ s/\[sub4\]/\)/g; #bracket )
                        $new[1] =~ s/\[sub5\]/:/g; #column
                        $new[1] =~ s/\[sub6\]/;/g; #semi column
                        $new[1] =~ s/\[sub7\]/\*/g; #semi column
                        $new[1] =~ s/\[sub8\]/</g; #lower
                        $new[1] =~ s/\[sub9\]/>/g; #higher
                        $new[1] =~ s/\[sub10\]/-/g; #minus
                        $new[1] =~ s/\[sub11\]/\+/g; #plus
                        $new[1] =~ s/\[sub12\]/\`/g; #hyphen`
                        $new[1] =~ s/\[sub13\]/\#/g; #
                        $new[1] =~ s/\[sub14\]/&/g; #&
                        $new[1] =~ s/\[sub15\]/\^/g; #&
                        $new[1] =~ s/\[sub16\]/\//g; #/
                        $new[1] =~ s/\[sub17\]/_/g; #_
                        $userTable_order{"<tr><td>" . $new[0] . "</td><td>" . $new[1] . "</td></tr>"} = '';
                    }
                }
            } else { #middle of the tree
                if ($all =~ /^([0-9\.]+),/) { #when the branch length is folled by , -> accession number
                    if ($all =~ /\(/) { #,((
                        my ($acc) = $all =~ /[\(]+(.*)$/;
                        push @acc_mysql, $acc; #save all accession
                        if ($acc =~ / -userHeader- /) { #get new and old name for user table
                            my @new = split(/ -userHeader- /, $acc);
                            #substitute accession special characters
                            $new[1] =~ s/_/ /g;#all spaces to _
                            $new[1] =~ s/^\s+|\s+$//g; #taxon name
                            $new[1] =~ s/\[sub1\]/\'/g; #single quote
                            $new[1] =~ s/\[sub2\]/,/g; #comma
                            $new[1] =~ s/\[sub3\]/\(/g; #bracket (
                            $new[1] =~ s/\[sub4\]/\)/g; #bracket )
                            $new[1] =~ s/\[sub5\]/:/g; #column
                            $new[1] =~ s/\[sub6\]/;/g; #semi column
                            $new[1] =~ s/\[sub7\]/\*/g; #semi column
                            $new[1] =~ s/\[sub8\]/</g; #lower
                            $new[1] =~ s/\[sub9\]/>/g; #higher
                            $new[1] =~ s/\[sub10\]/-/g; #minus
                            $new[1] =~ s/\[sub11\]/\+/g; #plus
                            $new[1] =~ s/\[sub12\]/\`/g; #hyphen`
                            $new[1] =~ s/\[sub13\]/\#/g; #
                            $new[1] =~ s/\[sub14\]/&/g; #&
                            $new[1] =~ s/\[sub15\]/\^/g; #&
                            $new[1] =~ s/\[sub16\]/\//g; #/
                            $new[1] =~ s/\[sub17\]/_/g; #_
                            $userTable_order{"<tr><td>" . $new[0] . "</td><td>" . $new[1] . "</td></tr>"} = '';
                        }
                    } else { #,
                        my ($acc) = $all =~ /,(.*)/;
                        push @acc_mysql, $acc; #save all accession
                        if ($acc =~ / -userHeader- /) { #get new and old name for user table
                            my @new = split(/ -userHeader- /, $acc);
                            #substitute accession special characters
                            $new[1] =~ s/_/ /g;#all spaces to _
                            $new[1] =~ s/^\s+|\s+$//g; #taxon name
                            $new[1] =~ s/\[sub1\]/\'/g; #single quote
                            $new[1] =~ s/\[sub2\]/,/g; #comma
                            $new[1] =~ s/\[sub3\]/\(/g; #bracket (
                            $new[1] =~ s/\[sub4\]/\)/g; #bracket )
                            $new[1] =~ s/\[sub5\]/:/g; #column
                            $new[1] =~ s/\[sub6\]/;/g; #semi column
                            $new[1] =~ s/\[sub7\]/\*/g; #semi column
                            $new[1] =~ s/\[sub8\]/</g; #lower
                            $new[1] =~ s/\[sub9\]/>/g; #higher
                            $new[1] =~ s/\[sub10\]/-/g; #minus
                            $new[1] =~ s/\[sub11\]/\+/g; #plus
                            $new[1] =~ s/\[sub12\]/\`/g; #hyphen`
                            $new[1] =~ s/\[sub13\]/\#/g; #
                            $new[1] =~ s/\[sub14\]/&/g; #&
                            $new[1] =~ s/\[sub15\]/\^/g; #&
                            $new[1] =~ s/\[sub16\]/\//g; #/
                            $new[1] =~ s/\[sub17\]/_/g; #_
                            $userTable_order{"<tr><td>" . $new[0] . "</td><td>" . $new[1] . "</td></tr>"} = '';
                        }
                    }
                }
            }
        }
        
        $newick =~ s/ -userHeader- [A-Za-z0-9_\-\.]+//g;
        
        foreach my $u (sort keys %userTable_order) {
            $userTable .= $u;
        }
        
        ####format newick tree
        my %gene;
        my %taxid_mysql;
        my $proIN = 0;
        
        my $filename = $upload_dir . "/" . (substr($folder, 0, 8)) . ".pro.txt";
        if (-e $filename) { ###if file is in the folder
            $sub += 1;
            open(IN_PRO, "<$filename");
            while (defined(my $input = <IN_PRO>)) {
                chomp($input);
                my @info = split(/\t/, $input);
                $gene{$info[0]} = $info[1];
                $proIN = 1;
            }
            close(IN_PRO);
        }
        
        $filename = $upload_dir . "/" . (substr($folder, 0, 8)) . ".tax.txt";
        if (-e $filename) { #if file is in the folder
            $sub += 2;
            open(IN_TAX, "<$filename");
            while (defined(my $input = <IN_TAX>)) {
                chomp($input);
                my @info = split(/\t/, $input);
                $taxid_mysql{$info[1]}{$info[0]} = '';
            }
            close(IN_TAX);
        }
        
        my $entry = 0;
        my $ask= '';
        foreach my $taxid (keys %taxid_mysql) {
            if ($entry > 0) {
                $ask .= " OR ";
            }
            $ask .= "(taxid='" . $taxid . "')";
            $entry++;
        }
        
        my $newick_domain = $newick;
        my $newick_phylum = $newick;
        my $newick_class = $newick;
        my $newick_order = $newick;
        my $newick_family = $newick;
        my $newick_genus = $newick;
        my $newick_species = $newick;
        
        
        if ($ask ne '') {
            #format tree
            #connect to mysql database
            my $dbh;
            my $sth;
            
            $dbh = DBI->connect ($dsn, $user_name, $password, { RaiseError => 1 });
            
            $sth = $dbh->prepare("SELECT * FROM taxid_taxonomy_pp WHERE ($ask)");
            
            $sth->execute(); #execute the prepared statement handle:
            
            my $composed;
            while (my @ary = $sth->fetchrow_array()) { #read results of a query, then clean up
                foreach my $acc (keys %{$taxid_mysql{$ary[0]}}) { #taxid - acc
                    if ($newick_domain =~ /$acc:/) {
                        #domain
                        $composed = $acc . ";" . $ary[1] . ";" . $gene{$acc} . ":";
                        $newick_domain =~ s/$acc:/$composed/g;
                        #phylum
                        $composed = $acc . ";" . $ary[2] . ";" . $gene{$acc} . ":";
                        $newick_phylum =~ s/$acc:/$composed/g;
                        #class
                        $composed = $acc . ";" . $ary[3] . ";" . $gene{$acc} . ":";
                        $newick_class =~ s/$acc:/$composed/g;
                        #order
                        $composed = $acc . ";" . $ary[4] . ";" . $gene{$acc} . ":";
                        $newick_order =~ s/$acc:/$composed/g;
                        #family
                        $composed = $acc . ";" . $ary[5] . ";" . $gene{$acc} . ":";
                        $newick_family =~ s/$acc:/$composed/g;
                        #genus
                        $composed = $acc . ";" . $ary[6] . ";" . $gene{$acc} . ":";
                        $newick_genus =~ s/$acc:/$composed/g;
                        #species
                        $composed = $acc . ";" . $ary[7] . ";" . $gene{$acc} . ":";
                        $newick_species =~ s/$acc:/$composed/g;
                    }
                    
                }
            }
            $sth->finish;
        } elsif ($proIN == 1) {
            my $composed;
            foreach my $acc (keys %gene) { #taxid - acc
                #domain
                $composed = $acc . ";;" . $gene{$acc} . ":";
                $newick_domain =~ s/$acc:/$composed/g;
                #phylum
                $newick_phylum =~ s/$acc:/$composed/g;
                #class
                $newick_class =~ s/$acc:/$composed/g;
                #order
                $newick_order =~ s/$acc:/$composed/g;
                #family
                $newick_family =~ s/$acc:/$composed/g;
                #genus
                $newick_genus =~ s/$acc:/$composed/g;
                #species
                $newick_species =~ s/$acc:/$composed/g;
                
            }
        }
        
        #pass number of accession numbers
        $count_acc = scalar @acc_mysql;
        
        #pass all phylogetic trees to html
        $newhtml = $newick_domain . "&" . $newick_phylum . "&" . $newick_class . "&" . $newick_order . "&" . $newick_family . "&" . $newick_genus . "&" . $newick_species;
        
    } else {
        ###ERROR
        $newhtmlSen = "<h3>Unexpected error. Please reload the page.</h3>";
        $newhtml = "ERROR";
    }
} elsif ($badEntry == 1) { #if error with BLAST parameters
    $newhtml = "ERROR";
    $newhtmlSen = "<h3>Please check the BLAST parameters. Only numeric values are allowed and the identity and coverage percentage must be set higher than 60%.</h3>";
    $count_acc = 0;
} else {
    $newhtml = "ERROR";
    $newhtmlSen = "<h3>Unexpected error. Please reload the page.</h3>";
}

if (($badEntry == 0) && ($noBLAST == 0)) {
    if ($input_tree eq 'gene')  {
        if ($clusterTable =~ /<td>/) {
            $newhtmlSen = "<h3>The BLAST search retrieved $count_acc_tot DNA sequences with an identity percentage higher than $identity_blast\%, a coverage higher than $coverage_blast\% and an e-value lower than $evalue_blast. However, the phylogenetic tree was constructed with $count_acc unique sequences representing $count_representative GenBank entries. The tree shows the GenBank accession numbers and the taxonomy of each entry (which rank can be changed with the designated buttons). Your sequences have been renamed and marked in <strong style='color: red;'>red</strong>: plase find the correspondence with the new assigned names in the table below. You can also hover on the nodes to get information about the GenBank sequence and some BLAST specifics for that alignment. You can select and deselect different sequences by clicking on the tree nodes. Once the sequences of interest have been selected, please click on the Create Consensus button and go to the next page.</h3>";
            $taxonomy = 'TRUE';
        } else {
            $newhtmlSen = "<h3>The BLAST search retrieved $count_acc_tot DNA sequences with an identity percentage higher than $identity_blast\%, a coverage higher than $coverage_blast\% and an e-value lower than $evalue_blast. The tree shows the GenBank accession numbers and the taxonomy of each entry (which rank can be changed with the designated buttons). Your sequences have been renamed and marked in <strong style='color: red;'>red</strong>: please find the correspondence with the new assigned names in the table below.You can also hover on the nodes to get information about the GenBank sequence and some BLAST specifics for that alignment. You can select and deselect different sequences by clicking on the tree nodes. Once the sequences of interest have been selected, please click on the Create Consensus button and go to the next page.</h3>";
            $taxonomy = 'TRUE';
        }
        
    } elsif ($input_tree eq 'sequence')  {
        if ($sub == 1) { #only prot
            if ($userTable =~ /<td>/) {
                $newhtmlSen = "<h3>The tree reports all the $count_acc entries you uploaded. The tree shows the sequence accessions. The accessions with a length of more than 20 characters were substituted by shorter sequence names assigned by PhyloPrimer and marked in <strong style='color: red;'>red</strong> in the tree, please find the correspondence with the new assigned names in the table below.You can hover on the nodes to get gene/protein for each entry. You can select and deselect different sequences by clicking on the tree nodes. Once the sequences of interest have been selected, please click on the Create Consensus button and go to the next page.</h3>";
            } else {
                $newhtmlSen = "<h3>The tree reports all the $count_acc entries you uploaded. The tree shows the sequence accessions. You can hover on the nodes to get gene/protein for each entry. You can select and deselect different sequences by clicking on the tree nodes. Once the sequences of interest have been selected, please click on the Create Consensus button and go to the next page.</h3>";
            }
            $taxonomy = 'FALSE';
        } elsif ($sub == 2) { #only tax
            if ($userTable =~ /<td>/) {
                $newhtmlSen = "<h3>The tree reports all the $count_acc entries you uploaded. The tree shows the sequence accessions and the taxonomy of each entry (which rank can be changed with the designated buttons). The accessions with a length of more than 20 characters were substituted by shorter sequence names assigned by PhyloPrimer and marked in <strong style='color: red;'>red</strong> in the tree, please find the correspondence with the new assigned names in the table below. You can select and deselect different sequences by clicking on the tree nodes. Once the sequences of interest have been selected, please click on the Create Consensus button and go to the next page.</h3>";
            } else {
                $newhtmlSen = "<h3>The tree reports all the $count_acc entries you uploaded. The tree shows the sequence accessions and the taxonomy of each entry (which rank can be changed with the designated buttons). You can select and deselect different sequences by clicking on the tree nodes. Once the sequences of interest have been selected, please click on the Create Consensus button and go to the next page.</h3>";
            }
            $taxonomy = 'TRUE';
        } elsif ($sub == 3) { #both pro and tax
            if ($userTable =~ /<td>/) {
                $newhtmlSen = "<h3>The tree reports all the $count_acc entries you uploaded. The tree shows the sequence accessions and the taxonomy of each entry (which rank can be changed with the designated buttons). The accessions with a length of more than 20 characters were substituted by shorter sequence names assigned by PhyloPrimer and marked in <strong style='color: red;'>red</strong> in the tree, please find the correspondence with the new assigned names in the table below. You can also hover on the nodes to get gene/protein for each entry. You can select and deselect different sequences by clicking on the tree nodes. Once the sequences of interest have been selected, please click on the Create Consensus button and go to the next page. </h3>";
            } else {
                $newhtmlSen = "<h3>The tree reports all the $count_acc entries you uploaded. The tree shows the sequence accessions and the taxonomy of each entry (which rank can be changed with the designated buttons). You can also hover on the nodes to get gene/protein for each entry. You can select and deselect different sequences by clicking on the tree nodes. Once the sequences of interest have been selected, please click on the Create Consensus button and go to the next page. </h3>";
            }
            $taxonomy = 'TRUE';
        } else { #none
            if ($userTable =~ /<td>/) {
                $newhtmlSen = "<h3>The tree reports all the $count_acc entries you uploaded. The tree shows the sequence accessions. The accessions with a length of more than 20 characters were substituted by shorter sequence names assigned by PhyloPrimer and marked in <strong style='color: red;'>red</strong> in the tree, please find the correspondence with the new assigned names in the table below. You can select and deselect different sequences by clicking on the tree nodes. Once the sequences of interest have been selected, please click on the Create Consensus button and go to the next page.</h3>";
            } else {
                $newhtmlSen = "<h3>The tree reports all the $count_acc entries you uploaded. The tree shows the sequence accessions. You can select and deselect different sequences by clicking on the tree nodes. Once the sequences of interest have been selected, please click on the Create Consensus button and go to the next page.</h3>";
            }
            $taxonomy = 'FALSE';
        }
    } elsif ($input_tree eq 'newick')  {
        if ($sub == 1) { #only prot
            if ($userTable =~ /<td>/) {
                $newhtmlSen = "<h3>The tree reports all the $count_acc entries you uploaded. The tree shows the sequence accessions. The accessions with a length of more than 20 characters were substituted by shorter sequence names assigned by PhyloPrimer and marked in <strong style='color: red;'>red</strong> in the tree, please find the correspondence with the new assigned names in the table below. You can hover on the nodes to get gene/protein for each entry. You can select and deselect different sequences by clicking on the tree nodes. Once the sequences of interest have been selected, please click on the Create Consensus button and go to the next page.</h3>";
            } else {
                $newhtmlSen = "<h3>The tree reports all the $count_acc entries you uploaded. The tree shows the sequence accessions. You can hover on the nodes to get gene/protein for each entry. You can select and deselect different sequences by clicking on the tree nodes. Once the sequences of interest have been selected, please click on the Create Consensus button and go to the next page.</h3>";
            }
            $taxonomy = 'FALSE';
        } elsif ($sub == 2) { #only tax
            if ($userTable =~ /<td>/) {
                $newhtmlSen = "<h3>The tree reports all the $count_acc entries you uploaded. The tree shows the sequence accessions and the taxonomy of each entry (which rank can be changed with the designated buttons). The accessions with a length of more than 20 characters were substituted by shorter sequence names assigned by PhyloPrimer and marked in <strong style='color: red;'>red</strong> in the tree, please find the correspondence with the new assigned names in the table below. You can select and deselect different sequences by clicking on the tree nodes. Once the sequences of interest have been selected, please click on the Create Consensus button and go to the next page.</h3>";
            } else {
                $newhtmlSen = "<h3>The tree reports all the $count_acc entries you uploaded. The tree shows the sequence accessions and the taxonomy of each entry (which rank can be changed with the designated buttons). You can select and deselect different sequences by clicking on the tree nodes. Once the sequences of interest have been selected, please click on the Create Consensus button and go to the next page.</h3>";
            }
            $taxonomy = 'TRUE';
        } elsif ($sub == 3) { #both pro and tax
            if ($userTable =~ /<td>/) {
                $newhtmlSen = "<h3>The tree reports all the $count_acc entries you uploaded. The tree shows the sequence accessions and the taxonomy of each entry (which rank can be changed with the designated buttons). The accessions with a length of more than 20 characters were substituted by shorter sequence names assigned by PhyloPrimer and marked in <strong style='color: red;'>red</strong> in the tree, please find the correspondence with the new assigned names in the table below. You can also hover on the nodes to get gene/protein for each entry. You can select and deselect different sequences by clicking on the tree nodes. Once the sequences of interest have been selected, please click on the Create Consensus button and go to the next page. </h3>";
            } else {
                $newhtmlSen = "<h3>The tree reports all the $count_acc entries you uploaded. The tree shows the sequence accessions and the taxonomy of each entry (which rank can be changed with the designated buttons). You can also hover on the nodes to get gene/protein for each entry. You can select and deselect different sequences by clicking on the tree nodes. Once the sequences of interest have been selected, please click on the Create Consensus button and go to the next page. </h3>";
            }
            $taxonomy = 'TRUE';
        } else { #none
            if ($userTable =~ /<td>/) {
                $newhtmlSen = "<h3>The tree reports all the $count_acc entries you uploaded. The tree shows the sequence accessions. The accessions with a length of more than 20 characters were substituted by shorter sequence names assigned by PhyloPrimer and marked in <strong style='color: red;'>red</strong> in the tree, please find the correspondence with the new assigned names in the table below. You can select and deselect different sequences by clicking on the tree nodes. Once the sequences of interest have been selected, please click on the Create Consensus button and go to the next page.</h3>";
            } else {
                $newhtmlSen = "<h3>The tree reports all the $count_acc entries you uploaded. The tree shows the sequence accessions. You can select and deselect different sequences by clicking on the tree nodes. Once the sequences of interest have been selected, please click on the Create Consensus button and go to the next page.</h3>";
            }
            $taxonomy = 'FALSE';
        }
    }
}

$userTable = $userTable . "</table>";
$clusterTable = $clusterTable . "</table>";

#result array
my @result;

if ($badEntry != 0) {
    push @result, "intro";
    push @result, $newhtmlSen;
    
    push @result, "result";
    push @result, $newhtml;
    
} else {
    
    my $out_file = $upload_dir . "/" . (substr($folder, 0, 8)) . ".treeInfo";
    open(my $fileI, ">${out_file}") or die "Could not open file '${out_file}' $!";
    
    #add description to user sequences
    foreach my $accession (keys %user) {
        my $composed1 = "<<" . $accession;
        my $composed2 = $composed1 . ";;" . $user{$accession};
        $newhtml =~ s/$composed1/$composed2/g;
    }

    push @result, "intro";
    push @result, $newhtmlSen;
    print $fileI "intro\t$newhtmlSen\n";
    
    push @result, "result";
    push @result, $newhtml;
    print $fileI "result\t$newhtml\n";
    
    push @result, "number";
    push @result, $count_acc;
    print $fileI "number\t$count_acc\n";
    
    push @result, "userTable";
    push @result, $userTable;
    print $fileI "userTable\t$userTable\n";
    
    push @result, "clusterTable";
    push @result, $clusterTable;
    print $fileI "clusterTable\t$clusterTable\n";
    
    #if PhyloPrimer needs to show taxonomy buttons
    push @result, "taxonomy";
    push @result, $taxonomy;
    print $fileI "taxonomy\t$taxonomy\n";
    
    close($fileI);
    
}

my $json = $op -> encode({
    @result
});
print $json;
