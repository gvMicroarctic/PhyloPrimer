#!/usr/bin/perl

#use strict;
use warnings;

use JSON; #if not already installed, just run "cpan JSON"
use CGI;
use File::Basename;

#connected to primer_input.html : check parameters set in Primer design page

#set path to folders:
my $path_html = 'path_to_html_folder';

#Set maximum size of uploaded file
$CGI::POST_MAX = 1024 * 50000; #50 MegaBytes

my $cgi = CGI->new;

print $cgi->header('application/json;charset=UTF-8');

##Section 0
my $project = $cgi->param('project'); #word
my $email = $cgi->param('email'); #word - needs to have
my $mailing_list = $cgi->param('mailing_list'); #user subscribed to mailing list

##Section 1
my $use = $cgi->param('use'); #how to use the script
my $folder = $cgi->param('folder'); #how to use the script
my $consensus = $cgi->param('consensus');
my $upload_dir = $path_html . "/uploadsPhyloprimer/" . $folder; #path to the folder
my $ST_input = (substr($folder, -2, 1));

##Section 2
my $plength_min = $cgi->param('plength_min'); #n - 10- 50
my $plength_max = $cgi->param('plength_max'); #n - 10 -50
my $tm_min = $cgi->param('tm_min'); #n
my $tm_max = $cgi->param('tm_max'); #n
my $maxrun = $cgi->param('maxrun'); #n
my $maxrepeat = $cgi->param('maxrepeat'); #n
my $gcclamp_min = $cgi->param('gcclamp_min'); #n
my $gcclamp_max = $cgi->param('gcclamp_max'); #n
my $mingc = $cgi->param('mingc'); #n
my $maxgc = $cgi->param('maxgc'); #n
my $numberwild_max = $cgi->param('numberwild_max');
my $hide_wild2 = $cgi->param('hide_wild2'); #true/""
my $hide_wild3 = $cgi->param('hide_wild3'); #true/""
my $hide_wild4 = $cgi->param('hide_wild4'); #true/""
my $end5_deg = $cgi->param('end5_deg'); #n
my $end3_deg = $cgi->param('end3_deg'); #n

##Section 3
my $amplicon_len_min = $cgi->param('amplicon_len_min');  #n - 0-1000
my $amplicon_len_max = $cgi->param('amplicon_len_max'); #n - 0-1000
my $tm_difference = $cgi->param('tm_difference'); #n
my $ta_min = $cgi->param('ta_min'); #n
my $ta_max = $cgi->param('ta_max'); #n

##Section 4
my $conserved = $cgi->param('conserved'); #yes/no
my $highF = $cgi->param('highF');
my $highR = $cgi->param('highR');
my $diffPos = $cgi->param('diffPos');

##Section 5
my $fasta_negative = $cgi->param('fasta_negative'); #fasta file
my $spec_mis = $cgi->param('spec_mis'); #positive integer
my $spec3_mis = $cgi->param('spec3_mis'); #positive integer
my $spec3_length = $cgi->param('spec3_length'); #positive integer

##Section 6
my $oli_hairpin = $cgi->param('oli_hairpin'); #n - negative
my $oli_self = $cgi->param('oli_self'); #n - negative
my $oli_cross = $cgi->param('oli_cross'); #n - negative
my $t_dG = $cgi->param('t_dG'); #n

##Section 7
my $mon_dG = $cgi->param('mon_dG'); #n
my $mg_dG = $cgi->param('mg_dG'); #n
my $oligo_dG = $cgi->param('oligo_dG'); #n
my $dNTP_dG = $cgi->param('dNTP_dG'); #n

##Section 8
my $species_sel = $cgi->param('species_sel'); #yes/no
if ($ST_input eq 'S') {
    $species_sel = '';
}
my $genus_sel = $cgi->param('genus_sel'); #yes/no
my $family_sel = $cgi->param('family_sel'); #yes/no
my $order_sel = $cgi->param('order_sel'); #yes/no
my $class_sel = $cgi->param('class_sel'); #yes/no
my $phylum_sel = $cgi->param('phylum_sel'); #yes/no
my $domain_sel = $cgi->param('domain_sel'); #yes/no

my $maximize_sel1 = $cgi->param('maximize_sel1'); #yes/no
my $maximize_sel2 = $cgi->param('maximize_sel2'); #yes/no
my $tm_sel = $cgi->param('tm_sel'); #yes/no
my $dG_sel = $cgi->param('dG_sel'); #yes/no
my $deg_sel = $cgi->param('deg_sel'); #yes/no

#color list for clusters
my %color;
$color{'1'} = 'DarkGreen';
$color{'2'} = 'DarkBlue';
$color{'3'} = 'DarkOrchid';
$color{'4'} = 'Orange';
$color{'5'} = 'SeaGreen';
$color{'6'} = 'OrangeRed';
$color{'7'} = 'Salmon';
$color{'8'} = 'YellowGreen';
$color{'9'} = 'Indigo';
$color{'10'} = 'Coral';
$color{'11'} = 'CadetBlue';
$color{'12'} = 'DarkGoldenRod';
$color{'13'} = 'DarkKhaki';
$color{'14'} = 'DeepPink';
$color{'15'} = 'DodgerBlue';
$color{'16'} = 'LightGreen';
$color{'17'} = 'MediumSlateBlue';
$color{'18'} = 'MediumVioletRed';
$color{'19'} = 'Peru';
$color{'20'} = 'SpringGreen';

#convert  data to JSON
my $op = JSON -> new -> utf8 -> pretty(1);

my $newhtml;
my $badEntry = 0;

my $highF_start;
my $highF_end;
my $highR_start;
my $highR_end;

#############################################
###############CHECK PARAMETERS##############
#############################################

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

##section 1
#re-check consensus

if ($consensus !~ /\A[acgtnryswkmbdhvn\-\.]+\z/i) {
    $badEntry = 1;
}

##Section 2
$newhtml = $newhtml . "<h4><b>PRIMERS</b></h4>";
#primer length
#only positive integers - ok
#larger than 50??? - to do
if (($plength_min !~ /\A[0-9]+\z/) or ($plength_max !~ /\A[0-9]+\z/)) {
    $newhtml = $newhtml . "<li><b style='color:red'\;>Primer length:</b> please check. Both values must be positive integers.</li>";
    $badEntry = 1;
} else {
    my $max_allow = 0; #?
    if ($plength_min < $plength_max) {
        $newhtml = $newhtml . "<li>Primer length: $plength_min - $plength_max bases.</li>";
        $max_allow = $plength_max;
    } elsif ($plength_min > $plength_max) {
        $newhtml = $newhtml . "<li><b style='color:red'\;>Primer length:</b> please check. Maximum primer size can not be lower than minimum oligo size.</li>";
        $badEntry = 1;
        $max_allow = $plength_min;
    } else {
        $newhtml = $newhtml . "<li>Primer length: $plength_min bases.</li>";
    }
}

$newhtml = $newhtml . "<br>";

#Melting temperature
#only positive integers - ok
if (($tm_min !~ /\A[0-9.]+\z/i) or ($tm_max !~ /\A[0-9.]+\z/i)) {
    $newhtml = $newhtml . "<li><b style='color:red'\;>T<sub>m</sub>:</b> please check. Both values must be positive numbers.</li>";
    $badEntry = 1;
} else {
    if ($tm_min < $tm_max) {
        $newhtml = $newhtml . "<li>T<sub>m</sub>: $tm_min - $tm_max &deg;C.</li>";
    } elsif ($tm_min > $tm_max) {
        $newhtml = $newhtml . "<li><b style='color:red'\;>T<sub>m</sub>:</b> please check. Maximum T<sub>m can not be lower than minimum T<sub>m.</li>";
        $badEntry = 1;
    } else {
        $newhtml = $newhtml . "<li>T<sub>m</sub>: $tm_min&deg;C.</li>";
    }
}

$newhtml = $newhtml . "<br>";

#Homopolymer length
#only positive integers - ok
if (($maxrun !~ /\A[0-9]+\z/) or ($maxrun < 2)) {
    $newhtml = $newhtml . "<li><b style='color:red'\;>Maximum homopolymer length:</b> please check. This value must be a positive integer higher than 1.</li>";
    $badEntry = 1;
} else {
    $newhtml = $newhtml . "<li>Maximum homopolymer length: $maxrun bases.</li>";
}

$newhtml = $newhtml . "<br>";

#Dinucleotide repeat length
#only positive integers - ok
if (($maxrepeat !~ /\A[0-9]+\z/) or ($maxrepeat < 4) or (($maxrepeat % 2) != 0)) {
    $newhtml = $newhtml . "<li><b style='color:red'\;>Maximum dinucleotide repeat length:</b> please check. This value must be a positive even integer higher than 3.</li>";
    $badEntry = 1;
} else {
    $newhtml = $newhtml . "<li>Maximum dinucleotide repeat length: $maxrepeat bases.</li>";
}

$newhtml = $newhtml . "<br>";

#GC clamp
#integer between 0-5
if (($gcclamp_min  !~ /\A[0-5]\z/i) or ($gcclamp_max  !~ /\A[0-5]\z/i)) {
    $newhtml = $newhtml . "<li><b style='color:red'\;>GC clamp:</b> please check. Both values must be positive intergers ranging between 0 and 5.</li>";
    $badEntry = 1;
} else {
    if ($gcclamp_min < $gcclamp_max) {
        $newhtml = $newhtml . "<li>GC clamp: $gcclamp_min - $gcclamp_max bases.</li>";
    } elsif ($gcclamp_min > $gcclamp_max) {
        $newhtml = $newhtml . "<li><b style='color:red'\;>GC clamp:</b> please check. Maximum GC clamp can not be lower than minimum GC clamp.</li>";
        $badEntry = 1;
    } else {
        if ($gcclamp_min == 1) {
            $newhtml = $newhtml . "<li>GC clamp: $gcclamp_min base.</li>";
        } else {
            $newhtml = $newhtml . "<li>GC clamp: $gcclamp_min bases.</li>";
        }
        
    }
}

$newhtml = $newhtml . "<br>";

#GC content
#integer between 0-100
if (($mingc  !~ /\A[0-9.]+\z/) or ($maxgc  !~ /\A[0-9.]+\z/) or ($mingc > 100) or ($maxgc > 100)) {
    $newhtml = $newhtml . "<li><b style='color:red'\;>GC content:</b> please check. Both values must be positive numbers ranging between 0 and 100.</li>";
    $badEntry = 1;
} else {
    if ($mingc < $maxgc) {
        $newhtml = $newhtml . "<li>GC content: $mingc - $maxgc%.</li>";
    } elsif ($mingc > $maxgc) {
        $newhtml = $newhtml . "<li><b style='color:red'\;>GC content:</b> please check. Maximum GC% can not be lower than minimum GC%.</li>";
        $badEntry = 1;
    } else {
        $newhtml = $newhtml . "<li>GC content: $mingc %.</li>";
    }
}

$newhtml = $newhtml . "<br>";

#wildcards
#Degenerate base numbers
#interger number
if ($numberwild_max !~ /\A[0-9]+\z/) {
    $newhtml = $newhtml . "<li><b style='color:red'\;>Maximum number of allowed degenerate bases:</b> please check. This value must be a positive integer.</li>";
    $badEntry = 1;
} else {
    
    if ($numberwild_max == 0) {
        $newhtml = $newhtml . "<li>Maximum number of allowed degenerate bases: no degenerate bases will be introduced.</li>";
    } else {
        if ($numberwild_max == 1) {
           $newhtml = $newhtml . "<li>Maximum number of allowed degenerate bases: $numberwild_max base.</li><ul>";
        } else {
           $newhtml = $newhtml . "<li>Maximum number of allowed degenerate bases: $numberwild_max bases.</li><ul>";
        }
        
        #Substitutions
        #true/false
        my @allow = ();
        if (($hide_wild2 eq "true") or ($hide_wild2 eq "") && ($hide_wild3 eq "true") or ($hide_wild3 eq "") && ($hide_wild4 eq "true") or ($hide_wild4 eq "")) {
            if ($hide_wild2 eq "true") {
                push @allow, ('R', 'Y', 'S', 'W', 'K', 'M');
            }
            if ($hide_wild3 eq "true") {
                push @allow, ('B', 'D', 'H', 'V');
            }
            if ($hide_wild4 eq "true") {
                push @allow, 'N';
            }
            $newhtml = $newhtml . "<li>Substitutions: " . join(',', @allow) . ".</li>";
        } else {
            $newhtml = $newhtml . "<li><b style='color:red'\;>Substitutions:</b> please check.</li>";
            $badEntry = 1;
        }
        
        #degenerate base at 5'end
        #integer between 0-5
        if ($end5_deg  !~ /\A[0-5]\z/i) {
            $newhtml = $newhtml . "<li><b style='color:red'\;>Avoid degenerate bases at the 5&rsquo; end for the first:</b> please check. This value must a positive integer ranging between 0-5.</li>";
        } else {
            $newhtml = $newhtml . "<li>Avoid degenerate bases at the 5&rsquo; end for the first: " . $end5_deg . " bases.</li>";
        }
        
        #degenerate base at 3'end
        #integer between 0-5
        if ($end3_deg  !~ /\A[0-5]\z/i) {
            $newhtml = $newhtml . "<li><b style='color:red'\;>Avoid degenerate bases at the 3&rsquo; end for the first:</b> please check. This value must a positive integer ranging between 0-5.</li>";
        } else {
            $newhtml = $newhtml . "<li>Avoid degenerate bases at the 3&rsquo; end for the first: " . $end3_deg . " bases.</li>";
        }
        
        $newhtml = $newhtml . "</ul>";
    }
}

$newhtml = $newhtml . "<br>";

##Section 3
$newhtml = $newhtml . "<h4><b>PRIMER PAIR</b></h4>";
#amplicon size
#only positive integers - ok
#cannot be longer than consensus - to do
if (($amplicon_len_min !~ /\A[0-9]+\z/) or ($amplicon_len_max !~ /\A[0-9]+\z/)) {
    $newhtml = $newhtml . "<li><b style='color:red'\;>Amplicon length:</b> please check. Both values need to be positive integers.</li>";
    $badEntry = 1;
} else {
    if ($amplicon_len_max > $amplicon_len_min) {
        $newhtml = $newhtml . "<li>Amplicon length: $amplicon_len_min - $amplicon_len_max bases.</li>";
    } elsif ($amplicon_len_min > $amplicon_len_max) {
        $newhtml = $newhtml . "<li><b style='color:red'\;>Amplicon length:</b> please check. Maximum amplicon size can not be lower than minimum amplicon size.</li>";
        $badEntry = 1;
    } elsif ($amplicon_len_min == $amplicon_len_max) {
        $newhtml = $newhtml . "<li>Amplicon length: $amplicon_len_min bases.</li>";
    }
}

$newhtml = $newhtml . "<br>";

#tm difference
#positive integers
if ($tm_difference !~ /\A[0-9.]+\z/i) {
    $newhtml = $newhtml . "<li><b style='color:red'\;>T<sub>m</sub> difference between forward and reverse primers:</b> please check. This value must be a positive number.</li>";
    $badEntry = 1;
} else {
    $newhtml = $newhtml . "<li>T<sub>m</sub> difference between forward and reverse primers: $tm_difference &deg;C.</li>";
}

$newhtml = $newhtml . "<br>";

#ta
if (($ta_min !~ /\A[0-9.]+\z/i) or ($ta_max !~ /\A[0-9.]+\z/i)) {
    $newhtml = $newhtml . "<li><b style='color:red'\;>T<sub>a</sub>:</b> please check. Both values must be positive numbers.</li>";
    $badEntry = 1;
} else {
    if ($ta_min < $ta_max) {
        $newhtml = $newhtml . "<li>T<sub>a</sub>: $ta_min - $ta_max &deg;C.</li>";
    } elsif ($ta_min > $ta_max) {
        $newhtml = $newhtml . "<li><b style='color:red'\;>T<sub>a</sub>:</b> please check. Maximum T<sub>m can not be lower than minimum T<sub>m</sub>.</li>";
        $badEntry = 1;
    } else {
        $newhtml = $newhtml . "<li>T<sub>a</sub>: $ta_min &deg;C.</li>";
    }
}

$newhtml = $newhtml . "<br>";


##section 4
$newhtml = $newhtml . "<h4><b>AREA SELECTION</b></h4>";
my $errorHigh = 0;
if ($conserved eq "yes") { #still to check
    #I need to check where $highF falls into consensus
    #it needs to be unique!
    my $consensus_len = length($consensus);
    my $inside = 0;
    $newhtml = $newhtml . "<li>Select specific sequence area where to design the primers: yes.</li><ul>";
    #forward area
    my $warningLenF = 0;
    my $warningLenF1 = 0;
    my $errorLenF1 = 0;
    my $warningWildF = 0;
    my $short = 0;
    #    my $minReq;
    #    my $minReqWarning;
    
    if ($highF ne "") {
        $inside += 1;
        if ($highF =~ /\A[acgtnryswkmbdhvn]+\z/i) {
            my $highF_len = length($highF);
            if ($highF_len < $plength_min) {
                $newhtml = $newhtml . "<li><b style='color:red'\;>Forward primer:</b> please check. The selected area is shorter than the minimum primer lenght.</li>";
                $short = 1;
            }
            if ($consensus =~ /$highF/) { #if highF is in consensus
                if ($consensus =~ /^$highF/) {
                    $highF_start = 1;
                    $highF_end = $highF_len;
                } elsif ($consensus =~ /$highF$/) {
                    $highF_start = $consensus_len - $highF_len + 1;
                    $highF_end = $consensus_len;
                } else {
                    my @checkF = split(/$highF/,$consensus);
                    my $split = scalar(@checkF);
                    if ($split == 2) {
                        $highF_start = length($checkF[0]) + 1;
                        $highF_end = length($checkF[0]) + $highF_len;
                    } else {
                        $newhtml = $newhtml . "<li><b style='color:red'\;>Forward primer:</b> please check. The selected area is present more than once in the consensus sequence, please change.</li>";
                    }
                }
            } else {
                $newhtml = $newhtml . "<li><b style='color:red'\;>Forward primer:</b> please check.</li>";
            }
            
            #warnings
            if ($short == 0) {
                if ($highF_len < 100) {
                    $warningLenF = 1;
                }
            }
            
            my $minReq = $consensus_len - $plength_min - $plength_min - $amplicon_len_min;
            my $minReqWarning = $consensus_len - $plength_min - $plength_min - $amplicon_len_min - 100;
            
            if ($highF_start > $minReq) {
                $errorLenF1 = 1;
            } elsif ($highF_start > $minReqWarning) {
                $warningLenF1 = 1;
            }
            
            
            if ($numberwild_max == 0) {
                if ($highF =~ /[nryswkmbdhvn]/i) {
                    my $wild = () = $highF =~ /[nryswkmbdhvn]/gi;
                    my $ratio = $wild/$highF_len;
                    if ($ratio >= 0.10) {
                        $warningWildF = 1;
                    }
                }
            }
        } else {
            $errorHigh += 1;
        }
    }
    #reverse area
    my $warningWildR = 0;
    my $warningLenR1 = 0;
    my $errorLenR1 = 0;
    my $warningLenR = 0;
    my $short = 0;
    if ($highR ne "") {
        $inside += 2;
        if ($highR =~ /\A[acgtnryswkmbdhvn]+\z/i) {
            my $highR_len = length($highR);
            if ($highR_len < $plength_min) {
                $newhtml = $newhtml . "<li><b style='color:red'\;>Reverse primer:</b> please check. The selected area is shorter than the minimum primer lenght.</li>";
                $short = 1;
            }
            if ($consensus =~ /$highR/) { #if highF is in consensus
                if ($consensus =~ /^$highR/) {
                    $highR_start = 1;
                    $highR_end = $highR_len;
                } elsif ($consensus =~ /$highR$/) {
                    $highR_start = $consensus_len - $highR_len + 1;
                    $highR_end = $consensus_len;
                } else {
                    my @checkR = split(/$highR/,$consensus);
                    my $split = scalar(@checkR);
                    if ($split == 2) {
                        $highR_start = length($checkR[0]) +1;
                        $highR_end = length($checkR[0]) + $highR_len;
                    } else {
                        $newhtml = $newhtml . "<li><b style='color:red'\;>Reverse primer:</b> please check. The selected area is present more than once in the consensus sequence, please change.</li>";
                    }
                }
            } else {
                $newhtml = $newhtml . "<li><b style='color:red'\;>Reverse primer:</b> please check.</li>";
            }
            
            #warnings
            if ($short == 0) {
                if ($highR_len < 100) {
                    $warningLenR = 1;
                }
            }
            
            my $minReq = $plength_min + $plength_min + $amplicon_len_min;
            my $minReqWarning = $plength_min + $plength_min + $amplicon_len_min + 100;
            
            if ($highR_start < $minReq) {
                $errorLenR1 = 1;
            } elsif ($highR_start < $minReqWarning) {
                $warningLenR1 = 1;
            }
            if ($numberwild_max == 0) {
                if ($highR =~ /[nryswkmbdhvn]/i) {
                    my $wild = () = $highR =~ /[nryswkmbdhvn]/gi;
                    my $ratio = $wild/$highR_len;
                    if ($ratio >= 0.10) {
                        $warningWildR = 1;
                    }
                }
            }
            
        } else {
            $errorHigh += 2;
        }
    }
    
    if ($inside == 0) { #no areas
        $newhtml = $newhtml . "<li><b style='color:red'\;>Select specific sequence area where to design the primers:</b> please check. You have not specified any area.</li>";
    } elsif ($inside == 1) { #both specified
        if ($errorHigh == 0) {
            if ($errorLenF1 == 0) {
                $newhtml = $newhtml . "<li>Forward primer: $highF_start-$highF_end base.</li>";
            } else {
                $newhtml = $newhtml . "<li><b style='color:red'\;>Forward primer:</b> please check. The forward selected area is too close to the 3' end of the consensus.</li>";
            }
            $newhtml = $newhtml . "<li>Reverse primer: no specified area.</li>";
        } else {
            $newhtml = $newhtml . "<li><b style='color:red'\;>Forward primer:</b> please check. The selection of only one area for the forward primer design is allowed.</li>";
            $newhtml = $newhtml . "<li>Reverse primer: no specified area.</li>";
        }
    } elsif ($inside == 2) {
        if ($errorHigh == 0) {
            $newhtml = $newhtml . "<li>Forward primer: no specified area.</li>";
            if ($errorLenR1 == 0) {
                $newhtml = $newhtml . "<li>Reverse primer: $highR_start-$highR_end base.</li>";
            } else {
                $newhtml = $newhtml . "<li><b style='color:red'\;>Reverse primer:</b> please check. The reverse selected area is too close to the 5' end of the consensus.</li>";
            }
        } else {
            $newhtml = $newhtml . "<li>Forward primer: no specified area.</li>";
            $newhtml = $newhtml . "<li><b style='color:red'\;>Reverse primer:</b> please check. The selection of only one area for the reverse primer design is allowed.</li>";
        }
    } else {
        if ($errorHigh == 0) {
            if ($highF_end > $highR_start) {
                $newhtml = $newhtml . "<li><b style='color:red'\;>Select specific sequence area where to design the primers:</b> please check. Forward primer needs to be before reverse primer on the consensus.</li>";
            } else {
                if ($amplicon_len_min > ($highR_end - $highF_start)) {
                    $newhtml = $newhtml . "<li><b style='color:red'\;>Select specific sequence area where to design the primers:</b> please check. The area between the two primers is smaller than min amplicon area.</li>";
                } elsif ($amplicon_len_max < (($highR_start + $plength_min) - ($highF_end -  $plength_min))) {
                    $newhtml = $newhtml . "<li><b style='color:red'\;>Select specific sequence area where to design the primers:</b> please check. The area between the two primers is longer than max amplicon area.</li>";
                } else {
                    $newhtml = $newhtml . "<li>Forward primer: $highF_start-$highF_end base.</li>";
                    $newhtml = $newhtml . "<li>Reverse primer: $highR_start-$highR_end base.</li>";
                }
            }
        } elsif ($errorHigh == 1) {
            $newhtml = $newhtml . "<li><b style='color:red'\;>Forward primer:</b> please check. The selection of only one area for the forward primer design is allowed.</li>";
            $newhtml = $newhtml . "<li>Reverse primer: $highR_start-$highR_end base.</li>";
        } elsif ($errorHigh == 2) {
            $newhtml = $newhtml . "<li>Forward primer: $highF_start-$highF_end base.</li>";
            $newhtml = $newhtml . "<li><b style='color:red'\;>Reverse primer:</b> please check. The selection of only one area for the reverse primer design is allowed.</li>";
        } elsif ($errorHigh == 3) {
            $newhtml = $newhtml . "<li><b style='color:red'\;>Forward primer:</b> please check. The selection of only one area for the forward primer design is allowed.</li>";
            $newhtml = $newhtml . "<li><b style='color:red'\;>Reverse primer:</b> please check. The selection of only one area for the reverse primer design is allowed.</li>";
        }
    }
    
    if ($warningLenF == 1) {
        $newhtml = $newhtml . "<li><b style='color:blue'\;>WARNING: Forward primer:</b> The forward selected area is smaller than 100 bp. Increasing the area would give your more chances to get suitable primers. </li>";
    }
    if ($warningLenF1 == 1) {
        $newhtml = $newhtml . "<li><b style='color:blue'\;>WARNING: Forward primer:</b> Beacuse of where you selected your reverse area, the forward area is smaller than 100 bp. Increasing the area would give your more chances to get suitable primers. </li>";
    }
    if ($warningWildF == 1) {
        $newhtml = $newhtml . "<li><b style='color:blue'\;>WARNING: Forward primer:</b> More than the 10% of the bases in the selected consensus are degenerate bases. Allowing the inclusion of degenerate bases would give your more chances to get suitable primers.</li>";
    }
    if ($warningLenR == 1) {
        $newhtml = $newhtml . "<li><b style='color:blue'\;>WARNING: Reverse primer:</b> The forward selected area is smaller than 100 bp. Increasing the area would give your more chances to get suitable primers. </li>";
    }
    if ($warningWildR == 1) {
        $newhtml = $newhtml . "<li><b style='color:blue'\;>WARNING: Reverse primer:</b> More than the 10% of the bases in the selected consensus are degenerate bases. Allowing the inclusion of degenerate bases would give your more chances to get suitable primers.</li>";
    }
    
    $newhtml = $newhtml . "</ul>";
    
} elsif ($conserved eq "no") {
    $newhtml = $newhtml . "<li>Select specific sequence area where to design the primers: no.</li>";
    $highF = '';
    $highR = '';
} else {
    $newhtml = $newhtml . "<li><b style='color:red'\;>Select specific sequence area where to design the primers:</b> please check.</li>";
    $highF = '';
    $highR = '';
}

if ($diffPos !~ /\A[0-9,\|]+\z/i) {
    $diffPos = 'no';
}

if ($conserved eq "no") {
    $newhtml = $newhtml . "<br>";
}

##Section 5
$newhtml = $newhtml . "<h4><b>BLAST CHECK</b></h4>";
#fasta negative
#check and load file
my $negative = "no";
if ($fasta_negative ne "") {  #if user uploaded a file
    if (!(defined($fasta_negative))) { #if the file is too big
        $newhtml = $newhtml . "<li><b style='color:red'\;>Negative fasta:</b> please check. The uploaded file is too big.</li>";
    } else { #if the file size is lower than the threshold
        my ( $name, $path, $extension ) = fileparse ( $fasta_negative, '..*' );
        $fasta_negative = $name . $extension;
        if ($fasta_negative =~ /^([a-zA-Z0-9_\.\-]+)$/) { #safe characters for file name
            my $upload_filehandle = $cgi->upload('fasta_negative');
            my $header = 0;
            my $count = 0;
            my $sequence = 0;
            my $badEntry = 0;
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
                    $accession =~ s/\s/_/g; #all spaces to _
                    $seq = "";
                } else {
                    if ($input ne "") {
                        $last = "sequence";
                        $count++;
                        if ($input !~ /\A[acgtnryswkmbdhvn]+\z/i) { #safe characters for sequence
                            $badEntry = 1;
                        } else {
                            $seq .= $input;
                        }
                        if ($count == 1) { #how many sequences
                            $sequence++;
                        }
                    }
                    if ($badEntry == 0) {
                        $print{$accession} = $seq; #accession = sequence
                    } else {
                        last;
                    }
                }
            }
            #tolerate if no header so only one sequence
            #BUT if more than one header, same number between headers and sequences
            my $filenameNew = (substr($folder, 0, 8)) . ".negativefasta";
            if (($badEntry == 0) && ($last eq "sequence") && (($header == $sequence) or (($header == 0) && ($sequence == 1)))) {
                $negative = "yes";
                if ($header > 1) {
                    $newhtml = $newhtml . "<li>Secondary BLAST check: you uploaded $header sequences.</li>";
                } else {
                    $newhtml = $newhtml . "<li>Secondary BLAST check: you uploaded $header sequence.</li>";
                }
                open (UPLOADFILE, ">$upload_dir/$filenameNew" ) or die "$!"; #print only if entry is okay
                binmode UPLOADFILE;
                foreach my $head (keys %print) {
                    print UPLOADFILE ">$head\n$print{$head}\n";
                }
                close(UPLOADFILE);
            } else {
                $newhtml = $newhtml . "<li><b style='color:red'\;>Secondary BLAST check:</b> please check. The file must be in fasta format. Please check the correct format <a id='format' href='https://www.cerealsdb.uk.net/cerealgenomics/phyloprimer/format_page.html#DNA_seq'>here</a>.</li>";
            }
        } else { #uncorrect file name
            $newhtml = $newhtml . "<li><b style='color:red'\;>Secondary BLAST check:</b> please check. The file name is uncorrect.</li>";
        }
    }
} else {
    $newhtml = $newhtml . "<li>Secondary BLAST check: you have not uploaded any fasta file for a secondary BLAST search.</li>";
}

$newhtml = $newhtml . "<br>";

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

##Section 6
$newhtml = $newhtml . "<h4><b>SECONDARY STRUCTURES</b></h4>";
#self dimer
#integer
if ($oli_self !~ /\A[0-9.-]+\z/) {
    $newhtml = $newhtml . "<li><b style='color:red'\;>Minimum self dimer &Delta;G:</b> please check. This value must be a number.</li>";
    $badEntry = 1;
} else {
    $newhtml = $newhtml . "<li>Minimum self dimer &Delta;G: $oli_self cal/mol.</li>";
}

$newhtml = $newhtml . "<br>";

#hairpin
#integer
if ($oli_hairpin  !~ /\A[0-9.-]+\z/) {
    $newhtml = $newhtml . "<li><b style='color:red'\;>Minimum hairpin &Delta;G:</b> please check. This value must be a number.</li>";
    $badEntry = 1;
} else {
    $newhtml = $newhtml . "<li>Minimum hairpin &Delta;G: $oli_hairpin cal/mol.</li>";
}

$newhtml = $newhtml . "<br>";

#primer dimer
if ($oli_cross !~ /\A[0-9.-]+\z/i) {
    $newhtml = $newhtml . "<li><b style='color:red'\;>Minimum cross dimer &Delta;G:</b> please check. This value must be a number.</li>";
    $badEntry = 1;
} else {
    $newhtml = $newhtml . "<li>Minimum cross dimer &Delta;G: $oli_cross cal/mol.</li>";
}

$newhtml = $newhtml . "<br>";

#temperature
if ($t_dG !~ /\A[0-9.]+\z/i) {
    $newhtml = $newhtml . "<li><b style='color:red'\;>Temperature:</b> please check. This value must be a positive number.</li>";
    $badEntry = 1;
} else {
    $newhtml = $newhtml . "<li>Temperature for secondary structure calculation: $t_dG &deg;C.</li>";
}

$newhtml = $newhtml . "<br>";

##Section 7
$newhtml = $newhtml . "<h4><b>PCR CONDITIONS</b></h4>";
#monovalent ion content
if ($mon_dG !~ /\A[0-9.]+\z/i) {
    $newhtml = $newhtml . "<li><b style='color:red'\;>Monovalent ion concentration:</b> please check. This value must be a positive number.</li>";
    $badEntry = 1;
} else {
    $newhtml = $newhtml . "<li>Monovalent ion concentration: $mon_dG mM.</li>";
}

$newhtml = $newhtml . "<br>";

#mg content
if ($mg_dG !~ /\A[0-9.]+\z/i) {
    $newhtml = $newhtml . "<li><b style='color:red'\;>Mg<sup>2+</sup> concentration:</b> please check. This value must be a positive number.</li>";
} else {
    $newhtml = $newhtml . "<li>Mg<sup>2+</sup> concentration: $mg_dG mM.</li>";
}

$newhtml = $newhtml . "<br>";

#oligo content
if ($oligo_dG !~ /\A[0-9.]+\z/i) {
    $newhtml = $newhtml . "<li><b style='color:red'\;>Oligo concentration:</b> please check. This value must be a positive number.</li>";
    $badEntry = 1;
} else {
    $newhtml = $newhtml . "<li>Oligo concentration: $oligo_dG uM.</li>";
}

$newhtml = $newhtml . "<br>";

#dNTP content
if ($dNTP_dG !~ /\A[0-9.]+\z/i) {
    $newhtml = $newhtml . "<li><b style='color:red'\;>dNTP concentration:</b> please check. This value must me a positive number.</li>";
    $badEntry = 1;
} else {
    $newhtml = $newhtml . "<li>dNTP concentration: $dNTP_dG uM.</li>";
}

$newhtml = $newhtml . "<br>";

##Section 8
$newhtml = $newhtml . "<h4><b>VISUALIZATION CRITERIA</b></h4>";

my $whichSel = 0;
my $whichSel_error = 0;
my @criteria;

if ($maximize_sel1 eq 'yes') {
    $newhtml = $newhtml . "<li>PhyloPrimer will assign scoring points if differing bases (between positive and negative consensus) fall into the primer sequence inside the 3' end (last 5 bases).</li><br>";
    push @criteria, "PhyloPrimer will assign scoring points if differing bases (between positive and negative consensus) fall into the primer sequence inside the 3' end (last 5 bases).";
    $whichSel++;
} elsif ($maximize_sel1 ne '') {
    $whichSel_error = 1;
}

if ($maximize_sel2 eq 'yes') {
    $newhtml = $newhtml . "<li>PhyloPrimer will assign scoring points if differing bases (between positive and negative consensus) fall into the primer sequence.</li><br>";
    push @criteria, "PhyloPrimer will assign scoring points if differing bases (between positive and negative consensus) fall into the primer sequence.";
    $whichSel++;
} elsif ($maximize_sel2 ne '') {
    $whichSel_error = 1;
}

if ($tm_sel eq 'yes') {
    $whichSel += 4;
    $newhtml = $newhtml . "<li>PhyloPrimer will assign scoring points to primer pairs that have a T<sub>m</sub> difference between the forward and the reverse primer lower than 1 &deg;C.</li><br>";
    push @criteria, "PhyloPrimer will assign scoring points to primer pairs that have a Tm difference between the forward and the reverse primer lower than 1 &deg;C.";
    $whichSel++;
} elsif ($tm_sel ne '') {
    $whichSel_error = 1;
}

if ($dG_sel eq 'yes') {
    $newhtml = $newhtml . "<li>PhyloPrimer will assign scoring points to primers and primer pairs that have &Delta;G values higher than -1 kcal/mol.</li><br>";
    push @criteria, "PhyloPrimer will assign scoring points to primers and primer pairs that have &Delta;G values higher than -1 kcal/mol.";
    $whichSel++;
} elsif ($dG_sel ne '') {
    $whichSel_error = 1;
}

if ($deg_sel eq 'yes') {
    if ($numberwild_max != 0) {
        $newhtml = $newhtml . "<li>PhyloPrimer will assign scoring points to primers that do not have degenerate bases in their sequences.</li><br>";
        push @criteria, "PhyloPrimer will assign scoring points to primers that do not have degenerate bases in their sequences.";
        $whichSel++;
    }
} elsif ($deg_sel ne '') {
    $whichSel_error = 1;
}

if ($species_sel eq 'yes') {
    $newhtml = $newhtml . "<li>PhyloPrimer will assign scoring points is the primers aligns to DNA sequences belonging to species that were selected from the dynamic tree.</li><br>";
    push @criteria, "PhyloPrimer will assign scoring points is the primers aligns to DNA sequences belonging to species that were selected from the dynamic tree.";
    $whichSel++;
} elsif ($species_sel ne '') {
    $whichSel_error = 1;
}

if ($genus_sel eq 'yes') {
    $newhtml = $newhtml . "<li>PhyloPrimer will assign scoring points is the primers aligns to DNA sequences belonging to genera that were selected from the dynamic tree.</li><br>";
    push @criteria, "PhyloPrimer will assign scoring points is the primers aligns to DNA sequences belonging to genera that were selected from the dynamic tree.";
    $whichSel++;
} elsif ($genus_sel ne '') {
    $whichSel_error = 1;
}

if ($family_sel eq 'yes') {
    $newhtml = $newhtml . "<li>PhyloPrimer will assign scoring points is the primers aligns to DNA sequences belonging to families that were selected from the dynamic tree.</li><br>";
    push @criteria, "PhyloPrimer will assign scoring points is the primers aligns to DNA sequences belonging to families that were selected from the dynamic tree.";
    $whichSel++;
} elsif ($family_sel ne '') {
    $whichSel_error = 1;
}

if ($order_sel eq 'yes') {
    $newhtml = $newhtml . "<li>PhyloPrimer will assign scoring points is the primers aligns to DNA sequences belonging to orders that were selected from the dynamic tree.</li><br>";
    push @criteria, "PhyloPrimer will assign scoring points is the primers aligns to DNA sequences belonging to orders that were selected from the dynamic tree.";
    $whichSel++;
} elsif ($order_sel ne '') {
    $whichSel_error = 1;
}

if ($class_sel eq 'yes') {
    $newhtml = $newhtml . "<li>PhyloPrimer will assign scoring points is the primers aligns to DNA sequences belonging to classes that were selected from the dynamic tree.</li><br>";
    push @criteria, "PhyloPrimer will assign scoring points is the primers aligns to DNA sequences belonging to classes that were selected from the dynamic tree.";
    $whichSel++;
} elsif ($class_sel ne '') {
    $whichSel_error = 1;
}

if ($phylum_sel eq 'yes') {
    $newhtml = $newhtml . "<li>PhyloPrimer will assign scoring points is the primers aligns to DNA sequences belonging to phyla that were selected from the dynamic tree.</li><br>";
    push @criteria, "PhyloPrimer will assign scoring points is the primers aligns to DNA sequences belonging to phyla that were selected from the dynamic tree.";
    $whichSel++;
} elsif ($phylum_sel ne '') {
    $whichSel_error = 1;
}

if ($domain_sel eq 'yes') {
    $newhtml = $newhtml . "<li>PhyloPrimer will assign scoring points is the primers aligns to DNA sequences belonging to domains that were selected from the dynamic tree.</li><br>";
    push @criteria, "PhyloPrimer will assign scoring points is the primers aligns to DNA sequences belonging to domains that were selected from the dynamic tree.";
    $whichSel++;
} elsif ($domain_sel ne '') {
    $whichSel_error = 1;
}

if ($whichSel == 0) {
    $newhtml = $newhtml . "<li>The primer pairs will be randomly selected for the result visualization.</li><br>";
    push @criteria, "The primer pairs will be randomly selected for the result visualization.";
}

if ($whichSel_error == 1) {
    $newhtml = $newhtml . "<li><b style='color:red'\;Visualization criteria:</b> there was an unexpected error in the visualization criteria. Please re-load the page.</li><br>";
}

$newhtml = $newhtml . "<br>";

if ($use eq "check") {
    $newhtml = $newhtml . "<b>If all the inputs are correct and you are happy with the inputs, please click on the submit button. If some input needs to be modified, please do it and then click the check button again.</b><br><br>";
} elsif ($use eq "oligo") {

    if ($badEntry == 0) {
        
        
        #Files to go back to Dynamic Tree Page
        my $upload_dirNew = $upload_dir;
        my $input_kind = chop($upload_dirNew); #get last letter and understand if a, b or c
        my $input_STC = chop($upload_dirNew); #get last letter and understand if T or S or C
        my $index = 1;
       
        if ($input_STC eq "T") {
            if ($upload_dirNew =~ /_/) {
                my $upload_dirNew0 = $upload_dirNew;
                ($upload_dirNew) = $upload_dirNew0 =~ /(.*)_/;
                ($index) = $upload_dirNew0 =~ /_(.*)/;
                $index++;
                my $subfix = (substr($folder, 0, 8));
                `cp ${upload_dir}/${subfix}.loadAcc ${upload_dir}/${subfix}.loadAccOld`;
                open(my $fileCl, ">${upload_dir}/${subfix}.loadAcc") or die;
                open(CLUSTER, "<${upload_dir}/${subfix}.allAccession") or die;
                my $accCount = 0;
                while (defined(my $input=<CLUSTER>)) {
                    chomp($input);
                    print $fileCl ",${input}($color{$index})";
                    $accCount++;
                }
                close(CLUSTER);
                open(PREV, "<${upload_dir}/${subfix}.loadAccOld") or die;
                while (defined(my $input=<PREV>)) {
                    chomp($input);
                    print $fileCl "$input\n";
                }
                close(PREV);
                print $fileCl "CLUSTER\t$index\t$project\t$color{$index}\t$accCount\n";
                close($fileCl);
            } else { #if no previous cluster file
                my $subfix = (substr($folder, 0, 8));
                open(CLUSTER, "<${upload_dir}/${subfix}.allAccession") or die;
                open(my $fileCl, ">${upload_dir}/${subfix}.loadAcc") or die;
                my $accCount = 0;
                while (defined(my $input=<CLUSTER>)) {
                    chomp($input);
                    print $fileCl ",${input}($color{$index})";
                    $accCount++;
                }
                print $fileCl "\nCLUSTER\t1\t$project\t$color{$index}\t$accCount";
                close(CLUSTER);
                close($fileCl);
                #`rm ${upload_dir}/${subfix}.allAccession`;
            }
            my $upload_dirNew1 = $upload_dirNew . "_" . $index . $input_STC . $input_kind;
            `cp -r $upload_dir $upload_dirNew1`;
        }
        
        #open folder from previous pages
        my $in_file = (substr($folder, 0, 8)) . ".info1";
        open(my $file, ">${upload_dir}/${in_file}");
        
        ##Section 0
        print $file "PROJECT\t$project\n";
        print $file "EMAIL\t$email\n";
        print $file "MAILING_LIST\t$mailing_list\n";
        
        ##Section 1
        print $file "CONSENSUS\t$consensus\n";
        
        ##Section 2
        print $file "LEN_MIN\t$plength_min\n";
        print $file "LEN_MAX\t$plength_max\n";
        print $file "TM_MIN\t$tm_min\n";
        print $file "TM_MAX\t$tm_max\n";
        print $file "MAXRUN\t$maxrun\n";
        print $file "MAXREPEAT\t$maxrepeat\n";
        print $file "GCCLAMP_MIN\t$gcclamp_min\n";
        print $file "GCCLAMP_MAX\t$gcclamp_max\n";
        print $file "GC_CONTENT_MIN\t$mingc\n";
        print $file "GC_CONTENT_MAX\t$maxgc\n";
        print $file "NUMBERWILD_MAX\t$numberwild_max\n";
        my @allow;
        if ($hide_wild2 eq "true") {
            push @allow, '2';
        }
        if ($hide_wild3 eq "true") {
            push @allow, '3';
        }
        if ($hide_wild4 eq "true") {
            push @allow, '4';
        }
        print $file "WHICHWILD\t", join(',', @allow), "\n";
        print $file "END5WILD\t$end5_deg\n";
        print $file "END3WILD\t$end3_deg\n";
        
        ##Section 3
        print $file "AMPL_SIZE_MIN\t$amplicon_len_min\n";
        print $file "AMPL_SIZE_MAX\t$amplicon_len_max\n";
        print $file "TM_DIFF\t$tm_difference\n";
        print $file "TA_MIN\t$ta_min\n";
        print $file "TA_MAX\t$ta_max\n";
        
        ##Section 4
        print $file "CONSERVED\t$conserved\n";
        print $file "HIGHF\t$highF\n";
        print $file "HIGHR\t$highR\n";
        print $file "HIGHF_START\t$highF_start\n";
        print $file "HIGHF_END\t$highF_end\n";
        print $file "HIGHR_START\t$highR_start\n";
        print $file "HIGHR_END\t$highR_end\n";
        print $file "DIFFERENT_POS\t$diffPos\n";
        
        ##Section 5
        print $file "NEGATIVE_FILE\t$negative\n";
        #print $file "SPEC_MIS\t$spec_mis\n";
        #print $file "SPEC3_MIS\t$spec3_mis\n";
        #print $file "SPEC3_LENGTH\t$spec3_length\n";
        
        ##Section 6
        print $file "OLI_SELF\t$oli_self\n";
        print $file "OLI_HAIRPIN\t$oli_hairpin\n";
        print $file "OLI_CROSS\t$oli_cross\n";
        print $file "T_DG\t$t_dG\n";
        
        ##Section 7
        print $file "MON_DG\t$mon_dG\n";
        print $file "MG_DG\t$mg_dG\n";
        print $file "OLIGO_DG\t$oligo_dG\n";
        print $file "DNTP_DG\t$dNTP_dG\n";
        
        ##Section 8
        print $file "SPECIES_SEL\t$species_sel\n";
        print $file "GENUS_SEL\t$genus_sel\n";
        print $file "FAMILY_SEL\t$family_sel\n";
        print $file "ORDER_SEL\t$order_sel\n";
        print $file "CLASS_SEL\t$class_sel\n";
        print $file "PHYLUM_SEL\t$phylum_sel\n";
        print $file "DOMAIN_SEL\t$domain_sel\n";
        
        print $file "MAXIMIZE_SEL1\t$maximize_sel1\n";
        print $file "MAXIMIZE_SEL2\t$maximize_sel2\n";
        print $file "TM_SEL\t$tm_sel\n";
        print $file "DG_SEL\t$dG_sel\n";
        print $file "DEG_SEL\t$deg_sel\n";
        
        close($file);
        
        #open folder from previous pages - user
        my $in_file = (substr($folder, 0, 8)) . ".info";
        
        open(my $file, ">${upload_dir}/${in_file}");
        
        ##Section 0
        print $file "\nUSER DATA\n";
        print $file "Job ID: $project\n";
        print $file "Email: $email\n";
        print $file "Mailing list: $mailing_list\n";
        
        ##Section 1
        print $file "Consensus: $consensus\n";
        
        print $file "\nPhyloPrimer mode: Primer pair design\n";
        
        ##Section 2
        print $file "\nPRIMERS\n";
        print $file "Primer length - minimum: $plength_min bases\n";
        print $file "Primer length - maximum: $plength_max bases\n";
        print $file "Melting temperature - minimum: $tm_min °C\n";
        print $file "Melting temperature - maximum: $tm_max °C\n";
        print $file "Homopolymer length - maximum: $maxrun bases\n";
        print $file "Dinucleotide repeat length - maximum: $maxrepeat bases\n";
        print $file "GC clamp - minimum: $gcclamp_min bases\n";
        print $file "GC clamp - maximum: $gcclamp_max bases\n";
        print $file "GC content - minimum: $mingc %\n";
        print $file "GC content - maximum: $maxgc %\n";
        print $file "Number of degenerate bases - maximum: $numberwild_max bases\n";
        my @allow;
        if ($hide_wild2 eq "true") {
            push @allow, '2';
        }
        if ($hide_wild3 eq "true") {
            push @allow, '3';
        }
        if ($hide_wild4 eq "true") {
            push @allow, '4';
        }
        print $file "Substitutions: ", join(',', @allow), "\n";
        print $file "Avoid degenerate bases at the 5' end for the first: $end5_deg\n";
        print $file "Avoid degenerate bases at the 3' end for the first: $end3_deg\n";
        
        ##Section 3
        print $file "\nPRIMER PAIR\n";
        print $file "Amplicon length - minimum: $amplicon_len_min bases\n";
        print $file "Amplicon length - maximum: $amplicon_len_max bases\n";
        print $file "Melting temperature difference between forward and reverse primers - maximum: $tm_difference °C\n";
        print $file "Annealing temperature - minimum: $ta_min °C\n";
        print $file "Annealing temperature - maximum: $ta_max °C\n";
        
        ##Section 4
        print $file "\nAREA SELECTION\n";
        print $file "Select specific sequence area where to design the primers: $conserved\n";
        print $file "Forward area: $highF\n";
        print $file "Reverse area: $highR\n";
        print $file "Forward area start position: $highF_start\n";
        print $file "Forward area end position: $highF_end\n";
        print $file "Reverse area start position: $highR_start\n";
        print $file "Reverse area end position: $highR_end\n";
        print $file "Differing positions between positive and negative consensus: $diffPos\n";
        
        ##Section 5
        print $file "\nBLAST CHECK\n";
        print $file "Secondary BLAST check: $negative\n";
        #print $file "Maximum number of mismatches in the entire oligo sequence: $spec_mis\n";
        #print $file "Maximum number of mismatches at the 3' oligo end: $spec3_mis\n";
        #print $file "Number of oligo bases intended as the 3' oligo end: $spec3_length\n";
        
        ##Section 6
        print $file "\nSECONDARY STRUCTURES\n";
        print $file "Self dimer ΔG - minimum: $oli_self kcal/mol\n";
        print $file "Hairpin ΔG - minimum: $oli_hairpin kcal/mol\n";
        print $file "Cross dimer ΔG - minimum: $oli_cross kcal/mol\n";
        print $file "Temperature for secondary structure calculation: $t_dG °C\n";
        
        ##Section 7
        print $file "\nPCR CONDITIONS\n";
        print $file "Monovalent ion concentration: $mon_dG mM\n";
        print $file "Mg2+ concentration: $mg_dG mM\n";
        print $file "Oligo concentration: $oligo_dG uM\n";
        print $file "dNTP concentration: $dNTP_dG uM\n";
        
        ##Section 8
        print $file "\nVISUALIZATION CRITERIA\n";
        foreach my $c (@criteria) {
            print $file "$c\n";
        }
    
        close($file);
               
    } else {
        $newhtml = "Something went wrong. Please reload the page.<br><br>";
    }
} else {
    $newhtml = $newhtml . "Something went wrong. Please reload the page.<br><br>";
}

my $json = $op -> encode({
    result => $newhtml
});
print $json;
