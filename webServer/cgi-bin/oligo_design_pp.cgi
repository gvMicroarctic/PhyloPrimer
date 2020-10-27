#! /usr/bin/perl -w

#use strict;
use warnings;
use CGI;

#connected to tree_input.html, sequence_input.html and tree_input.cgi : display Oligo pages

#set path to folders:
my $path_html = 'path_to_html_folder';

my $cgi = CGI->new;

##get name of the folder
my $defSetDeg = $cgi->param('defSet');

my ($defSet, $defDeg) = split(/-/, $defSetDeg);

#get if S (sequence) or T (tree)
my $input_ST = chop($defSet); #get last letter and understand if T or S

my $input_kind = chop($defSet); #get last letter and understand if a, b or c

#re-construct folder name
my $folder = $defSet . $input_ST . $input_kind;

#path to the folder
my $upload_dir = $path_html . "/uploadsPhyloprimer/" . $folder; #path to the folder

my $count;
my %wildcard;

#wildcards
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

my %pos_count;
my %neg_count;
my %pos_countNew;
my %neg_countNew;
my %bold;

#name of the consensus file
my $in_file = (substr($folder, 0, 8)) . ".consensus.positive";
#read consensus file
my $positive;
open(CONSENSUS, "<${upload_dir}/${in_file}");
while (defined(my $input = <CONSENSUS>)) {
    if ($input !~ /^>/) {
        $positive = $input;
    }
}
close(CONSENSUS);

my @pos = ($positive =~ m/./g); #create hashes with position and base
$count = 0;
foreach my $p (@pos) {
    $count++;
    $pos_count{$count} = $p;
}

my $positiveClean = $positive;
$positiveClean =~ s/-//g;

#name of the consensus file
$in_file = $upload_dir . "/" . (substr($folder, 0, 8)) . ".consensus.negative";
my $insideF = 0;
if (-e $in_file) { ###if file is in the folder
    $insideF = 1;
    #read consensus file
    my $negative;
    open(CONSENSUS, "<${in_file}");
    while (defined(my $input = <CONSENSUS>)) {
        if ($input !~ /^>/) {
            $negative = $input;
        }
    }
    close(CONSENSUS);
    my @neg = ($negative =~ m/./g); #create hashes with position and base
    $count = 0;
    foreach my $n (@neg) {
        $count++;
        $neg_count{$count} = $n;
    }
}

if ($insideF == 1) {
    my $cNew = 0;
    my $inside = 0;
    foreach my $c (sort {$a <=> $b} keys %pos_count) {
        if (($pos_count{$c} eq "-") && ($neg_count{$c} ne "-")) { #positive:gap - negative:no gap
            if ($inside == 0) {
                $bold{$cNew} = '';
            }
            $inside++;
        } elsif (($pos_count{$c} ne "-") && ($neg_count{$c} eq "-")) { #positive:no gap - negative:gap
            $cNew++;
            $pos_countNew{$cNew} = $pos_count{$c};
            $neg_countNew{$cNew} = $neg_count{$c};
            if ($inside > 0) {
                $bold{$cNew} = '';
            }
            $inside=0;
        } elsif (($pos_count{$c} ne "-") && ($neg_count{$c} ne "-")) {
            $cNew++;
            $pos_countNew{$cNew} = $pos_count{$c};
            $neg_countNew{$cNew} = $neg_count{$c};
            if ($inside > 0) {
                $bold{$cNew} = '';
            }
            $inside=0;
        }
    }
} else {
    my $cNew = 0;
    foreach my $c (sort {$a <=> $b} keys %pos_count) {
        if ($pos_count{$c} ne "-") {
            $cNew++;
            $pos_countNew{$cNew} = $pos_count{$c};
        }
    }
}

my $consensus;
my $diff = 0;
my $diffPos;

if ($insideF == 1) {
    my $pos_temp;
    my $neg_temp;
    ##if the bases are different I print them in bold
    ##if there are wildcards, I print the bases as different only if the wildcards do not overlap with any base!
    foreach my $c (sort {$a <=> $b} keys %pos_countNew) {
        if (($c % 50) == 0) { #every 50 bases (50 because I am dividing the sequences evry 50 bases in the visualisation)
            if ($pos_countNew{$c} ne $neg_countNew{$c}) { #if bases are different
                if ($neg_countNew{$c} eq '-') {
                    if (defined($bold{$c})) {
                        $pos_temp .= "<b id='d'>" . $pos_countNew{$c} . "</b>"; #dark blue
                        $diffPosD .= $c . ",";
                        $diff += 14;
                        $nogreen{$c} = '';
                    } else {
                        $pos_temp .= "<b id='c'>" . $pos_countNew{$c} . "</b>"; #blue
                        $diffPosC .= $c . ",";
                        $diff += 14;
                    }
                } else {
                    if ((!(defined($wildcard{$pos_countNew{$c}}))) && (!(defined($wildcard{$neg_countNew{$c}})))) { #if not wildcards
                        if (defined($bold{$c})) {
                            $pos_temp .= "<b id='b'>" . $pos_countNew{$c} . "</b>"; #dark red
                            $diffPosB .= $c . ",";
                            $diff += 14;
                            $nogreen{$c} = '';
                        } else {
                            $pos_temp .= "<b id='a'>" . $pos_countNew{$c} . "</b>"; #red
                            $diffPosA .= $c . ",";
                            $diff += 14;
                        }
                    } elsif ((defined($wildcard{$pos_countNew{$c}})) && (defined($wildcard{$neg_countNew{$c}}))) { #if both wildcards
                        if ((($pos_countNew{$c} eq 'R') && ($neg_countNew{$c} eq 'Y')) or (($pos_countNew{$c} eq 'S') && ($neg_countNew{$c} eq 'W')) or (($pos_countNew{$c} eq 'K') && ($neg_countNew{$c} eq 'M')) or (($neg_countNew{$c} eq 'Y') && ($pos_countNew{$c} eq 'R')) or (($neg_countNew{$c} eq 'W') && ($pos_countNew{$c} eq 'S')) or (($neg_countNew{$c} eq 'M') && ($pos_countNew{$c} eq 'K'))){ #if the two wildcards do not overlap
                            if (defined($bold{$c})) {
                                $pos_temp .= "<b id='b'>" . $pos_countNew{$c} . "</b>"; #dark red
                                $diffPosB .= $c . ",";
                                $diff += 14;
                                $nogreen{$c} = '';
                            } else {
                                $pos_temp .= "<b id='a'>" . $pos_countNew{$c} . "</b>"; #red
                                $diffPosA .= $c . ",";
                                $diff += 14;
                            }
                        } else {
                            $pos_temp .= $pos_countNew{$c};
                        }
                    } elsif (defined($wildcard{$pos_countNew{$c}})) { #if wildcard only if positive consensus
                        if (defined($wildcard{$pos_countNew{$c}}{$neg_countNew{$c}})) {
                            $pos_temp .= $pos_countNew{$c};
                        } else {
                            if (defined($bold{$c})) {
                                $pos_temp .= "<b id='b'>" . $pos_countNew{$c} . "</b>"; #dark red
                                $diffPosB .= $c . ",";
                                $diff += 14;
                                $nogreen{$c} = '';
                            } else {
                                $pos_temp .= "<b id='a'>" . $pos_countNew{$c} . "</b>"; #red
                                $diffPosA .= $c . ",";
                                $diff += 14;
                            }
                        }
                    } elsif (defined($wildcard{$neg_countNew{$c}})) { #if wildcard only if negative consensus
                        if (defined($wildcard{$neg_countNew{$c}}{$pos_countNew{$c}})) {
                            $pos_temp .= $pos_countNew{$c};
                        } else {
                            if (defined($bold{$c})) {
                                $pos_temp .= "<b id='b'>" . $pos_countNew{$c} . "</b>"; #dark red
                                $diffPosB .= $c . ",";
                                $diff += 14;
                                $nogreen{$c} = '';
                            } else {
                                $pos_temp .= "<b id='a'>" . $pos_countNew{$c} . "</b>"; #red
                                $diffPosA .= $c . ",";
                                $diff += 14;
                            }
                        }
                        
                    }
                }
                
            } else { #if bases are the same
                $pos_temp .= $pos_countNew{$c};
                $neg_temp .= $neg_countNew{$c};
            }
            $consensus .= $pos_temp . "<br>";
            $pos_temp = '';
            $diff += 4;
        } else {
            if ($pos_countNew{$c} ne $neg_countNew{$c}) {
                if ($neg_countNew{$c} eq '-') {
                    if (defined($bold{$c})) {
                        $pos_temp .= "<b id='d'>" . $pos_countNew{$c} . "</b>"; #dark blue
                        $diffPosD .= $c . ",";
                        $diff += 14;
                        $nogreen{$c} = '';
                    } else {
                        $pos_temp .= "<b id='c'>" . $pos_countNew{$c} . "</b>"; #blue
                        $diffPosC .= $c . ",";
                        $diff += 14;
                    }
                } else {
                    if ((!(defined($wildcard{$pos_countNew{$c}}))) && (!(defined($wildcard{$neg_countNew{$c}})))) {
                        if (defined($bold{$c})) {
                            $pos_temp .= "<b id='b'>" . $pos_countNew{$c} . "</b>"; #dark red
                            $diffPosB .= $c . ",";
                            $diff += 14;
                            $nogreen{$c} = '';
                        } else {
                            $pos_temp .= "<b id='a'>" . $pos_countNew{$c} . "</b>"; #red
                            $diffPosA .= $c . ",";
                            $diff += 14;
                        }
                    } elsif ((defined($wildcard{$pos_countNew{$c}})) && (defined($wildcard{$neg_countNew{$c}}))) {
                        if ((($pos_countNew{$c} eq 'R') && ($neg_countNew{$c} eq 'Y')) or (($pos_countNew{$c} eq 'S') && ($neg_countNew{$c} eq 'W')) or (($pos_countNew{$c} eq 'K') && ($neg_countNew{$c} eq 'M')) or (($neg_countNew{$c} eq 'Y') && ($pos_countNew{$c} eq 'R')) or (($neg_countNew{$c} eq 'W') && ($pos_countNew{$c} eq 'S')) or (($neg_countNew{$c} eq 'M') && ($pos_countNew{$c} eq 'K'))){
                            if (defined($bold{$c})) {
                                $pos_temp .= "<b id='b'>" . $pos_countNew{$c} . "</b>"; #dark red
                                $diffPosB .= $c . ",";
                                $diff += 14;
                                $nogreen{$c} = '';
                            } else {
                                $pos_temp .= "<b id='a'>" . $pos_countNew{$c} . "</b>"; #red
                                $diffPosA .= $c . ",";
                                $diff += 14;
                            }
                        } else {
                            $pos_temp .= $pos_countNew{$c};
                        }
                    } elsif (defined($wildcard{$pos_countNew{$c}})) {
                        if (defined($wildcard{$pos_countNew{$c}}{$neg_countNew{$c}})) {
                            $pos_temp .= $pos_countNew{$c};
                        } else {
                            if (defined($bold{$c})) {
                                $pos_temp .= "<b id='b'>" . $pos_countNew{$c} . "</b>"; #dark red
                                $diffPosB .= $c . ",";
                                $diff += 14;
                                $nogreen{$c} = '';
                            } else {
                                $pos_temp .= "<b id='a'>" . $pos_countNew{$c} . "</b>"; #red
                                $diffPosA .= $c . ",";
                                $diff += 14;
                            }
                        }
                        
                    } elsif (defined($wildcard{$neg_countNew{$c}})) {
                        if (defined($wildcard{$neg_countNew{$c}}{$pos_countNew{$c}})) {
                            $pos_temp .= $pos_countNew{$c};
                        } else {
                            if (defined($bold{$c})) {
                                $pos_temp .= "<b id='b'>" . $pos_countNew{$c} . "</b>"; #dark red
                                $diffPosB .= $c . ",";
                                $diff += 14;
                                $nogreen{$c} = '';
                            } else {
                                $pos_temp .= "<b id='a'>" . $pos_countNew{$c} . "</b>"; #red
                                $diffPosA .= $c . ",";
                                $diff += 14;
                            }
                        }
                    }
                }
            } else {
                $pos_temp .= $pos_countNew{$c};
            }
        }
    }
    $consensus .= $pos_temp;
    
    #Add green flag to consensus
    my $len = length($consensus);
    my @consensus = split(//, $consensus);
    my $consensusNew;
    $count = 0;
    foreach my $c (0..$len) {
        if ($consensus[$c] =~ /[A-Z]/) {
            $count++;
            if ((defined($bold{$count})) && (!(defined($nogreen{$count})))) {
                $consensusNew .= "<b id='e'>" . $consensus[$c] . "</b>"; #green
                $diffPosE .= $count . ",";
                $diff += 14;
            } else {
                $consensusNew .= $consensus[$c];
            }
        } else {
            $consensusNew .= $consensus[$c];
        }
    }
    $diffPos = $diffPosA . "|" . $diffPosB . "|" . $diffPosC . "|" . $diffPosD . "|" . $diffPosE;
    $consensus = $consensusNew;
} else { #if no negative consensus
    my $pos_temp;
    foreach my $c (sort {$a <=> $b} keys %pos_countNew) {
        if (($c % 50) == 0) { #every 50 bases (50 because I am dividing the sequences evry 50 bases in the visualisation)
            $pos_temp .= $pos_countNew{$c};
            $consensus .= $pos_temp . "<br>";
            $diffPos .= $c . ",";
            $diff += 4;
            $pos_temp = '';
        } else {
            $pos_temp .= $pos_countNew{$c};
        }
    }
    $consensus .= $pos_temp;
}

my $degOption = "'1'";
#deg bases
if ($defDeg > 20) {
    $degOption = "'3'";
} elsif ($defDeg > 10) {
    $degOption = "'2'";
}

#.taxonomy file exists
my $out_file_species = $upload_dir . "/" . (substr($folder, 0, 8)) . ".taxonomy";
my $taxonomyVis = "'NO'";
if (-e $out_file_species) {
    $taxonomyVis = "'YES'";
}

$consensus = "\"" . $consensus . "\"";
$consensusClean = "\"" . $positiveClean . "\"";
$consensus_diff = "'" . $diff . "'";
$folder = "'" . $folder . "'";
$diffPos = "'" . $diffPos . "'";
#$input_ST = "'" . $input_ST . "'";



#load consensus, save and send to html printing
#check folder
#check consensus

my $html = <<"END HTML";
Content-Type: text/html

<HTML>
<link rel='stylesheet' href='../phyloprimer/phyloprimer.css'>
<HEAD>
<meta http-equiv='Content-Type' content='text/html; charset=utf-8'/>
<TITLE>Oligo Design</TITLE>
<style>

#title {
    display: block;
    margin-left: auto;
    margin-right: auto;
    width: 540px;
    height: auto;
}

h1 {
    font-size: 50px;
}

#navbar {
    list-style-type: none;
    margin: 0;
    padding: 0;
    overflow: hidden;
    height: 60px;
    background-color: #008CBA;
}

li a {
    display: block;
    color: white;
    text-align: center;
    padding-top: 20px;
    text-decoration: none;
}

li a:hover {
    background-color: #00008B !important;
}

a {
	height: 60px;
}

body {
margin:30px;
}

#primerLi, #probeLi, #bothLi {
	float: left;
	width: 33.333%;
	height: 60px;
{

</style>
<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.4.0/jquery.min.js"></script>
<SCRIPT TYPE="TEXT/JAVASCRIPT">
\$(document).ready(function() {
    \$("#includedContent").load("../phyloprimer/primer_input.html");
    \$("#primerB").css('background-color', '#00008B');
    \$("#probeB").css('background-color', '#008CBA');
    \$("#bothB").css('background-color', '#008CBA');
    \$("#primerB").click(function(){
        \$("#includedContent").load("../phyloprimer/primer_input.html");
        \$("#primerB").css('background-color', '#00008B');
        \$("#probeB").css('background-color', '#008CBA');
        \$("#bothB").css('background-color', '#008CBA');
    })
    \$("#probeB").click(function(){
	\$("#includedContent").load("../phyloprimer/oligo_input.html");
        \$("#primerB").css('background-color', '#008CBA');
        \$("#probeB").css('background-color', '#00008B');
        \$("#bothB").css('background-color', '#008CBA');
    })
    \$("#bothB").click(function(){
	\$("#includedContent").load("../phyloprimer/probe_input.html");
        \$("#primerB").css('background-color', '#008CBA');
        \$("#probeB").css('background-color', '#008CBA');
        \$("#bothB").css('background-color', '#00008B');
    })
})
    
consensus = $consensus;
consensusClean = $consensusClean;
folder = $folder;
consensus_diff = $consensus_diff;
diffPos = $diffPos;
degOption = $degOption;
taxonomyVis = $taxonomyVis;

//cookies
function getCookie(cname) {
    var name = cname + "=";
    var decodedCookie = decodeURIComponent(document.cookie);
    var ca = decodedCookie.split(';');
    for (var i = 0; i < ca.length; i++) {
        var c = ca[i];
        while (c.charAt(0) == ' ') {
            c = c.substring(1);
        }
        if (c.indexOf(name) == 0) {
            return c.substring(name.length, c.length);
        }
    }
    return "";
}
    
var cookieBack = '';

function checkCookie() {
    if (getCookie("email") !== "") {
        
        cookieBack = {
            
            email : getCookie("email"),
            
            plength_min:getCookie("plength_min"),
            plength_max:getCookie("plength_max"),
            tm_min:getCookie("tm_min"),
            tm_max:getCookie("tm_max"),
            maxrun:getCookie("maxrun"),
            maxrepeat:getCookie("maxrepeat"),
            gcclamp_min:getCookie("gcclamp_min"),
            gcclamp_max:getCookie("gcclamp_max"),
            mingc:getCookie("mingc"),
            maxgc:getCookie("maxgc"),
            hide_wild2:getCookie("hide_wild2"),
            hide_wild3:getCookie("hide_wild3"),
            hide_wild4:getCookie("hide_wild4"),
            end5_deg:getCookie("end5_deg"),
            end3_deg:getCookie("end3_deg"),
            
            //only oligo
            anti:getCookie("anti"), //probe
            sense:getCookie("sense"), //probe
            both:getCookie("both"), //probe
            probe:getCookie("probe"),
            primer:getCookie("primer"),
            
            //only probe
            plength_min_1:getCookie("plength_min_1"),
            plength_max_1:getCookie("plength_max_1"),
            tm_min_1:getCookie("tm_min_1"),
            tm_max_1:getCookie("tm_max_1"),
            maxrun_1:getCookie("maxrun_1"),
            maxrepeat_1:getCookie("maxrepeat_1"),
            gcclamp_min_1:getCookie("gcclamp_min_1"),
            gcclamp_max_1:getCookie("gcclamp_max_1"),
            mingc_1:getCookie("mingc_1"),
            maxgc_1:getCookie("maxgc_1"),
            hide_wild2_1:getCookie("hide_wild2_1"),
            hide_wild3_1:getCookie("hide_wild3_1"),
            hide_wild4_1:getCookie("hide_wild4_1"),
            end5_deg_1:getCookie("end5_deg_1"),
            end3_deg_1:getCookie("end3_deg_1"),
            tm_difference_1:getCookie("tm_difference_1"),
            
            amplicon_len_min:getCookie("amplicon_len_min"),
            amplicon_len_max:getCookie("amplicon_len_max"),
            tm_difference:getCookie("tm_difference"),
            ta_min:getCookie("ta_min"),
            ta_max:getCookie("ta_max"),
            
            spec_mis:getCookie("spec_mis"),
            spec3_mis:getCookie("spec3_mis"),
            spec3_length:getCookie("spec3_length"),
            
            oli_self:getCookie("oli_self"),
            oli_hairpin:getCookie("oli_hairpin"),
            oli_cross:getCookie("oli_cross"),
            t_dG:getCookie("t_dG"),
            
            mon_dG:getCookie("mon_dG"),
            mg_dG:getCookie("mg_dG"),
            oligo_dG:getCookie("oligo_dG"),
            dNTP_dG:getCookie("dNTP_dG"),
            
            tm_sel:getCookie("tm_sel"),
            dG_sel:getCookie("dG_sel"),
            deg_sel:getCookie("deg_sel"),
            maximize_sel1:getCookie("maximize_sel1"),
            maximize_sel2:getCookie("maximize_sel2"),
            
            species_sel:getCookie("species_sel"),
            genus_sel:getCookie("genus_sel"),
            family_sel:getCookie("family_sel"),
            order_sel:getCookie("order_sel"),
            class_sel:getCookie("class_sel"),
            phylum_sel:getCookie("phylum_sel"),
            domain_sel:getCookie("domain_sel"),
            
            type:getCookie("type")
            
        }
    }
}

</SCRIPT>
</HEAD>
<body onload="checkCookie()">
<h1 align="center"><img src="../phyloprimer/figures/OligoDesign.png" alt="title" id="title"></h1>
<nav>
<ul id ="navbar">
<li id="primerLi"><a id="primerB" selected="selected">PRIMER DESIGN</a></li>
<li id="bothLi"><a id="bothB">PRIMER AND PROBE DESIGN</a></li>
<li id="probeLi"><a id="probeB">SINGLE OLIGO DESIGN</a></li>
</ul>
</nav>

<div id="includedContent"></div>

</body>
END HTML

print $html;
