#! /usr/bin/perl -w

#use strict;
use warnings;
use CGI;

#connected to Oligo_check_pp.pl : visualise oligo check results

#set path to folders:
my $path_html = 'path_to_html_folder';

my $cgi = CGI->new;

###get name of the folder
my $defSet = $cgi->param('defSet');

#archive name
$archive = $defSet . ".tar.gz";

my $input_STC = chop($defSet); #get last letter and understand if T, S or C

my $input_kind = chop($defSet); #get last letter and understand if a, b or c

#re-construct folder name
my $folder = $defSet . $input_STC . $input_kind;

#path to the folder
my $upload_dir = $path_html . "/analysesPhyloprimer/" . $folder; #path to the folder

my %hairpin;
my %self;

my $primer;
my $seq;

#Self-dimer data
open(INFO, "<${upload_dir}/selfDimer.txt");
while (defined(my $input = <INFO>)) {
    chomp($input);
    if ($input =~ /^>/) {
        my @info = split(/\t/, $input);
        $primer = $info[1];
    } elsif ($input =~ /^@/) {
        my @info = split(/\t/, $input);
        if ($info[1] == 0) {
            $seq .= "<pre>PhyloPrimer did not find any self dimer with a &Delta;G lower than 0.</pre>";
            $self{$primer} = $seq;
            $seq = '';
        } else {
            $self{$primer} = $seq;
            $seq = '';
        }
    } elsif ($input eq "") {
        $seq .= "<br>";
    } elsif ($input =~ /^dG:/) {
        $input =~ s/dG/\&Delta\;G/g;
        $seq .= "<pre>" . $input . "</pre>";
    } else {
        $seq .= "<pre>" . $input . "</pre>";
    }
}
close(INFO);

#Hairpin data
open(INFO, "<${upload_dir}/hairpin.txt");
while (defined(my $input = <INFO>)) {
    chomp($input);
    if ($input =~ /^>/) {
        my @info = split(/\t/, $input);
        $primer = $info[1];
    } elsif ($input =~ /^@/) {
        my @info = split(/\t/, $input);
        if ($info[1] == 0) {
            $seq .= "<pre>PhyloPrimer did not find any hairpin with a &Delta;G lower than 0.</pre>";
            $hairpin{$primer} = $seq;
            $seq = '';
        } else {
            $hairpin{$primer} = $seq;
            $seq = '';
        }
        
    } elsif ($input eq "") {
        $seq .= "<br>";
    } else {
        if ($input =~ /\\/) {
            $input =~ s/\\//g;
            $seq .= "<pre>" . $input . "&#92</pre>";
        } elsif ($input =~ /^dG:/) {
            $input =~ s/dG/\&Delta\;G/g;
            $seq .= "<pre>" . $input . "</pre>";
        } else {
            $seq .= "<pre>" . $input . "</pre>";
        }
    }
}
close(INFO);


my $pair = 0;
my $consensus;
my $tableBlast;
my $difference;
my %pairSeq;
my %pairInfo;
my $project;

#file with all primer information
my $last;
my %primerInfo;
open(INFO, "<${upload_dir}/info.primer");
while (defined(my $input = <INFO>)) {
    chomp($input);
    if ($input =~ /^BLAST_TABLE/) {
        my @info = split(/\t/, $input);
        $tableBlast = $info[1];
    } elsif ($input =~ /^PROJECT/) {
        my @info = split(/\t/, $input);
        $project = $info[1];
    } elsif ($input =~ /^PIECHART_F/) {
        my @info = split(/\t/, $input);
        $tablePieChartF = $info[1];
    } elsif ($input =~ /^USER_TABLE/) {
        my @info = split(/\t/, $input);
        $tableUser = $info[1];
    } elsif ($input =~ /^SPECIES/) {
        my @info = split(/\t/, $input);
        $sp = $info[1];
    } elsif ($input =~ /^OLIGO/) {
        $pair++;
        my @info = split(/\t/, $input);
        $pairSeq{$pair}{'F'} = $info[1];
        $pairInfo{$pair} .= $info[1] . "!";
        $last = $info[1];
    } elsif ($input =~ /^START/) {
        my @info = split(/\t/, $input);
        $pairInfo{$pair} .= $info[1] . "!";
        $primerInfo{$last}{'START'} = $info[1];
    } elsif ($input =~ /^LENGTH/) {
        my @info = split(/\t/, $input);
        $pairInfo{$pair} .= $info[1] . "!";
        $primerInfo{$last}{'LEN'} = $info[1];
    } elsif ($input =~ /^GC/) {
        my @info = split(/\t/, $input);
        $pairInfo{$pair} .= $info[1] . "!";
    } elsif ($input =~ /^TM/) {
        my @info = split(/\t/, $input);
        $pairInfo{$pair} .= $info[1] . "!";
    } elsif ($input =~ /^SELF/) {
        my @info = split(/\t/, $input);
        if ($info[1] =~ />=/) {
            $pairInfo{$pair} .= "&ge; 0" . "!";
        } else {
            $pairInfo{$pair} .= $info[1] . "!";
        }
    } elsif ($input =~ /^HAIR/) {
        my @info = split(/\t/, $input);
        if ($info[1] =~ />=/) {
            $pairInfo{$pair} .= "&ge; 0" . "!";
        } else {
            $pairInfo{$pair} .= $info[1] . "!";
        }
    } elsif ($input =~ /^\\\\/) {
        my $for = $pairSeq{$pair}{'F'};
        if ($for =~ /-/) {
            ($for) = $pairSeq{$pair}{'F'} =~ /-(.*)/;
        }
        $pairInfo{$pair} .= $self{$for} . "!" . $hairpin{$for} . "!" . $sp; #add info for dG
    }
}
close(INFO);

$firstValue = "'" . $pairInfo{1} . "'"; #send first value when the page opens
my $folderLink = $folder; #folder for download the archive
my $nameFile = substr($folder, 0, 8);
$folder = "'" . $folder . "'"; #folder
$tableBlast = "\"" . $tableBlast . "\""; #Blast table
$tablePieChartF = "'" . $tablePieChartF . "'"; #tablePieChartF
$tableUser = "\"" . $tableUser . "\""; #Blast user table

my $radio;

foreach my $pair (sort {$a <=> $b} keys %pairInfo) { #print only 50 combinations
    if ($pair == 1) {
        $radio .= "<input type='radio' name='pair' class='pair' value='$pairInfo{$pair}' checked='checked'>";
        $radio .= "<span class='prefixMain'>$pairSeq{$pair}{'F'}<br></span>";
    } elsif ($pair <= 100) {
        $radio .= "<input type='radio' name='pair' class='pair' value='$pairInfo{$pair}'>";
        $radio .= "<span class='prefixMain'>$pairSeq{$pair}{'F'}<br></span>";
    } else {
        last;
    }
}

my $html = <<"END HTML";
Content-Type: text/html

<html>
<head>
<meta http-equiv='Content-Type' content='text/html; charset=utf-8'/>
<link rel='stylesheet' href='../phyloprimer/phyloprimer.css'>
<title>PhyloPrimer Results</title>
    <style type="text/css">
        h2 {
            text-align: center;
        }

        table {
            font-family: courier;
        }

        .prefix {
            font-size: 15px;
        }

        .prefix1 {
            font-size: 20px;
        }

        /* space radioButtons */
        input[type="radio"] {
            margin: 10px 5px 0 0;
        }

        #wrapper {
            width: 100%;
            margin: 30px;
            background-color: #fff;
        }

        #top-nav {
            position: fixed;
            left: 0;
            right: 0;
            top: 0;
            height: 120px;
            width: 100%;
        }

        #title {
            display: block;
            margin-left: auto;
            margin-right: auto;
            margin-bottom: 50px;
            width: 800px;
            height: auto;
        }


        #oligoList {
            position: fixed;
            width: 28.5%;
            height: calc(100% - 120px);
            left: 0;
            right: 0;
            overflow-y: scroll;
            top: 200px;
            margin-left: 1%;
            padding-right: 1%;
            padding-left: 1%;
            /* border-style: dotted; */
        }

        #oligoEmpty {
            position: fixed;
            margin: 200px 0 0 31.85%;
            width: 0.5%;
            left: 0;
            top: 0;
            right: 0;
            overflow-y: scroll;
            height: calc(100% - 120px);
            background-color: #00008B;
        }

        #oligoInfo {
            margin: 200px 0 0 33%;
            width: calc(65.5% - 30px);
            padding-right: 1%;
            padding-left: 1%;
            position: fixed;
            left: 0;
            top: 0;
            right: 0;
            overflow-y: scroll;
            height: calc(100% - 120px);
            /*border-style: dotted;*/
            
            /* to give priority at dropdown menu*/
        }


        /*#oligoInfo_int {
            width: 100%;
            z-index: 10000;
            height: calc(100% - 120px);
            position: flex;
            overflow: scroll;
        }*/


        /* NAVBAR - start */
        #oligoTitle {
            position: fixed;
            left: 0;
            right: 0;
            top: 125px;
            height: 60px;
            width: 98%;
            margin: 0;
            margin: 0% 1% 0% 1%;
        }

        ul {
            padding: 0;
            background: #008CBA;
        }

        #navbar {
            margin: 0;
            height: 60px;
            background: #008CBA;
            display: block;
        }

        .maindG {
            height: 40px;
        }

        .mainTax {
            height: 40px;
        }

        ul li ul.dropdown {
            min-width: 100%;
            /* Set width of the dropdown */
            background: #008CBA;
            display: none;
            position: absolute;
            z-index: 9999;
            left: 0;
        }

        ul li:hover ul.dropdown {
            /* Display the dropdown */
            display: block;
        }

        ul li ul.dropdown li {
            display: block;
        }


        ul li ul.dropdown1 {
            min-width: 100%;
            /* Set width of the dropdown */
            background: #008CBA;
            display: none;
            position: absolute;
            z-index: 9999;
            left: 0;
        }

        ul li:hover ul.dropdown1 {
            /* Display the dropdown */
            display: block;

        }

        ul li ul.dropdown1 li {
            display: block;
        }

        a {
            height: 40px;
            text-align: center;
            padding-top: 20px;
        }

        #self_navbar,
        #hair_navbar,
        #domain_navbar,
        #phylum_navbar,
        #class_navbar,
        #order_navbar,
        #family_navbar,
        #genus_navbar,
        #species_navbar {
            height: 40px;
            text-align: center;
            padding-top: 9px;

        }

        ul li {
            display: inline-block;
            position: relative !important;
            line-height: 21px;
            text-align: left;
        }

        ul li a {
            display: block;
            color: white;
            text-decoration: none;
            justify-content: center;
        }

        li a:hover {
            background-color: #00008B !important;
        }

        #oligo_navbarLi,
        #empty_navbarLi,
        #consensus_navbarLi,
        #dG_navbarLi,
        #taxonomy_navbarLi,
        #userCheck_navbarLi,
        #download_navbarLi {
            float: left;
        }

        #consensus_navbarLi:hover,
        #dG_navbarLi:hover,
        #taxonomy_navbarLi:hover,
        #userCheck_navbarLi:hover,
        #download_navbarLi:hover {
            color: white;
            height: 60px;
        }

        .maindG:hover,
        .mainTax:hover {
            color: white;
            height: 40px;
        }

        #oligo_navbarLi {
            width: 31.5%;
            height: 99%;
            background-color: #00008B;
        }

        #empty_navbarLi {
            width: 0.5%;
            height: 100%;
            background-color: white;
        }

        #consensus_navbarLi,
        #dG_navbarLi,
        #taxonomy_navbarLi,
        #userCheck_navbarLi,
        #download_navbarLi {
            height: 100%;
        }

        /* NAVBAR - end */

        /* PAGES - start*/
        #consensus_info,
        #dG_info,
        #taxonomy_info,
        #userCheck_info,
        #download_info {
            width: 100%;
            height: 100%;
            justify-content: center;
            top: 0;
            left: 0;
            margin: 0 auto;
            display: none;
            /*border-style: solid;*/
        }

        /* PAGES - end*/

        /* TABLE - start*/

        table {
            width: 100%;
            border-collapse: collapse;
        }

        td,
        th {
            border-bottom: 1px solid #ddd;
            border-left: 1px solid #ddd;
            padding: 15px;
            text-align: center;
        }

        .canvasjs-chart-credit {
            display: none;
        }

        /*divide different sections tableBlast*/
        #dividingCell {
            all: unset;
            border-top: 3px solid #00008B;
            border-bottom: 3px solid #00008B;
            height: 10;
        }

        /*scroll horizontally div table*/
        #taxonomy_table, #userCheck_table, #pair_info {
            width: 100%;
            overflow-x: auto;
            white-space: nowrap;
            border-collapse: collapse;
            border-right: 3px solid #00008B;
            border-left: 3px solid #00008B;
        }

        #taxonomy_contentF {
            height: 500px;
            width: 100%;
            margin: 0px auto;
            text-align: center;
        }

        /* TABLE - end */


	/* User list entry */
	#userList {
	   display: none;
	}

/* Search section */
#search {
background-color: lightblue;
padding-top: 2px;
padding-bottom: 10px;
padding-right: 10px;
padding-left: 10px;
}

/* Taxonomy dropdown */

/*the container must be positioned relative:*/
.autocomplete {
position: relative;
display: inline-block;
}
#dropTaxa {
background-color: white;
border: 1px solid transparent;
padding: 10px;
width: 100%;
font-size: 10px;
}

#dropTaxa:hover {
background-color: #f5f5f0;
}

.autocomplete-items {
position: absolute;
border: 1px solid #d4d4d4;
    border-bottom: none;
    border-top: none;
    z-index: 99;
    /*position the autocomplete items to be the same width as the container:*/
top: 100%;
left: 0;
right: 0;
}

.autocomplete-items div {
cursor: pointer;
    background-color: white;
    border-bottom: 1px solid #d4d4d4;
}

/*when hovering an item:*/
.autocomplete-items div:hover {
    background-color: #f5f5f0;
}

/*when navigating through the items using the arrow keys:*/
.autocomplete-active {
    background-color: DodgerBlue !important;
color: #ffffff;
}

[id^="clear_taxon"], #add_taxon {
width: 20px;
border: none;
text-align: center;
background-color: #008CBA;
}

#search_button {
width: 80px;
border: none;
text-align: center;
background-color: #008CBA;
}

    </style>
</head>

<body>

    <div id="wrapper">
        <div id="top-nav">
            <br>
            <img src="../phyloprimer/figures/PhyloPrimerResult.png" alt="title" id="title">
        </div>

        <div id="oligoTitle" align="center">
            <ul id="navbar">
                <li class="main" id="oligo_navbarLi"><a id="oligo_navbar">OLIGO LIST</a></li>
                <li class="main" id="empty_navbarLi"><a id="empty_navbar">&nbsp;</a></li>
                <li class="main" id="consensus_navbarLi"><a id="consensus_navbar">SPECIFICS</a>
                </li>
                <li class="main" id="dG_navbarLi"><a id="dG_navbar">&Delta;G &#9662;</a>
                    <ul class="dropdown">
                        <li class="maindG" id="self_navbarLi"><a id="self_navbar">Self dimers</a></li>
                        <li class="maindG" id="hair_navbarLi"><a id="hair_navbar">Hairpins</a></li>
                    </ul>
                </li>
                <li class="main" id="taxonomy_navbarLi"><a id="taxonomy_navbar">GENBLANK BLAST &#9662;</a>
                    <ul class="dropdown1">
                        <li class="mainTax" id="domain_navbarLi"><a id="domain_navbar">Domain</a></li>
                        <li class="mainTax" id="phylum_navbarLi"><a id="phylum_navbar">Phylum</a>
                        </li>
                        <li class="mainTax" id="class_navbarLi"><a id="class_navbar">Class</a></li>
                        <li class="mainTax" id="order_navbarLi"><a id="order_navbar">Order</a></li>
                        <li class="mainTax" id="family_navbarLi"><a id="family_navbar">Family</a>
                        </li>
                        <li class="mainTax" id="genus_navbarLi"><a id="genus_navbar">Genus</a></li>
                        <li class="mainTax" id="species_navbarLi"><a id="species_navbar">Species</a>
                        </li>
                    </ul>
                </li>
                <li class="main" id="userCheck_navbarLi"><a id="userCheck_navbar">YOUR BLAST</a></li>
                <li class="main" id="download_navbarLi"><a id="download_navbar">ARCHIVE</a></li>

            </ul>
        </div>

        <div id="oligoList">

            <h2>$project</h2>
            <div id="oligoList_intro"><h3>Please select the oligo you want to explore from the list below. All the information can be also found in the archive which is downloadable from the Archive window. Use the following research section for highlighting the oligos that reflect certain criteria:</h3></div>

<div id = "search">


<h3><b>Filtering criteria</b></h3>
Taxonomy <br>
<div class="input_fields_wrap">

<div class="autocomplete" style="width:300px;"><input style="width:300px;" type="text" id="searchTax" class="searchTaxon"></div>
<button id="add_taxon" onclick="addTax()">+</button>
</div>

<br>Melting temperature (<sup>o</sup>C) <br>
<input type="text" id="searchTm_min" placeholder="minimum"> -
<input type="text" id="searchTm_max" placeholder="maximum"></br>

<br>Oligo length (bases) <br>
<input type="text" id="searchLen_min" placeholder="minimum"> -
<input type="text" id="searchLen_max" placeholder="maximum"></br>
<br>
<button type="button" id = "search_button" onclick="searchEng()">SEARCH</button></br>

</div>


</br>

<div id='pair'>
$radio
</div>

            <br>
            <br>
            <br>
            <br>
            <br>
            <br>
        </div>

        <div id="oligoEmpty"></div>

        <div id="oligoInfo">
            <div id="consensus_info">
                <h2>Specifics</h2>
            <div id="consensus_intro"></div>
                <h3>The table reports the specifics for selected oligo. More information about the &Delta;G and T<sub>m</sub> can be found in the manual. All the reported position at the oligo position at the 5’ end.</h3>
                <div id="pair_info"></div>
                <br>
                <br>
                <br>
                <br>
                <br>
                <br>
            </div>


            <div id="dG_info">
                <div id="dG_content"></div>
                <br>
                <br>
                <br>
                <br>
                <br>
                <br>
            </div>

            <div id="taxonomy_info">
                <div id="taxonomy_intro"></div>
                <h3>The oligos were blasted against the nt database which is a taxonomical comprehensive subset of the GenBank database. The pie chart shows the number of hits for each taxon that were matched by the selected oligo. Please keep in mind that the abundance of each taxon may be skewed towards the taxa that are more represented in the BLAST database.</h3>
                <div id="taxonomy_contentF"></div>
                <br>
                
                <div id="taxonomy_table_mess"></div>
                <div id="taxonomy_table"></div>
                <br>
                <br>
                <br>
                <br>
                <br>
                <br>
            </div>
            <div id="userCheck_info">
                <h2>Your BLAST</h2>
                <div id="userCheck_content"></div>
                <br>
                <div id="userCheck_table"></div>
                <br>
                <br>
                <br>
                <br>
                <br>
                <br>
            </div>
            <div id="download_info">
                <h2>Archive</h2>
                <h3>The archive with all the files you need for a deeper exploration of the data can be downloaded
                <a href="https://www.cerealsdb.uk.net/cerealgenomics/phyloprimer/analysesPhyloprimer/${folderLink}/PhyloPrimer_${nameFile}.zip"
                    download>here</a>.</h3>

<h3>The archive contains two folders:</h3>
<ol>
<li><b>results</b> which contains files reporting:</li>
<br>
<ol type="A">
<li>the list of all the oligos</li>
<br>
<ol type="a">
<li>oligoList.txt</li>
</ol>
<br>
<li>the secondary structures</li>
<br>
<ol type="a">
<li>selfDimer.txt</li>
<li>hairpin.txt</li>
</ol>
<br>
<li>the BLAST search outputs</li>
<br>
<ol type="a">
<li>genbankBLAST.m8 - BLAST search vs the GenBank database</li>
<li id="userList">userBLAST.m8 - BLAST search vs the your fasta file</li>
</ol>
</ol>
<br>
<li><b>inputs</b> which contains your input files and files PhyloPrimer used for the primer design (more details on
this in README file).</li>
</ol>
            </div>
        </div>
        <br>
        <br>
        <br>
        <br>
        <br>
        <br>
    </div>




</body>

<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.4.0/jquery.min.js"></script>
<script src="https://www.cerealsdb.uk.net/cerealgenomics/phyloprimer/canvasjs.min.js"></script>
<script TYPE="TEXT/JAVASCRIPT">

//search engine
function searchEng() {
    
    var filterTax = [];
    \$(".searchTaxon").each(function(){ //all taxonomy for all the input fields
        filterTax.push(\$(this).val().toLowerCase());
    });
    
    var filterTm_min = document.getElementById("searchTm_min").value;
    var filterTm_max = document.getElementById("searchTm_max").value;
    if (filterTm_min == "") {
        filterTm_min = 0;
    }
    if (filterTm_max == "") {
        filterTm_max = 1000;
    }
    
    var filterLen_min = document.getElementById("searchLen_min").value;
    var filterLen_max = document.getElementById("searchLen_max").value;
    if (filterLen_min == "") {
        filterLen_min = 0;
    }
    if (filterLen_max == "") {
        filterLen_max = 1000;
    }
    
    var names = document.getElementsByClassName("prefixMain");
    var nodes = document.getElementsByClassName("pair");
    
    for (i = 0; i < nodes.length; i++) {
        
        res = nodes[i].value.split("!");
                
        if ((res[3] >= filterTm_min) && (res[3] <= filterTm_max) && (res[1] >= filterLen_min) && (res[1] <= filterLen_max)) {
            var inside = 1;
            if (filterTax[0] !== '') {
                for (a=0; a < filterTax.length; a++) {
                    if (nodes[i].value.toLowerCase().includes(filterTax[a])) {
                        
                    } else {
                        inside = 0;
                    }
                }
                if (inside == 1) {
                    names[i].style.backgroundColor = "lightblue";
                } else {
                    names[i].style.backgroundColor = "white";
                }
            } else {
                names[i].style.backgroundColor = "lightblue";
            }
        } else {
            names[i].style.backgroundColor = "white";
        }
        
    }
    
}

    //oligoInfo has z-index of 1 when the user hovers on the navbar, so than the dropdown menu can be visualized. And it has an high z-index when not on the the navbar so then the scrolling works!!
    \$("#oligoTitle").hover(
        function() {
            \$("#oligoInfo").css('z-index', -1);
    },
        function() {
            \$("#oligoInfo").css('z-index', 900);
    }
    );



    //GET DATA FROM CGI AND FIRST VALUE
    \$(document).ready(function () {
        //values from cgi
        folder = $folder;
        firstValue = $firstValue;
        tableBlast = $tableBlast;
        tablePieChartF = $tablePieChartF;
        tableUser = $tableUser;

        tableBlast = tableBlast.replace(/\\[sub1\\]/g, "'");
        tableBlast = tableBlast.replace(/\\[sub2\\]/g, ",");
        tableBlast = tableBlast.replace(/\\[sub3\\]/g, "(");
        tableBlast = tableBlast.replace(/\\[sub4\\]/g, ")");
        tableBlast = tableBlast.replace(/\\[sub5\\]/g, ":");
        tableBlast = tableBlast.replace(/\\[sub7\\]/g, "*");
        tableBlast = tableBlast.replace(/\\[sub8\\]/g, "<");
        tableBlast = tableBlast.replace(/\\[sub9\\]/g, ">");
        tableBlast = tableBlast.replace(/\\[sub10\\]/g, "-");
        tableBlast = tableBlast.replace(/\\[sub11\\]/g, "+");
        tableBlast = tableBlast.replace(/\\[sub12\\]/g, "'");
        tableBlast = tableBlast.replace(/\\[sub13\\]/g, "#");
        tableBlast = tableBlast.replace(/\\[sub14\\]/g, "&");
        tableBlast = tableBlast.replace(/\\[sub15\\]/g, "^");
        tableBlast = tableBlast.replace(/\\[sub16\\]/g, "/");

        tableUser = tableUser.replace(/\\[sub1\\]/g, "'");
        tableUser = tableUser.replace(/\\[sub2\\]/g, ",");
        tableUser = tableUser.replace(/\\[sub3\\]/g, "(");
        tableUser = tableUser.replace(/\\[sub4\\]/g, ")");
        tableUser = tableUser.replace(/\\[sub5\\]/g, ":");
        tableUser = tableUser.replace(/\\[sub7\\]/g, "*");
        tableUser = tableUser.replace(/\\[sub8\\]/g, "<");
        tableUser = tableUser.replace(/\\[sub9\\]/g, ">");
        tableUser = tableUser.replace(/\\[sub10\\]/g, "-");
        tableUser = tableUser.replace(/\\[sub11\\]/g, "+");
        tableUser = tableUser.replace(/\\[sub12\\]/g, "'");
        tableUser = tableUser.replace(/\\[sub13\\]/g, "#");
        tableUser = tableUser.replace(/\\[sub14\\]/g, "&");
        tableUser = tableUser.replace(/\\[sub15\\]/g, "^");
        tableUser = tableUser.replace(/\\[sub16\\]/g, "/");
        tableUser = tableUser.replace(/_/g, " ");
        tableUser = tableUser.replace(/\\[sub17\\]/g, "_");

        //correct nav bar if user table
        if (tableUser != "no") {
            document.getElementById("consensus_navbarLi").style.width = "13.6%";
            document.getElementById("dG_navbarLi").style.width = "13.6%";
            document.getElementById("taxonomy_navbarLi").style.width = "13.6%";
            document.getElementById("userCheck_navbarLi").style.width = "13.6%";
            document.getElementById("download_navbarLi").style.width = "13.6%";
            \$("#consensus_info").css('display', 'block');
            \$("#userList").css('display', 'list-item');
            place = 'consensus';
        } else {
            document.getElementById("consensus_navbarLi").style.width = "17%";
            document.getElementById("dG_navbarLi").style.width = "17%";
            document.getElementById("taxonomy_navbarLi").style.width = "17%";
            document.getElementById("userCheck_navbarLi").style.width = "0%";
            document.getElementById("download_navbarLi").style.width = "17%";
            \$("#consensus_info").css('display', 'block');
            place = 'consensus';
        }
        
        
        
        //OPEN ON CONSENSUS PAGE WITH FIRST VALUE
        navbarColors('#00008B', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA');
        navbarDisplay('block', 'none', 'none', 'none', 'none');
        res = firstValue.split("!");
        var infoTable = "<table><tr><th>Sequence</th><th>Length</th><th>GC%</th><th>T<sub>m (<sup>o</sup>C)</sub></th><th>Self Dimer &Delta;G (kcal/mol)</th><th>Hairpin &Delta;G (kcal/mol)</th></tr><tr><td>" + res[0] + "</td><td>" + res[1] + "</td><td>" + res[2] + "</td><td>" + res[3] + "</td><td>" + res[4] + "</td><td>" + res[5] + "</td></tr></table>";
        document.getElementById("pair_info").innerHTML = infoTable;
        
        combinedTrim = res[0];
        if (tableBlast === "none") {
            document.getElementById("taxonomy_intro").innerHTML = "<h2>Nt info</h2>";
            document.getElementById("taxonomy_contentF").style.height = "30px";
            document.getElementById("taxonomy_contentF").innerHTML = "<h3>The BLAST search against the nucleotide BLAST database did not produce any result for any of the PhyloPrimer primers.</h3>";
            tablePieChartSubF = 'zero';
            tableBlastSub = 'zero';
        } else {
            tablePieChartSubF = tablePieChartF.match("<" + combinedTrim + ">(.*?)" + "<>"); //primer pair pie chart
            tableBlastSub = tableBlast.match(">" + combinedTrim + "<table>" + "(.*?)" + "</table>"); //select table Blast
        }

        //dG texts
        textSelf = "<h3>Only the self dimers that showed a &Delta;G lower than 0 <sup>o</sup>C are reported here.</h3>";
        textHair = "<h3>Only the chairpins that showed a &Delta;G lower than 0 <sup>o</sup>C reported here.</h3>";
    });


    \$(document).ready(function () {
        //CONSENSUS TAB
        \$("#consensus_navbar").click(function () {
            navbarColors('#00008B', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA');
            navbarDisplay('block', 'none', 'none', 'none', 'none');
            var infoTable = "<table><tr><th>Sequence</th><th>Length</th><th>GC%</th><th>T<sub>m (<sup>o</sup>C)</sub></th><th>Self Dimer &Delta;G (kcal/mol)</th><th>Hairpin &Delta;G (kcal/mol)</th></tr><tr><td>" + res[0] + "</td><td>" + res[1] + "</td><td>" + res[2] + "</td><td>" + res[3] + "</td><td>" + res[4] + "</td><td>" + res[5] + "</td></tr></table>";
            document.getElementById("pair_info").innerHTML = infoTable;
            place = 'consensus';
        })

        //DELTA G TABS
        \$("#self_navbar").click(function () {
            navbarColors('#008CBA', '#00008B', '#00008B', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA');
            navbarDisplay('none', 'block', 'none', 'none', 'none');
            infodG = "<h2>&Delta;G - Self Dimers</h2>" + textSelf + res[6] + "</br>";
            document.getElementById("dG_content").innerHTML = infodG;
            place = 'self';
        })
        \$("#hair_navbar").click(function () {
            navbarColors('#008CBA', '#00008B', '#008CBA', '#008CBA', '#00008B', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA');
            navbarDisplay('none', 'block', 'none', 'none', 'none');
            infodG = "<h2>&Delta;G - Hairpins</h2>" + textHair + res[7] + "</br>";
            document.getElementById("dG_content").innerHTML = infodG;
            place = 'hair';
        })

        //NT TAXONOMY TABS
        \$("#domain_navbar").click(function () {
            navbarColors('#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#00008B', '#00008B', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA');
            navbarDisplay('none', 'none', 'block', 'none', 'none');
            if ((tablePieChartSubF !== 'zero') && (tableBlastSub !== 'zero')) {
                taxonomyPage(tablePieChartSubF, tableBlastSub, 'Domain');
            }
            place = 'domain';
        })
        \$("#phylum_navbar").click(function () {
            navbarColors('#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#00008B', '#008CBA', '#00008B', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA');
            navbarDisplay('none', 'none', 'block', 'none', 'none');
            if ((tablePieChartSubF !== 'zero') && (tableBlastSub !== 'zero')) {
                taxonomyPage(tablePieChartSubF, tableBlastSub, 'Phylum');
            }
            place = 'phylum';
        })
        \$("#class_navbar").click(function () {
            navbarColors('#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#00008B', '#008CBA', '#008CBA', '#00008B', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA');
            navbarDisplay('none', 'none', 'block', 'none', 'none');
            if ((tablePieChartSubF !== 'zero') && (tableBlastSub !== 'zero')) {
                taxonomyPage(tablePieChartSubF, tableBlastSub, 'Class');
            }
            place = 'class';
        })
        \$("#order_navbar").click(function () {
            navbarColors('#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#00008B', '#008CBA', '#008CBA', '#008CBA', '#00008B', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA');
            navbarDisplay('none', 'none', 'block', 'none', 'none');
            if ((tablePieChartSubF !== 'zero') && (tableBlastSub !== 'zero')) {
                taxonomyPage(tablePieChartSubF, tableBlastSub, 'Order');
            }
            place = 'order';
        })
        \$("#family_navbar").click(function () {
            navbarColors('#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#00008B', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#00008B', '#008CBA', '#008CBA', '#008CBA', '#008CBA');
            navbarDisplay('none', 'none', 'block', 'none', 'none');
            if ((tablePieChartSubF !== 'zero') && (tableBlastSub !== 'zero')) {
                taxonomyPage(tablePieChartSubF, tableBlastSub, 'Family');
            }
            place = 'family';
        })
        \$("#genus_navbar").click(function () {
            navbarColors('#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#00008B', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#00008B', '#008CBA', '#008CBA', '#008CBA');
            navbarDisplay('none', 'none', 'block', 'none', 'none');
            if ((tablePieChartSubF !== 'zero') && (tableBlastSub !== 'zero')) {
                taxonomyPage(tablePieChartSubF, tableBlastSub, 'Genus');
            }
            place = 'genus';
        })
        \$("#species_navbar").click(function () {
            navbarColors('#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#00008B', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#00008B', '#008CBA', '#008CBA');
            navbarDisplay('none', 'none', 'block', 'none', 'none');
            if ((tablePieChartSubF !== 'zero') && (tableBlastSub !== 'zero')) {
                taxonomyPage(tablePieChartSubF, tableBlastSub, 'Species');
            }
            place = 'species';
        })

        //USER BLAST TAB
        \$("#userCheck_navbar").click(function () {
            navbarColors('#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#00008B', '#008CBA');
            navbarDisplay('none', 'none', 'none', 'block', 'none');
            if (tableUser === "none") {
                    document.getElementById("taxonomy_contentF").style.height = "30px";
                    document.getElementById("taxonomy_contentF").innerHTML = "<h3>The BLAST search against the your fasta sequences did not produce any result for any of the PhyloPrimer primers.</h3>";
                } else {
                    infoUserSub = tableUser.match(">" + combinedTrim + "<" + "(.*?)" + "><table>"); //select user info Blast
                    if (infoUserSub === null) {
                        document.getElementById("userCheck_content").innerHTML = "<h3>PhyloPrimer did not find any match between your fasta file and these primers.</h3>";
                        document.getElementById("userCheck_table").innerHTML = "";
                    } else {
                        infoUserData = infoUserSub[1].split(";");
                        if (infoUserData[1] == 0) {
                            document.getElementById("userCheck_content").innerHTML = "<h3>You submitted " + infoUserData[0] + " DNA sequences which PhyloPrimer blasted against these oligos. Of all the sequences, none was hit from the selected oligo. Please find the details of the matches in the following table.</h3>";
                        } else if (infoUserData[1] == 1) {
                            document.getElementById("userCheck_content").innerHTML = "<h3>You submitted " + infoUserData[0] + " DNA sequences which PhyloPrimer blasted against these oligos. Of all the sequences, " + infoUserData[1] + " was hit from the selected oligo. Please find the details of the matches in the following table.</h3>";
                        } else {
                            document.getElementById("userCheck_content").innerHTML = "<h3>You submitted " + infoUserData[0] + " DNA sequences which PhyloPrimer blasted against these oligos. Of all the sequences, " + infoUserData[1] + " were hit from the selected oligo. Please find the details of the matches in the following table.</h3>";
                        }
                        tableUserSub = tableUser.match(">" + combinedTrim + "<" + infoUserSub[1] + "><table>(.*?)" + "</table>"); //select user table Blast
                        document.getElementById("userCheck_table").innerHTML = "<table><tr><th>Sequence</th><th>Alignment</th><th>Alignment start</th><th>Alignment end</th></tr>" + tableUserSub[1] + "</table>";
                    }
                }
            place = 'user';
        })

        //DOWNLOAD TAB
        \$("#download_navbar").click(function () {
            navbarColors('#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#008CBA', '#00008B');
            navbarDisplay('none', 'none', 'none', 'none', 'block');
            place = 'download';
        })
    });

        \$('input:radio[name=pair]').change(function () {

            res = this.value.split("!");
            combinedTrim = res[0];
            
            if (tableBlast === "none") {
                document.getElementById("taxonomy_intro").innerHTML = "<h2>Nt info</h2>";
                document.getElementById("taxonomy_contentF").style.height = "30px";
                document.getElementById("taxonomy_contentF").innerHTML = "<h3>The BLAST search against the nucleotide BLAST database did not produce any result for any of the PhyloPrimer primers.</h3>";
                tablePieChartSubF = 'zero';
                tableBlastSub = 'zero';
            } else {
                tablePieChartSubF = tablePieChartF.match("<" + combinedTrim + ">(.*?)" + "<>"); //primer pair pie chart
                tableBlastSub = tableBlast.match(">" + combinedTrim + "<table>" + "(.*?)" + "</table>"); //select table Blast
            }

            if (place === 'consensus') {
                var infoTable = "<table><tr><th>Sequence</th><th>Length</th><th>GC%</th><th>T<sub>m (<sup>o</sup>C)</sub></th><th>Self Dimer &Delta;G (kcal/mol)</th><th>Hairpin &Delta;G (kcal/mol)</th></tr><tr><td>" + res[0] + "</td><td>" + res[1] + "</td><td>" + res[2] + "</td><td>" + res[3] + "</td><td>" + res[4] + "</td><td>" + res[5] + "</td></tr></table>";
                document.getElementById("pair_info").innerHTML = infoTable;
            } else if (place === 'self') {
                infodG = "<h2>&Delta;G - Self Dimers</h2>" + textSelf + res[6] + "</br>";
                document.getElementById("dG_content").innerHTML = infodG;
            } else if (place === 'hair') {
                infodG = "<h2>&Delta;G - Hairpins</h2>" + textHair + res[7] + "</br>";
                document.getElementById("dG_content").innerHTML = infodG;
            } else if (place === 'domain') {
                if ((tablePieChartSubF !== 'zero') && (tableBlastSub !== 'zero')) {
                    taxonomyPage(tablePieChartSubF, tableBlastSub, 'Domain');
                }
            } else if (place === 'phylum') {
                if ((tablePieChartSubF !== 'zero') && (tableBlastSub !== 'zero')) {
                    taxonomyPage(tablePieChartSubF, tableBlastSub, 'Phylum');
                }
            } else if (place === 'class') {
                if ((tablePieChartSubF !== 'zero') && (tableBlastSub !== 'zero')) {
                    taxonomyPage(tablePieChartSubF, tableBlastSub, 'Class');
                }
            } else if (place === 'order') {
                if ((tablePieChartSubF !== 'zero') && (tableBlastSub !== 'zero')) {
                    taxonomyPage(tablePieChartSubF, tableBlastSub, 'Order');
                }
            } else if (place === 'family') {
                if ((tablePieChartSubF !== 'zero') && (tableBlastSub !== 'zero')) {
                    taxonomyPage(tablePieChartSubF, tableBlastSub, 'Family');
                }
            } else if (place === 'genus') {
                if ((tablePieChartSubF !== 'zero') && (tableBlastSub !== 'zero')) {
                    taxonomyPage(tablePieChartSubF, tableBlastSub, 'Genus');
                }
            } else if (place === 'species') {
                if ((tablePieChartSubF !== 'zero') && (tableBlastSub !== 'zero')) {
                    taxonomyPage(tablePieChartSubF, tableBlastSub, 'Species');
                }
            } else if (place === 'user') {
                
                if (tableUser === "none") {
                    document.getElementById("taxonomy_contentF").style.height = "30px";
                    document.getElementById("taxonomy_contentF").innerHTML = "<h3>The BLAST search against the your fasta sequences did not produce any result for any of the PhyloPrimer primers.</h3>";
                } else {
                    infoUserSub = tableUser.match(">" + combinedTrim + "<" + "(.*?)" + "><table>"); //select user info Blast
                    if (infoUserSub === null) {
                        document.getElementById("userCheck_content").innerHTML = "<h3>PhyloPrimer did not find any match between your fasta file and these primers.</h3>";
                        document.getElementById("userCheck_table").innerHTML = "";
                    } else {
                        infoUserData = infoUserSub[1].split(";");
                        if (infoUserData[1] == 0) {
                            document.getElementById("userCheck_content").innerHTML = "<h3>You submitted " + infoUserData[0] + " DNA sequences which PhyloPrimer blasted against these oligos. Of all the sequences, none was hit from the selected oligo. Please find the details of the matches in the following table.</h3>";
                        } else if (infoUserData[1] == 1) {
                            document.getElementById("userCheck_content").innerHTML = "<h3>You submitted " + infoUserData[0] + " DNA sequences which PhyloPrimer blasted against these oligos. Of all the sequences, " + infoUserData[1] + " was hit from the selected oligo. Please find the details of the matches in the following table.</h3>";
                        } else {
                            document.getElementById("userCheck_content").innerHTML = "<h3>You submitted " + infoUserData[0] + " DNA sequences which PhyloPrimer blasted against these oligos. Of all the sequences, " + infoUserData[1] + " were hit from the selected oligo. Please find the details of the matches in the following table.</h3>";
                        }
                        tableUserSub = tableUser.match(">" + combinedTrim + "<" + infoUserSub[1] + "><table>(.*?)" + "</table>"); //select user table Blast
                        document.getElementById("userCheck_table").innerHTML = "<table><tr><th>Sequence</th><th>Alignment</th><th>Alignment start</th><th>Alignment end</th></tr>" + tableUserSub[1] + "</table>";
                    }
                }
            }
            //no download because always the same
    });

    //FUNCTIONS:

    //CHANGE VISUALISED PAGE
    function navbarDisplay(N1, N2, N3, N4, N5) {
        \$("#consensus_info").css('display', N1);
        \$("#dG_info").css('display', N2);
        \$("#taxonomy_info").css('display', N3);
        \$("#userCheck_info").css('display', N4);
        \$("#download_info").css('display', N5);
    }

    //CHANGE NAVBAR COLOR
    function navbarColors(N1, N2, N3, N4, N5, N6, N7, N8, N9, N10, N11, N12, N13, N14, N15) {
        \$("#consensus_navbar").css('background-color', N1);
        \$("#dG_navbar").css('background-color', N2);
        \$("#self_navbar").css('background-color', N3);
        \$("#hair_navbar").css('background-color', N5);
        \$("#taxonomy_navbar").css('background-color', N6);
        \$("#domain_navbar").css('background-color', N7);
        \$("#phylum_navbar").css('background-color', N8);
        \$("#class_navbar").css('background-color', N9);
        \$("#order_navbar").css('background-color', N10);
        \$("#family_navbar").css('background-color', N11);
        \$("#genus_navbar").css('background-color', N12);
        \$("#species_navbar").css('background-color', N13);
        \$("#userCheck_navbar").css('background-color', N14);
        \$("#download_navbar").css('background-color', N15);
    }

    //RENDER PIE CHARTS
    function renderChart(env, oligo, tablePieChartDiv, rank) {
        if (rank == "Domain") {
            var data = tablePieChartDiv[0].split(",");
        } else if (rank == "Phylum") {
            var data = tablePieChartDiv[1].split(",");
        } else if (rank == "Class") {
            var data = tablePieChartDiv[2].split(",");
        } else if (rank == "Order") {
            var data = tablePieChartDiv[3].split(",");
        } else if (rank == "Family") {
            var data = tablePieChartDiv[4].split(",");
        } else if (rank == "Genus") {
            var data = tablePieChartDiv[5].split(",");
        } else if (rank == "Species") {
            var data = tablePieChartDiv[6].split(",");
        }
        var dataPoints = new Array();
        var arrayLength = data.length;
        for (var i = 0; i < arrayLength; i++) {
            var dataNew = data[i].split("-");
            dataNew[0] = dataNew[0].replace(/\\[sub1\\]/g, "'");
            dataNew[0] = dataNew[0].replace(/\\[sub2\\]/g, ",");
            dataNew[0] = dataNew[0].replace(/\\[sub3\\]/g, "(");
            dataNew[0] = dataNew[0].replace(/\\[sub4\\]/g, ")");
            dataNew[0] = dataNew[0].replace(/\\[sub5\\]/g, ":");
            dataNew[0] = dataNew[0].replace(/\\[sub7\\]/g, "*");
            dataNew[0] = dataNew[0].replace(/\\[sub8\\]/g, "<");
            dataNew[0] = dataNew[0].replace(/\\[sub9\\]/g, ">");
            dataNew[0] = dataNew[0].replace(/\\[sub10\\]/g, "-");
            dataNew[0] = dataNew[0].replace(/\\[sub11\\]/g, "+");
            dataNew[0] = dataNew[0].replace(/\\[sub12\\]/g, "'");
            dataNew[0] = dataNew[0].replace(/\\[sub13\\]/g, "#");
            dataNew[0] = dataNew[0].replace(/\\[sub14\\]/g, "&");
            dataNew[0] = dataNew[0].replace(/\\[sub15\\]/g, "^");
            dataNew[0] = dataNew[0].replace(/\\[sub16\\]/g, "/");
            var objArray = {
                y: dataNew[1],
                label: dataNew[0],
                color: dataNew[2]
            };
            dataPoints.push(objArray);
        }
        var chart = new CanvasJS.Chart(env, {
            animationEnabled: false,
            title: {
                text: oligo, // primer name
                fontSize: 20,
                fontFamily: "tahoma"
            },
            data: [{
                type: "pie",
                startAngle: 340,
                yValueFormatString: "##0.00\"%\"",
                indexLabel: "{label}",
                dataPoints
            }]
        });
        chart.render();
    }

    //TAXONOMY PAGE

    function taxonomyPage(tablePieChartSubF, tableBlastSub, rank) {
        rankPass = rank;
        if ((tablePieChartSubF === null) && (tableBlastSub === null)) {
            document.getElementById("taxonomy_intro").innerHTML = "<h2>GenBank BLAST</h2>";
            document.getElementById("taxonomy_contentF").style.height = "30px";
            document.getElementById("taxonomy_contentF").innerHTML = "<h3>The BLAST search against the nucleotide BLAST database did not match any entry for either the forward or the reverse primer.</h3>";
            document.getElementById("taxonomy_table_mess").innerHTML = "";
            document.getElementById("taxonomy_table").innerHTML = "";
        } else if ((tablePieChartSubF === null) && (tableBlastSub !== null)) {
            document.getElementById("taxonomy_intro").innerHTML = "<h2>GenBank BLAST - " + rank + "</h2>";
            document.getElementById("taxonomy_contentF").style.height = "30px";
            document.getElementById("taxonomy_contentF").innerHTML = "<h3>None of the entries in the nucleotide BLAST database matched with both the forward and the reverse primers.</h3>";
            tableBlastSubNew = tableBlastSub[1] + "</table>";
            document.getElementById("taxonomy_table_mess").innerHTML = "<h3>In the following table you can find the details of the BLAST matches.</h3>";
            document.getElementById("taxonomy_table").innerHTML = "<table><tr><th>Genbank sequence</th><th>Alignment</th><th>Alignment start</th><th>Alignment end</th><th class ='D'>Domain</th><th class ='P'>Phylum</th><th class ='C'>Class</th><th class ='O'>Order</th><th class ='F'>Family</th><th class ='G'>Genus</th><th class ='S'>Species</th></tr>" + tableBlastSubNew;
            if (rank === 'Domain') {
                \$(".D").show();
            } else {
                \$(".D").hide();
            }
            if (rank === 'Phylum') {
                \$(".P").show();
            } else {
                \$(".P").hide();
            }
            if (rank === 'Class') {
                \$(".C").show();
            } else {
                \$(".C").hide();
            }
            if (rank === 'Order') {
                \$(".O").show();
            } else {
                \$(".O").hide();
            }
            if (rank === 'Family') {
                \$(".F").show();
            } else {
                \$(".F").hide();
            }
            if (rank === 'Genus') {
                \$(".G").show();
            } else {
                \$(".G").hide();
            }
            if (rank === 'Species') {
                \$(".S").show();
            } else {
                \$(".S").hide();
            }
        } else {
            document.getElementById("taxonomy_intro").innerHTML = "<h2>GenBank BLAST - " + rank + "</h2>";
            document.getElementById("taxonomy_contentF").style.height = "500px";
            document.getElementById("taxonomy_contentF").style.width =
            document.getElementById("taxonomy_info").offsetWidth;
            tablePieChartDivF = tablePieChartSubF[1].split(":");
            renderChart("taxonomy_contentF", "", tablePieChartDivF, rank);
            var c = document.getElementsByClassName("canvasjs-chart-canvas")[0];
            var ctx = c.getContext("2d");
            ctx.fillStyle = 'white';
            ctx.fillRect(0, 489, (document.getElementById("taxonomy_info").offsetWidth/5), 40);
            ctx.stroke();
            tableBlastSubNew = tableBlastSub[1] + "</table>";
            document.getElementById("taxonomy_table_mess").innerHTML = "<h3>In the following table you can find the details of the BLAST matches.</h3>";
            document.getElementById("taxonomy_table").innerHTML = "<table><tr><th>Genbank sequence</th><th>Alignment</th><th>Alignment start</th><th>Alignment end</th><th class ='D'>Domain</th><th class ='P'>Phylum</th><th class ='C'>Class</th><th class ='O'>Order</th><th class ='F'>Family</th><th class ='G'>Genus</th><th class ='S'>Species</th></tr>" + tableBlastSubNew;
            if (rank === 'Domain') {
                \$(".D").show();
            } else {
                \$(".D").hide();
            }
            if (rank === 'Phylum') {
                \$(".P").show();
            } else {
                \$(".P").hide();
            }
            if (rank === 'Class') {
                \$(".C").show();
            } else {
                \$(".C").hide();
            }
            if (rank === 'Order') {
                \$(".O").show();
            } else {
                \$(".O").hide();
            }
            if (rank === 'Family') {
                \$(".F").show();
            } else {
                \$(".F").hide();
            }
            if (rank === 'Genus') {
                \$(".G").show();
            } else {
                \$(".G").hide();
            }
            if (rank === 'Species') {
                \$(".S").show();
            } else {
                \$(".S").hide();
            }
           
        }
    }

//resize of the pie chart
\$(document).ready(function () {
    \$(window).resize(function() {
        document.getElementById("taxonomy_contentF").style.width = document.getElementById("taxonomy_info").offsetWidth;
        renderChart("taxonomy_contentF", "", tablePieChartDivF, rankPass);
        var c = document.getElementsByClassName("canvasjs-chart-canvas")[0];
        var ctx = c.getContext("2d");
        ctx.fillStyle = 'white';
        ctx.fillRect(0, 489, (document.getElementById("taxonomy_info").offsetWidth/5), 40);
        ctx.stroke();
    });
});

//taxonomy dropdown
function autocompleteDrop(id) {
    var text = id.value;
    \$.ajax({
    url: 'https://www.cerealsdb.uk.net/cerealgenomics/cgi-bin/dynamicDropdown_pp.pl',
    data: {
        'text': text,
    },
    type: 'POST',
    success: function (resp) {
        
        if (resp.result == "empty") {
            \$('.autocomplete-items').remove(); //remove previous dropdown when no taxa were found
        } else {
            arr = resp.result.split("!");
            var a, b, i, val = text;
            \$('.autocomplete-items').remove(); //remove previous dropdown
            a = document.createElement("DIV");
            a.setAttribute("id", text + "autocomplete-list");
            a.setAttribute("class", "autocomplete-items");
            id.parentNode.appendChild(a);
            
            for (i = 0; i < (arr.length-1); i++) {
                /*check if the item starts with the same letters as the text field value:*/
                    if (arr[i].substr(0, val.length).toUpperCase() == val.toUpperCase()) {
                        /*create a DIV element for each matching element:*/
                            b = document.createElement("DIV");
                        b.setAttribute("class", "autocomplete-itemsSec");
                        b.innerHTML += "<input id='dropTaxa' value='" + arr[i] + "'>";
                        b.addEventListener("click", function(e) {
                            /*insert the value for the autocomplete text field:*/
                                id.value = this.getElementsByTagName("input")[0].value;
                            /*close the list of autocompleted values,
                            (or any other open lists of autocompleted values:*/
                            //\$('.autocomplete-items').remove(); //remove previous dropdown when no taxa were found
                            //closeAllLists();
                        });
                        a.appendChild(b);
                    }
            }
            
        }
        
    },
    error: function () { alert("An error occoured"); }
    });
    
    //remove dropdown when click
    /*execute a function when someone clicks in the document:*/
    document.addEventListener("click", function (e) {
        console.log("click");
        \$('.autocomplete-items').remove(); //remove previous dropdown when no taxa were found
    });
    
    \$('#pair,#search,#oligoEmpty').mouseenter(function() {
        \$('.autocomplete-itemsSec').remove(); //remove previous dropdown when no taxa were found
    });
    
    
    
}

var num = 0;
var max = 0;

//add taxonomy
function addTax() {
    num++;
    max++;
    console.log(num + " and " + max);
    if (max < 10) {
        \$(".input_fields_wrap").append("<div class='autocomplete_wrap" + num + "'><div class='autocomplete" + num + "' style='width:300px;display:inline-block;position:relative;'><input style='width:300px;' type='text' onkeyup='dropdown(" + num +")' id='searchTax" + num + "' class='searchTaxon'></div>&nbsp;<button id='clear_taxon" + num + "' onclick='removeTax(" + num + ")'> - </button></div>"); //add input box
    } else {
        num--;
        max--;
    }
    
}

\$( "#searchTax" ).keyup(function() {
    autocompleteDrop(document.getElementById("searchTax"));
});

function dropdown(seq) {
    el = "searchTax" + seq;
    autocompleteDrop(document.getElementById(el));
}

function removeTax(seq) {
    max--;
    el = ".autocomplete_wrap" + seq;
    //remove the object - the input filed
    \$(el).remove(); //remove previous dropdown when no taxa were found
}

</script>

</body>
END HTML

print $html;




