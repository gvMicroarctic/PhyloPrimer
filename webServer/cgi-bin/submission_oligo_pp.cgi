#! /usr/bin/perl -w

use strict;
use warnings;
use CGI;

# #connected to primer_input.html, probe_input.html, oligo_input.html and check_input.html : display submission page

my $cgi = CGI->new;

##get name of the folder
my $defSet = $cgi->param('defSet');

#get if S (sequence) or T (tree)
my $input_kind = chop($defSet); #get last letter and understand if a, b or c
my $input_STC = chop($defSet); #get last letter and understand if T or S or C
my $index = 1;
my $defSet0;

if ($defSet =~ /_/) {
    ($defSet0) = $defSet =~ /(.*)_/;
    ($index) = $defSet =~ /_(.*)/;
    $index++;
    $defSet = "\"" . $defSet0 . "_" . $index . $input_STC . $input_kind . "\"";
} else {
    $defSet = "\"" .  $defSet . "_" . $index . $input_STC . $input_kind . "\"";
}

$input_STC = "\"" . $input_STC . "\"";

my $html = <<"END HTML";
Content-Type: text/html

<HTML>
<link rel='stylesheet' href='../phyloprimer/phyloprimer.css'>
<HEAD>
<meta http-equiv='Content-Type' content='text/html; charset=utf-8'/>
<TITLE>PhyloPrimer submission</TITLE>
<style>
.wrapper1 {
    width:100%;
    height:70%;
    text-align: center;
    position: relative;
}
.wrapper2 {
    width:100%;
    height:30%;
    text-align: center;
    position: relative;
}

img {
    width:40%;
    height:auto;
    top: 0;
    left: 0;
    right: 0;
    bottom: 0;
    display: inline-block;
    position: absolute;
    top: 50%;
    left: 50%;
    transform: translate(-50%, -50%);
}
.button {
    background-color: #008CBA;
    border: none;
    color: white;
    padding: 15px 32px;
    text-align: center;
    text-decoration: none;
    font-size: 16px;
    cursor: pointer;
    display: inline-block;
    position: absolute;
    top: 10%;
    left: 50%;
    transform: translate(-50%, -50%);
}
.button:hover {
    background-color: #00008B;
}
</style>
</HEAD>
<body>
<div class="wrapper1">
<img src="../phyloprimer/figures/finalResult.png" alt="submission">
</div>
<div class="wrapper2">
<input type="button" name="back_tree" class="button" value="Go back to Dynamic Selection page" id="back_tree">
</div>
<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.4.0/jquery.min.js"></script>
<script type="application/javascript">
folder = $defSet;
STC = $input_STC;
if (STC !== "T") {
    \$('.button').hide();
}
\$("#back_tree").click(function () {
    var newUrl = "https://www.cerealsdb.uk.net/cerealgenomics/cgi-bin/tree_input.cgi?defSet=" + folder;
    document.location = newUrl;
})
</script>
</body>

END HTML
print $html;
