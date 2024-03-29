#! /usr/bin/perl -w

#use strict;
use warnings;
use CGI;

#connected to submission_primer_pp.cgi : display Dynamic Selction page after one or more clusters have been previously selected

#set path to folders:
my $path_html = 'path_to_html_folder';

my $cgi = CGI->new;

##get name of the folder
my $defSet = $cgi->param('defSet');
my $defSetNew = $defSet;
#thre will be a hyphen
#get if S (sequence) or T (tree)
my $input_kind = chop($defSetNew); #get last letter and understand if a, b or c
my $input_ST = chop($defSetNew); #get last letter and understand if T or S
#re-construct folder name
my $folder = $defSetNew . $input_kind;

#open file with all information
my $fileIN = (substr($defSet, 0, 8)) . ".treeInfo";
open(IN, "<${path_html}/uploadsPhyloprimer/${defSet}/${fileIN}") or die;
my $intro;
my $result;
my $number;
my $userTable;
my $clusterTable;
my $taxonomy;
while (defined(my $input = <IN>)) {
    chomp($input);
    my @info = split(/\t/, $input); #intro paragraph
    if ($info[0] eq "intro") {
        $intro = $info[1];
    } elsif ($info[0] eq "result") { #phylogenetic tree
        $result = $info[1];
    } elsif ($info[0] eq "number") { #number of nodes
        $number = $info[1];
    } elsif ($info[0] eq "userTable") { #user table
        $userTable = $info[1];
    } elsif ($info[0] eq "clusterTable") { #cluster table
        $clusterTable = $info[1];
    } elsif ($info[0] eq "taxonomy") { #taxonomy
        $taxonomy = $info[1];
    }
}
close(IN);

#open file with pre-selected sequences
my $fileACC = (substr($defSet, 0, 8)) . ".loadAcc";
open(IN, "<${path_html}/uploadsPhyloprimer/${defSet}/${fileACC}") or die;

#my $accessionPre = ",CP012749.1(red),NOYX01000015.1(blue),CP015235.1(green)";
my $accessionPre;

my $tableCluster = "<table><tr><th>Cluster</th><th>Sequences</th></tr>";

while (defined(my $input = <IN>)) {
    chomp($input);
    if ($input =~ /^CLUSTER/) {
        my @info = split(/\t/,$input);
        $tableCluster .= "<tr><td><strong style='color: " . $info[3] . ";'>" . $info[2] . "</strong></td><td>" . $info[4] . "</td></tr>";
    } else {
        $accessionPre = "\"" . $input . "\"";
    }
}
close(IN);

$tableCluster .= "</table>";

$intro = "\"" . $intro . "\"";
$result = "\"" . $result . "\"";
$number = "\"" . $number . "\"";
$userTable = "\"" . $userTable . "\"";
$clusterTable = "\"" . $clusterTable . "\"";
$taxonomy = "\"" . $taxonomy . "\"";
$folder = "\"" . $folder . "\"";

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
            width: 740px;
            height: auto;
        }


        #wrap_phylocanvas {
            /*margin-left: 20px;
            margin-right: 20px;*/
            clear: both;
            vertical-align: bottom;
            /*margin-top: 50px;*/
            /*outline: 1px solid red;*/
            overflow: hidden;
            width: 100%;
        }

        #wrap_phylocanvas_cont {
            margin: 3px;
            width: 99%;
            /*outline: 1px solid black;*/
        }

        table {
            border-collapse: collapse;
            border-left: 3px solid #00008B;
            border-right: 3px solid #00008B;
        }

        td,
        th {
            border-bottom: 1px solid #ddd;
            border-left: 1px solid #ddd;
            padding: 15px;
        }

        textarea {
            vertical-align: top;
            resize: none;
        }

        #in_sequence {
            width: 800px;
            height: 200px;
            overflow: scroll;
            margin-left: 25px;
        }

        .button {
            background-color: #008CBA;
            border: none;
            color: white;
            padding: 15px 32px;
            text-align: center;
            text-decoration: none;
            display: inline-block;
            font-size: 16px;
            margin: 4px 2px;
            cursor: pointer;
        }

        .button:hover {
            background-color: #00008B !important;
        }

        #species_tree {
            /*background-color: #00008B;*/

        }

        /* Layout different sections*/
        #container {
            display: flex;
        }

        #container_title {
            color: white;
            background-color: #00008B;
            margin: 0.15em 0.75em 0 0;
            -webkit-writing-mode: vertical-lr;
            -ms-writing-mode: tb-lr;
            writing-mode: vertical-lr;
            text-orientation: upright;
            -webkit-font-feature-settings: "vkrn", "vpal";
            text-align: center;
        }

        #container1 {
            display: none;
        }

        #container_title1 {
            color: white;
            background-color: #00008B;
            margin: 0.15em 0.75em 0 0;
            -webkit-writing-mode: vertical-lr;
            -ms-writing-mode: tb-lr;
            writing-mode: vertical-lr;
            text-orientation: upright;
            -webkit-font-feature-settings: "vkrn", "vpal";
            text-align: center;
        }

        #container_title2 {
            color: white;
            background-color: #00008B;
            margin: 0.15em 0.75em 0 0;
            -webkit-writing-mode: vertical-lr;
            -ms-writing-mode: tb-lr;
            writing-mode: vertical-lr;
            text-orientation: upright;
            -webkit-font-feature-settings: "vkrn", "vpal";
            text-align: center;
        }

        #container_consensus {
            display: none;
            width: 100%;
            /* -ms-word-break: break-all;
            word-break: break-all; */
        }

        .break {
            -ms-word-break: break-all;
            word-break: break-all;
            font-size: 15px;
        }

        #container_title3 {
            color: white;
            background-color: #00008B;
            margin: 0.15em 0.75em 0 0;
            -webkit-writing-mode: vertical-lr;
            -ms-writing-mode: tb-lr;
            writing-mode: vertical-lr;
            text-orientation: upright;
            -webkit-font-feature-settings: "vkrn", "vpal";
            text-align: center;
        }

        /*top part above phylogenetic tree*/

        #top_phylo {
            /* background: red; */
            height: 70px;
            width: 95%;
            margin-left: 3px;

        }

        #taxonomy_button {
            width: 80%;
            float: left;
        }

        #container_phylo {
            display: flex;
        }

        /* slider */
        #slider {
            width: 20%;
            height: 20px;
            text-align: center;
            background: #008CBA;
            margin: 50 auto;
            display: inline-block;
            margin-top: 20px;
        }

        #zoomMinus {
            height: 100%;
            float: left;
            width: 6%;
            background: #008CBA;
            font-size: 30px;
            display: -webkit-flexbox;
            display: -ms-flexbox;
            display: -webkit-flex;
            display: flex;
            -webkit-flex-align: center;
            -ms-flex-align: center;
            -webkit-align-items: center;
            align-items: center;
            justify-content: center;
        }

        #zoomPlus {
            height: 100%;
            float: right;
            width: 6%;
            background: #008CBA;
            font-size: 20px;
            display: -webkit-flexbox;
            display: -ms-flexbox;
            display: -webkit-flex;
            display: flex;
            -webkit-flex-align: center;
            -ms-flex-align: center;
            -webkit-align-items: center;
            align-items: center;
            justify-content: center;
        }

        #zoomBody {
            height: 15px;
            margin: 0 auto;
            background: #008CBA;
            display: inline-block;
            width: 88%;
        }

        #zoom {
            -webkit-appearance: none;
            width: 100%;
            height: 15px;
            background: #008CBA;
            outline: none;
            opacity: 1;
            -webkit-transition: .2s;
            transition: opacity .2s;
        }

        #phylocanvas {
            /*background: #008CBA;*/
            height: 10px;
            width: 100%;
            padding-top: 0px;
            padding-bottom: 60px;
            overflow: hidden;
            float: left;
            /*outline: 5px solid green;*/
        }

        #myCanvas {
            /*outline: 1px dashed black;*/
        }

        #zoom:hover {
            /*opacity: 1;*/
        }

        #zoom::-webkit-slider-thumb {
            -webkit-appearance: none;
            appearance: none;
            width: 25px;
            height: 25px;
            background: #00008B;
            opacity: 1;
            cursor: pointer;
        }

        #zoom::-moz-range-thumb {
            width: 25px;
            height: 25px;
            background: #00008B;
            opacity: 1;
            cursor: pointer;
        }

        /*remove hover when scrolling*/
        /* .disable-hover,
        .disable-hover * {
            pointer-events: none !important;
        } */

        /* clear button */
        #clear_button {
            display: flex;
        }

        /* disable red border on firefox */
        input:invalid {
            box-shadow: none;
        }

        a:link {
            text-decoration: none;
            color: #00008B;
        }

        a:visited {
            text-decoration: none;
            color: #00008B;
        }

        a:hover {
            text-decoration: underline;
            color: #00008B;
        }

        a:active {
            text-decoration: underline;
            color: #00008B;
        }
</style>
<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.4.0/jquery.min.js"></script>
<script src="https://www.cerealsdb.uk.net/cerealgenomics/phyloprimer/phylocanvas.js"></script>
<script src="https://unpkg.com/infinite-scroll@3/dist/infinite-scroll.pkgd.min.js"></script>
<script type="application/javascript">

let resp = new Object();
resp.intro = $intro;
resp.result = $result;
resp.number = $number;
resp.userTable = $userTable;
resp.clusterTable = $clusterTable;
resp.taxonomy = $taxonomy;

accessionPre = $accessionPre; //pre-selected sequences
folder = $folder; //folder

\$(document).ready(function () {
        \$("#species_tree").css('background-color', '#00008B');
                   
        //data to phylogenetic design
        global_x = "XX";
        //print taxonomy buttons only if taxonomy information
        if (resp.taxonomy == "TRUE") {
            document.getElementById('taxonomy_button').style.display = 'block'; //show taxonomy button
        } else {
            document.getElementById('taxonomy_button').style.display = 'none'; //hide taxonomy button
        }

        //document.getElementById('phylo_error').style.display = 'none'; //hide the error message
        document.getElementById('consensus_button').style.display = 'block'; //show consensus_button
        var newick_all = resp.result;
        newick = newick_all.split("&");
        document.getElementById('phylocanvas_intro').innerHTML = resp.intro;
        //Load and draw tree
        document.getElementById('overlay').style.display = 'none'; //hide the loading wheel
        document.getElementById('wrap_phylocanvas').style.display = 'block';
        document.getElementById("container_title2").style.background = '#00008B';
        document.getElementById("phylocanvas").innerHTML = "";

        if (resp.userTable.includes("<td>")) {
            document.getElementById("userSeq_table").innerHTML = resp.userTable; //visualize this table only if there is data inside.
        }

        if (resp.clusterTable.includes("<td>")) {
            document.getElementById("clusterSeq_table").innerHTML = resp.clusterTable; //visualize this table only if there is data inside.
        }

        var rem = document.getElementById('phylocanvas').offsetHeight; //rem

        acc_number = resp.number;
        acc_height = (acc_number * 40) + 100;

        \$("#phylocanvas").css('height', acc_height);
        \$("#phylocanvas").css('width', '100%');
        phylocanvas = new PhyloCanvas.Tree('phylocanvas');
        phylocanvas.passmaxx('INITIAL'); ///
        phylocanvas.clusterSelected(accessionPre); //highligth accession numbers
        phylocanvas.setNodeSize(4);
        phylocanvas.setTreeType('rectangular');
        phylocanvas.nodeAlign = true;

        phylocanvas.load(newick[6]);
        phylocanvas.draw();

        var width = document.getElementById('phylocanvas').offsetWidth;
        var height = document.getElementById('phylocanvas').offsetHeight;

        new_height = (height * global_zoom); //global_zoom is the initial zoom
        \$("#phylocanvas").css('height', new_height);

        //value for slider
        slider_value = global_zoom * 100;
        slider_max = slider_value * 2;
        document.getElementById("slider").innerHTML = "<div id='zoomMinus' class='zoomAll'><img src='../phyloprimer/minus.png' width='100%' height='100%'></div><div id='zoomBody' class='zoomAll'><input id='zoom' type='range' min='0' max='" + slider_max + "' value='" + slider_value + "' onclick='getZoom()'></input></div><div id='zoomPlus' class='zoomAll'><img src='../phyloprimer/plus.png' width='100%' height='100%'></div>";

        global_x = global_maxx;
        count = 0 //translating highlighted to all resized tree
        //retrieved = phylocanvas.selectedNodes;
})

</script>
</HEAD>
<body>

<title>Dynamic selection</title>

<div id="overlay" class="hide">
    <div id="loader"></div>
</div>
<img src="../phyloprimer/figures/DynamicSelectionTitle.png" alt="title" id="title">
</br>
<hr>

<div id="loader_wrap">
    <div id="container_phylo">
        <h3 id="container_title2">&nbsp;&nbsp;&nbsp;DYNAMIC TREE&nbsp;&nbsp;&nbsp;</h3>
        <div id='phylo_error'></div>
        <div id="wrap_phylocanvas">
            <div id="wrap_phylocanvas_cont">
                <div id="phylocanvas_intro"></div>
                <h3>The pre-selected clusters are:</h3>
                $tableCluster
                </br>
                <div id="top_phylo">
                    <div id="taxonomy_button">
                        <input type="button" name="phylum_tree" class="button" value="Phylum" id="phylum_tree">
                        <input type="button" name="tree_deselect" class="button" value="Class" id="class_tree">
                        <input type="button" name="order_tree" class="button" value="Order" id="order_tree">
                        <input type="button" name="family_tree" class="button" value="Family" id="family_tree">
                        <input type="button" name="genus_tree" class="button" value="Genus" id="genus_tree">
                        <input type="button" name="species_tree" class="button" value="Species" id="species_tree">
                    </div>
                    <div id="slider"></div>
                </div>
                <div id="phylocanvas" width="300" height="30"></div>
                </br>
                </br>
                <canvas id="myCanvas" width="120" height="30"></canvas>

                </br>
                <input type="button" name="tree_select" class="button" value="Select all the nodes" id="tree_select"
                    onclick="select_all()">
                <input type="button" name="tree_deselect" class="button" value="Deselect all the nodes"
                    id="tree_deselect" onclick="deselect()">
                </br>
                </br>
            </div>
            <div id="userSeq_table"></div>
            <div id="clusterSeq_table"></div>
        </div>
    </div>
</div>

</br>

<div id="consensus_button">
    <input type="button" name="tree_acc" class="button" value="Create Consensus" id="tree_acc" onclick="show()">
</div>

<!-- consensus error -->
</br>
<div id="container_consensus" class="hide">
    <h3 id="container_title3">&nbsp;CONSENSUS&nbsp;</h3>
    <div id="consesus_message"></div>
</div>

</br>
<div id="oligo_button" class="hide">
    <input type="button" name="oligo_design" class="button" value="Oligo Design" id="oligo_design">
</div>


<script type="application/javascript">

function show() {
    document.getElementById('overlay').style.display = 'block';
}

//get zoom value from slider bar and re-draw the phylogenetic tree
function getZoom() {
    var zoom = (\$("#zoom").val() / 100);
    phylocanvas.imposeZoom(zoom); //send zoom value to phylocanvas
    phylocanvas.draw();
}

function select_all() {
    phylocanvas.allSelect();
}
function deselect() {
    phylocanvas.clearSelect();
}
\$("#phylum_tree").click(function () {
    retrieved = phylocanvas.selectedNodes;
    \$("#phylum_tree").css('background-color', '#00008B');
    \$("#class_tree").css('background-color', '#008CBA');
    \$("#order_tree").css('background-color', '#008CBA');
    \$("#family_tree").css('background-color', '#008CBA');
    \$("#genus_tree").css('background-color', '#008CBA');
    \$("#species_tree").css('background-color', '#008CBA');
    phylocanvas.passmaxx(global_x);
    phylocanvas.prevSelected(retrieved);
    phylocanvas.setNodeSize(4);
    phylocanvas.setTreeType('rectangular');
    phylocanvas.nodeAlign = true;
    phylocanvas.load(newick[1]); //speficy attributes before loading the tree
    phylocanvas.draw();
    document.getElementById("slider").innerHTML = "<div id='zoomMinus' class='zoomAll'><span>-</span></div><div id='zoomBody' class='zoomAll'><input id='zoom' type='range' min='0' max='" + slider_max + "' value='" + slider_value + "' onclick='getZoom()'></input></div><div id='zoomPlus' class='zoomAll'>+</div>";
});

\$("#class_tree").click(function () {
    retrieved = phylocanvas.selectedNodes;
    \$("#phylum_tree").css('background-color', '#008CBA');
    \$("#class_tree").css('background-color', '#00008B');
    \$("#order_tree").css('background-color', '#008CBA');
    \$("#family_tree").css('background-color', '#008CBA');
    \$("#genus_tree").css('background-color', '#008CBA');
    \$("#species_tree").css('background-color', '#008CBA');
    phylocanvas.passmaxx(global_x);
    phylocanvas.prevSelected(retrieved);
    phylocanvas.setNodeSize(4);
    phylocanvas.setTreeType('rectangular');
    phylocanvas.nodeAlign = true;
    phylocanvas.load(newick[2]); //speficy attributes before loading the tree
    phylocanvas.draw();
    document.getElementById("slider").innerHTML = "<div id='zoomMinus' class='zoomAll'><span>-</span></div><div id='zoomBody' class='zoomAll'><input id='zoom' type='range' min='0' max='" + slider_max + "' value='" + slider_value + "' onclick='getZoom()'></input></div><div id='zoomPlus' class='zoomAll'>+</div>";
});
\$("#order_tree").click(function () {
    count = 0 //translating highlighted to all resized tree
    retrieved = phylocanvas.selectedNodes;
    \$("#phylum_tree").css('background-color', '#008CBA');
    \$("#class_tree").css('background-color', '#008CBA');
    \$("#order_tree").css('background-color', '#00008B');
    \$("#family_tree").css('background-color', '#008CBA');
    \$("#genus_tree").css('background-color', '#008CBA');
    \$("#species_tree").css('background-color', '#008CBA');
    phylocanvas.passmaxx(global_x);
    phylocanvas.prevSelected(retrieved);
    phylocanvas.setNodeSize(4);
    phylocanvas.setTreeType('rectangular');
    phylocanvas.nodeAlign = true;
    phylocanvas.load(newick[3]); //speficy attributes before loading the tree
    phylocanvas.draw();
    document.getElementById("slider").innerHTML = "<div id='zoomMinus' class='zoomAll'><span>-</span></div><div id='zoomBody' class='zoomAll'><input id='zoom' type='range' min='0' max='" + slider_max + "' value='" + slider_value + "' onclick='getZoom()'></input></div><div id='zoomPlus' class='zoomAll'>+</div>";
});
\$("#family_tree").click(function () {
    count = 0 //translating highlighted to all resized tree
    retrieved = phylocanvas.selectedNodes;
    \$("#phylum_tree").css('background-color', '#008CBA');
    \$("#class_tree").css('background-color', '#008CBA');
    \$("#order_tree").css('background-color', '#008CBA');
    \$("#family_tree").css('background-color', '#00008B');
    \$("#genus_tree").css('background-color', '#008CBA');
    \$("#species_tree").css('background-color', '#008CBA');
    phylocanvas.passmaxx(global_x);
    phylocanvas.prevSelected(retrieved);
    phylocanvas.setNodeSize(4);
    phylocanvas.setTreeType('rectangular');
    phylocanvas.nodeAlign = true;
    phylocanvas.load(newick[4]); //speficy attributes before loading the tree
    phylocanvas.draw();
    document.getElementById("slider").innerHTML = "<div id='zoomMinus' class='zoomAll'><span>-</span></div><div id='zoomBody' class='zoomAll'><input id='zoom' type='range' min='0' max='" + slider_max + "' value='" + slider_value + "' onclick='getZoom()'></input></div><div id='zoomPlus' class='zoomAll'>+</div>";
});
\$("#genus_tree").click(function () {
    count = 0 //translating highlighted to all resized tree
    retrieved = phylocanvas.selectedNodes;
    \$("#phylum_tree").css('background-color', '#008CBA');
    \$("#class_tree").css('background-color', '#008CBA');
    \$("#order_tree").css('background-color', '#008CBA');
    \$("#family_tree").css('background-color', '#008CBA');
    \$("#genus_tree").css('background-color', '#00008B');
    \$("#species_tree").css('background-color', '#008CBA');
    phylocanvas.passmaxx(global_x);
    phylocanvas.prevSelected(retrieved);
    phylocanvas.setNodeSize(4);
    phylocanvas.setTreeType('rectangular');
    phylocanvas.nodeAlign = true;
    phylocanvas.load(newick[5]); //speficy attributes before loading the tree
    phylocanvas.draw();
    document.getElementById("slider").innerHTML = "<div id='zoomMinus' class='zoomAll'><span>-</span></div><div id='zoomBody' class='zoomAll'><input id='zoom' type='range' min='0' max='" + slider_max + "' value='" + slider_value + "' onclick='getZoom()'></input></div><div id='zoomPlus' class='zoomAll'>+</div>";
});
\$("#species_tree").click(function () {
    count = 0 //translating highlighted to all resized tree
    retrieved = phylocanvas.selectedNodes;
    \$("#phylum_tree").css('background-color', '#008CBA');
    \$("#class_tree").css('background-color', '#008CBA');
    \$("#order_tree").css('background-color', '#008CBA');
    \$("#family_tree").css('background-color', '#008CBA');
    \$("#genus_tree").css('background-color', '#008CBA');
    \$("#species_tree").css('background-color', '#00008B');
    phylocanvas.passmaxx(global_x);
    phylocanvas.prevSelected(retrieved);
    phylocanvas.setNodeSize(4);
    phylocanvas.setTreeType('rectangular');
    phylocanvas.nodeAlign = true;
    phylocanvas.load(newick[6]); //speficy attributes before loading the tree
    phylocanvas.draw();
    document.getElementById("slider").innerHTML = "<div id='zoomMinus' class='zoomAll'><span>-</span></div><div id='zoomBody' class='zoomAll'><input id='zoom' type='range' min='0' max='" + slider_max + "' value='" + slider_value + "' onclick='getZoom()'></input></div><div id='zoomPlus' class='zoomAll'>+</div>";
});

//pass data to consensus
//create file and submit confirmation
\$(document).ready(function () {
    \$("#tree_acc").click(function () {
        //merge all accession numbers together
        var accession = phylocanvas.selectedNodes.join();
        //data to perl script
        defSet = folder;
        \$.ajax({
            url: 'tree_consensus_pp.pl',
            data: {
                'defSet': defSet,
                'accession': accession
            },
            type: 'POST',
            success: function (resp) {
                //loader
                document.getElementById('overlay').style.display = 'none'; //hide the loading wheel
                if (resp.result == "SUCCESS") {
                    document.getElementById("container_title3").style.background = '#00008B';
                    if (resp.neg !== 'X') {
                        document.getElementById("consesus_message").innerHTML = "<h3>PhyloPrimer calculated two consensus sequences: one out of the selected sequences (positive consensus) and the other one out of the sequences that were not selected (negative consensus)." +
                            " The oligos will be designed on the consensus sequences and by default PhyloPrimer will report in the Result visualization mainly the oligos that mazimixe the difference between the two consensus sequences.</h3><ul><li class='break'>positive consensus: " +
                            resp.pos + "</li></br><li class='break'>negative consensus: " + resp.neg + "</li></ul><h3>If you are happy about the result, please click on the Oligo Design button and go to the next page.</h3>";
                    } else {
                        document.getElementById("consesus_message").innerHTML = "<h3>PhyloPrimer calculated one consensus sequence out of the selected sequences.</h3><ul><li class='break'>positive consensus: " +
                            resp.pos + "</li></br></ul><h3>If you are happy about the result, please click on the Oligo Design button and go to the next page.</h3>";
                    }
                    document.getElementById('container_consensus').style.display = 'flex';
                    document.getElementById('oligo_button').style.display = 'block'; // Proceed to the next page anyway
                } else if (resp.result == "ERROR_POS") {
                    document.getElementById("container_title3").style.background = '#EA2845';
                    document.getElementById('container_consensus').style.display = 'flex';
                    document.getElementById("consesus_message").innerHTML = "<h3>Please select the DNA sequences you want PhyloPrimer to use for the oligo design. Once selected from the dynamic tree, click on the Create Consensus button.</h3>";
                    document.getElementById('oligo_button').style.display = 'none'; // Proceed to the next page anyway
                } else if (resp.result == "ERROR_PERC") {
                    document.getElementById("container_title3").style.background = '#00008B';

                    if (resp.perc <= 20) {
                        if (resp.neg !== 'X') {
                            document.getElementById("consesus_message").innerHTML = "<h3>PhyloPrimer calculated two consensus sequences: one out of the selected sequences (positive consensus) and the other one out of the sequences that were not selected (negative consensus)." +
                                " The oligos will be designed on the consensus sequences and by default PhyloPrimer will report in the Result visualization mainly the oligos that mazimixe the difference between the two consensus sequences.<ul><li class='break'>positive consensus: " +
                                resp.pos + "</li><li class='break'>negative consensus: " + resp.neg + "</li></ul></h3>" + "<h3>The " + resp.perc + "% of the bases in the consensus sequence are represented by degenerated bases.</h3><h3>If you are happy about the result, please click on the Oligo Design button and go to the next page.</h3>";
                        } else {
                            document.getElementById("consesus_message").innerHTML = "<h3>PhyloPrimer calculated one consensus sequence out of the selected sequences.<ul><li class='break'>positive consensus: " +
                                resp.pos + "</li></ul></h3>" + "<h3>The " + resp.perc + "% of the bases in the consensus sequence are represented by degenerated bases.</h3><h3>If you are happy about the result, please click on the Oligo Design button and go to the next page.</h3>";
                        }
                    } else {
                        if (resp.neg !== 'X') {
                            document.getElementById("consesus_message").innerHTML = "<h3>PhyloPrimer calculated two consensus sequences: one out of the selected sequences (positive consensus) and the other one out of the sequences that were not selected (negative consensus)." +
                                " The oligos will be designed on the consensus sequences and by default PhyloPrimer will report in the Result visualization mainly the oligos that mazimixe the difference between the two consensus sequences.<ul><li class='break'>positive consensus: " +
                                resp.pos + "</li><li class='break'>negative consensus: " + resp.neg + "</li></ul></h3>" + "<h3><b style='color:blue';>WARNING:</b> The " + resp.perc + "% of the bases in the consensus sequence are represented by degenerated bases. You may have to create more then one oligo assay (e.g. primer pair) in order to target all the selected sequences.</h3><h3>If you are happy about the result, please click on the Oligo Design button and go to the next page.</h3>";
                        } else {
                            document.getElementById("consesus_message").innerHTML = "<h3>PhyloPrimer calculated one consensus sequence out of the selected sequences.<ul><li class='break'>positive consensus: " +
                                resp.pos + "</li></ul></h3>" + "<h3><b style='color:blue';>WARNING:</b> The " + resp.perc + "% of the bases in the consensus sequence are represented by degenerated bases. You may have to create more then one oligo assay (e.g. primer pair) in order to target all the selected sequences.</h3><h3>If you are happy about the result, please click on the Oligo Design button and go to the next page.</h3>";
                        }
                    }

                    document.getElementById('container_consensus').style.display = 'flex';
                    document.getElementById('oligo_button').style.display = 'block'; // Proceed to the next page anyway
                } else {
                    alert(resp.result);
                    document.getElementById("container_title3").style.background = '#EA2845';
                    document.getElementById('container_consensus').style.display = 'flex';
                    document.getElementById("consesus_message").innerHTML = "<h3>Unexpected error. Please reload the page.</h3>";
                }
                percDeg = resp.perc;
            },
            error: function () { alert("An error occoured"); }
        });
    })
    \$("#oligo_design").click(function () {
        defSet1 = defSet + "T"; //I'll have to remove a
        var newUrl = "oligo_design_pp.cgi?defSet=" + defSet1 + "-" + percDeg;
        document.location = newUrl;
    });

})
</script>

</body>
END HTML

print $html;
