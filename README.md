# PhyloPrimer

PhyloPrimer was born from a collaboration between the University of Bristol and ENOVEO (https://enoveo.com), a company involved in bioremediation and microbial monitoring projects that looks into a wide variety of environmental samples with the aim of identifying specific genes and microorganisms. The idea was to create a user- friendly platform to i) easily retrieve DNA sequences similar to a gene of interest, ii) design oligos (e.g. primer and probes) suitable for different PCR applications, iii) check the oligos for the presence of secondary structures and iv) study each oligo with *in silico* tests for determining which organisms they can potentially target. The strength of PhyloPrimer is the possibility of retrieving the DNA sequences, used for the oligo design, through a dynamic tree. The tool was designed with the study of microbial communities in mind, so we do not encourage the use of PhyloPrimer for the study of eukaryotic organisms: the dynamic tree is built only with microorganism sequences and the tool does not take in consideration the presence of intron regions. Also, PhyloPrimer uses a consensus approach and will introduce degenerate primers only if necessary and if specified by the user, but we do not suggest the use of it if the aim is specifically to design degenerate primers. PhyloPrimer database is a modified version of the GenBank database (release 239) and it was last updated in October 2020.

Link to PhyloPrimer homepage: https://www.cerealsdb.uk.net/cerealgenomics/phyloprimer/

Look at the manual to have move information.

#### PhyloPrimer scripts

PhyloPrimer runs on a remote server provided from the University of Bristol. The current server has 48 CPUs (64-bit Intel(R) Xeon(R) CPU E5-2680 v3 at 2.50GHz). Only 4 PhyloPrimer processes at one time are allowed on the server, the excess processes enter a queue. On average, the oligo design requires 40/50 minutes whereas the oligo check requires 5/10 minutes. The web interface was implemented in html and JavaScript. PhyloPrimer is coded in Perl, JavaScript, HTML, CSS and MySQL. Two JavaScript packages were used: a modified version of PhyloCanvas v 1.7.3 (phylocanvas.org) and CanvasJS v 2.3.2 (canvasjs.com).

All the scripts can be found in the folder [webServer](webServer).

#### PhyloPrimer databases

Here are the links for the two PhyloPrimer sequence databases:

- DB1 : https://www.cerealsdb.uk.net/database_pp/DB1_pp/DB1.fasta.gz
- DB2 : https://www.cerealsdb.uk.net/database_pp/DB2_pp/DB2.fasta.gz

And the links to the PhyloPrimer taxonomy databases:

- DB1 accession-taxid correspondances : https://www.cerealsdb.uk.net/database_pp/taxonomy/DB1_accession_taxid.txt.gz
- DB2 accession-taxid correspondances : https://www.cerealsdb.uk.net/database_pp/taxonomy/DB2_accession_taxid.txt.gz
- taxid-taxonomy correspondances : https://www.cerealsdb.uk.net/database_pp/taxonomy/taxid_taxonomy.txt.gz
