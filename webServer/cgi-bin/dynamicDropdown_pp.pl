#!/usr/bin/perl -w

use DBI;
use strict;
use CGI qw(:standard);
use CGI::Carp qw(fatalsToBrowser);
use JSON;

#connected to phyloprimerResultsPrimer.cgi, phyloprimerResultsProbe.cgi, phyloprimerResultsOligo.cgi, phyloprimerResultsCheckPrimer.cgi, phyloprimerResultsCheckProbe.cgi and phyloprimerResultsCheckOligo.cgi : render taxa for result page dropdown

#set mySQL parameters
my $dsn = "mysql_database";
my $user_name = "mysql_user";
my $password = "mysql_password";

my $dbh;
my $sth;

my (@ary);
my (@tot);
my ($tot_j);
my ($ary_s);

my $cgi = CGI->new;
print $cgi->header(-type => "application/json", -charset => "utf-8");

#convert  data to JSON
my $op = JSON -> new -> utf8 -> pretty(1);


my $company = $cgi->param('text');

#connect to database
$dbh = DBI->connect ($dsn, $user_name, $password, { RaiseError => 1 });

#issue query
$sth = $dbh->prepare("SELECT taxa FROM taxa_pp WHERE taxa LIKE '$company%' LIMIT 20");

#execute the prepared statement handle:
$sth->execute ();
#read results of a query, then clean up

my $tot;
my $inside = 0;
while (@ary = $sth->fetchrow_array ()) {
    $inside=1;
    $tot .= $ary[0] . "!";
}

$sth->finish ();

#result array
my @result;

if ($inside == 0) {
    $tot = 'empty';
}

push @result, "result";
push @result, $tot;

my $json = $op -> encode({
    @result
});
print $json;

