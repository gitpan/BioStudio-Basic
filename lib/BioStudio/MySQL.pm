package BioStudio::MySQL;
require Exporter;

use BioStudio::Basic qw($VERNAME);
use DBI;
use Carp;
use Bio::DB::SeqFeature::Store;

use strict;
use warnings;

use vars qw($VERSION @ISA @EXPORT_OK %EXPORT_TAGS);
$VERSION = '1.00';

@ISA = qw(Exporter);
@EXPORT_OK = qw(  
  drop_database
  list_databases
  create_database
  load_database
  fetch_database
);
%EXPORT_TAGS = (all => [qw(drop_database list_databases create_database 
  load_database fetch_database)]);

################################################################################
##################### MySQL Database Interaction Functions #####################
################################################################################

sub list_databases
{
	my ($BS) = @_;
	my $dbh = DBI->connect('dbi:mysql:mysql', $BS->{mysql_user}, $BS->{mysql_pass}, 
	  { RaiseError => 1, AutoCommit => 1});
	my %dblist;
	my $sth = $dbh->prepare(q{SHOW DATABASES}) 
	  or croak "Unable to prepare show databases: ". $dbh->errstr."\n";
	$sth->execute or croak "Unable to exec show databases: ". $dbh->errstr."\n";
	my $aref;
	while ($aref = $sth->fetchrow_arrayref) 
	{	
		$dblist{$aref->[0]}++;
	} 
	$sth->finish;
	$dbh->disconnect();
	return \%dblist;
}

sub create_database
{
	my ($chrname, $BS) = @_;
	my $dbh = DBI->connect('dbi:mysql:mysql', $BS->{mysql_user}, $BS->{mysql_pass}, 
	  { RaiseError => 1, AutoCommit => 1});
	$dbh->do("create database $chrname;");
	$dbh->do("grant select on $chrname.* to nobody\@localhost;");
	$dbh->do("flush privileges;");
	$dbh->disconnect();
	return;
}

sub load_database
{
	my ($chrname, $BS, $altchrname) = @_;
  my $gffget = $altchrname  ? $altchrname : $chrname;
	my ($species, $chromname) = ($1, $2) if ($gffget =~ $VERNAME);
	my $fileloc = "$BS->{genome_repository}/$species/chr$chromname/$gffget.gff";
	my @args = ($BS->{bioperl_bin} . "/bp_seqfeature_load.pl", "--noverbose");
	push @args, "-c", "-d", $chrname, $fileloc, "-f";
	push @args, "--user", $BS->{mysql_user}, "-p", "$BS->{mysql_pass}";
	$SIG{CHLD} = 'DEFAULT';
  system(@args) == 0 or croak "system @args failed: $!";
	return;
}

sub fetch_database
{
  my ($chrname, $BS, $write) = @_;
  my $writeflag = $write  ? 1 : 0;
  my $db = Bio::DB::SeqFeature::Store->new(
        -adaptor => "DBI::mysql",
        -dsn => "dbi:mysql:$chrname",
        -user => $BS->{mysql_user},
        -pass => $BS->{mysql_pass},
        -write => $writeflag
  );
  return $db;
}

sub drop_database
{
	my ($chrname, $BS) = @_;
	my $dbh = DBI->connect('dbi:mysql:mysql', $BS->{mysql_user}, $BS->{mysql_pass}, 
	  { RaiseError => 1, AutoCommit => 1});
	$dbh->do("drop database $chrname;");
	$dbh->do("flush privileges;");
	$dbh->disconnect();
	return;
}


1;
__END__

=head1 NAME

BioStudio::MySQL

=head1 VERSION

Version 1.00

=head1 DESCRIPTION

BioStudio functions for MySQL interaction.

=head1 FUNCTIONS

=head2 list_databases()
  Given the BioStudio config hashref, return a hashref where the keys are all
  chromosomes available in the MySQL database.

=head2 create_database()
  Given a chromosome name and the BioStudio config hashref, create a mysql 
  database to store the chromosome in MySQL. This function does not load the
  database!

=head2 load_database()
  Given a chromosome name, the BioStudio config hashref, and optionally, the
  name of an alternate chromosome, this function loads a MySQL database (which
  must have been previously created, see create_database()) with a GFF file. The
  file is the one corresponding to the first chromosome provided unless the
  alternate is defined, in which case the alternate is loaded into a database
  named after the first chromosome.

=head2 fetch_database()
  Given a chromosome name and the BioStudio config hashref, returns a 
  Bio::DB::SeqFeature::Store interface for the MySQL database containing the 
  chromosome. An optional write flag sets whether or not the interface will 
  support adding, deleting, or modifying features

=head2 drop_database()
  Given a chromosome name and the BioStudio config hashref, deletes the MySQL
  database containing the chromosome.

=head1 AUTHOR

Sarah Richardson <notadoctor@jhu.edu>.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2011, BioStudio developers
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the Johns Hopkins nor the
      names of the developers may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE DEVELOPERS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=cut
