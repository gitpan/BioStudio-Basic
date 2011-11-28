package BioStudio::BLAST;
require Exporter;

use BioStudio::Basic qw($VERNAME);
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::SeqFeature::Store;

use strict;
use warnings;

use vars qw($VERSION @ISA @EXPORT_OK %EXPORT_TAGS);
$VERSION = '1.00';

@ISA = qw(Exporter);
@EXPORT_OK = qw(
  make_BLAST_db
  make_megaBLAST_index
);
%EXPORT_TAGS = (all => [qw(make_BLAST_db make_megaBLAST_index)]);

################################################################################
############################### BLAST functions ################################
################################################################################

sub _make_FASTA
{
  my ($gff_hsh, $BS, $label) = @_;
  my $FASTAfile = $BS->{blast_directory} . "/" . $label . ".fa";
  print "seeking $FASTAfile\n";
  unless (-e $FASTAfile)
  {
    my $out = Bio::SeqIO->new(-file => ">$FASTAfile", -format => 'Fasta')
      || die "can't make fasta file $FASTAfile for output ($!)";
    foreach my $chr (keys %$gff_hsh)
    {
      print "Loading $chr...\n";
      my $db = Bio::DB::SeqFeature::Store->new( 
          -adaptor => 'memory',
          -dir     => $gff_hsh->{$chr} );
      my $chrname = $2 if ($chr =~ $VERNAME);
      my $seqid = "chr$chrname";
      my $bases = $db->fetch_sequence($seqid);
      unless ($bases)
      {
        print "\t CAN'T FIND BASES for $seqid!!!\n";
        next;
      }
      my $seq = Bio::Seq->new( -id => $seqid, -seq => $bases);
      $out->write_seq($seq);
    }
  }
  return;
}

sub make_BLAST_db
{
  my ($gff_hsh, $BS, $label) = @_;
  my $BLASTdb = $BS->{blast_directory} . "/" . $label . ".nsq";
  my $FASTAfile = $BS->{blast_directory} . "/" . $label. ".fa";
  if (! -e $BLASTdb)
  {
    if (! -e $FASTAfile)
    {
      _make_FASTA($gff_hsh, $BS, $label);
      system("wait");
    }
    my @args = ("$BS->{makeblastdb} -input_type fasta -in $FASTAfile -dbtype nucl -parse_seqids -out $BS->{blast_directory}/$label");
    $SIG{CHLD} = 'DEFAULT';
    system(@args) == 0 or die "system @args failed: $!";
  }
  return $BS->{blast_directory} . "/" . $label;
}

sub make_megaBLAST_index
{
  my ($BLASTdb, $BS, $label) = @_;
  my $megaBLASTchk = $BS->{blast_directory} . "/mb" . $label . ".00.idx";
  my $megaBLASTidx = $BS->{blast_directory} . "/mb" . $label;
  if (! -e $megaBLASTchk)
  {
    my @args = ("$BS->{makembindex}", "-input", $BLASTdb);
    push @args, "-output", $megaBLASTidx, "-iformat", "blastdb";
    $SIG{CHLD} = 'DEFAULT';
    system(@args) == 0 or die "system @args failed: $!";
  }
  return $megaBLASTidx;
}

1;

__END__

=head1 NAME

BioStudio::BLAST

=head1 VERSION

Version 1.00

=head1 DESCRIPTION

BioStudio functions for BLAST functionality

=head1 FUNCTIONS

=head2 _make_FASTA()
  given a hashref of chromosomes (where the key is the name and the value is the 
  path, see BioStudio::Basic::gather_versions()), the BioStudio config hashref,
  and a label, creates a FASTA file that contains all of their sequences for 
  BLAST database creation

=head2 make_BLAST_db()
  given a hashref of chromosomes (where the key is the name and the value is the 
  path, see BioStudio::Basic::gather_versions()), the BioStudio config hashref,
  and a label, creates a BLAST database containing all of the chromosome
  sequences in the hashref

=head2 make_megaBLAST_index()
  given the name of a BLAST database, a BioStudio config hashref, and a label,
  creates a megaBLAST index to speed BLASTing

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