package BioStudio::GFF3;
require Exporter;

use BioStudio::Basic qw($VERNAME &print_as_fasta);
use URI::Escape;
use Perl6::Slurp;
use Carp;

use strict;
use warnings;

use vars qw($VERSION @ISA @EXPORT_OK %EXPORT_TAGS);
$VERSION = '1.00';

@ISA = qw(Exporter);
@EXPORT_OK = qw(
  get_GFF_comments 
  make_GFF3
  gff3_string
  
);
%EXPORT_TAGS = (all => [qw(get_GFF_comments make_GFF3 gff3_string)]);
  
my %PHASES = (0 => 1, 1 => 1, 2 => 1);
my $COMMENT			=	qr/^\#/;

################################################################################
################################ GFF3 Functions ################################
################################################################################

sub get_GFF_comments
{
	my ($chrname, $BS) = @_;
	my ($species, $chromname) = ($1, $2) if ($chrname =~ $VERNAME);
	my $fileloc = "$BS->{genome_repository}/$species/chr$chromname/$chrname.gff";
	carp "can't find $chrname for parsing GFF comments" unless( -e $fileloc);
  my @comments = slurp $fileloc;
  @comments = grep {$_ =~ $COMMENT} @comments;
  @comments = grep {$_ !~ /^\#\#/} @comments;
	return \@comments;
}

sub gff3_string
{
  my ($feat) = @_;
  my ($seqid, $source, $type, $start, $end, $score, $strand, $phase) = 
  ($feat->seq_id(), $feat->source_tag(), $feat->primary_tag(), $feat->start(), 
   $feat->end(), $feat->score(), $feat->strand(), $feat->phase());
  $score = $score ? $score  : ".";
  $phase = $phase && exists($PHASES{$phase}) ? $phase  : ".";
  $strand = $strand == "-1" ? "-" : $strand == "1"  ? "+" : ".";
  my $string = "$seqid\t$source\t$type\t$start\t$end\t$score\t$strand\t$phase\t";
  if ($feat->has_tag("load_id"))
  {
    $string .= "ID=" . $feat->Tag_load_id . ";";
    $string .= "Name=" . $feat->Tag_load_id . ";";
  }
  if ($feat->has_tag("parent_id"))
  {
    $string .= "Parent=" . $feat->Tag_parent_id . ";";
  }
  foreach my $tag (sort{$a cmp $b} $feat->get_all_tags())
  {
    next if ($tag eq "load_id" || $tag eq "parent_id");
    my @vals = $feat->get_tag_values($tag);
    if (scalar(@vals))
    {
      my $attstr = join(",", @vals);
      if ($tag eq "Note")
      {
        $string .= "$tag=" . uri_escape($attstr, '^\s^\d^\w^\-^\.^\_^\,^\(^\)') . ";";
      }
      else
      {
        $string .= "$tag=" . $attstr . ";";
      }
    }
  }
  return $string;
}

sub make_GFF3
{
  my ($pa, $BS, $newcomments) = @_;
  my @NEWGFF;
  
  #comments
  my $gffheader = shift @{$pa->{COMMENTS}};
  push @NEWGFF, $gffheader, @$newcomments, @{$pa->{COMMENTS}};

  #features
  my $seq_stream = $pa->{DB}->get_seq_stream()  or die "failed to get_seq_stream()";
  my @featarr;
  while (my $seq = $seq_stream->next_seq) 
  {
    push @featarr, $seq;
  }
  @featarr = sort {$a->start <=> $b->start || (($b->end - $b->start) <=> ($a->end - $a->start))} @featarr;
  push @NEWGFF, gff3_string($_) . "\n" foreach (@featarr);

  #sequence
  push @NEWGFF, @{print_as_fasta($pa->{EDITCHR}, $pa->{SEQID})};
  
	local $| = 1;
  open (OUT, ">$pa->{NEWFILE}") || die "can't open test output, $pa->{NEWFILE} $!\n";
  print OUT @NEWGFF;
	close OUT; 
  return;
}

1;

__END__

=head1 NAME

BioStudio::GFF3

=head1 VERSION

Version 1.00

=head1 DESCRIPTION

BioStudio functions for parsing and creating GFF3 files

=head1 Functions

=head2 get_GFF_comments()
  Given a chromosome name and the BioStudio config hashref, returns an arrayref
  containing all of the comments and directives from a GFF3 file.

=head2 gff3_string()

  Stupid Bioperl stupid Bio::DB::SeqFeature stupid no gff3_out

=head2 make_GFF3

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
