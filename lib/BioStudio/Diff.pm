package BioStudio::Diff;
require Exporter;

use BioStudio::Basic qw(flatten_subfeats);
use GeneDesign::Basic;
use Text::Diff;
use Time::Format qw(%time);
use Perl6::Slurp;
use Digest::MD5;
#use Bio::AlignIO;

use strict;
use warnings;

use vars qw($VERSION @ISA @EXPORT_OK %EXPORT_TAGS);
$VERSION = '1.00';

@ISA = qw(Exporter);
@EXPORT_OK = qw(
  compare_features
  compare_dbs
  compare_feature_sequences
  compare_feature_translations
  compare_comments
  %CHANGES
);
%EXPORT_TAGS = (all => [qw(compare_features compare_dbs %CHANGES
  compare_feature_sequences compare_feature_translations compare_comments)]);
     
our %CHANGES = (
    0  => [ -1, "deleted feature",           "LOST FEATURES" ],
    1  => [ 1,  "added feature",             "ADDED FEATURES" ],
    2  => [ -1, "lost subfeature",           "LOST SUBFEATURES" ],
    3  => [ 1,  "gained subfeature",         "GAINED SUBFEATURES" ],
    4  => [ -2, "lost sequence",             "LOST SEQUENCE" ],
    5  => [ 2,  "gained sequence",           "GAINED SEQUENCES" ],
    6  => [ 0,  "change in translation",     "CHANGES IN TRANSLATION" ],
    7  => [ 0,  "change in sequence",        "CHANGES IN SEQUENCE" ],
    8  => [ -2, "lost annotation",           "LOST ANNOTATIONS" ],
    9  => [ 2,  "gained annotation",         "GAINED ANNOTATIONS" ],
    10 => [ 0,  "changed annotation",        "CHANGES IN ANNOTATIONS" ],
    11 => [ 0,  "changed subfeature",        "CHANGES TO SUBFEATURES" ],
    -1 => [ 0,  "change in region sequence", "" ],
    -2 => [ 0,  "change in region length",   "" ],
    -3 => [ -1, "lost region",               "" ],
    -4 => [ 1,  "gained region",             "" ]
);

################################################################################
################################ Diff Functions ################################
################################################################################

sub compare_dbs
{
  my ($pa, $BS, $DB1, $DB2) = @_;
  my (%feats1, %feats2, @DIFFS) = ((), (), ());
  
  my $iterator1 = $DB1->get_seq_stream();
  while (my $feature = $iterator1->next_seq)
  {
    $feats1{$feature->Tag_load_id} = $feature;
  }
  
  my $iterator2 = $DB2->get_seq_stream();
  while (my $feature = $iterator2->next_seq)
  {
    $feats2{$feature->Tag_load_id} = $feature;
  }

  foreach my $featid (grep {exists ($feats2{$_}) } keys %feats1)
  {
    push @DIFFS, compare_features($pa, $BS, $feats1{$featid}, $feats2{$featid});
  }
  foreach my $featid (grep {! exists ($feats2{$_}) } keys %feats1)
  {
    my $feat1 = $feats1{$featid};
    next if ($feat1->has_tag("parent_id"));
    push @DIFFS, [0, $feat1, $feat1->primary_tag() . " " . $featid];
  }
  foreach my $featid (grep {! exists ($feats1{$_}) } keys %feats2)
  {
    my $feat2 = $feats2{$featid};
    next if ($feat2->has_tag("parent_id"));  
    push @DIFFS, [1, $feat2, $feat2->primary_tag() . " " . $featid];
  }
  return @DIFFS;
}

sub compare_features
{
  my ($pa, $BS, $feat1, $feat2) = @_;
  my $pfeat1 = $feat1->primary_tag() . " " . $feat1->Tag_load_id;
  my $pfeat2 = $feat2->primary_tag() . " " . $feat2->Tag_load_id;
  my $f1score = $feat1->score ? $feat1->score : 0;
  my $f2score = $feat2->score ? $feat2->score : 0;
  my $f1phase = $feat1->phase ? $feat1->phase : 0;
  my $f2phase = $feat2->phase ? $feat2->phase : 0;
  my @Changes = ();
  my $subseqchangeflag = 0;
  ##Check for existence subfeatures and compare
  if (scalar($feat1->get_SeqFeatures) || scalar($feat2->get_SeqFeatures))
  {
    my @allsubs1 = flatten_subfeats($feat1);
    my @allsubs2 = flatten_subfeats($feat2);
    my %types1 = map {$_->primary_tag() => 1} @allsubs1; 
    my %types2 = map {$_->primary_tag() => 1} @allsubs2;
    foreach my $type ( grep {exists $types2{$_}} keys %types1)
    {
      my %subs1 = map {$_->Tag_load_id => $_} grep {$_->primary_tag eq $type} @allsubs1;
      my %subs2 = map {$_->Tag_load_id => $_} grep {$_->primary_tag eq $type} @allsubs2;
      foreach my $sub1id (grep {exists $subs2{$_}} keys %subs1)
      {
        foreach my $subchange(compare_features($pa, $BS, $subs1{$sub1id}, $subs2{$sub1id}))
        {
          if ($$subchange[0] > 7)
          {     
            push @Changes, [11, $feat1, "$pfeat1 " . $CHANGES{11}->[1] . " $type $sub1id (" . $CHANGES{$$subchange[0]}->[1] . ")"];
          }   
          else
          {
            $subseqchangeflag++;
          }      
        }
        delete $subs1{$sub1id};
        delete $subs2{$sub1id};
      }
      push @Changes, [2, $feat1, "$pfeat1 " . $CHANGES{2}->[1] . " $type $_"] foreach (keys %subs1);
      push @Changes, [3, $feat2, "$pfeat2 " . $CHANGES{3}->[1] . " $type $_"] foreach (keys %subs2);
    }
    foreach my $type (grep { ! exists $types2{$_}} keys %types1)
    {
      foreach my $subfeat (grep {$_->primary_tag eq $type} flatten_subfeats($feat1))
      {
        push @Changes, [2, $feat1, "$pfeat1 " . $CHANGES{2}->[1] . " $type " . $subfeat->Tag_load_id];
      }
    }
    foreach my $type ( grep {! exists $types1{$_}} keys %types2)
    {
      foreach my $subfeat(grep {$_->primary_tag eq $type} flatten_subfeats($feat2))
      {
        push @Changes, [3, $feat2, "$pfeat2 " . $CHANGES{3}->[1] . " $type " . $subfeat->Tag_load_id];
      }
    }
  }
  ##Check to see if the feature length has changed between 2 and 1, and if it has gotten shorter or longer
  if (($feat1->end - $feat1->start) != ($feat2->end - $feat2->start))# && $feat1->primary_tag() ne "chromosome")	
  {
    my $len1 = $feat1->end - $feat1->start + 1;
    my $len2 = $feat2->end - $feat2->start + 1;
    my $bpdiff = abs($len1 - $len2);
    if ($len1 > $len2)
    {
      push @Changes, [4, $feat1, "$pfeat1 " . $CHANGES{4}->[1] . " (" . $bpdiff . "bp)", $bpdiff];
    }
    if ($len1 < $len2)
    {
      push @Changes, [5, $feat2, "$pfeat2 " . $CHANGES{5}->[1] . " (" . $bpdiff . "bp)", $bpdiff];
    }
  }
  ##Check to see if source, score, orientation, or phase have changed.  
  if ($feat1->strand() ne $feat2->strand())
  {
    my $old = $feat1->strand();
    my $new = $feat2->strand();
    push @Changes, [10, $feat1, "$pfeat1 " . $CHANGES{10}->[1] . " from strand = $old to strand = $new"];
  }
  if ($f1score ne $f2score)
  {
    my $old = $feat1->score();
    my $new = $feat2->score();
    push @Changes, [10, $feat1, "$pfeat1 " . $CHANGES{10}->[1] . " from score = $old to score = $new"];
  }
  if ($f1phase ne $f2phase)
  {
    my $old = $feat1->phase();
    my $new = $feat2->phase();
    push @Changes, [10, $feat1, "$pfeat1 " . $CHANGES{10}->[1] . " from phase = $old to phase = $new"];
  } 
  if ($feat1->source_tag() ne $feat2->source_tag())
  {
    my $old = $feat1->source_tag();
    my $new = $feat2->source_tag();
    push @Changes, [10, $feat1, "$pfeat1 " . $CHANGES{10}->[1] . " from source = $old to source = $new"];
  }  
  
  ##Check to see if annotations have changed
  my @tags1 = $feat1->get_all_tags();
  my @tags2 = $feat2->get_all_tags();
  foreach my $tag (grep {$feat2->has_tag($_)} @tags1)
  {
    my @vals1 = sort $feat1->get_tag_values($tag);
    my @vals2 = sort $feat2->get_tag_values($tag);
    while (scalar(@vals1) && scalar(@vals2))
    {
      my ($val1, $val2) = (shift @vals1, shift @vals2);
      next if ($val1 eq $val2);
    }
  }
  foreach my $losttag (grep {! $feat2->has_tag($_)} @tags1)
  {
    my @lostvals = $feat1->get_tag_values($losttag);
    push @Changes, [8, $feat1, "$pfeat1 " . $CHANGES{8}->[1] . " $losttag = @lostvals"];
  }
  foreach my $addtag (grep {! $feat1->has_tag($_)} @tags2)
  {
    my @addvals = $feat2->get_tag_values($addtag);
    push @Changes, [9, $feat2, "$pfeat2 " . $CHANGES{9}->[1] . " $addtag = @addvals"];
  }
  
  ##If feature is a CDS and the translation switch is on, check to see if the translation of the feature has changed
  if ($feat1->primary_tag() eq "CDS" && ($pa->{TRX} || $pa->{TRXLN}))
  {
  	my $orf1 = $feat1->seq->seq;
  	my $orf2 = $feat2->seq->seq;
  	my $old = translate($orf1, $f1phase + 1, $pa->{CODON_TABLE});
  	my $var = translate($orf2, $f2phase + 1, $pa->{CODON_TABLE});
    if ((!$old && length($orf1) >=3) || (!$var && length($orf2) >=3))
    {
      print "original $feat1 " if (!$old);
      print "variant $feat2 " if (!$var);
      print "has no translation, strangely!!\n";
    }
  	elsif ($old ne $var)
  	{ 
  	  print "$feat1\n$orf1\n$orf2\n\n";
      my $alignment = "";
      if ($pa->{TRXLN})
      {
        my $id = $pa->{ID} ? $pa->{ID} : Digest::MD5::md5_hex(time().{}.rand().$$);
        $alignment = "\n";
        $alignment .= compare_feature_translations($BS, $pa->{CODON_TABLE}, $feat1, $feat2, $id);
        $alignment = "$alignment\n";
      }
  		push @Changes, [6, $feat2, "$pfeat2 " . $CHANGES{6}->[1] . "$alignment"];
  	}
  }
  ##If feature is not a chromosome and the check sequence switch is on, check to see if the sequence of the feature has changed
  if (($pa->{SEQ} || $pa->{ALN}) && $subseqchangeflag == 0 ) #&& $feat1->primary_tag() ne "chromosome")
  {
    my $old = $feat1->seq->seq;
    my $var = $feat2->seq->seq;
    if ($old ne $var)
    {
      my $alignment = "";
      if ($pa->{ALN} && $feat1->primary_tag() ne "chromosome")
      {
        my $id = $pa->{ID} ? $pa->{ID} : Digest::MD5::md5_hex(time().{}.rand().$$);
        $alignment = "\n";
        $alignment .= compare_feature_sequences($BS, $feat1, $feat2, $id);
        $alignment = "$alignment\n";
      }
      push @Changes, [7, $feat2, "$pfeat2 " . $CHANGES{7}->[1] . "$alignment"];
    }
  }
  return @Changes;
}

sub compare_feature_translations
{
  my ($BS, $CODON_TABLE, $feat1, $feat2, $id) = @_;
  my $seq1 = $BS->{blast_directory} . "/pold$id";
  my $seq2 = $BS->{blast_directory} . "/pnew$id";
  my $f1phase = $feat1->phase ? $feat1->phase : 0;
  my $f2phase = $feat2->phase ? $feat2->phase : 0;
  my $outfile = $BS->{blast_directory} . "/pbl2seq.out";
  open (SEQ1, ">$seq1");
  print SEQ1 ">pold\n", translate($feat1->seq->seq, $f1phase + 1, $CODON_TABLE), "\n";
  close SEQ1;
  open (SEQ2, ">$seq2");
  print SEQ2 ">pnew\n", translate($feat2->seq->seq, $f2phase + 1, $CODON_TABLE). "\n";
  close SEQ2;  
  my @args = ("$BS->{legacyblast} bl2seq -p blastp -i $seq1 -j $seq2 -o $outfile");
  $SIG{CHLD} = 'DEFAULT';
  system(@args) == 0 or die "system @args failed: $!";
  my $alignment = slurp($outfile);
  return $alignment;
}

sub compare_feature_sequences
{
  my ($BS, $feat1, $feat2, $id) = @_;
  my $seq1 = $BS->{blast_directory} . "/nold$id";
  my $seq2 = $BS->{blast_directory} . "/nnew$id";
  my $outfile = $BS->{blast_directory} . "/bl2seq.out";
  open (SEQ1, ">$seq1");
  print SEQ1 ">old\n", $feat1->seq->seq;
  close SEQ1;
  open (SEQ2, ">$seq2");
  print SEQ2 ">new\n", $feat2->seq->seq;
  close SEQ2;
  my @args = ("$BS->{legacyblast} bl2seq -p blastn -i $seq1 -j $seq2 -o $outfile -G 11 -E 2");
  $SIG{CHLD} = 'DEFAULT';
  system(@args) == 0 or die "system @args failed: $!";
#  my $factory = Bio::Tools::Run::StandAloneBlastPlus->new();
#  $factory->bl2seq(-method => 'blastn', -query=> $feat2->seq, -subject => $feat1->seq, -outfile => $outfile);
#  my $alnstream = Bio::AlignIO->new(-file=>$outfile, -format=> 'bl2seq') or die "$!";
#  my $aln = $alnstream->next_aln();
  my $alignment = slurp($outfile);
  return $alignment;
}

sub compare_comments
{
	my ($ref1, $ref2) = @_;
	my $diff = diff $ref1, $ref2;
	return $diff;
}

1;

__END__

=head1 NAME

BioStudio::Basic

=head1 VERSION

Version 1.00

=head1 DESCRIPTION

Basic BioStudio functions

=head1 FUNCTIONS

=head2 compare_dbs()

=head2 compare_features()

=head2 compare_feature_translations()

=head2 compare_feature_sequences()

=head2 compare_comments()

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
