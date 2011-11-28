package BioStudio::Basic;
require Exporter;

use GeneDesign::Basic;
use File::Find;
use Text::Wrap qw($columns &wrap);
use Perl6::Slurp;
use Config::Auto;
use Bio::DB::SeqFeature::Store;

use strict;
use warnings;

use vars qw($VERSION @ISA @EXPORT_OK %EXPORT_TAGS);
$VERSION = '1.00';

@ISA = qw(Exporter);
@EXPORT_OK = qw(
  configure_BioStudio
  fetch_custom_features
  fetch_custom_markers
  fetch_enzyme_lists
  make_mask
  mask_combine
  mask_filter
  get_src_path
  get_genome_list
  gather_versions
  ORF_compile
  get_feature_sequence
  check_new_sequence
  flatten_subfeats
  gene_names
  allowable_codon_changes
  print_as_fasta
  rollback
  $VERNAME
  %SPECIES
);
%EXPORT_TAGS = (all => [qw(configure_BioStudio fetch_custom_features 
  fetch_custom_markers fetch_enzyme_lists make_mask mask_combine mask_filter 
  get_src_path get_genome_list gather_versions ORF_compile get_feature_sequence 
  check_new_sequence flatten_subfeats gene_names allowable_codon_changes 
  print_as_fasta rollback $VERNAME %SPECIES)]);

our $VERNAME		=	qr/([\w]+)_chr([\w\d]+)_(\d+)_(\d+)/;

our %SPECIES		= ("yeast" => 1, 
                   "ecoli" => 2, 
                   "bsubtilis" => 6, 
                   "dradiodurans" => 7, 
                   "mgenitalium" => 8);         

################################################################################
########################### Customization  Functions ###########################
################################################################################
sub configure_BioStudio
{
	my ($conf_dir) = @_;
	my $BS = Config::Auto::parse("BioStudio.conf", ("path"=>$conf_dir));
  my @featkeys = grep {$_ =~ /\w+\.colors/} keys %$BS;
  $BS->{COLORS} = {};
  foreach my $prefeat (@featkeys)
  {
    my @temp = split(/\./, $prefeat);
    my $feat = $temp[0];
    $BS->{COLORS}->{$feat} = {};
    foreach my $entry (@{$BS->{$prefeat}})
    {
      my @entries = split(":", $entry);
      $BS->{COLORS}->{$feat}->{$entries[0]} = $entries[-1];
    }
  }
  my @featkey = grep {$_ =~ /\w+\.color\b/} keys %$BS;
  $BS->{COLOR} = {};
  foreach my $prefeat (@featkey)
  {
    my @temp = split(/\./, $prefeat);
    my $feat = $temp[0];
    $BS->{COLOR}->{$feat} = $BS->{$prefeat};
  }    
  return $BS;
}

sub fetch_custom_features
{
  my ($BS) = @_;
  my %BS_FEATURES;
  my @feats = slurp($BS->{feature_file});
  foreach my $featline (@feats)
  {
    my @featatts = split(/[\t\n]+/, $featline);
    my $feature = {};
    $feature->{NAME} = $featatts[0];
    $feature->{KIND} = $featatts[1];
    $feature->{SOURCE} = $featatts[2];
    $feature->{SEQ} = uc $featatts[3];
    $BS_FEATURES{$featatts[0]} = $feature;
  }
  return \%BS_FEATURES;
}

sub fetch_custom_markers
{
  my ($BS) = @_;
  opendir(DIR, $BS->{marker_dir});
  my @markers = readdir(DIR);
  closedir DIR;
  @markers = grep {$_ =~ /\.gff\Z/} @markers;
  my %markers;
  foreach my $marker (@markers)
  {
    my $name = $1 if ($marker =~/([\w\d]+)\.gff\Z/);
    my $markerdb = Bio::DB::SeqFeature::Store->new(
      -adaptor => 'memory',
      -gff     => "$BS->{marker_dir}/$marker",
      -index_subfeatures => "true");
    my @alls = $markerdb->get_features_by_type("region");
    my $all = $alls[0];
    my @CDSes = $markerdb->get_features_by_type("CDS");
    my $seqid = $CDSes[0]->seq_id;
    my $region = $markerdb->fetch_sequence($seqid);
    my $newmarker = {};
    $newmarker->{NAME} = $name;
    $newmarker->{SEQ} = $region;
    $newmarker->{DB} = $markerdb;
    if ($all->has_tag("color"))
    {
      $newmarker->{COLOR} = "#" . join("", $all->get_tag_values("color"));
    }
    $markers{$name} = $newmarker;
  }
  return \%markers;
}

sub fetch_enzyme_lists
{
  my ($BS) = @_;
  opendir(DIR, $BS->{enzyme_dir});
  my @lists = readdir(DIR);
  closedir DIR;
  @lists = grep {$_ =~ /\w+/} @lists;
  return @lists;
}


################################################################################
############################## Masking  Functions ##############################
################################################################################
sub make_mask
{
	my ($length, $featlist, $offset) = @_;
	my @MASK = (0) x $length;
	$offset = 0 if (! $offset);
	foreach my $q (@{$featlist})
	{
		for my $x ($q->start - $offset - 1 .. $q->end - $offset - 1)
		{
            next if ($x < 0);
			$MASK[$x] = $MASK[$x]+1;	#could always mask essential genes by making them negative numbers
			last if ($x == $length - 1);
		}
	}
	return join("", @MASK);
}

sub mask_combine
{
	my ($mask1, $mask2) = @_;
	my @MASK = split("", $mask1);
	my @MASK2 = split("", $mask2);
	for my $q (0..scalar(@MASK)-1)
	{
		$MASK[$q] += $MASK2[$q];
	}
	return join("", @MASK);
}

sub mask_filter
{
	my ($mask) = @_;
	my $coords = substr($mask, 0, 1) eq "0"	?	[0]	:	[];
	while ($mask =~ /(?=[^0]0|0[^0])/ig) {	push @$coords, (pos $mask)+1;	}
	push @$coords, length($mask) if (substr($mask, -1, 1) eq "0");
	return $coords;
}


################################################################################
######################### Genome Repository  Functions #########################
################################################################################
sub get_src_path
{
	my ($chrname, $BS) = @_;
	my ($species, $chromname) = ($1, $2) if ($chrname =~ $VERNAME);
	my $fileloc = "$BS->{genome_repository}/$species/chr$chromname/$chrname.gff";
  return $fileloc;
}

sub get_genome_list
{
  my ($BS) = @_;
  my @list;
  my @genomes;
	find sub { push @list, $File::Find::name}, $BS->{genome_repository};
	@list = grep {	$_ =~ /\.gff\Z/} @list;
  foreach (@list)
  {
    push @genomes, $1 if ($_ =~ /($VERNAME)/);
  }
	return @genomes; 
}

sub gather_versions
{
	my ($species, $target, $bsconfig) = @_;
	my $verhsh;
	my @gffs;
	find sub {  push @gffs, $File::Find::name, -d && "/" }, $$bsconfig{genome_repository};
	if ($target == 0)
	{
		foreach my $name (@gffs)
		{
			$$verhsh{$1} = $name if ($name =~ /($species\_chr\d+\_0\_00)\.gff/);
		}
	}
	elsif ($target == -1)
	{
		foreach my $name (@gffs)
		{
			$$verhsh{$1} = $name if ($name =~ /($species\_chr\d+\_\d+\_\d+)\.gff/);
		}
		my @oldnames;
		my %chrseen;
		foreach my $name (keys %$verhsh)
		{
			my ($chr, $genv, $chrv) = ($1, $2, $3) if ($name =~ $VERNAME);
			if (exists $chrseen{$chr})
			{
				my ($agenv, $achrv) = ($2, $3) if ($chrseen{$chr} =~ $VERNAME);
				if ($achrv < $chrv || $agenv < $genv)
				{
					push @oldnames, $chrseen{$chr};
					$chrseen{$chr} = $name;
				}
				else
				{
					push @oldnames, $name;
				}
			}
			else
			{
				$chrseen{$chr} = $name;
			}
		}
		delete @{$verhsh}{@oldnames};
	}
	else
	{
		foreach my $name (@gffs)
		{
			$$verhsh{$1} = $name if ($name =~ /($species\_chr\d+\_\d+\_$target)\.gff/);
		}	
	}
	return $verhsh;
}

sub rollback
{
  my ($chrname, $BS) = @_;
	print "Deleting gff file...";
	my $path = get_src_path($chrname, $BS);
	system "rm $path" == 0 || print "can't delete gff file!";
	print "Dropping mysql database...";
	drop_database($chrname, $BS);
  return;
}


################################################################################
######################### Editing and Markup Functions #########################
################################################################################
sub ORF_compile
{
	my ($genearr) = @_;
	my $final = {};
	foreach my $gene (@$genearr)
	{
    my $genename = $gene->Tag_load_id;
    my @subs = flatten_subfeats($gene);
    my @CDSes = grep {$_->primary_tag eq "CDS"} @subs;
    @CDSes = sort {$b->start <=> $a->start} @CDSes if ($gene->strand == -1);
    @CDSes = sort {$a->start <=> $b->start} @CDSes if ($gene->strand == 1);
    $$final{$genename} .= $_->seq->seq foreach (@CDSes);
	}	
	return $final;
}

sub get_feature_sequence
{
	my ($feat, $seq) = @_;
	return substr($seq, $feat->start-1, $feat->end - $feat->start + 1);
}

sub flatten_subfeats
{
  my ($feature) = @_;
  my @subs = $feature->get_SeqFeatures();
  push @subs, $_->get_SeqFeatures foreach (@subs);
  return @subs;
}

sub gene_names
{
  my ($genelist, $BS) = @_;
  my $DISPLAY = {};
  foreach my $gene (@$genelist)
  {    
    my $essstatus = $gene->Tag_essential_status;
    my $orfstatus = $gene->Tag_orf_classification;
    my $displayname = $gene->has_tag("gene")  
                  ?  $gene->Tag_load_id . " (" . $gene->Tag_gene . ")"  
                  :  $gene->Tag_load_id;
    $DISPLAY->{$gene->Tag_load_id} = $displayname;
  }  
  return $DISPLAY;
}

sub allowable_codon_changes
{
	my ($cod1, $cod2, $CODON_TABLE) = @_;
	my %result;
	for my $orient (0..1)
	{
		$result{$orient} = {};
		for my $offset (0..2)
		{
			my %union = my %isect = ();
			my $seq1 = $offset != 0	
				?	"N" x $offset . $cod1 . "N" x (3 - $offset)
				:	$cod1;
			my $seq2 = $offset != 0
				?	"N" x $offset . $cod2 . "N" x (3 - $offset)
				:	$cod2;
			my $qes1 = $orient	
				?	complement($seq1, $orient)
				:	$seq1;
			my $qes2 = $orient	
				?	complement($seq2, $orient)
				:	$seq2;
			my @set1 = amb_translation($qes1, $CODON_TABLE, 1);
			my @set2 = amb_translation($qes2, $CODON_TABLE, 1);
			foreach my $pep (@set1, @set2) 
			{ 
				$union{$pep}++ && $isect{$pep}++;
			}
			$result{$orient}->{$_}++ foreach keys %isect;
		#	print "$offset and $orient:\n";
		#	print "\t$seq1 - $qes1: @set1\n\t$seq2 - $qes2: @set2\n";
		}
	}
	return \%result;
}

sub check_new_sequence
{
  my ($feat) = @_;
  return 1 unless ($feat->has_tag("newseq"));
  my $newseq = $feat->Tag_newseq;
  my $featseq = $feat->seq->seq;
  return 1 if ($newseq eq $featseq || $newseq eq complement($featseq, 1));
  return 0;  
}

################################################################################
########################### Miscellaneous  Functions ###########################
################################################################################
sub print_as_fasta
{
  my ($sequence, $seqid) = @_;
  my $FASTArr = ["##FASTA\n"];
  $columns = 81;
  push @$FASTArr, ">" . $seqid . "\n";
  push @$FASTArr, wrap("","", $sequence), "\n";
  push @$FASTArr, "\n";
  return $FASTArr;
}

1;

__END__

=head1 NAME

BioStudio::Basic - basic functions for the BioStudio synthetic biology framework
 
=head1 VERSION

Version 1.00

=head1 DESCRIPTION

Basic BioStudio functions

=head1 FUNCTIONS

=head2 configure_BioStudio()
  This function loads the configuration file into a hash ref. You must pass it
  the path to the directory containing the configuration file; it will use 
  Config::Auto to ``magically'' parse the file.  

=head2 fetch_custom_features()
  Pass the config hashref, receive a hashref of the custom features defined in
  the BioStudio configuration directory. Each feature has four attributes: NAME,
  KIND, SOURCE, and SEQ

=head2 fetch_custom_markers()
  Pass the config hashref, receive a hashref of the custom markers defined in
  the BioStudio configuration directory. Each marker is a GFF file that gets
  read into the attributes NAME, SEQ, DB (a Bio::DB::SeqFeature::Store), and 
  COLOR (if a color is defined in the GFF file).

=head2 fetch_enzyme_lists()
  Pass the config hashref, receive an array that contains the names of the 
  enzyme lists in the BioStudio configuration directory. Each list is a 
  GeneDesign compatible list of restriction enzyme recognition sites.

=head2 make_mask()
  Given a length, a reference to a list full of Bio::SeqFeatures, and optionally
  an offset, returns a string of integers where each positon corresponds to a 
  base of sequence, and the integer represents the number of features that
  overlap that base. Obviously limited to ten overlapping features before a
  serious bug sets in :(

=head2 mask_combine()
  Takes two string masks (see make_mask()) and adds them. Returns the merged
  mask.

=head2 mask_filter()
  Takes a string mask (see make_mask()) and returns a listref of break 
  coordinates; that is, where does feature sequence end and interfeature 
  sequence begin, and where does interfeature sequence end and feature sequence
  begin? For example, if the mask is "0001100033221100", the resulting list 
  would be [0 3 5 8 14 16], meaning that features exist from 4 to 5 and 9 to 14.
  Intergenic sequence coordinates can thus be pulled out by hashing the array,
  %inter = @{mask_filter($mask)} where each key + 1 is the left coordinate, 
  and the value is the right coordinate.

=head2 get_src_path()
  Given a chromosome name and the config hashref, returns the absolute path to 
  that chromosome in the BioStudio genome repository.

=head2 get_genome_list()
  Given the config hashref, returns a list of all chromosomes in the BioStudio
  genome repository.

=head2 gather_versions()
  Given a species, a target, and the config hashref, returns a hashref of all 
  chromosomes in the species in the BioStudio genome repository that match the 
  target.  The target is an integer that represents a version.  
  If target is set to 0, we will return every wildtype version.
  If target is set to -1, we will return every latest version.
  For any other target (1, 3, 5) we will return that particular version.

=head2 rollback()
  Given a chromosome name and the BioStudio config hashref, removes that
  chromosome from the BioStudio genome repository.

=head2 ORF_compile()
  given a reference to an array full of Bio::SeqFeature gene objects, returns a 
  reference to a hash with gene ids as keys and concatenated 5' to 3' coding 
  sequences as values

=head2 get_feature_sequence()
  For when you can't use the Bio::SeqFeature seq function. Given a 
  Bio::SeqFeature compliant feature and a sequence, returns the sequence that
  the coordinates of the feature indicate.

=head2 check_new_sequence()
  Best when used as a confirmation that your edits went as expected. Given a 
  Bio::SeqFeature compliant feature that has a``newseq'' attribute, checks if
  the newseq and the actual sequence occupied by the feature are the same

=head2 flatten_subfeats()
  Given a seqfeature, iterate through its subfeatures and add all their subs to
  one big array. Mainly need this when CDSes are hidden behind mRNAs in genes.

=head2 gene_names()
  Given a list of Bio::SeqFeature gene objects and the BioStudio config hashref,
  returns a hash where each gene id is the key to a display friendly string.

=head2 allowable_codon_changes()
  Given two codons (a from, and a to) and a GeneDesign codon table hashref, this
  function generates every possible peptide pair that could contain the from 
  codon and checks to see if the peptide sequence can be maintained when the
  from codon is replaced by the to codon.  This function is of particular use 
  when codons are being changed in genes that overlap one another.

=head2 print_as_fasta()
  takes a sequence as a string and a sequence id and returns an 80 column FASTA
  formatted sequence block as an array reference

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