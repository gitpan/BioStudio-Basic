package BioStudio::GBrowse;
require Exporter;

use BioStudio::Basic qw($VERNAME &configure_BioStudio);
use File::Find;
use Time::Format qw(%time);
use Config::Auto;
use Perl6::Slurp;

use strict;
use warnings;

use vars qw($VERSION @ISA @EXPORT_OK %EXPORT_TAGS);
$VERSION = '1.00';

@ISA = qw(Exporter);
@EXPORT_OK = qw(
  gene_color
  write_conf_file
  rollback_gbrowse
  update_gbrowse
  get_gbrowse_src_list
  make_link
  gbrowse_gene_names
);
%EXPORT_TAGS = (all => [qw(gene_color write_conf_file rollback_gbrowse 
  update_gbrowse get_gbrowse_src_list make_link gbrowse_gene_names)]);
  
################################################################################
############################## GBrowse  Functions ##############################
################################################################################
sub gbrowse_gene_names
{
  my ($genelist, $pa, $BS) = @_;
  my $DISPLAY = {};
  my ($pre, $prep, $post) = ("<code style=\"color:", ";\">", "</code>");
  foreach my $gene (@$genelist)
  { 
    my $essstatus = $gene->Tag_essential_status;
    my $orfstatus = $gene->Tag_orf_classification;
    my $displayname = $gene->has_tag("gene")  
                  ?  $gene->Tag_load_id . " (" . $gene->Tag_gene . ")"  
                  :  $gene->Tag_load_id;
    my $fulldisplay;
    if ($essstatus ne "Nonessential")
    {
      $fulldisplay = $pre . $BS->{COLORS}->{gene}->{$essstatus} . $prep . $displayname . $post;
    }
    else
    {
      $fulldisplay = $pre . $BS->{COLORS}->{gene}->{$orfstatus} . $prep . $displayname . $post;
    }
    $DISPLAY->{$gene->Tag_load_id} = make_link($gene, $fulldisplay, $pa, $BS);
  }
  return $DISPLAY;
}

sub make_link
{
  my ($feat, $desc, $pa, $BS) = @_;
  my $chr = $pa->{CHROMOSOME} ? $pa->{CHROMOSOME} : $pa->{OLDCHROMOSOME}  ? $pa->{OLDCHROMOSOME}  : $pa->{NEWCHROMOSOME};
  my $href = "http://$BS->{this_server}/cgi-bin/gb2/gbrowse/$chr/?start=" . $feat->start . ";stop=" . $feat->end . ";ref=" . $pa->{SEQID} . ";";
  return "<a href=\"$href\" target=\"_blank\" style=\"text-decoration:none\">$desc</a>";
}

sub get_gbrowse_src_list
{
	my ($BS) = @_;
	my @srcs;
	find sub { push @srcs, $File::Find::name}, $BS->{conf_repository};
	@srcs = grep {	$_ =~ /\.conf\Z/} @srcs;
  my @sources;
  foreach (@srcs)
  {
    push @sources, $1 if ($_ =~ /($VERNAME)/);
  }
	return @sources;
}

sub write_conf_file
{
	my ($BS, $newgffname, $note) = @_;
	my ($SPECIES, $CHRNAME) = ($1, $2) if ($newgffname =~ $VERNAME);
  my $SEQID    = "chr" . $CHRNAME;
	my $ref = slurp($BS->{reference_conf_file});
	$ref =~ s/\*LANDMARK\*/$CHRNAME/g;
	$ref =~ s/\*VERSION\*/$newgffname/g;
	$ref =~ s/\*VERSIONNOTE\*/$note/g;
  mkdir $BS->{conf_repository} unless (-e $BS->{conf_repository});
  my $spedir = $BS->{conf_repository} . "/" . $SPECIES;
  mkdir $spedir unless (-e $spedir);
  my $confpath = $spedir . "/" . $SEQID . "/";
  mkdir $confpath unless (-e $confpath);
  $confpath .= $newgffname . ".conf";
	open (my $confout, ">", $confpath) || die "can't open new conf file $confpath : $!";
	print $confout $ref;
	close $confout;
  my $confstream = slurp($BS->{GBrowse_conf_file});
  if ($confstream !~ /$newgffname/)
  {
    open (my $gbout, ">>", $BS->{GBrowse_conf_file}) || die "can't open GBrowse conf file $BS->{GBrowse_conf_file} : $!";
    print $gbout "\n[$newgffname]\n";
    print $gbout "description     = $newgffname\n";
    print $gbout "path            = $confpath\n";
    close $gbout;
  }
	return;	
}

sub update_gbrowse
{
	my ($BS, $pa) = @_;
	
	print "Reloading mysql database...\n\n\n";
  BioStudio::MySQL::load_database($pa->{NEWCHROMOSOME}, $BS);
	
	print "Creating new conf file...\n\n\n";
  write_conf_file($BS, $pa->{NEWCHROMOSOME}, "$pa->{NEWCHROMOSOME} created from $pa->{OLDCHROMOSOME} $time{'yymmdd'}", $pa->{SEQID});
    
  my $newlink = "http://" . $BS->{"this_server"} . "/cgi-bin/gb2/gbrowse/$pa->{NEWCHROMOSOME}/";
	print "The new chromosome is <a href=\"$newlink\">$pa->{NEWCHROMOSOME}</a>.\n\n\n";
	return;
}

sub rollback_gbrowse
{
  my ($BS, $pa) = @_;
		
	print "Deleting conf file...";
	system "rm $BS->{conf_repository}/$pa->{DROPPER}.conf" == 0 
	  || print "can't delete conf file!";
  my $GBconf = $BS->{GBrowse_conf_file};
  my $GBconftmp = $BS->{GBrowse_conf_file} . ".tmp";
  my @args = ("sed -e \"/$pa->{DROPPER}/d\" $GBconf >$GBconftmp");
  $SIG{CHLD} = 'DEFAULT';
  system (@args) == 0 || die ("oh no, can't edit $GBconf? $!");
  system "mv $GBconftmp $GBconf";
  return;
}

sub gene_color 
{
  my $feat = shift;
  my $confpath = shift;
  my $BS = configure_BioStudio("/Users/Shared/BioStudio/");
  return "darkblue" unless ($BS);
  return $BS->{COLOR}->{$feat->primary_tag} if (exists $BS->{COLOR}->{$feat->primary_tag});
  my $essstat = $feat->has_tag('essential_status') ? $feat->Tag_essential_status  : "";
  my $orfstat = $feat->has_tag('orf_classification') ? $feat->Tag_orf_classification : "";
  return $BS->{COLORS}->{gene}->{$essstat} if ($essstat eq "Essential");
  return $BS->{COLORS}->{gene}->{$essstat} if ($essstat eq "fast_growth");
  
  return $BS->{COLORS}->{gene}->{"transposable_element"} if ($orfstat =~ /transpos/ig);
  return $BS->{COLORS}->{gene}->{$orfstat} if (exists $BS->{COLORS}->{gene}->{$orfstat});
  return $BS->{COLORS}->{gene}->{$essstat} if (exists $BS->{COLORS}->{gene}->{$essstat});
  return "yellow";
}

1;
__END__

=head1 NAME

BioStudio::GBrowse

=head1 VERSION

Version 1.00

=head1 DESCRIPTION

BioStudio functions for GBrowse

=head1 Functions

=head2 make_link()
  Given a list of Bio::SeqFeature gene objects, a BioStudio parameter hashref,
  and the BioStudio config hashref, generate a hashref such that the key is the
  gene's id and the value is a display friendly link to GBrowse.

=head2 make_link()
  Given a Bio::SeqFeature compliant object, a string description, a BioStudio
  parameter hashref, and a BioStudio config hashref, generate a link for that 
  object in GBrowse.

=head2 get_gbrowse_src_list()
  Given the BioStudio config hashref, return a list of all of the chromosomes
  that are available through GBrowse

=head2 write_conf_file()
  Given the BioStudio config hashref, the name of the chromosome, and a note, 
  create a configuration file by replacing values in the BioStudio template
  of a GBrowse conf file.

=head2 update_gbrowse()
  Given a BioStudio config hashref and a BioStudio parameter hashref, update 
  GBrowse - reload the database, write the configuration file, and pass a link.

=head2 rollback_gbrowse()
  Given a BioStudio config hashref and a BioStudio parameter hashref, remove
  a source from GBrowse. 

=head2 gene_color()
  Designed to be used from within GBrowse. Sets the color of a gene, based on
  its essential status or orf classification

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
