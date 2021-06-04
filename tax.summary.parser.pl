#!/usr/bin/perl -w
use warnings;
use strict;
use FindBin;                # where was script installed?
use lib $FindBin::Bin;      # use that dir for libs, too
use Getopt::Long;
use Data::Dumper;

#parses a mothur-style tax.summary into a table for plotting simple barplots in Excel/R
#script expects the taxonomic levels defined below

my $quiet=1;
my %taxlvl = qw ( 0 root 1 bacteria 2 phylum 3 class 4 order 5 family 6 genus 7 species);
my %rankID;

my $opt = {taxlevel => 'phylum', keep => '', maxtax => 15};
if (! @ARGV){usage()};
GetOptions ($opt, 'in=s', 'maxtax=i', 'taxlevel=s', 'keep=s', 'prepend_phylum', 'prepend_family', 'prepend_order', 'prepend_class', 'prepend_genus' );
if ($opt->{help}){usage()};

#load hash of user defined "keeps"
my %keep;
if ($opt->{'keep'} ne "")
{
    my @k=split(/\,/,$opt->{'keep'});
    foreach my $k (@k)
    {
	$keep{$k}=1;
    }
}
if (scalar(keys(%keep))>0)
  {
    print STDERR "Keeping top ".$opt->{'maxtax'}." tax AND ".join(" ",keys(%keep))."\n";
  }
else
  {
    print STDERR "Keeping top ".$opt->{'maxtax'}." tax\n";
  }

open(TAX,"<$opt->{'in'}") or die "can not open input file\n";
my $newlin=<TAX>;
chomp($newlin);
my @header=split(/\t/,$newlin);
my %select;
while (my $newlin=<TAX>)
  {
    chomp($newlin);
    $newlin=~s/\"/_/g;
    $newlin=~s/\'/_/g;

    my @dat=split(/\t/,$newlin);
    my ($taxlevel,$rankID,$taxon,$daughterlevels,$total,@samplecounts)=split(/\t/,$newlin);
    $rankID{$rankID}=$taxon;
    if ($taxlvl{$taxlevel} eq $opt->{taxlevel})
      {
	$select{$rankID}=[ @dat ];
      }
  }
close(TAX);

my @tax_list;
#sort by taxon with the largest total counts
print join("\t",@header)."\n";
my $t=0;
foreach my $rid (sort {$select{$b}[4]<=>$select{$a}[4]} keys %select)
{
    #This is a very kludgy way to pull the parent ranks for appending
    if ($opt->{'prepend_genus'})
    {
	if ($select{$rid}[1]=~/^(\d+\.\d+\.\d+\.\d+\.\d+\.\d+\.\d+)\./)
	{
	    $select{$rid}[2]=$rankID{$1}.":".$select{$rid}[2];
	}
    }

    if ($opt->{'prepend_family'})
    {
	if ($select{$rid}[1]=~/^(\d+\.\d+\.\d+\.\d+\.\d+\.\d+)\./)
	{
	    $select{$rid}[2]=$rankID{$1}.":".$select{$rid}[2];
	}
    }
    
    if ($opt->{'prepend_order'})
    {
	if ($select{$rid}[1]=~/^(\d+\.\d+\.\d+\.\d+\.\d+)\./)
	{
	    $select{$rid}[2]=$rankID{$1}.":".$select{$rid}[2];
	}
    }
    
    if ($opt->{'prepend_class'})
    {
	if ($select{$rid}[1]=~/^(\d+\.\d+\.\d+\.\d+)\./)
	{
	    $select{$rid}[2]=$rankID{$1}.":".$select{$rid}[2];
	}
    }
    
    if ($opt->{'prepend_phylum'})
    {
	if ($select{$rid}[1]=~/^(\d+\.\d+\.\d+)\./)
	{
	    $select{$rid}[2]=$rankID{$1}.":".$select{$rid}[2];
	}
    }
    
    if ($select{$rid}[2] eq 'unclassified')
    {
	my $parent=$select{$rid}[1];
	$parent=~s/\.\d+$//;
	$select{$rid}[2]=$rankID{$parent}.":".$select{$rid}[2];
    }
    
    #write values
    if (exists($opt->{'maxtax'}) && $t >= $opt->{'maxtax'} && ! exists($keep{$rid}) )
    {
	#if not one of the maxtax of use selected...
	$select{'_other'}[0]=0;
	$select{'_other'}[1]=0;
	$select{'_other'}[2]='_other';
	$select{'_other'}[3]=0;
	for (my $i=4;$i<@{$select{$rid}};$i++)
	{
	    if ($select{$rid}[$i] !~ /^\d+$/)
	    {
		print STDERR "$rid (".$rankID{$rid}.") has non-numeric value (".$select{$rid}[$i]."), skipping cell\n";
		next;
	    }
	    $select{'_other'}[$i]+=$select{$rid}[$i];
	}
    }
    else
    {
	print join("\t",@{$select{$rid}})."\n";
	push(@tax_list,$select{$rid}[2]); 
    }
    $t++;
}
if (exists($select{'_other'})){print join("\t",@{$select{'_other'}})."\n"};

if ($quiet==0){
  print STDERR Dumper(\@tax_list);
  print STDERR scalar(keys(@tax_list))." taxa in table\n";
}

###############################
# SUBROUTINES

sub usage
  {
    print STDERR "tax.summary.parser.pl -in mothur.tax.summary [-taxlevel][-maxtax][-prepend_phylum|order|class|family]\n";
    exit(0);
  }
