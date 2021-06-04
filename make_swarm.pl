#!/usr/bin/perl
use warnings;
use strict;
use FindBin;                # where was script installed?
use lib $FindBin::Bin;      # use that dir for libs, too
use File::Temp qw/ tempfile tempdir /;

my $debug=0;
my $exec="now";



if (@ARGV<3)
{
    print STDERR "creates a swarm file from a command with [replacements] and a file list\n";
    print STDERR "make_swarm.pl command [FOO] [BAR] [BAZ] -- (file1 file2 ...|my.filelist)\n";
    print STDERR " - remember to \\escape quotes\n";
    print STDERR " - FOO will be replaced with the filename\n";
    print STDERR " - BAR will be replaced with root filename (no path, last .ext removed)\n";
    print STDERR " - BAZ will be replaced with bare file (assuming FOO=file.R1.name)\n";
    print STDERR " - RAZ will be replaced the filename R2 replacing R1\n";
    print STDERR "ex: make_swarm.pl now cp FOO BAR.bak -- *.pl #make backup copies of pl files\n";
    print STDERR "swarm file will get a unique name swarmXXXX\n";
    exit(0);
}

#get command and files
my @files;
my $command;
my @opts=@ARGV;
my $opt="";

while ($opt ne '--')
{
    $opt=shift(@opts);
    if ($opt ne '--'){$command.=$opt.' '};
}

#detect if a list file was supplied instead of files
if (scalar(@opts)==1 && $opts[0]=~/\.filelist/)
  {
     #a list file was supplied instead of file names
     open(TMP,"<$opts[0]") or die "can't open list file\n";
     while (my $newlin=<TMP>)
       {
         chomp($newlin);
         if ($newlin=~/^\#/ || $newlin=~/^\s*$/){next};
         push(@files,$newlin);
       }
     close(TMP);
  }
else
  {
     @files=@opts;
  }

print STDERR "$command\n".join(",",@files)."\n" if ($debug);

my ($ofh, $ofile) = tempfile('swarmXXXX', UNLINK => 0);
foreach my $f (@files)
{
    my $com=$command;
    my $root=$f;
    my $noext=$f;
    my $readsub=$f;
    $root=~s/.+\///g;
    $root=~s/\.[^\.]+$//;
    #$noext=~s/.+\///g;  #path
    $noext=~s/\..+$//g; #extensions, this would remove long things that aren't really extensions
    $noext=~s/(\.\w{2,5}){1,3}$//g; #extensions, extensions should be 2-5 characters gz,fastq, up to 3 of them
    $readsub=~s/\.R1\./.R2./; #it would be better if this could handle _R1_ and _R1. etc...
    $com=~s/FOO/$f/g;
    $com=~s/BAR/$root/g;
    $com=~s/BAZ/$noext/g;
    $com=~s/RAZ/$readsub/g;

    print $ofh "$com\n";
}
print $ofile."\n";
close($ofh);
