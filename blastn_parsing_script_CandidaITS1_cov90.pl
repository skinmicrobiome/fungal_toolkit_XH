#!/usr/bin/perl
use strict;
use warnings;

#in: blastx output file, -outfmt 6
my $in = shift @ARGV or die $!;
#fa: fasta file for the query
my $fa = shift @ARGV or die $!;
my $aln_len_cutoff = shift @ARGV or die $!;
my $pid_cutoff = shift @ARGV or die $!;
#out: PID>=%95 and alignment length>=50bp for passing matches
my $out = shift @ARGV or die $!;

my %fa_hash = ();
open (my $FA, $fa) or die "can't open $fa\n";
my $def;
while (my $fa_line = <$FA>){
    chomp $fa_line;
    if ($fa_line =~ m/^>(\S+)/) {
	$def = $1;
	$fa_hash{$def} = 1;
    }
}
close $FA;

my %hash = ();
open ( my $IN , $in ) or die "can't open $in\n";
while ( my $line = <$IN> ){
    chomp $line;
    
#ICX1UM101DZNY2  AB054035.1_Candida_tropicalis_  95.714  140     1       4       84      220     138     1       3.89e-58        212
#ICX1UM101DZNY2  LT837795.1_Candida_tropicalis_K 95.000  140     2       4       84      220     138     1       4.74e-57        208
#ICX1UM101DZNY2  AB437070.1_Candida_tropicalis_J 95.000  140     2       4       84      220     138     1       4.74e-57        208
#ICX1UM101DZNY2  KC408969.1_Candida_tropicalis_  94.964  139     2       4       84      219     137     1       1.66e-56        206
#ICX1UM101DZNY2  KC408966.1_Candida_tropicalis_  94.964  139     2       4       84      219     137     1       1.66e-56        206
#ICX1UM101DZNY2  KC408971.1_Candida_tropicalis_  94.286  140     2       4       84      219     138     1       5.78e-56        205
#ICX1UM101DZNY2  AB437077.1_Candida_tropicalis_J 94.286  140     3       4       84      220     138     1       2.02e-55        203
#ICX1UM101DZNY2  HM771639.1_Candida_tropicalis_C 94.286  140     2       5       84      220     137     1       2.46e-54        199
#ICX1UM101DZNY2  AB467291.1_Candida_tropicalis_T 94.286  140     1       6       84      220     136     1       2.99e-53        196
#ICX1UM101DZNY2  AB369915.1_Candida_albicans_B   75.524  143     27      6       83      220     140     1       1.66e-18        80.6

#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

    my @array = split /\s+/, $line;
    my $query = $array[0];
    my $hit = $array[1];
    my $pid = $array[2];
    my $aln_len = $array[3];
    my $bitscore = $array[-1];
    if ($aln_len >= $aln_len_cutoff && $pid >= $pid_cutoff) {
	$hash{$query} -> {$bitscore} -> {$hit} = "$query\t$hit\t$bitscore\t$pid\t$aln_len";
    }
}
close $IN;
#my $hash_key_num = scalar keys %hash;
#print "$hash_key_num\n";
#exit;

open ( my $OUT, ">$out" ) or die "can't open $out\n";
foreach my $query_key ( sort keys %hash ){ 
    my $score_counter;
    my $species;
    foreach my $score_key (sort {$b<=>$a} keys %{$hash{$query_key}}) {
	$score_counter++;
	if ($score_counter == 1) {
	    foreach my $hit_key (sort keys %{$hash{$query_key}{$score_key}}) {
		my @species_array = split /_/,$hit_key;
		    $species = "$species_array[1]"."_"."$species_array[2]";
	    }
	}
    }
    print $OUT "$query_key\tFungi;Ascomycota;Saccharomycetes;Saccharomycetales;Saccharomycetales_family_incertae_sedis;Candida;$species;\n";
}
foreach my $def_key (sort keys %fa_hash) {
    unless (defined $hash{$def_key}) {
	print $OUT "$def_key\tunknown;unknown_unclassified(0);unknown_unclassified(0);unknown_unclassified(0);unknown_unclassified(0);unknown_unclassified(0);unknown_unclassified(0);\n";
    }
}
close $OUT;
exit;
