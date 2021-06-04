#!/usr/bin/perl
use strict;
use warnings;

#in: blastx output file, -outfmt 6
my $in = shift @ARGV or die $!;
#fasta file for the query sequences to get trimmed sequence
my $fa = shift @ARGV or die $!;
#out: PID and coverage for matches, separated for intra- and inter-specific matches
#my $out1 = shift @ARGV or die $!;
#my $out2 = shift @ARGV or die $!;

my %fa_hash = ();
open (my $FA, $fa) or die "can't open $fa\n";
my $def;
while (my $fa_line = <$FA>){
    chomp $fa_line;
    if ($fa_line =~ /^>(\S+)/){
	$def = $1;
    }
    else{
	$fa_hash{$def} .= $fa_line;
    }
}
close $FA;

my %hash = ();
#open ( my $OUT1, ">$out1" ) or die "can't open $out1\n";
#open ( my $OUT2, ">$out2" ) or die "can't open $out2\n";
open ( my $IN , $in ) or die "can't open $in\n";
while ( my $line = <$IN> ){
    chomp $line;
    
#AB040215.1_Candida_albicans_	AB040215.1_Candida_albicans_	100.000	139	0	0	1	139	1	139	3.77e-72	257
#AB040215.1_Candida_albicans_	AB365317.1_Candida_albicans_E	99.286	140	0	1	1	139	1	140	1.76e-70	252
#AB040215.1_Candida_albicans_	HE860438.1_Candida_albicans_	98.561	139	2	0	1	139	1	139	8.17e-69	246
#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

    my @array = split /\s+/, $line;
    my $query = $array[0];
    my $hit = $array[1];
    #my $pid = $array[2];
    my $bitscore = $array[-1];
    #my $aln_len = $array[3]; 
    $hash{$query} -> {$bitscore} -> {$hit} = "$query\t$hit";
}
close $IN;
#my $hash_key_num = scalar keys %hash;
#print "$hash_key_num\n";
#exit;

#my $cov;
my $best_hit_counter_intra = 0;
my $best_hit_counter_inter = 0;
foreach my $query_key ( sort keys %hash ){
    my @query_array = split /_/, $query_key;
    my $query_species = $query_array[2];
    my $counter = 0;
    foreach my $score_key (sort {$b<=>$a} keys %{$hash{$query_key}}) {
	#next if $query_key eq $hit_key;
	$counter++;
	if ($counter == 1) {
	    foreach my $hit_key (sort keys %{$hash{$query_key}{$score_key}}) {
		my @hit_array = split /_/, $hit_key;
		my $hit_species = $hit_array[2];
		#my @array = split /\t/, $hash{$query_key}{$hit_key};
		$best_hit_counter_intra++ if ($query_species eq $hit_species);
		$best_hit_counter_inter++ if ($query_species ne $hit_species); 
	    }
	}
    }
}
print "The number of best hits that are within the same species: $best_hit_counter_intra\nThe number of best hits that are not within the same species: $best_hit_counter_inter\n";
#close $OUT1;
#close $OUT2;
exit;
