#!/usr/bin/perl
use strict;
use warnings;

#in: blastx output file, -outfmt 6
my $in = shift @ARGV or die $!;
#fa: fasta file for the query
my $fa = shift @ARGV or die $!;
my $cov_cutoff = shift @ARGV or die $!;
my $pid_cutoff = shift @ARGV or die $!;
#out: PID>=%70 and coverage>=60% for passing matches
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
    
#MS4335685-600V3_1_1119_16002_17873/1	73912456_Root;Fungi;Ascomycota;Eurotiomycetes;Onygenales;Onygenales_family_incertae_sedis;Coccidioides	99.500	200	1	0	1	200	356	157	1.32e-98	356	100
#MS4335685-600V3_1_1119_16002_17873/1	157678735_Root;Fungi;Ascomycota;Eurotiomycetes;Onygenales;Onygenaceae;Uncinocarpus	88.500	200	14	4	1	198	356	164	4.94e-66	248	99
#MS4335685-600V3_1_1119_16002_17873/1	28975603_Root;Fungi;Ascomycota;Eurotiomycetes;Onygenales;Onygenaceae;Aphanoascus	88.325	197	13	4	1	197	319	133	2.10e-64	242	99
#MS4335685-600V3_1_1119_16002_17873/1	28975607_Root;Fungi;Ascomycota;Eurotiomycetes;Onygenales;Onygenales_family_incertae_sedis;Chrysosporium	87.879	198	14	4	1	197	322	134	7.33e-64	241	99
#MS4335685-600V3_1_1119_16002_17873/1	213876640_Root;Fungi;Ascomycota;Eurotiomycetes;Onygenales;unclassified_Onygenales_family;unclassified_Onygenales_genus	87.310	197	15	4	1	197	316	130	1.09e-61	233	99
#MS4335685-600V3_1_1119_16002_17873/1	28975594_Root;Fungi;Ascomycota;Eurotiomycetes;Onygenales;Onygenaceae;Aphanoascus	87.310	197	14	5	1	197	319	134	1.32e-60	230	99
#MS4335685-600V3_1_1119_16002_17873/1	157678744_Root;Fungi;Ascomycota;Eurotiomycetes;Onygenales;Onygenales_family_incertae_sedis;Chrysosporium	86.802	197	16	4	1	197	344	158	1.32e-60	230	99
#MS4335685-600V3_1_1119_16002_17873/1	28975598_Root;Fungi;Ascomycota;Eurotiomycetes;Onygenales;Onygenaceae;Aphanoascus	86.869	198	15	5	1	197	321	134	4.62e-60	228	99
#MS4335685-600V3_1_1119_16002_17873/1	7208636_Root;Fungi;Ascomycota;Eurotiomycetes;Onygenales;Onygenaceae;Aphanoascus	85.787	197	18	4	1	197	306	120	6.86e-58	221	99
#MS4335685-600V3_1_1119_16002_17873/1	157678743_Root;Fungi;Ascomycota;Eurotiomycetes;Onygenales;Onygenaceae;Uncinocarpus	84.343	198	25	2	1	198	353	162	1.02e-55	214	99

#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs

    my @array = split /\s+/, $line;
    my $query = $array[0];
    my $hit = $array[1];
    my $pid = $array[2];
    my $bitscore = $array[-2];
    my $qcov = $array[-1];
    if ($qcov >= $cov_cutoff && $pid >= $pid_cutoff) {
	$hash{$query} -> {$bitscore} -> {$pid} -> {$qcov} -> {$hit} = "$query\t$hit\t$pid\t$bitscore\t$qcov";
    }
}
close $IN;
#my $hash_key_num = scalar keys %hash;
#print "$hash_key_num\n";
#exit;

open ( my $OUT, ">$out" ) or die "can't open $out\n";
foreach my $query_key ( sort keys %hash ){ 
    my $score_counter;
    my $genus;
    foreach my $score_key (sort {$b<=>$a} keys %{$hash{$query_key}}) {
	$score_counter++;
	my $pid_counter;
	if ($score_counter == 1) {
	    #my $hash_key_num = scalar keys %{$hash{$query_key}{$score_key}};
		#if ($hash_key_num > 1) {
		    #foreach my $hit_key (sort keys %{$hash{$query_key}{$score_key}}) {
			#print "$hash{$query_key}{$score_key}{$hit_key}\n";
		    #}
	    #}
	    foreach my $pid_key (sort {$b<=>$a} keys %{$hash{$query_key}{$score_key}}) {
		$pid_counter++;
		my $qcov_counter;
		if ($pid_counter == 1) {
		    foreach my $qcov_key (sort {$b<=>$a} keys %{$hash{$query_key}{$score_key}{$pid_key}}) {
			$qcov_counter++;
			my $hit_counter;
			if ($qcov_counter == 1) {
			    foreach my $hit_key (sort keys %{$hash{$query_key}{$score_key}{$pid_key}{$qcov_key}}) {
				$hit_counter++;
				if ($hit_counter == 1 && $hit_key =~ m/_Root;(\S+)/) { 
				    $genus = $1;
				}
			    }
			}
		    }
		}
	    }
	}
    }
    print $OUT "$query_key\t$genus;\n";
}
foreach my $def_key (sort keys %fa_hash) {
    unless (defined $hash{$def_key}) {
	print $OUT "$def_key\tunknown;unknown_unclassified(0);unknown_unclassified(0);unknown_unclassified(0);unknown_unclassified(0);unknown_unclassified(0);\n";
    }
}
close $OUT;
exit;
