#!/usr/bin/perl
use strict;
use warnings;
use POSIX;

my $in = shift @ARGV or die $!;
my $out = shift @ARGV or die $!;

my %fa_hash = ();
open (my $FA, $in) or die "can't open $in\n";
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

open (my $OUT, ">$out") or die "can't open $out\n";
my $seq_len;
my $errors;
my $random_position;
my $old_base;
my $random_base;
my $new_seq;
my @base_array = ('A','T','C','G');
foreach my $key (keys %fa_hash) {
    $seq_len = length($fa_hash{$key});
    #print "$seq_len\n";
    $errors = ceil($seq_len * 0.01);
    #print "$errors\n";
    for (1..$errors) {
	$random_position = int(rand($seq_len));
	#print "$random_position\n";
	$old_base = substr($fa_hash{$key},$random_position,1);
	#print "$old_base\n";
      MARK: $random_base = $base_array[rand(@base_array)];
	#print "$random_base\n";
	if ($random_base ne $old_base) {
	    $new_seq = substr($fa_hash{$key},$random_position,1,$random_base);
	}
	else {
	    goto MARK;
	}
    }
    print $OUT ">$key\n$fa_hash{$key}\n";
}
close $OUT;
exit;

