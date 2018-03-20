#! /usr/bin/env perl
use strict;
use warnings;
use Test::Script::Run;
use Test::More;

my ($return,$stdout,$stderr)=run_script('../gff2bed12.pl',["test.gff"]);

my $ref = "test.bed12";
open(REF,"<",$ref)||die;
my $ref_string;
while(<REF>){
	$ref_string .= $_;
}
close(REF)||die;

my $got = &parser($stdout);
my $expected = &parser($ref_string);

is_deeply($got,$expected,'BED12 output as expected');

done_testing();

sub parser{
	my $p		= $_[0];
	my @p_array 	= ();
	my %p_hash;
	foreach my $bed_line (split("\n",$p)){
		my @bed_fields 		= split("\t",$bed_line);
		my $bed_string		= join("",@bed_fields);
		$p_hash{$bed_fields[6]} = $bed_line;
	}
	foreach my $key (sort {$a <=> $b} keys %p_hash){
		push(@p_array,$p_hash{$key});
	}
	return(\@p_array);
}
