#! /usr/bin/perl
use strict;
use warnings;

# Author Daniel.Amsel@ime.fraunhofer.de

# BED12
# For example: chr1 11873 14409 uc001aaa.3 0 + 11873 11873 0 3 354,109,1189, 0,739,1347,
# chrom         := chromosome
# start         := zero-based [
# end           := one-based excluded )
# name          := feature name
# score         := .
# strand        := [+|-]
# thickStart    := start codon position (CDS start)
# thickEnd      := stop codon position  (CDS stop)
# itemRgb       := .
# blockCount    := number of exons 
# blockSizes    := list of block sizes
# blockStarts   := list of block starts ; starting with the first nucleotide of the transcript 



my $annot_file	= $ARGV[0]; # gff format
my %annotation;	# {rna1234}  from Parent=rna1234 - second column in description
open (FILE,"<",$annot_file) || die "Can not open $annot_file\n";
while(<FILE>){
	chomp;
	next if (/^#/);
	my ($chr,undef,$type,$start,$stop,undef,$strand,undef,$description)    = split("\t",$_);
	my @description_list	= split(";",$description);
	my @parent_id_split	= split("=",$description_list[1]);
	my $parent_id		= $parent_id_split[1];
	next unless ($parent_id =~/^rna/);
	push(@{$annotation{$parent_id}{$type}},$start);
        push(@{$annotation{$parent_id}{$type}},$stop);

	# use the exon description to add an info key for strand and transcript_id 
	if((lc($type) eq "exon") and (not exists $annotation{$parent_id}{"info"}) and (index($description, "transcript_id") != -1)){
	        my @description_tail    = split("transcript_id",$description);
	        my $transcript_id;
	        ($transcript_id = $description_tail[1]) =~/\w{2}_\d+\.\d+/;
        	$transcript_id          =~ s/^\s+|\s+$|\"|=|;//g; 
		my @info	= ($chr,$strand,$transcript_id);
		$annotation{$parent_id}{"info"} = \@info;
	}
	
}
close(FILE) || die "Can not close $annot_file\n";


foreach ( keys(%annotation) ){
	my $pid		= $_;
	my %pid_hash	= %{$annotation{$pid}};

	next if (not exists $pid_hash{"exon"});
	next if (not exists $pid_hash{"CDS"});
	next if (not exists $pid_hash{"info"});

	my @exons		= sort {$a <=> $b} @{$pid_hash{"exon"}};
	my @cds			= sort {$a <=> $b} @{$pid_hash{"CDS"}};
	my @info		= @{$pid_hash{"info"}};
	

	my $bed_start		= $exons[0] - 1;
	my $bed_thick_start 	= $cds[0]   - 1;
	my $block_count		= scalar(@exons)/2;
	

        my @block_lengths;              # exon1_stop - exon1_start + 1
        my @block_start_positions;
	push(@block_start_positions,0);
        for (my $i = 0; $i < scalar(@exons)-2;$i++){
                my $tmp_block_start     = $exons[$i];
                $i++;
                my $tmp_block_stop      = $exons[$i];
                my $tmp_block_len       = $tmp_block_stop - $tmp_block_start + 1;
                push(@block_lengths,$tmp_block_len);

                my $tmp_block_start_2   = $exons[$i+1];
                my $tmp_block_start_positions = $tmp_block_start_2 - $tmp_block_start;
                push(@block_start_positions,$tmp_block_start_positions);
        }
	push(@block_lengths,($exons[-1]-$exons[-2]+1));

        my $block_sizes         = join(",",@block_lengths);
        my $block_starts        = join(",",@block_start_positions);

	print "$info[0]\t$bed_start\t$exons[-1]\t$info[2]\t0\t$info[1]\t$bed_thick_start\t$cds[-1]\t0\t$block_count\t$block_sizes\t$block_starts\n";

}



