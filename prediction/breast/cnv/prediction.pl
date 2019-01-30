#!/usr/bin/perl

if(scalar(@ARGV)<1) {
	print "perl cv_partition.pl 0\n";
	die;
}

$path0="../../../data/trimmed_set/";
$path1="../../../data/normalization/rna_imputation_normalize_by_sample_avg_sd/"; 
$path2="../../../data/cv_set/";
$path3="./"; 

%cnv=();
open PROTEOME, "${path1}trimmed_breast_proteome.txt.scaled.avg_imputation" or die; 
$line=<PROTEOME>;
chomp $line;
@ids=split "\t", $line;
shift @ids;
while($line=<PROTEOME>){
        chomp $line;
        @vals=split "\t", $line;
        $gene=shift @vals;
        $i=0;
        foreach $sample (@ids){
                $cnv{$sample}{$gene}=$vals[$i];
                $i++;
        }
}
close PROTEOME;

#%cnv=();
open PROTEOME, "${path1}trimmed_breast_cnv.txt.scaled.avg_imputation" or die; 
$line=<PROTEOME>;
chomp $line;
@ids=split "\t", $line;
shift @ids;
while($line=<PROTEOME>){
        chomp $line;
        @vals=split "\t", $line;
        $gene=shift @vals;
        $i=0;
        foreach $sample (@ids){
                $cnv{$sample}{$gene}=$vals[$i];
                $i++;
        }
}
close PROTEOME;

foreach $num (@ARGV){
	open TEST_P, "${path2}test_breast_phospho_${num}.txt" or die;	
	$line=<TEST_P>;
	chomp $line;
	@test_ids=split "\t", $line;
	shift @test_ids;

	open OUTPUT_PRED, ">${path3}prediction_breast_phospho_${num}.txt" or die;
	print OUTPUT_PRED "$line\n";

	while($line=<TEST_P>){
		chomp $line;
		@vals=split "\t", $line;
		$gene=shift @vals;

		print OUTPUT_PRED "$gene";
		foreach $id (@test_ids){
			$tmp=$cnv{$id}{$gene};
			print OUTPUT_PRED "\t$tmp";
		}
		print OUTPUT_PRED "\n";
	}
	close TEST_P;
	close OUTPUT_PRED;
}

