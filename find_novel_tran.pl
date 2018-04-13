#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;

my $file=$ARGV[0]; ### ref gtf file
my $gtf=$ARGV[1];  ### gtf file assembly by cufflinks or cuffmerge
my %seqinfor;
my $outdir;
my $gtf_basename=basename($gtf);
GetOptions
(
	"o:s"=>\$outdir,
);
`mkdir $outdir` if(!-d "$outdir");

open OUT1,">$outdir/$gtf_basename.novel_transcripts.gtf\n";

open IN,$file || die "$!";
while(<IN>){
	chomp;
	next if /^#/;
	my @c=split(/\t/,$_);
	@c[3,4]=@c[4,3] if($c[3]>$c[4]);
	if($file=~/gtf$/){
		next unless($c[2] =~ /exon/);
		my $gene_id=$1 if($c[8]=~/gene_id\s+"(\S+)";/);
		my $tran_id=$1 if($c[8]=~/transcript_id\s+"(\S+)";/);
		push @{$seqinfor{$c[0]}{$gene_id}{$tran_id}{exon}},[@c];   ### Here, mRNA means transcriptome	
	}
	if($file=~/gff3$/ or $file=~/gff$/){
		my ($gene_id,$tran_id);
		if($c[2]=~/gene/){
			$gene_id=$1 if($c[8]=~/ID=([^;]+);/);
		}if($c[2]=~/CDS/ or $c[2]=~/exon/){
			$tran_id=$1 if($c[8]=~/Parent=([^;]+);/ or $c[8]=~/PARENT=([^;]+);/);
			push @{$seqinfor{$c[0]}{$gene_id}{$tran_id}{$c[2]}},[@c];
		}
	}
}
close IN;
my %pos_gene;
my %flag;
foreach my $scaf(sort keys %seqinfor){
	foreach my $gene(sort keys %{$seqinfor{$scaf}}){
		foreach my $id(sort keys %{$seqinfor{$scaf}{$gene}}){
			if(!$seqinfor{$scaf}{$gene}{$id}{exon}){
				@{$seqinfor{$scaf}{$gene}{$id}{exon}}=@{$seqinfor{$scaf}{$gene}{$id}{CDS}};
			}
			@{$seqinfor{$scaf}{$gene}{$id}{exon}}=sort {$a->[3]<=>$b->[3]} @{$seqinfor{$scaf}{$gene}{$id}{exon}};
			@{$pos_gene{$scaf}{$seqinfor{$scaf}{$gene}{$id}{'exon'}[0][6]}{$id}}=($seqinfor{$scaf}{$gene}{$id}{'exon'}[0][3],$seqinfor{$scaf}{$gene}{$id}{'exon'}[-1][4]);
		}
	}
}
my (%know,%novel,%total);
open IN,$gtf;
while(<IN>){
	chomp;
	my @c=split("\t",$_);
	next if($c[2] eq "gene");
	my $id=$1 if($c[8]=~/transcript_id "(\S+)";/);
	my $gene=$1 if($c[8]=~/gene_id "(\S+)";/);
	$total{'tran'}{$id}=1;
	$total{'gene'}{$gene}=1;
	if($c[2] eq "transcript" or $c[2] eq "mRNA"){
		$flag{$id}=0;
		foreach my $ref_id(sort keys %{$pos_gene{$c[0]}{$c[6]}}){
			if($c[3]<$pos_gene{$c[0]}{$c[6]}{$ref_id}[1] && $c[4]>$pos_gene{$c[0]}{$c[6]}{$ref_id}[0]){
				$flag{$id}=1;
				$know{'tran'}{$id}=1;
				$know{'gene'}{$gene}=1;
			}	
		}
		if($flag{$id}==0){
			print OUT1"$_\n";
			print "$gene\t$c[0]\t$c[3]\t$c[4]\n";
			$novel{'tran'}{$id}=1;
			$novel{'gene'}{$gene}=1;
		}
	}else{
		print OUT1"$_\n" if($flag{$id}==0);
	}
}
close IN;
close OUT1;

open OUT,">$outdir/$gtf_basename.known_and_novel_stat";
my $a=keys %{$know{'tran'}};
my $b=keys %{$know{'gene'}};
my $c=keys %{$novel{'tran'}};
my $d=keys %{$novel{'gene'}};
my $e=keys %{$total{'tran'}};
my $f=keys %{$total{'gene'}};
print OUT"Total Assembled Transcripts\t$e\nTotal Assembled Loci\t$f\n";
print OUT"In Known Transcripts\t$a\nCorresponding Known Loci\t$b\n";
print OUT"In Novel Transcripts\t$c\nCorresponding Novel Loci\t$d\n";
close OUT;

