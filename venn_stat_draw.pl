#!/usr/bin/perl
#!/usr/bin/perl
use strict;
use warnings;
use FileHandle;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw(basename dirname);
use Cwd;
use FindBin qw($Bin $Script);
use List::Util qw(max min);
chomp $Bin;
#####时间和版本#######################################################################################################
my $BEGIN_TIME=time();
my $version="1.0";
########################################

####使用程序需要输入的参数############################################################################################
my ($Sample_room_with_special,$in_file,@sample_headname,$lower_limit,$Diff_characteristic,@sample_file,$sample_file_head,$od);
GetOptions(		
	"help|?" =>\&USAGE,
	"SRWS!"=>\$Sample_room_with_special,
	"in_file:s"=>\$in_file,
	"sample_headname:s{,}"=>\@sample_headname,
	"lower_limit:s"=>\$lower_limit,
	"DC!"=>\$Diff_characteristic,
	"sample_file:s{,}"=>\@sample_file,
	"sample_file_head!"=>\$sample_file_head,
	"od:s"=>\$od,
) or &USAGE;
&USAGE unless (($Sample_room_with_special and $in_file and @sample_headname and $lower_limit) or ($Diff_characteristic and @sample_file) and $od );
&MKDIR($od);
$od=ABSOLUTE_DIR($od);

if ($Sample_room_with_special) {
	unless ($in_file and @sample_headname and $lower_limit) {print "please input -in_file and -sample_headname and -lower_limit\n";die;	}
	my %sample_head;
	foreach my $sample_head (@sample_headname) {
		my ($sample,$head)=split /,/,$sample_head;
		$sample_head{$sample}=$head;
	}
	my $sample_num_temp=scalar @sample_headname;
	$in_file=ABSOLUTE_DIR($in_file);
	open (IN,$in_file) or die $!;
	my $file_head=<IN>;
	chomp $file_head;
	my @file_head=split /\t/,$file_head;
	my %head_col;	
	foreach my $col (0..$#file_head) {
		$head_col{$file_head[$col]}=$col;
	}
	close IN;
	&MKDIR("$od/sample_list");
	foreach my $sample (keys %sample_head) {
		open (IN,$in_file) or die $!;
		open OUT,">$od/sample_list/$sample.list";
		my $file_temp_head=<IN>;
		while (<IN>) {
			chomp;
			my ($id,$num)=(split /\t/,$_)[0,$head_col{$sample_head{$sample}}];
			if ($num>=$lower_limit) {
				print OUT"$id\n";
			}
		}
		close OUT;
	}
	my $all_sample_file;
	foreach my $sample (keys %sample_head) {
		$all_sample_file .= "$sample,$od/sample_list/$sample.list ";
	}
	my $cmd1="perl $Bin/veen_statistic_temp.pl -f $all_sample_file -od $od >$od/veen_statistic.out 2>$od/veen_statistic.err";
	`$cmd1`;
	if($sample_num_temp <= 5){
		my $cmd2="perl $Bin/drawVeenv1.2.pl -in $od/list.txt -od $od -p venn >$od/veen_draw.out 2>$od/veen_draw.err";
		`$cmd2`;
	}
	if (-s "$od/stat.txt") {
		open (IN,"$od/stat.txt") or die $!;
		my %num_gene;
		while (<IN>) {
			chomp;
			my ($ids,$genes)=(split /\s+/,$_)[0,2];

			my @ids=split /,/,$ids;
			my $id_num=scalar @ids;
			$num_gene{$id_num}=$genes;
		}
		close IN;

	}


}elsif ($Diff_characteristic) {
	unless ($Diff_characteristic and @sample_file) {print "please input -sample_file\n";die;}
	my $all_sample_file=join " ",@sample_file;
	my $sample_num_temp=scalar @sample_file;
	my $real_head;
	if ($sample_file_head) {
		$real_head=" -head ";
	}else{
		$real_head="";	
	}
	my $cmd1="perl $Bin/veen_statistic_temp.pl -f $all_sample_file -od $od $real_head >$od/veen_statistic.out 2>$od/veen_statistic.err";
	`$cmd1`;
	if($sample_num_temp<=5){
		my $cmd2="perl $Bin/drawVeenv1.2.pl -in $od/list.txt -od $od -p venn >$od/veen_draw.out 2>$od/veen_draw.err";
		`$cmd2`;
	}
	if (-s "$od/stat.txt") {
		open (IN,"$od/stat.txt") or die $!;
		my %num_gene;
		while (<IN>) {
			chomp;
			my ($ids,$genes)=(split /\s+/,$_)[0,2];

			my @ids=split /,/,$ids;
			my $id_num=scalar @ids;
			$num_gene{$id_num}=$genes;
		}
		close IN;
	}

}else{
	print "please choose the SRWS or DC\n\n";
}


###################################################################################################################

sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	#rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

#####输入说明子程序####################################################################################################
sub USAGE
{
my $usage=<<"USAGE";
Program: 
Version: $version
Contact: lv ran <lvran18931993760\@163.com>
Description:

Usage:
if the input file(infile) is like :
	ID	T1	T2	T3	T4	......
	gene1	1	2	1	4	......
	gene2	5	6	9	1	......
	.......
	perl $0 -SRWS -in_file infile -sample_headname sample1,sample1 sample2,sample2 sample3,sample3 -lower_limit 0 -od outputdir
	-SRWS	the input file format is like fpkm result;
	=in_file	the input file;
	-sample_headname	sample name,sample header in input file(sample name should not start be number);
	-lower_limit	the lower limit if gene is expression (should be bigger than 0)
	-od output dir
if the input file id like:
	T1.list
	ID ......
	gene1	......
	gene2	......
	.....	......

	T2.list
	ID ......
	gene3	......
	gene2	......
	.....	......
	perl $0 -DC -sample_file T1,T1.list T2,T2.list -sample_file_head -od outputdir
	-DC		the input file are many list that have geneId at first col;
	-sample_file	sample name,sample file(sample name should not start be number)
	-sample_file_head	if sample file have header then give the parameter.else do not give this parameter; 
	-od	output dir
#############################################
USAGE
	print $usage;
	exit;
}

#####获取文件绝对路径####################################################################################################
sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	
	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;chomp $return;
	}
	else
	{
		warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}
