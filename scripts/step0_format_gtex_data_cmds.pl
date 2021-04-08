#!/usr/local/perl/5.24.0/x86_64/gcc-4.9.3/nonet/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Find;
use File::Path;
#use Thread::Semaphore;
use threads;
use FindBin;
use File::Basename;
use Cwd 'abs_path';



#
my $Rscript = "/accre/arch/easybuild/software/MPI/GCC/8.2.0/OpenMPI/3.1.4/R/3.6.0/bin/Rscript";

############

my $script = "$FindBin::RealBin/lib"; ##
my $user = `echo \$USER`;
my $pwd = `pwd`;
my $jobname_pre = $$;

############
##Parameter
############
my $opt_help;
my $get_config = "$FindBin::RealBin/calling.config";
my $opt_out_dir = "";
my $where = ""; #local or qsub or sbatch
my $email = "";
my $nthreads = 20;
##
my $sleep = 1 ;

#########
my $input1 = "";
my $input2 = "";
my $genotype = "";

#############
my $genelist;
my $generateMatrix;
my $geneID = "../ref_data/gene.500k.id";
my $rsid = "../ref_data/hg38.vcf";

#############

###################

GetOptions("help|h"  => \$opt_help,
	"get_config|get_config=s" => \$get_config,
	"out_dir|output_dir=s" => \$opt_out_dir,
	"where|where=s" => \$where,
	"email|email=s" => \$email,
	"nthreads|nthreads=f" => \$nthreads,


	#####
	"input1|input1=s" => \$input1,
	"input2|input2=s" => \$input2,
	"rsid|rsid=s" => \$rsid,
	"genotype|genotype=s" => \$genotype,
	
####

	"sleep|sleep=f" => \$sleep,
	"jobname_pre|jobname_pre=f" => \$jobname_pre,
	
######
	"genelist|genelist" => \$genelist,	
	"generateMatrix|generateMatrix" => \$generateMatrix,	

	);
	
#########################################
## how to use
#######################################
&usage() if $opt_help;

die "-get_config is not provided!\n" if !$get_config;

#die "-output_dir is not provided!\n" if !$opt_out_dir;
die "-where is not provided. -local, -qsub or -sbatch!\n" if ($where ne "local" && $where ne "qsub" && $where ne "sbatch");
die "-email is not provided when running in $where !\n" if ($where eq "qsub" || $where eq "sbatch") && !$email;
die "-nthreads is not provided when running in local !\n" if $where eq "local" && !$nthreads;	
	
#########################################
## show options
#######################################

&options();

############################
# sleep time in different where.
######################
if (($where eq "qsub" || $where eq "sbatch") && $sleep == 1200)
{
	$sleep = 1200;
}
elsif($where eq "local" && $sleep <= 1)
{
	$sleep = 1;
}

#########################################
## config files and check all softwares and databases
#######################################

my %hash_config = &get_config($get_config);

#####################################################
##Make dir
#####################################################

chop($opt_out_dir) if $opt_out_dir =~ /\/$/;

eval{mkpath($opt_out_dir,0,0755)};
if($@)
{
   warn("Make path [$opt_out_dir] failed:\n$@");
}

#my $data_dir = "$opt_out_dir/rawdata"; #fastq files

my $log_dir = "$opt_out_dir/logs"; #log files

mkdir ($log_dir,0755) unless (-d $log_dir);


######################
my $stats_log = "$opt_out_dir/$jobname_pre.stats.log"; #if stop because somethingwrong, delete the step* running to rerun
my $process = "$opt_out_dir/process.log";
my $error = "$opt_out_dir/$jobname_pre.error.log";
my $commands = "$opt_out_dir/$jobname_pre.commands.log";


#######################################################################

################
##
#################

################
##
#################

###########
my $input;

if ($genelist)
{
	$input = &genelist($nthreads, $opt_out_dir, \%hash_config,$input1,$input2,$genotype);	
}

if ($generateMatrix)
{
	$input = &generateMatrix($nthreads, $opt_out_dir, \%hash_config, $input,$input1,$input2,$genotype);	
}


##################
#################
sub generateMatrix()
{
	my ($nthreads, $opt_out_dir, $hash_config, $input, $input1,$input2, $genotype) = @_;

	if ($nthreads>30 && !$nthreads) 
	{
		$nthreads = 30;
	}

	my $nthreads_perjob = 1;
	my $job_run = int ($nthreads/$nthreads_perjob);

	#####
	open(TO,">>$process") or die();
	print TO &nowtime."\t generateMatrix starts running\n";

	print "\n###################\n generateMatrix starts running; -->\n";
	print "\n***********check config files\n";
	print "***********\n";
	
	my @commandline;
	
	my $i = 0;

	open(OUT, ">$input/jobs.txt") or die "there is not exist $input/jobs.txt \n";
	
	open(GENE, "$input/geneid.output") or die "there is not exist $input/geneid.output \n";
	while(my $line=<GENE>)
	{	
		chomp($line);
		next if $line =~ /^#/;
		
		my @f=split(/\s+/,$line);

		my $gene = $f[0];
		my $chr = $f[1];
			$chr =~ s/chr//;
		my $start = $f[2];
		my $end = $f[3];
		
		my @gene2 = split("\\.",$f[7]);
			my $gene2 = $gene2[0];
		my $start2 = $f[4];
		my $end2 = $f[5];
		
		$commandline[$i] = qq/perl $script\/step1_genMatrix_gtex.pl -input "$gene:$chr:$start:$end:$gene2:$start2:$end2" -input1 $input1 -input2 $input2 -genotype $genotype -out_dir $opt_out_dir -samples $input\/samples.output -rsid $rsid /;
		
		#$commandline[$i] .= qq/ && $Rscript $script\/omics.R --x $opt_out_dir\/$gene.$gene2.x_all --y $opt_out_dir\/$gene.$gene2.y_all   --output $opt_out_dir\/$gene.$gene2 /;
		
		
		print OUT $commandline[$i]."\n";

		#print $i."\n";

		$i++;

		#if($i>200){last;}
	}

	close GENE;
	close OUT;
	close TO;	
	
	my $whichSubroutine = \&stand_command;
	#&multifork($job_run,$whichSubroutine,@commandline);
}



sub genelist()
{
	my ($nthreads, $opt_out_dir, $hash_config, $input1,$input2,$genotype) = @_;

	if ($nthreads>30 && !$nthreads) 
	{
		$nthreads = 30;
	}

	my $nthreads_perjob = 1;
	my $job_run = int ($nthreads/$nthreads_perjob);

	#####
	open(TO,">>$process") or die();
	print TO &nowtime."\t genelist starts running\n";

	print "\n###################\n genelist starts running; -->\n";
	print "\n***********check config files\n";
	print "***********\n";

###########
	my %hashExpr;
	my %hashSplicng;	
	my %hashTotal;
	
	my %hashGene;
	
	my %hashSamples;
	
	##############
	open(GENE, $geneID) or die "there is not exist $geneID \n";
	while(my $line=<GENE>)
	{	
		chomp($line);
		next if $line =~ /^#/;
		
		my @f=split(/\s+/,$line);
					
		my $gene = $f[0];

		$hashGene{$gene}=join("\t",@f[1..@f-1]);
	}

	close GENE;
		
	##############
	open(EXPRESS,"gunzip -c $input1 | ") or die "there is not exist $input1 \n";
	while(my $line=<EXPRESS>)
	{	
		chomp($line);
		next if $line =~ /^##/;
		my @f=split(/\s+/,$line);
		
		if($line =~ /^#/){
			for(my $k=4;$k<@f;$k++)
			{
				$hashSamples{$f[$k]}++;				
			}
			next;
		}

		my $chr = $f[0];
		my $start = $f[1];
		my $end = $f[2];
		my @gene = split("\\.",$f[3]);

		$hashExpr{$gene[0]}++;
		$hashTotal{$gene[0]}++;
	}

	close EXPRESS;
	################
	
	open(EXPRESS,"gunzip -c $input2 | ") or die "there is not exist $input2 \n";
	while(my $line=<EXPRESS>)
	{	
		chomp($line);
		next if $line =~ /^##/;
		my @f=split(/\s+/,$line);
		
		if($line =~ /^#/){
			for(my $k=4;$k<@f;$k++)
			{
				$hashSamples{$f[$k]}++;				
			}
			next;
		}

		my $chr = $f[0];
		my $start = $f[1];
		my $end = $f[2];
		my @gene = split(":",$f[3]);
		my @gene2 = split("\\.",$gene[4]);

		$hashSplicng{$gene2[0]}++;
		$hashTotal{$gene2[0]}++;
	}
	close EXPRESS;
	
	################
	open(GENOTYPE,"gunzip -c $genotype | ") or die "there is not exist $genotype \n";
	while(my $line=<GENOTYPE>)
	{	
		chomp($line);
		next if $line =~ /^##/;
		my @f=split(/\s+/,$line);
		
		if($line =~ /^#/){
			for(my $k=9;$k<@f;$k++)
			{
				$hashSamples{$f[$k]}++;
			}
			last;
		}
	}
	close GENOTYPE;	
	
	
	open OUT, ">$opt_out_dir/geneid.output" or die "something wrong in include normalization wig files"; 

	foreach my $key (keys %hashTotal)
	{
		my $a = 0;
		my $b = 0;
		if (exists $hashExpr{$key}) {$a=$hashExpr{$key};}
		if (exists $hashSplicng{$key}) {$b=$hashSplicng{$key};}
		
		if(exists $hashGene{$key})
		{
			print OUT join("\t",$key, $hashGene{$key}, $a,$b )."\n";
		}
		# else {print join("\t",$key, 0, $a,$b )."\n";}
		
		
	}
	close OUT;
	
	open OUT2, ">$opt_out_dir/samples.output" or die "something wrong in include normalization wig files"; 
	
	foreach my $key (keys %hashSamples)
	{		
		if($hashSamples{$key}==3)
		{
			print OUT2 $key."\n";
		}		
	}
		
	close OUT2;

	close TO;
	
	return $opt_out_dir;
}



##############################################################################################################
### basic subprogress
##############################################################################################################

sub nowtime()
{
	my @time=(localtime)[5,4,3,2,1,0];   
	$time[0]+=1900;   
	$time[1]+=1;   
	my $nowtime=sprintf("%04u/%02u/%02u\t%02u:%02u:%02u",@time);   
	return $nowtime;
}

sub cluster_bsub()
{
	my ($command_line,$bsub_options) = @_;

	my $jobname_pre = @$bsub_options[2];
	my $signal_total = @$command_line;
	my $opt_working_dir = `pwd`; chomp($opt_working_dir);


	my $cluster_submit;

	if ($where eq "qsub") 
	{
		$cluster_submit = $hash_config{"bsub"};
		open(CLUSTER, "|$cluster_submit -pd `pwd` -pm @$bsub_options[0] -ph @$bsub_options[1] -jobn @$bsub_options[2] -email @$bsub_options[3] -e");
	}
	elsif($where eq "sbatch")
	{
		$cluster_submit = $hash_config{"ssub_array2"};
		open(CLUSTER, "|$cluster_submit -pd `pwd` -pm @$bsub_options[0] -ph @$bsub_options[1] -cpt @$bsub_options[4] -jobn @$bsub_options[2] -email @$bsub_options[3] -e");
	}

	print CLUSTER join("\n",@$command_line)."\n";

	close CLUSTER;

	while()
	{
		my $signal = 0;	

		if(-e "$opt_working_dir/sbatch.$jobname_pre.jobs.signal")
		{
			open IN, "<$opt_working_dir/sbatch.$jobname_pre.jobs.signal";

			while(my $line=<IN>)
			{
				chomp($line);

				next if $line !~ /\S+/;
				next if $line =~ /^#/;

				my @f = split(/\./,$line);

				if ($f[0] eq $jobname_pre) 
				{
					$signal++;
				}		 
			}
			close IN;

			if ($signal == $signal_total) 
			{
				last;
			}
		}

		sleep($sleep);
	}

	my $out_file = ` less $opt_working_dir/tmp/$jobname_pre*.out | grep error `;

	 if ($out_file) 
	 {
	 	die "error: something wrong\n";
	 }
}


sub get_config()
{
	(my $get_config) =@_;

	my %hash_config;

	open(FROM,"$get_config") or die "there is not exist $get_config\n";
	while(my $line=<FROM>)
	{	
		chomp($line);

		next if $line !~ /\S+/ || $line =~ /^#/;
		my @tmp=split(/\s+/,$line);

		if($tmp[0]){$hash_config{$tmp[0]} = $tmp[1];}
	}
	close FROM;

	return %hash_config;
}

sub multifork
{
	(my $num, my $whichSubroutine, my @array) = @_;

	my $j=0;
	my $thread;

	my $array_num = @array;

	while()
	{
		last if($j>=$array_num);
		while(scalar(threads->list())<$num && $j < $array_num)  
		{		
			threads->new(
						$whichSubroutine,$array[$j]
				);
			$j++;
		}

		foreach $thread(threads->list(threads::all))
		{
			if($thread->is_joinable())
			{
				$thread->join();
				#print scalar(threads->list()),"\t$j\t",localtime(time),"\n";
			}
			if (my $err = $thread->error()) 
			{
				warn("Thread error: $err\n");
			}

		}
		sleep($sleep);		
	}

	foreach $thread(threads->list(threads::all))
	{
		$thread->join();

		if (my $err = $thread->error()) 
		{
			warn("Thread error: $err\n");
		}
		#print scalar(threads->list()),"\t$j\t",localtime(time),"\n";       
	}
}

sub stand_command()
{
	(my $cmd) = @_;

	my $run_status = 0;

	#print join("\t",$j,@array-1+1,$cmd)."\n";
	
	open(TO,">>$commands") or die();
	print TO "$cmd\n";
	close TO;

	$run_status = system("$cmd");

	if ($run_status != 0)
	{
		&sendmail("$cmd");

		open(TO,">>$error") or die();
		print TO &nowtime."\n************\n$cmd\n is wrong\n************\n";
		close TO;

		die "\n************\n$cmd\n is wrong\n************\n";
	}
}


sub sendmail()
{
	(my $command_line) = @_;

	if ($email)
	{
		my $from="vmpsched\@vmpsched.vampire";
		my $to="$email";
		my $subject="An error";

		my $sendmailpath="/usr/sbin/sendmail";

		my $message = "An error has occurred processing your job, see below.\n$command_line\n\nfrom cgg lab\n";

		open (SENDMAIL, "| $sendmailpath -t") or die "Cannot open $sendmailpath: $!";

		print SENDMAIL "Subject: $subject\n";
		print SENDMAIL "From: $from\n";
		print SENDMAIL "To: $to\n\n";

		print SENDMAIL "$message";

		close (SENDMAIL);
	}
}

sub checkfile()
{
	my @file = @_;

	for (my $i = 0; $i < @file; $i++) 
	{		
		print STDERR "\tthe $file[$i] file exists -->";
		if (-e "$file[$i]") 
		{
			print STDERR "\tOK\n";
		}
		else 
		{
			die "\tERROR\n################\n$file[$i] file does not exist\nPlease check the original files\n#################\n";
		}
	}

}



sub options()
{
		
	print "\n*******************************************************\n";	
	print "*************************\n";

	print "\n*** The script was complied on Jan 9 2020 15:21:20 ***\n";
	print "\tUsage: [OPTIONS]\n";

	print "Options\n";
	
	print "\n*************************\n";
	print "*******************************************************\n";
	
}


sub usage()
{
		
	print "\n*******************************************************\n";	
	print "*************************\n";

	print "\n*** The script was complied on Jan 9 2020 15:21:20 ***\n";
	print "\tUsage: [OPTIONS]\n";

	

	print "\n*************************\n";
	######

	
	
	##################################
	#
    print "\n step0: perl step0_format_gtex_data_cmds.pl -genelist -generateMatrix  -input1 Brain_Cortex.v8.normalized_expression.bed.gz  -input2 Brain_Cortex.v8.leafcutter_phenotypes.bed.gz   -genotype GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz -out_dir matrix -rsid hg38.vcf -where local -nthreads 8 \n"
	
	#
    print "\n less matrix/jobs.txt |cut -f 1 -d '&' |bash \n"
		

	
	####
	print "\n*************************\n";
	print "*******************************************************\n\n";	

	exit(0);
}





	
	
	
	
	
	
	
	
	
	
	
	
	
