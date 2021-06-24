#!/usr/local/perl/5.24.0/x86_64/gcc-4.9.3/nonet/bin/perl
use strict;
use warnings;
use Getopt::Long;

############
##Parameter
############
my $opt_help;
my $opt_out_dir = "";

#########
my $input = "";
my $input1 = "";
my $input2 = "";
my $rsid = "";
my $genotype = "";
my $samples = "";
my $output = "";

my $AF = 0.01;
my $an = 1676;
#################
##software
######################
# the location of tabix
my $tabix = "bin/tabix";

###############

GetOptions("help|h"  => \$opt_help,
	"out_dir|output_dir=s" => \$opt_out_dir,
	#####
	"input|input=s" => \$input,
	"input1|input1=s" => \$input1,
	"input2|input2=s" => \$input2,
	"rsid|rsid=s" => \$rsid,
	"genotype|genotype=s" => \$genotype,
	"samples|samples=s" => \$samples,
	
	"AF|AF=f" => \$AF,
	);
	
################
## get samples	
##########################
my @gene =split(":",$input); #"$gene:$chr:$start:$end"
#mkdir ("$opt_out_dir/$gene[0]",0755) unless (-d "$opt_out_dir/$gene[0]");

#############
my @samples;
my $i=0;
open(SAMPLE, $samples) or die "there is not exist $samples \n";
while(my $line=<SAMPLE>)
{	
	chomp($line);
	next if $line =~ /^#/;
	
	my @f=split(/\s+/,$line);	
	$samples[$i]=$f[0];
	$i++;
}
close SAMPLE;

#############
## get y_all

my @header_index;

open(EXPRESS, "$tabix $input1 chr$gene[1]:$gene[2]-$gene[3] -h |") or die "there is not exist $input1 \n";

open OUT, ">$opt_out_dir/$gene[0].$gene[4].y_all" or die "there is not exist $opt_out_dir/$gene[0].$gene[4].y_all \n";

while(my $line=<EXPRESS>)
{	
	chomp($line);
	next if $line =~ /^##/;

	my @f=split(/\s+/,$line);
	
	if($line =~ /^#/){
		for(my $k2 = 0; $k2<@samples; $k2++)
		{
			for(my $k=0;$k<@f;$k++)
			{
				if($samples[$k2] eq $f[$k]){				
					push @header_index, $k;			
				}
			}		
		}
		next;
	}
	

	my $chr = $f[0];
		$chr =~ s/chr//;
	my $start = $f[1];
	my $end = $f[2];
	
	my @gene2 = split("\\.",$f[3]);

	if ($gene2[0] eq $gene[0])
	{	
		print OUT join("\t","Expression",@f[@header_index])."\n";
	
	}
}	

close EXPRESS;
close OUT;		
	
####################
@header_index = ();

open(EXPRESS, "$tabix $input2 chr$gene[1]:$gene[2]-$gene[3] -h |") or die "there is not exist $input2 \n";

open OUT, ">>$opt_out_dir/$gene[0].$gene[4].y_all" or die "there is not exist $opt_out_dir/$gene[0].$gene[4].y_all \n";

my $k3 =1;
while(my $line=<EXPRESS>)
{	
	chomp($line);
	next if $line =~ /^##/;

	my @f=split(/\s+/,$line);
	
	if($line =~ /^#/){
		for(my $k2 = 0; $k2<@samples; $k2++)
		{
			for(my $k=0;$k<@f;$k++)
			{
				if($samples[$k2] eq $f[$k]){				
					push @header_index, $k;			
				}
			}		
		}
		next;
	}
	
	my $chr = $f[0];
		$chr =~ s/chr//;
	my $start = $f[1];
	my $end = $f[2];
	my @gene2 = split(":",$f[3]);
	my @gene3 = split("\\.",$gene2[4]);

	if ($gene3[0] eq $gene[0])
	{	
		print OUT join("\t","Splice".$k3,@f[@header_index])."\n";
		
		$k3++;
	}
}	

close EXPRESS;
close OUT;

#######################
## get rs id

my %hashRSid;

open(INPUT, $rsid) or die "there is not exist $rsid \n";
while(my $line=<INPUT>)
{	
	chomp($line);
	next if $line =~ /^#/;
	
	my @f=split(/\s+/,$line);
	
	my $chr = $f[0];
			$chr =~ s/chr//;
	my $start = $f[1];	
	my $id = $f[2];
	
	$hashRSid{"$chr:$start"}=$id;
}
close INPUT;


##########################


@header_index = ();

open(GENOTYPE, "$tabix $genotype chr$gene[1]:$gene[2]-$gene[3] -h |") or die "there is not exist $genotype \n";
open OUT, ">$opt_out_dir/$gene[0].$gene[4].x_all" or die "there is not exist $opt_out_dir/$gene[0].$gene[4].x_all \n";

while(my $line=<GENOTYPE>)
{	
	chomp($line);
	next if $line =~ /^##/;
	
	my @f=split(/\s+/,$line);
	
	if($line =~ /^#/){
		for(my $k2 = 0; $k2<@samples; $k2++)
		{
			for(my $k=0;$k<@f;$k++)
			{				
				if($samples[$k2] eq $f[$k]){				
					push @header_index, $k;					
				}
			}		
		}
		next;	
	}
############
##info
	my $chr = $f[0];
		$chr =~ s/chr//;

	my $alt = 0;
	my @info=split(";",$f[7]);
	
	my @an=split("\=",$info[0]);
	my @alt=split("\=",$info[1]);	
		
	
	if($an[1]<$an)
	{next;}
	
	if($alt[1]<$AF  || $alt[1]>1-$AF)
	{next;}
	
##########	
	my @f2;
	for(my $k=9;$k<@f;$k++)
	{
		my @a=split(":",$f[$k]);
		
		if( $a[0] eq "0|0" ||  $a[0] eq "0/0")
		{
			$f2[$k] = 0;
		}elsif ($a[0] eq "0|1" || $a[0] eq "1|0" || $a[0] eq "0/1" )
		{
			$f2[$k] = 1;
		}elsif ($a[0] eq "1|1" || $a[0] eq "1/1" )
		{
			$f2[$k] = 2;
		}else{ $f2[$k] = "NA"}
	
	}
	
	
	if(exists $hashRSid{"$chr:$f[1]"})
	{
		print OUT join("\t","$chr",$hashRSid{"$chr:$f[1]"},0,$f[1],$f[3],$f[4],$alt[1],@f2[@header_index])."\n";
	}
	
	
			
}

close GENOTYPE;	
close OUT;




