#!/bin/env perl
## extract region of an alignment based on position, assume all sequences are of same aligned length
## Original script written by Shannon Soucy 2013-10-08
## Modified by Michael C. Nelson 2013-11-04

# INPUT<- File name of aligned fasta file of sequences
# OUTPUT-> INPUT_resistance.txt gives you a file detailing the codon of the 83 a.a. and it's resistance signficance.

$seq=$ARGV[0];

# these variables set the positions in the alignment that the user wants keep
$start=84;
$stop=87;
$length = $stop-$start;

#grabbing all files with extension '.fasta'
	@name=split(/\./,$seq);
	#split the filename to get out the prefix to use on the out file
	open (IN, "< $seq") or die "Cannot open file!!\n";
	#open in file
	open (OUT, "> $name[0]_resistance.txt") or die "Problem opening output file!!";
	#open out file
	while (defined($line=<IN>)){
	#loop through sequences
		#if the line begins with a carrot print it to the outfile
		if ($line=~m/>/){
			( $put = $line ) =~ s/>(denovo[0-9]*) .*\n/$1/g;
			print OUT "$put\t";
						}
		#else the line doesn't have a carrot so extract the substring and print that to the out file 
		else{
			$string=substr($line,"$start","$length");
            if ($string=~m/AT[ATC]/){ #REGEX for Isoleucine
				print OUT "CiproR\t$string\n";
			}
			else{
				print OUT "CiproS\t$string\n";
			}
		}
								}
