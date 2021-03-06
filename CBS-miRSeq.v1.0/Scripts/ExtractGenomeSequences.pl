
###################################################################################################################################
# ExtractGenomeSequences.pl
#
# This program uses samtools to extract sequences listed in gff format from a fasta formatted file. This is primarily made for 
# the miRspring pipeline but can be used for any number of uses.  
# For more information please visit http://mirspring.victorchang.edu.au.
#
#
# The following command will list the options available with this script:
#   perl ExtractGenomeSequences.pl
#
# Version 1.1
# Author: David Thomas Humphreys
# Copyright (C) 2015, Victor Chang Cardiac Research Institute
#
#
#     Version 1.1
#	- added ability to read gff3 mirbase formatted files
#	- fixed samtools detection bug
#	- Ensure miRbase defined stem loop sequence is upper case.
#
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    see <http://www.gnu.org/licenses/> for GNU license.

####################################################################################################################################


#!/usr/bin/perl
use strict;

my $DebugMode;
####------SETUP A------------------------------------------------------------------------------------------------------------------------- 
#### Before running the software it is recommended to set the following SIX global variables. 

## The three variables below should contain the DIRECTORY of the miRbase files. 
my $Precursor_GFF_Dir = " SET THE DIRECTORY WHERE mirbase gff files are located";
my $Precursor_Seq_Dir = " SET THE DIRECTORY WHERE mirbase precursor files (hairpin.fa) are located";
my $Mature_Seq_Dir = "SET THE DIRECTORY WHERE mirbase mature files (mature.fa) are located";
my $Samtools_Directory = '';	# ONLY Set this if the samtools directory if not set in the global path variable.

## This states the default mapping strategy you use to map data sets, either to a genome or a local version of miRbase. 
my $ReferenceFormat = 0;	# 0 = genome	1 = mirbase

## This states the default flanking sequence you provide in making a miRspring document. 
## If you are using unmodified miRbase files set this to 0 
my $Flank = 35;			# Flanking sequence on each side of precursor
####-------------------------------------------------------------------------------------------------------------------------------

my $parameters = {};

sub load_default_reference_files()
{
####------SETUP B------------------------------------------------------------------------------------------------------------------------- 
#### Before running the software it is recommended to set the default miRBase filenames


	$parameters->{gff} = "$Precursor_GFF_Dir/ ENTER FILENAME HERE, if using miRBase files use: $parameters->{species}.gff2", if (! defined $parameters->{gff});
	$parameters->{precursor_seq_file} = "$Precursor_Seq_Dir/ ENTER FILENAME HERE, if using miRBase files use: hairpin.fa", if (! defined $parameters->{precursor_seq_file});
####-------------------------------------------------------------------------------------------------------------------------------

	$parameters->{flank} = $Flank, if (! defined $parameters->{flank});
	$parameters->{reference_format} = $ReferenceFormat, if (! defined $parameters->{reference_format});
       $parameters->{colour_mm} = 6, if (! defined $parameters->{colour_mm});
}


my $miRspring;		# Global hash to save final output
my $ChromStats;         	# Global hash to save chromosome numbers
my $Flanking_Seq;

my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my  @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = gmtime();
my $year = 1900 + $yearOffset;
my $GMTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
print "Time stamp: $GMTime\n";


if ($DebugMode ne 'ON')
# Test to see if samtools is installed correctly
{	my $Samtools_Version = '';
	my $cmd = "$Samtools_Directory"."samtools";
	my @SamtoolTest = `$cmd 2>&1`;		# The "2>&1" redirects STDERR to STDOUT. 

	foreach (@SamtoolTest)
	{	my $Temp;
		if ($_ =~ m/^(Version: .*?)$/)
		{	$Samtools_Version = $1;	}
	}
	if ($Samtools_Version eq '')
	{ 	print "\nError with the following samtools path or installation: $cmd\nPlease check\n@SamtoolTest\n\n";	
		exit(0);
	}
	else {	print "\nDetected samtools installation, $Samtools_Version\n\n";		}
}	


# Test to see if enough parameters have been passed to script.
sub usage {
	print "\nUsage: $0 \n\n\t ";

    	print "REQUIRED \n\t";
    	print "-gff <Gene features file (gff), default is defined in script> \n\t";
	print "-fa <reference fasta file>\n\t";
    	print "-s <three letter species code>\n\t";
    #	print "-out <output file name with full path> \n\n\t\n";
    	print "Optional\n";
    	print "-flank <extra flanking sequence, default is 0>\n\n\n";

    	exit(1);
}
if(scalar(@ARGV) == 0){
    usage();
}

# Parse the Command Line
&parse_command_line($parameters, @ARGV);

sub parse_command_line {

    my($parameters, @ARGV) = @_;

    my $next_arg;

    while(scalar @ARGV > 0)
    {
        $next_arg = shift(@ARGV);
        if($next_arg eq "-s"){ $parameters->{species} = shift(@ARGV); }
        elsif($next_arg eq "-gff"){ $parameters->{gff} = shift(@ARGV); }
        elsif($next_arg eq "-fa"){ $parameters->{fasta} = shift(@ARGV); }
    #    elsif($next_arg eq "-out"){ $parameters->{output} = shift(@ARGV); }
	elsif($next_arg eq "-flank"){ $parameters->{flank} = shift(@ARGV); }

        else { print "Invalid argument: $next_arg"; usage(); }
    }

	# Check all essential options have been entered.
    	my $Error_Log = '';
    	$Error_Log .= "Species undefined (-s)\n\t", if (! defined $parameters->{species});
   # 	$Error_Log .= "No output file defined (-out)\n\t", if (! defined $parameters->{output});
    	$Error_Log .= "No input reference fasta file defined (-fa)\n\t", if (! defined $parameters->{fasta});
     	$Error_Log .= "mirbase gff file not defined (-gff)\n\t", if (! defined $parameters->{gff});
    	if ($Error_Log ne '')
    	{	print "\nThe following essential option(s) have not been defined:\n\t$Error_Log\n\n";
		usage();
	}

	print "\n\n**************    Parameters configures for fasta script ********************************\n";
	while (my ($key, $value) = each %$parameters)
	{	print "$key\t$value\n";	}
	print "********************************************************************************************\n";

}

# Read header of BAM file to obtain reference layout and other information
my $HeaderFields;	# Will contain string of relevant header information
my $RefFields;	# A hash of all reference contig names


# open mirbase precursor coordinate file
open PRECURSOR_GFF, $parameters->{gff} or die "Cannot open precursor file $parameters->{gff}";
#### Example of what this file should look like
# 1	.	miRNA	20669091	20669163	.	+	.	ACC="MI0000249"; ID="mmu-mir-206";
#### gff3 similar to gff2 but also has stem loop entries. 3rd column identified the entry:
# chr1	.	miRNA_primary_transcript	12425986	12426106	.	+	.	ID=MI0021869;Alias=MI0021869;Name=mmu-mir-6341
# chr1	.	miRNA	12426016	12426038	.	+	.	ID=MIMAT0025084;Alias=MIMAT0025084;Name=mmu-miR-6341;Derives_from=MI0021869

# Identify if file is gff2 or gff3 from file extension:
my $GFF_file_format = 2;	# Default is gff2 file type.
$GFF_file_format = 3, if ($parameters->{gff} =~ m/\.gff3$/);


my (@miR_Name, @miR_Start, @miR_End, @miR_Strand); 	# Storage for miR coordinates
my ($i, $Index) = (0,0,0);
my $Unique_miR_ID;					# Hash that stores miR name index ID
$RefFields->{"Max"} = 1;
while (<PRECURSOR_GFF>)
{	chomp();	# Remove return character
	if ($_ =~ /^chr(.*?)$/)
	{	$_ = $1;	}
	my @Record = split("\t",$_);
	next, if (scalar(@Record) < 7);	# gff records should have 8 fields\
	next, if (($GFF_file_format == 3) && ($Record[2] eq 'miRNA'));	# Only want stem loop entries, not mature miRNA entries.


	if ($RefFields->{$Record[0]} == '')
	{	# New chromosome entry.
		$RefFields->{$Record[0]} = $RefFields->{"Max"};			# Save the index
		$RefFields->{$RefFields->{"Max"}." label"} = $Record[0];	# Save the name/label
		$RefFields->{"Max"}++;

		$i = 0;
	}			

	my $Stem_Loop_ID='';
	if ($GFF_file_format == 2)
	{	$Record[8] =~ m/\"(\w\w\w-\w\w\w-.*?)\"/;	
	#	$Record[8] = m/($parameters->{species}-.*?)\"/;		# This is potentially more accurate but may be (rare) cases where it could hinder.
		$Stem_Loop_ID = $1;
	}
	elsif ($GFF_file_format == 3)
	{	$Record[8] =~ m/.*?\=(\w\w\w-\w\w\w-.*?)$/;	
		$Stem_Loop_ID = $1;
	}

	$miR_Name[$RefFields->{$Record[0]}][$i] = $Stem_Loop_ID;
	$Unique_miR_ID->{$Stem_Loop_ID} = $Index;
	$miR_Start[$RefFields->{$Record[0]}][$i]= $Record[3];
	$miR_End[$RefFields->{$Record[0]}][$i]= $Record[4];
	$miR_Strand[$RefFields->{$Record[0]}][$i] = $Record[6];
#	print "\nLoaded $Stem_Loop_ID at position $i.. chromosome $Record[0]";
	$i++;
	$Index++;
}
close PRECURSOR_GFF;



# samtools faix reference SAMFILE:chr x-y
for(my $i = 1; $i < $RefFields->{"Max"}; $i++)
{       my $j = 0;
	my $Chrom = $i;

	while($miR_Start[$i][$j] > 0)
	{	
						
		my $Start = $miR_Start[$i][$j] - $parameters->{flank};
		my $End = $miR_End[$i][$j] + $parameters->{flank};

		my $Command = "$Samtools_Directory"."samtools faidx $parameters->{fasta} chr".$RefFields->{"$i label"}.":$Start-$End";
		my @FastaRecord = `$Command`;
		my $HairpinSeq = FormatFasta($i, $Start, $End, $j, @FastaRecord);

		if ($miR_Strand[$i][$j] eq '-')
		{	$HairpinSeq = reverse($HairpinSeq);
			$HairpinSeq =~ tr/ATCG/TAGC/;
		}

		if ($parameters->{flank} > 0)
		{	my $TotalLength = length($HairpinSeq);
			my $FirstBit = lc(substr($HairpinSeq,0,$parameters->{flank}));
			my $Precursor = uc(substr($HairpinSeq,$parameters->{flank},$TotalLength-$parameters->{flank}-$parameters->{flank}));
			my $LastBit = lc(substr($HairpinSeq,$TotalLength-$parameters->{flank},$parameters->{flank}));
			#print "\n$HairpinSeq\n$FirstBit$Precursor$LastBit";
			$HairpinSeq= "$FirstBit$Precursor$LastBit";
		}

		print "\n>$miR_Name[$i][$j],chr".$RefFields->{"$i label"}.":$miR_Start[$i][$j]-$miR_End[$i][$j]\n$HairpinSeq";
		$j++;
	}
}



print "\n";

exit(0);









sub FormatFasta
{       my $Chr = shift;
	my $Start = shift;
	my $End = shift;
	my $Index = shift;
	my @Fasta_Input = @_;

	my $FinalOutput; 

	foreach(@Fasta_Input)
	{	chomp($_);
		if ($_ =~ m/^>.*?/)
		{	# Header, do not need to save this
			#$FinalOutput = $miR_Name[$Chr][$Index];
		}
		else
		{	# Sequence
			$FinalOutput .= $_;

	        }

	}
	return $FinalOutput;
}



