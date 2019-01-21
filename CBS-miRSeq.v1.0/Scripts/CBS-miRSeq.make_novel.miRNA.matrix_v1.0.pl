#!/usr/bin/perl
use strict;
use warnings;
use File::DosGlob qw(bsd_glob);
use File::Find;

############################ Create novel miRNA expr  matrix #######################
## Description:
## Using reference and overlapped coordinates produce novel miRNA genome coordinates 
## expression matrix for DE analysis across sample
## Author: Kesharwani Rupesh Kumar
## Date: 12 feb 2015
## Copyright (c) 2019 Kesharwani RK
## Parameters used:
## <files.dir> - input directory (*.refer and *.common)
####################################################################################

# # my $date = localtime();
# # print "Local Date and Time $date\n";

if ($#ARGV < 0) {
    die "\nUSAGE : ./CBS-miRSeq.make_novel.miRNA.matrix_v1.0.pl <files.dir> \n";
    exit 1;
}

# reads inputs dir or die
my ($file_dir) = @ARGV;
chdir $file_dir or die "Cannot find $file_dir: $!\n";

#read 'reference file' into a hash:
my $ref = glob ("*.refer");
my %ref;
open( my $ref_fh, "<", $ref );
while (<$ref_fh>) {
    my ( $first, $second, $third ) = split;

    #turn the first three fields into space delimited key.
    $ref{"$first\t$second\t$third"} = ();
}

# open each of the files.
my @files = glob ("*.common");
foreach my $input (@files) {
    open( my $input_fh, "<", $input );
    my %current;
    while (<$input_fh>) {

        #line by line, extract 'first 3 fields' to use as a key.
        #then 'value' which we store.
        my ( $first, $second, $third, $value ) = split;
        $current{"$first\t$second\t$third"} = $value;
    }

    #refer to 'reference file' and insert matching value or zero into
    #the array.
    foreach my $key ( keys %ref ) {
        push( @{ $ref{$key} }, $current{$key} ? $current{$key} : 0 );
    }
}

foreach my $key ( keys %ref ) {
    print join( "\t", $key, @{ $ref{$key} } );
	print "\n";
}