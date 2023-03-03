#!/usr/bin/perl
use Bio::SeqIO;
use Getopt::Long;

# use strict;
# use warnings;
# use Text::CSV;

=head1 DESCRIPTION

	Computes the number of N in a fasta sequence. List of sequence id and number of N is given in a tab separated output file.

=head1 USAGE

	perl Count_N_in_seq.pl -i file.fasta [-o output_file]

=head1 OPTIONS

	-i or -input:	fasta file	(REQUIRED)
	-o or -out:	Output File	(DEFAULT: [input_name].N_number.csv)

	-h or -help:	This Documentation

=head1 AUTHOR

	Damien Richard, 2023

=cut

my $fasta_file; my $out; my $hl;

GetOptions(
	"i|input=s" => \$fasta_file,
	"o|out=s" => \$out,
	"h|help" => \$hl,
);

die `pod2text $0` unless $fasta_file;
die `pod2text $0` if $hl;

if(!$out){$out = $fasta_file . ".N_number.csv";}
# End of user interface code, beginning of program: 

my $fasta_in ; my $nb_n;
open (OUT, ">", $out) ;


$fasta_in  = Bio::SeqIO->new( -file => $fasta_file, -format => 'fasta' );

print "Processing fasta file ...\n";
print $fasta_file . "\n" ;

while ( my $seq = $fasta_in->next_seq() ) {
my $id  = $seq->display_id();
my $seqstr = uc($seq->seq()); # actual sequence as a string
@matches = $seqstr =~ /N/g ; 
$nb_n = scalar @matches ;
print OUT $id . "\t" . $nb_n . "\n";

}

close OUT;





