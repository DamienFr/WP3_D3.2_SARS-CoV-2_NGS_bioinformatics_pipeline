#!/usr/bin/perl
use Bio::SeqIO;
use Getopt::Long;

# use strict;
# use warnings;
# use Text::CSV;

=head1 DESCRIPTION

	Computes the pairwise number of SNP between a set of fasta sequences that need to be of identical length. List of all pairwise SNP numbers is outputed to a file while the mean number is given in STDOUT (the terminal).

=head1 USAGE

	perl Pairwise_SNP_distance_from_fasta -i in.fasta [-o output_file]

=head1 OPTIONS

	-i or -input:	fasta file	(REQUIRED)
	-i2 or -input2:	fasta file	(REQUIRED)
	-r or -restrict:	To compute only certain pairwise comparisons, provide a tab delimited file with comparisons to make as rows	(OPTIONAL)
	-o or -out:	Output File	(DEFAULT: [input_name].pairwise_SNP_number.csv)
	-o2 or -out2:	Summarized Output File	(DEFAULT: no summarized output)
	-m or -min:	Min size of each sequence of each pair in order to consider the pair (DEFAULT: no min size of sequence applied)

	-h or -help:	This Documentation

=head1 AUTHOR

	Damien Richard, 2020

=cut

my $fasta_file; my $out; my $out2; my $hl; my $restrict; my $fasta_file2; my $min_size;  

GetOptions(
	"i|input=s" => \$fasta_file,
	"i2|input2=s" => \$fasta_file2,
	"r|restrict=s" => \$restrict,
	"o|out=s" => \$out,
	"o2|out2=s" => \$out2,
	"m|min=s" => \$min_size,
	"h|help" => \$hl,
);

die `pod2text $0` unless $fasta_file;
die `pod2text $0` if $hl;

if(!$out){$out = $fasta_file . ".pairwise_SNP_number.csv";}
# End of user interface code, beginning of program: 

my $fasta_in ; my %h; my %done; my $nb_comparison; my $sum; my $seq_length; my $seq_number;
my $fasta_in2 ;  my %h2; my %done2; my $nb_comparison; my $sum2; my $seq_length2; my $seq_number2;

$fasta_in  = Bio::SeqIO->new( -file => $fasta_file,       -format => 'fasta' );
$fasta_in2  = Bio::SeqIO->new( -file => $fasta_file2,       -format => 'fasta' );


if( $restrict ){

 open (GB, "<", $restrict) or die "could not open $restrict for read in openfile subroutine\n";
 while (my $line = <GB>) {
  $line =~ s/\r?\n//g;
 my @fields = split /\t/,$line;
# print STDERR $line . "\n" ;
  $comp_todo{$line} = 1;

 }
 close GB ; 
}

print "Processing fasta file ...\n";
print $fasta_file . "\n" ;
print $fasta_file2 . "\n" ;

while ( my $seq = $fasta_in->next_seq() ) {
my $id  = $seq->display_id();
my $seqstr = $seq->seq(); # actual sequence as a string
my $current_seq_length = $seq->length();
$seq_number++;
if( $seq_number == 1 ){ $seq_length = $current_seq_length }else{ if($current_seq_length != $seq_length ){ die "All sequences are not the same length ! Seq $id is $current_seq_length bp long whereas first seq of the dataset was $seq_length bp"} }
$h{$id} = $seqstr;
}

while ( my $seq2 = $fasta_in2->next_seq() ) {
my $id2  = $seq2->display_id();
my $seqstr2 = $seq2->seq(); # actual sequence as a string
my $current_seq_length2 = $seq2->length();
$seq_number2++;
if( $seq_number2 == 1 ){ $seq_length2 = $current_seq_length2 }else{ if($current_seq_length2 != $seq_length2 ){ die "All sequences are not the same length ! Seq $id2 is $current_seq_length2 bp long whereas first seq of the dataset was $seq_length2 bp"} }
$h2{$id2} = $seqstr2;
}


# check noms dupliqués
# check length
open (OUT, ">", $out) ;

my %exclusion ;

FIRST: foreach $key1 ( keys %h){
SECOND: foreach $key2 ( keys %h2){
if( $key1 ne $key2 && !exists($hdone{$key1 . "|" . $key2})  && !exists($hdone{$key2 . "|" . $key1})  && !exists($exclusion{$key1}) && !exists($exclusion{$key2})  ){

if( ! $restrict || (exists( $comp_todo{$key1 . "\t" . $key2} ) || exists( $comp_todo{$key2 . "\t" . $key1} ))){

$hdone{$key1 . "|" . $key2} = "done";

my $seq1 = lc $h{$key1} ;
my $seq2 = lc $h2{$key2} ;

$seq1 =~ s/-/n/g;
$seq2 =~ s/-/n/g;


my %hseq1 ;
my %hseq2 ;
my $offset = 0;
my $result = index($seq1, "n", $offset);

while ($result != -1) {
$hseq1{$result}++; 
$offset = $result + 1;
$result = index($seq1, "n", $offset);
}
if ((keys %hseq1) > 10000){
$exclusion{$key1} = 1; last FIRST;
} ;
# foreach $key (sort { $a <=> $b } keys %hseq1){print $hseq1{$key} .  "," . $key . "\n"};

$offset = 0;
$result = index($seq2, "n", $offset);

while ($result != -1) {
$hseq2{$result}++; 
$offset = $result + 1;
$result = index($seq2, "n", $offset);
}
if ((keys %hseq2) > 10000){
$exclusion{$key1} = 1; last SECOND;
} ;
# foreach $key (sort { $a <=> $b } keys %hseq2){print $hseq2{$key} .  "," . $key . "\n"};

# place N in seq2 according to positions of N in seq1
foreach my $k ( keys %hseq1){ $k = $k;
substr( $seq2, $k, 1 ) =~ s/./n/ };

# place N in seq1 according to positions of N in seq2
foreach my $k ( keys %hseq2){ $k = $k ;
substr( $seq1, $k, 1 ) =~ s/./n/ };

$seq1 =~ s/n//g;
$seq2 =~ s/n//g;

undef %hseq1;
undef %hseq2;

# print STDERR "!" . $seq1 . "!\n" ;
# print STDERR "!" . $seq2 . "!\n" ;


if( ! $min_size || ( length($seq1) > $min_size && length($seq2) > $min_size )){
my $count = ( $seq1 ^ $seq2 ) =~ tr/\0//c;

# my $count = ( $h{$key1} ^ $h2{$key2} ) =~ tr/\0//c;
$sum+=$count;
$nb_comparison++;
# print STDERR $nb_comparison . "\n";
print OUT $key1 .  "," . $key2 .  "," . $count . "\n";

# print STDERR $h{$key1} .  "," . $key1 ;
# print STDERR $h{$key2} .  "," . $key2 . " $count\n";
}
} # end of the conditional if for restricted comparisons to make

}
}
} 

close OUT;
if( $nb_comparison > 0){
my $mean_pairwise_diff = $sum / $nb_comparison;


print STDERR $mean_pairwise_diff . " mean pairwise SNPs\n";

if($out2){
open (OUT2, ">", $out2) ;
print OUT2 $fasta_file .  "\t" . $fasta_file2 .  "\t" . $mean_pairwise_diff . "\n";
close OUT2;
}

}






