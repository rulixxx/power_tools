#!/usr/bin/perl
# get the frequency of each barcode in an index file

use Getopt::Std;
use strict;

my %opt;
getopts('pa', \%opt);

my $narg = @ARGV;
die ("\nbarcodeFreq.pl index1.fastq [ index2.fastq ]\n 
      Compute the frequency of each barcode present in an index fastqs.
      By default only process the first 1000000 sequences\n
            -a  process all the fastq
            -p  include percenages in output\n\n") if ($narg == 0 );

my $command = 'cat';
$command = 'zcat' if ( $ARGV[0] =~ /\.gz$/ );

if (defined $opt{'a'} )
   {
   open(FILE1,"$command $ARGV[0] |") or die ("$ARGV[0] file not found\n");
   open(FILE2,"$command $ARGV[1] |") or die ("$ARGV[1] file not found\n") if ($narg == 2);
   }
else
   {
   open(FILE1,"$command $ARGV[0] | head -4000000 |") or die ("$ARGV[0] file not found\n");
   open(FILE2,"$command $ARGV[1] | head -4000000 |") or die ("$ARGV[1] file not found\n") if ($narg == 2);
   }

my $nlines = 0;
my (%hash, $f1, $f2);
while($f1 = <FILE1> )
   {
   $nlines++;
   $f2 = '';
   $f2 = ' ' . <FILE2> if ( $narg == 2);
   next if ( ($nlines % 4) != 2 );
   chomp($f1);
   chomp($f2);
   $hash{ $f1 . $f2} ++;
   }
my @order = sort { $hash{$a} <=> $hash{$b} } keys %hash;
$nlines /= 4;
my ($ii,$freq,$pct);
foreach $ii (@order)
   {
   $freq = $hash{$ii};
   if ( defined $opt{'p'} )
      {
      $pct = sprintf ("%.2f", $freq / $nlines * 100);
      print "$ii\t$freq\t$pct\n";
      }
   else
      {
      print "$ii\t$freq\n";
      }
   }
