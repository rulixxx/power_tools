#!/usr/bin/perl
#

use List::Util;
use Getopt::Std;
use POSIX;
use File::Temp qw/ tempfile /;
use strict;

#required installations of samtools and blastn + DB
my $samtools = "/apps/SAMTOOLS/1.3/AVX0/GCC/bin/samtools";
my $blastn = "/apps/BLAST+/2.2.28/bin/blastn";
my $blastnDB = "/scratch/devel/jcamps/blast/nt";

my %opt;
getopts('c:n:e:kfums', \%opt);
my $n=@ARGV;

die "\ndoBlast.pl [ -c cores ] [-n num] [ -e e-value ] [-m] [-u] [-k] [-s] file\n\n" .
   "Do a blast search from selected reads in a bam file ( unmapped default ) or from all reads in a fastq file.\n\n" .
   "   -c  secifiy number of cores (default 1 )\n".
   "   -n  specify number of sequences to use ( default 1000)\n".
   "   -e  specify an e-value for Blast ( default 1e-10 )\n".
   "   -m  do blast on multiple mapping reads for a BAM\n".
   "   -u  do blast on unique mapping reads for a BAM\n".   
   "   -k  keep fasta files and print their locations\n".
   "   -s  print output to stdout\n\n\n" if ($n == 0);

my $file = $ARGV[0];
my $fastqgz = 0;
my $fastq = 0;
my $bam = 0;
my $eval = '1e-10';
$eval = $opt{'e'} if ( defined $opt{'e'} );
#don't delete the temp files
$File::Temp::KEEP_ALL = 1 if ( defined $opt{'k'} );
my $nseq = 1000;
$nseq = $opt{'n'} if ( defined $opt{'n'} );

my ($command,$count);
#identify the file type from its extension
if ( $file =~ /fastq\.gz$/ ) #compressed fastq
   {
   $fastqgz = 1;
   $command = "zcat $file"; 
   $count = int(`$command | wc -l`);
   $count = int( $count / 4);
   open(IFILE,"$command |" ) or die "Cant open pipe to $file\n";
   }
elsif ( $file =~ /\.fastq$/ ) #uncompressed fastq
   {
   $fastq = 1;
   $command = "cat $file";
   $count = int(`$command | wc -l`);
   $count = int( $count / 4);
   open(IFILE,"< $file" ) or die "Cant open file $file\n";
   }
elsif ( $file =~ /\.bam$/ ) #BAM file
   {
   $bam = 1;
   if ( defined $opt{'u'} )
      {
      #unique maps ( pair mapped and has a map quality > 20 )
      $command = "$samtools view -f2 -q20 $file";
      $count = int(`$command | wc -l`);
      open(IFILE,"$command |" ) or die "Cant open pipe to $file\n";
      }
   elsif ( defined $opt{'m'} ) 
      {
      #multi map ( pair mapped in proper pair and has a map quality < 20 )
      $command = "$samtools view -f2  $file  | awk '{ FS = \"\t\"}; { if (\$5 < 20) print \$0 }'";
      $count = int(`$command | wc -l`);
      open(IFILE,"$command |" ) or die "Cant open pipe to $file\n";
      }
   else
      {
      #unmaped ( default)
      $command = "$samtools view -f68 $file";
      $count = int(`$command | wc -l`);
      open(IFILE,"$command  | " ) or die "Cant open pipe to $file\n";
      }
   }
else
   {
   exit("File type not identified\n");
   }

my $threads = 1;
$threads = $opt{'c'} if ( defined $opt{'c'} );
my @BOUT= ();
my @FASTA = ();
my ($i, $iFASTA, $ifname, $iBOUT, $ibname, $nmin );

#store the file names and handles for the input and output files
for ($i =0; $i < $threads; $i++)
   {
   ($iFASTA, $ifname) = tempfile( UNLINK => 1); #child process in file (delete default)
   ($iBOUT,$ibname) = tempfile( UNLINK => 1); #child process out file (delete default)
   push(@BOUT, {'handle' => $iBOUT ,'name' => $ibname } );
   push(@FASTA, {'handle' => $iFASTA ,'name' => $ifname } );
   }

$nmin = floor($nseq/$threads);
$nseq = $opt{'n'} if defined $opt{'n'};
my ( $rn, %rands, $l1, $l2, $l3, $l4 );
my ( $ifile, @vals, $str );
srand(); #initialize the random number generator

my $nsel = 0;
#select lines at random
while ( $nsel < $nseq  )
   {
   $rn = int(rand($count)) + 1;
   if ( not exists $rands{$rn} )
      {
      $rands{$rn} = 1;
      $nsel++;
      }
   }

my $tseq = keys %rands;
my $iline = 0;
my $mline = 0;
#create the input fasta files (= threads)
if ( $fastqgz || $fastq ) #we have a fastq file
   {
   while(1)
      {
      $iline++;
      $l1 = <IFILE>;
      $l2 = <IFILE>;
      $l3 = <IFILE>;
      $l4 = <IFILE>;
      next if (not defined $rands{"$iline"});
      $mline++;
      $l1 =~ s/^@/>/;
      $ifile = ($mline - 1 ) % $threads;
      print {$FASTA[$ifile]{'handle'}} $l1;
      print {$FASTA[$ifile]{'handle'}} $l2;
      last if ($mline == $tseq );
      }   
   }
elsif ($bam)  #we have a bam file
   {
   while(<IFILE>)
      {
      $iline++;
      next if (not defined $rands{"$iline"});
      $mline++;
      @vals = split('\t');
      $str = ">".$vals[0]."\n".$vals[9]."\n";
      $ifile = ($mline - 1 ) % $threads;
      print {$FASTA[$ifile]{'handle'}} $str;
      last if ($mline == $tseq);
      }
   }
close(IFILE);

#flush the fasta files
my ($val, $handle);
foreach $iFASTA (@FASTA)
   {
   $val = $iFASTA -> {'name'};
   $handle = $iFASTA -> {'handle'};
   print STDERR "fasta file: $val\n" if ( defined $opt{'k'} );
   $handle -> flush;
   }

#launch n instances of blastn
my ( $j, $blastCom, $fname, $bname, $bb, %kid);
for ($i=0;$i<$threads;$i++)
   {
   $j = $i + 1;
   #print "Launching: Blast # $j\n";
   $fname = $FASTA[$i]{'name'};
   $bname = $BOUT[$i]{'name'};
   $bb = "echo blast$i";
   $blastCom = "$blastn -db $blastnDB -outfmt \"7 score evalue stitle\" -num_threads 1 -max_target_seqs 1 -query $fname  -out $bname -evalue $eval -word_size 18"; #keep only top target, default word size 11 set to 18 to speed up.
   defined(my $cpid = fork) or warn $! and next;
   $kid{$cpid} = undef, next if $cpid;
   exec $blastCom;
   exit 0;
   }
delete $kid{wait()} while %kid;

my( @result);
#store the resulting hits
foreach $iBOUT (@BOUT)
   {
   $val = $iBOUT -> {'name'};
   $handle = $iBOUT -> {'handle'};
   while(<$handle>)
      {
      if ( $_ =~ /# (\d+) hits found/ )
         {
         next if ($1 == 0);
         $_ = <$handle>;
         push @result,$_;
         }
      }
   }
#sort output based on description to make it more readable
@result = sort { (split(' ',$a))[2] cmp (split(' ',$b))[2] } @result;

my ($fileName, $prefix, $fc, $lane, $index, $sql, $doQuery, $subproject, $spath, $newFileName, $outFile, $linkFile);
if ( defined $opt{'s'} )
   { 
   print @result;
   }
else
   {
   open(OUT,"> $file".".blast");
   foreach $i ( @result )
      {
      print OUT "$i";
      }
   close OUT;
   }

#clean up the tempoarary files
foreach $iBOUT (@BOUT)
   {
   $handle = $iBOUT -> {'handle'};
   close $handle;
   }
foreach $iFASTA (@FASTA)
   {
   $handle = $iFASTA -> {'handle'};
   close $handle;
   }
