#!/usr/bin/perl
#
use Getopt::Std;
use strict;

my %opts;
getopts('shaovf:t:q:', \%opts);
my $nargs = @ARGV;
if ( defined $opts{'h'} or $nargs == 0  )
   {
   print "\n";
   print "   grepFastq.pl [-avso] [-f pattern_file] [-q quality] [-t tiles] [ pattern ] fastq\n";
   print "\n";
   print "      Find the fastq register(s) that match the specified pattern.\n";
   print "      Default is to match sequence ids, stop at the first match\n";
   print "      and output to stdout.\n";
   print "\n";
   print "         -s match sequence instead of the sequence id\n";
   print "         -f specify a file with patterns\n";
   print "         -t remove the specified tiles (comma separated list)\n";
   print "         -a find all posible matches, otherwise find only first match\n";
   print "         -q return reads with mean quality greater/equal to this value\n";
   print "         -v invert match\n";
   print "         -o output to a file instead of stdout\n";
   print "\n";
   exit 0;
   }


die "cannot specify -t and -f switches at the same time\n" if ( defined $opts{'t'} and defined $opts{'f'} );
die "cannot specify -t and -q switches at the same time\n" if ( defined $opts{'t'} and defined $opts{'q'} );
die "cannot specify -f and -q switches at the same time\n" if ( defined $opts{'f'} and defined $opts{'q'} );
die "cannot specify -t,-q or -f and a pattern at the same time\n" if ( (defined $opts{'f'} or defined $opts{'t'} or defined $opts{'q'} ) and $nargs > 1 );

print $ARGV;
my $all = 0;
$all = 1 if ( defined $opts{'a'});
my $fseq = 0;
$fseq = 1 if ( defined $opts{'s'});
my $doTiles = 0;
$doTiles = 1 if ( defined $opts{'t'});
my $invert = 0;
$invert = 1 if ( defined $opts{'v'});
*OUT = \*STDOUT;

my $expressionFile;
my @tiles;
my $expression;
my $fastq = $ARGV[0];
my $qcutoff;

if (defined $opts{'f'} )
   {
   $expressionFile = $opts{'f'} ;
   }
elsif (defined $opts{'t'} )
   {
   @tiles = split(/,/,$opts{'t'});
   }
elsif (defined $opts{'q'} )
   {
   $qcutoff =  $opts{'q'} + 0.0;
   }
else
   {
   $expression = $ARGV[0];
   $fastq = $ARGV[1];
   }
die "$fastq file not found!\n" if ( not -e $fastq);
if ( $fastq =~ /\.gz$/ )
   {
   open(FASTQ,"zcat $fastq |");
   }
else
   {
   open(FASTQ,"<$fastq");
   }

my $nout;
if ( defined $opts{'o'} ) 
   {
   if ( $fastq =~ /(.+)\.gz$/ )
      {
      $nout = $1;
      }
   else
      {
      $nout = $fastq . ".tmp"; 
      }
   open(OUT, ">$nout" );
   }

#get patterns from a file
my ($head,$seq,$plus,$quality,$match,$comp);
if  (defined $opts{'f'} )
   {
   open(PATT,"< $expressionFile") or die "File $expressionFile not found\n";
   my @expressions = ();
   while(<PATT>)
      {
      chomp;
      push(@expressions,$_);
      }
   close(PATT);
   while(<FASTQ>)
      {
      $head = $_;
      $seq = <FASTQ>;
      $plus = <FASTQ>;
      $quality = <FASTQ>;
      $match = 0;
      foreach $expression (@expressions)
         {
         if ( $fseq )
            {
            $comp = $seq;
            }
         else
            {
            $comp = $head;
            }
         if ( $invert )
            {
            $match = 1 if ( $comp !~ /$expression/);
            }
         else
            {
            $match = 1 if ( $comp =~ /$expression/);
            }
         if ($match)
            {
            print OUT $head;
            print OUT $seq;
            print OUT $plus;
            print OUT $quality;
            last if ( not $all);
            }
         }
      }
   }
#filter by tiles
elsif (defined $opts{'t'} )
   {
   my (@vals,$tile);
   while(<FASTQ>)
      {
      $head = $_;
      $seq = <FASTQ>;
      $plus = <FASTQ>;
      $quality = <FASTQ>;
      $match = 0;
      @vals = split(/ /,$head);
      @vals = split(/:/,$vals[0]);
      $tile = $vals[4];
      die "Check tile field\n" if ( length($tile) != 4 );
      if ( $invert )
         {
         next if ( $tile ~~ @tiles );
         }
      else
         {
         next if ( not ( $tile ~~ @tiles) );
         }
      print OUT $head;
      print OUT $seq;
      print OUT $plus;
      print OUT $quality;
      }
   }
#filter by average quality
elsif ( defined $opts{'q'} )
   {
   my ($iq, $sumq,$aveq);
   while(<FASTQ>)
      {
      $head = $_;
      $seq = <FASTQ>;
      $plus = <FASTQ>;
      $quality = <FASTQ>;
      $sumq = 0;
      chomp($quality);
      foreach  $iq (split //, $quality )
         {
         $sumq += ord($iq) - 33;
         }
      $aveq = $sumq/length($quality);
      if ( $invert )
         {
         $match = 1 if ( $aveq < $qcutoff);
         }
      else
         {
         $match = 1 if ( $aveq >= $qcutoff);
         }
      if ($match)
         {
         print OUT $head;
         print OUT $seq;
         print OUT $plus;
         print OUT "$quality\n";
         last if ( not $all);
         }
      }
   }
#filter by header or sequence
else
   {
   while(<FASTQ>)
      {
      $head = $_;
      $seq = <FASTQ>;
      $plus = <FASTQ>;
      $quality = <FASTQ>;
      $match = 0;
      if ( $fseq )
         {
         $comp = $seq;
         }
      else
         {
         $comp = $head;
         }
      if ( $invert )
         {
         $match = 1 if ( $comp !~ /$expression/);
         }
      else
         {
         $match = 1 if ( $comp =~ /$expression/);
         }
      if ($match)
         {
         print OUT $head;
         print OUT $seq;
         print OUT $plus;
         print OUT $quality;
         last if ( not $all);
         }
      }
   }
close(OUT);
