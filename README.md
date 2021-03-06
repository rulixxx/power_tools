# FASTQ power tools
A series of scripts that I have developed for my day-to-day work with FASTQ files.

### barcodeFreq.pl
This simple script take one or two index fastq files and produces a frequency (and percentage) of the occurrence of each barcode. It is usually the first tool to use to examine experiments that had problems demultiplexing.

        barcodeFreq.pl index1.fastq [ index2.fastq ]
 
            Compute the frequency of each barcode present in an index fastqs.
            By default only process the first 1000000 sequences

                    -a  process all the fastq
                    -p  include percenages in output

Sample output (dual index, frequency , percentage):
        
        TGGCTAATCA CTACATGCCT	12059	1.21
        GGAATTGTAA AGGATGTGCT	12164	1.22
        GCATAAGCTT TGCGACGGAA	13088	1.31
        CCACCAGGCA ATTCCATAAG	13130	1.31
        CGTTATTCTA AACCTTATGG	13156	1.32
        GTGCAGACAG ACGGATGGTA	13367	1.34
        TTGGAATTCC GTATTGACGT	13531	1.35
        AACATCGCGC ACAAGTGGAC	14018	1.40
        AGATCCATTA TGGTAGAGAT	14518	1.45
        TGGCGCGAAC AGTATCAGTT	14846	1.48
        GGAGCGTGTA ACTTACGGAT	15595	1.56
        GATCAAGGCA ATTAACAAGG	15860	1.59
        TATTGCGCTC CCTAACACAG	16675	1.67
        GAAGCGGCAC AGTAGAGCCG	17017	1.70
        CATTAGTGCG TATGTAGTCA	17760	1.78
        AAAAAAAAAA AGATCTCGGT	36275	3.63

### grepFastq.pl
Efficient grep utility for fastq files, can search either read IDs or sequences. Very useful for filtering fastq files or quickly finding specific reads.

    grepFastq.pl [-avso] [-f pattern_file] [-q quality] [-t tiles] [ pattern ] fastq

      Find the fastq register(s) that match the specified pattern.
      Default is to match sequence ids, stop at the first match
      and output to stdout.

         -s match sequence instead of the sequence id
         -f specify a file with patterns
         -t remove the specified tiles (comma separated list)
         -a find all posible matches, otherwise find only first match
         -q return reads with mean quality greater/equal to this value
         -v invert match
         -o output to a file instead of stdout
         
 Sample usage:
        
        grepFastq.pl 1101:3671:1314 file.fastq
        
        @K00310:254:H5FLMBBXY:8:1101:3671:1314 1:N:0:CTTGGATG+TCGCTATC
        AAGGCCCAGTTAGATCCTCCCCTCTTAGCTGTCTTAAAGAATGGCAGACTGACAGTTGGTCAGAAGATTATTCTTCATGGAGCAGAACTGGTGGGCTCTCC
        +
        AA<AFJJJJJJAAJJJJFJJJJJJJFJFFFJAFFFFJJJJJA-7<F-F--7FFJ7----7<J7FJ<JJFFFAF--AAFJFJFFJJJJJFFFJJJJJJJFJJ
         

### doBlast.pl
Organism identification for reads in a BAM or fastq file using the blastn database. The script requires an istallation of samtools, blastn and Blast database. It is capable of running parallel instances to try to speed up the search.
Its main usefullness has been to look at unmapped fractions in BAM alignment files to try to detect possible contaminating sequences. I have also used it sucessfully to do quick estimates of percentages of ribosomal & mitochondrial sequences present in a library prior to alignment.

    doBlast.pl [ -c cores ] [-n num] [ -e e-value ] [-m] [-u] [-k] [-s] file

    Do a blast search from selected reads in a bam file ( unmapped default ) or from reads in a fastq file.
  
      -c  secifiy number of cores (default 1 )
      -n  specify number of sequences to use ( default 1000)
      -e  specify an e-value for Blast ( default 1e-10 )
      -m  do blast on multiple mapping reads for a BAM
      -u  do blast on unique mapping reads for a BAM
      -k  keep fasta files and print their locations
      -s  print output to stdout

Sample output ( score, evalue, subject title )

    32	2e-06	Kocuria rhizophila DC2201 DNA, complete genome
    28	4e-04	Mesorhizobium opportunistum WSM2075, complete genome
    86	2e-36	Metarhizium flavoviride strain ARSEF 2025 NADH dehydrogenase subunit 1 (nad1) and NADH dehydrogenase subunit 4 genes, partial cds; mitochondrial
    76	7e-31	Metarhizium anisopliae strain ME1 mitochondrion, complete genome
    101	9e-45	Microlaena stipoides chloroplast, partial genome
    44	5e-13	Mycobacterium abscessus chromosome, complete sequence
    34	2e-07	Mycobacterium smegmatis str. MC2 155, complete genome
    72	1e-28	Mycobacterium gilvum PYR-GCK, complete genome
    51	6e-17	Mycobacterium rhodesiae NBB3, complete genome
    55	4e-19	Nectria haematococca mpVI 77-13-4 predicted protein, mRNA
    84	3e-35	Neotyphodium sinofestucae strain Fnj4604 laccase gene, complete cds
    56	1e-19	Neurospora crassa OR74A DNA polymerase epsilon partial mRNA
    101	9e-45	Oryza rufipogon voucher AusTRCF 309313 chloroplast, complete genome
    74	1e-29	Oryza sativa Indica Group strain WA-CMS mitochondrion, complete genome
    73	3e-29	Oryza sativa Indica Group strain WA-CMS mitochondrion, complete genome 


