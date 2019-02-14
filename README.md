# power_tools
A series of scripts that I have developed for my day-to-day bioinformatic chores.

### doBlast.pl
Organism identification for reads in a BAM or fastq file using the blastn database. The script requires an istallation of samtools, blastn and blast database. It is capable of running parallel instances to try to speed up the search.

    doBlast.pl [ -c cores ] [-n num] [ -e e-value ] [-m] [-u] [-k] [-s] file

    Do a blast search from selected reads in a bam file ( unmapped default ) or from all reads in a fastq file.
  
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
