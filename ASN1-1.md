*Note: This assignment has a number of associated questions you'll find as you go through this that need to be submitted on ELMS. Some of these ask you to explain or guess why certain things may be the way that they are. I am not grading on correctness here, rather I am grading based on whether you have completed the assignment and put some thought into the answers. Feel free to discuss these with your group or with me in lab. We will likely go over some of the more important material during class as well.

# Submission

This assignment will be linked on ELMS and there should be a submission form for ASN1-1. You should submit this assignment individually. Feel free to work with and discuss the process and questions with classmates or peer mentors, but make please write your own answers and pick your own dataset for the last section!

# 1-i. Connecting to our computing resources.

## Connecting to the command line

Since the datasets we will be using are quite large (entire experiments can take more than a terabyte of space uncompressed) substantial parts of the analysis cannot easily be done on personal computers. We have a computer set up in the lab with the software we will need pre-installed with plenty of storage space. 

For the majority of the work, we will be connected to the command line interface remotely. We will be using a method called SSH or Secure SHell to do so. If you already have opinions or preferences about what software to use for this, feel free to use whichever method you like. 

For simplicity, I recommend downloading a free program called Termius https://termius.com/ to connect. 

To get started connecting to our server, after installing you should go to the "Hosts" menu and choose "+ NEW HOST". Once you do that you can fill in the Label, Address, and Username fields. Your username will be your university Directory ID. Unfortunately, we could not get integration with the university systems, so your password will not match. I've set your password to be your UID number. Since this is predictable, we'll have you change this immediately. This also happens to be the first grade point available for ASN2!

![Termius host menu](Images/Termius1.png)

Finally, click the save button and the connection will be stored in the list of Hosts. Double click on the TOD-Compute option and Termius should prompt you for your password (currently your UID). Enter that and you should be connected! If it shows a prompt mentioning that the host key is new you can accept this, this should only happen the first time you connect and the program will remember these details later.

![Termius connection](Images/Termius2.png)

Once you are connected, the first thing you should do is change your password. This is as simple as running one command:

```bash

passwd

```

Running `passwd` will prompt you for your current password, then ask you to enter your new password twice. For security, this program doesn't show your password or even asterisks as you are typing it, so don't worry when nothing shows up.

Successfully connecting to `tod-compute` and changing your account password so that I cannot login with your old password gets 1/4 points for ASN2. If you would like you can go back and edit the TOD-Compute host entry in Termius to save your password so that you do not need to type it in every time you connect.

At some points, we will also want to be able to transfer files to and from our personal computers and the server. Unfortunately, I am not a fan of any of the free programs that will work on all computers. I will suggest you downloading and install WinSCP if your personal computer runs Windows (https://winscp.net/eng/download.php) or Fetch if it runs macOS (https://fetchsoftworks.com/fetch/download/). WinSCP is free to use, but Fetch requires a license. Fortunately, this is free for Academic purposes, so once you get Fetch installed you should apply for the free license (https://fetchsoftworks.com/fetch/free).

### Connecting with WinSCP

Connecting with either program is similar to connecting with Termius. When you open WinSCP it will automatically open the Login window and default the connection method to "SFTP" which we will be using.

![WinSCP connection](Images/WinSCP1.png)

You can save this connection in the same way. Once you have this all set up, you can hit "Login" at the bottom. This again may warn you about the host key, go ahead and accept that again in this program.

This will open up a two-panel window similar to Windows Explorer, with your computer/files on the left side and the remote server files on the right side. Transferring files is as simple and clicking and dragging where you would like them to go.

### Connecting with Fetch

Connecting with Fetch should be nearly identical, on the first run the program will open up the New Connection window where you can put the hostname and username. Fetch defaults to the "FTP" connection mode which is insecure and not enable on our server. Simply use the dropdown and select "SFTP" and connect. Unlike WinSCP Fetch has only a single side showing the remote server files. To transfer files to or from the server you can click and drag from a Finder window.

![Fetch connection](Images/Fetch.png)

At this point, we should be set up with everything we need to use the computing resources (at least until we get to differential expression analysis and visualization in a few weeks!)

# 1-ii. Understanding our reference genome

## Introduction to shell usage for bioinformatics

The goal of the second part of this assignment is both to give you some very basic exposure to using computers via a command line interface instead of a more familar GUI or browser interface. This assignment will only utilize the bare minimum amount of tools and programs necessary to get started looking at the types of data and representations of biology we use in transcriptomics. Next week, you will use a DataCamp module to practice and gain exposure to a wider variety of useful tools and techniques.

If you have never used a linux/unix shell before, I strongly encourage you to work through the first free chapter of the DataCamp Introduction To Shell course, available at https://www.datacamp.com/courses/introduction-to-shell-for-data-science

## Reference genome sequences

To begin, we want to get more familiar with the genome and transcriptome we will be studying initially. In our case, this main genome is the human reference genome. We will be using the most recent version, version 38 from NCBI Ensembl (GRCh38). We will be interacting with the data in a few ways.

- Directly on the raw sequence data at the whole-genome, chromosome, or individual transcript levels
- Using annotated versions of the genomic data where genes, regulatory elements, or other features are listed by location
- With a web-based interacting database allowing us to search and explore a wide range of database resources from one location

To save some time and space during these steps, you'll be using only one chromosome at a time for this exercise, chromosome 22.

## Obtaining the reference genome

The web-based database for the genome assembly is available at the Ensembl website (https://useast.ensembl.org/Homo_sapiens/Info/Index). Another resource for viewing this information is available on the UCSC website and GRCh38/hg38 (https://genome.ucsc.edu/cgi-bin/hgGateway).

The most current human genome sequences are available on the EMBL-EBI FTP (ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/). These are compressed FASTA files containing raw assembled genome sequence. There are three variants of the files, one that is labeled "dna" and contains the full sequence. The other two are labeled "dna_rm" and "dna_sm". These are the same sequence, but with some repetitive or duplicated elements removed or "masked" out. We will be using the main "dna" files.

First, we should make a directory to do our work for this assignment. You can name this whatever you like or put it wherever you like.

```bash

cd ~
mkdir ASN2
cd ASN2

```

We can now download the sequence for chromosome 22 directly from the EMBL-EBI FTP.

```bash

wget ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
ls -lh

```

Since this file is compressed, we won't be able to work with it directly. The `.gz` suffix indicates that it is a "gzip" compressed file. We can decompress this file with the command `gzip -d [filename]`. 

Note: Whenever I give you a command with an element in `[square brackets]` that is just a placeholder. You should replace everything, including the brackets, with the text you need (in this case the filename of the file you just downloaded).

Q1) What does this do to the file and file name? (Check with `ls -lh` before and after)

Take a look inside the chr22 sequence file to look at the general format of the FASTA files:

```bash

head Homo_sapiens.GRCh38.dna.chromosome.22.fa

```

Then to look at a sort of summary view, use the `wc` command to count the number of lines, "words", and characters in the file:

```bash 

wc Homo_sapiens.GRCh38.dna.chromosome.22.fa

```

Q2) How many lines and characters are there in our sequence file? Use either the Ensembl or UCSC genome browsers to search for chr22 and find out how many DNA base pairs are present in chromosome 22. How well does this match up with what we find from the file and why might they be different?

The beginning of our sequence file may look a little funny. Use the line count you got from the `wc` command to take a look in the middle of the file. To easily pull out some lines in the middle, you can use a combination of the `head` and `tail` commands to get the front of the file up to a certain point, and then the end of that section. Aim for a number of lines roughly halfway through the file.

The `|` symbol is a special character that takes the output from one command, and sends it as input to the second command. In this case, instead of printing 425,000 lines of DNA to the screen, it sends those to the `tail` command, which instead will print the last 10 lines of its input.

```bash

head -n 425000 [filename] | tail

```

Q3) Do the beginning and middle of the file look similar? Why might one part be different than another? (Speculation is good here, I'm not judging based on correctness)

Q4) How many of each base (A, C, G, T, or N) are present in chromosome 22? There are multiple ways to answer this. Hint: using `grep -o ` with the `-o` flag outputs every match on a new line instead of entire lines. You can combine this with the `wc` command in a similar way to how we just combined `head` and `tail`.

Advanced question worth an extra half hour of lab time: AQ1) Restriction sites are DNA elements targeted by special enzymes known as "restriction enzymes" which cut DNA. How many EcoRI sites are present in chromosome 22? How did you determine this? (This need to be correct, think carefully about how one big contiguous sequence of chromosome is represented in this FASTA file)

### The purpose of the reference genome files

The FASTA reference genome sequence is essential for determining where any particular DNA/RNA sequence in our sequencing data derived from. Read alignment software like TopHat, BWA, or HISAT2 will require these FASTA files to match the short DNA sequence fragments to. After we know what part of the genome a particular fragment comes from, we can infer which transcripts this fragment may have been a part of and in aggregate how many of each transcript there were in the original sample.

## Reference transcriptome

The reference genome is a great resource and has every bit of genome information we know about at this point, but it represents the entire physical existence of the DNA chromosomes. We care about how that information gets transcribed into RNA, and that process involves selecting only certain bits of DNA, transcribing them into RNA, and then cutting and pasting bits of that RNA back together to generate the final transcript. The information that we need to predict how this happens is present in some annotation we will look at next, but the transcript sequences themselves are also available for us. These can be found on the same Ensembl FTP as the genome sequences in a different folder (ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/cdna/). 

These are labeled as cDNA, as opposed to RNA or transcripts. Since modern sequencing methods require DNA, to determine the sequence of RNA transcripts, a biology trick is necessary to "reverse transcribe" RNA into DNA. Since this generates a piece of DNA complementary to the RNA strand, this is known as "complementary DNA" or "cDNA". There are two types of sequence files available one labeled "abinitio" and one labeled "all". The *ab initio* sequences are a limited set of sequences predicted using simple models directly from the genome file. The *all* sequence files contain several times more transcripts, types of transcripts, and variant transcripts (isoforms) and represent all of the known transcripts compiled from experimental databases. 

Let's download each version:

```bash

wget ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.abinitio.fa.gz
wget ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gzip -d Homo_sapiens.GRCh38.cdna.abinitio.fa.gz
gzip -d Homo_sapiens.GRCh38.cdna.all.fa.gz

```

Since the FASTA format begins each sequence with a special line and follows with any number of lines of sequence, multiple sequences may be included in a single file. This is an example sequence file, the names and sequence will be different.

```

>Sequence1
CGGCATTATCGGAGCTAGCTAGCTGACTG
CGATGAGCATTACGGTA
>Sequence2
GCTAGCGATGCGGCGATTAGCGGGGAGAG
ATCGATGGGACGATGTATGCGTATG

```

The sequence name line always starts with the `>` character.

Q5) How many transcripts are in each version of the cDNA sequence files? (Hint: if you use the `>` symbol in a command line and you want to mean it is a text character, you should put it in quote `">"`, otherwise it will be used as a special character like the `|`.)

### The purpose of the reference transcript files

Not all methods for transcriptome analysis use direct read alignment to genomes. In recent years, tools such as Sailfish, Kallisto, or Salmon have become available to rapidly and directly estimate the original number of transcripts directly from the individual DNA fragment sequences. These approaches use lists of individual transcript sequences instead of the original genome sequence. Although these methods do not give you full alignments useful for things like variant calling or other techniques and cannot identify any new transcripts, they can be invaluable for efficiently quantifying the number of each transcript for expression analysis. Since this is the main focus of this FIRE stream, we will likely lean heavily on these transcript-focused approaches.

## Reference genome annotation

With only the raw DNA sequence, or even the raw cDNA sequences, there isn't a lot we can do. I certainly can't figure out much biologically meaningful information manually from one cDNA sequence, much less tens of thousands. Fortunately, since we began sequencing and studying the human genome, thousands of researchers have spent decades investigating every bit of the human genome. Since this is science, the community has published, collected, and organized this information in a variety of ways to let use associate bits of the human genome with relevant biological information. In general, this process is called *annotation*. 

There are several common file formats for containing gene-level annotation of genomes, but perhaps the most common is GTF/GFF. Both file types are very similar, with some subtle differences (https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/). Both formats are simple text files containing one feature annotation per line. Lines that start with `#` are "comment" lines and don't represent features but do include information that might be useful for orienting yourself with the file. Each normal line contains nine columns separated by tabs. The first eight columns are simple names, numbers, or symbols representing the name, type, and location of the feature. The final column can be much more complicated and contains a list of attributes in name=value pairs. The details of the GFF3 format may be found at https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

Since these features are referenced to specific numeric locations in a sequence file, any annotations much match up to a specific reference genome sequence. The Ensembl FTP has files available both for individual chromosomes as well as the entire chromosome or transcript cDNA files. 

```bash

wget ftp://ftp.ensembl.org/pub/release-86/gff3/homo_sapiens/Homo_sapiens.GRCh38.86.chromosome.22.gff3.gz
gzip -d Homo_sapiens.GRCh38.86.chromosome.22.gff3.gz
head -n 20 Homo_sapiens.GRCh38.86.chromosome.22.gff3

```

Q6) The third column of the gff3 files contain the type of feature. How many different types of features are present on chromosome 22? How many coding sequences (CDS) are there? How many standard "gene" annotations are there? What could this tell you about the relationship between "genes" and coding sequences? You may want to refer to other resources like in the introduction to this section.

###  The purpose of the reference sequence annotation

When running alignment-based approaches, the exon and transcript level annotations are essential for determining what RNA/cDNA sequences might exist. Individual strings of DNA letters may show up in our RNA sequencing data that do not exist together in the original genome. These might be due to events like splicing which fuse two distant sequences together, or from non-templated addition like in the poly-A tails at the end of mRNAs.  Aligners like TopHat or HISAT2 will use these annotations to properly map these types of reads.

Beyond the process of alignment, these annotation files contain the information that lets us know what gene or gene variant is associated with any particular transcript. Knowing what sequence is what gene is essential for interpreting the results, understanding what changes and what doesn't change and therefore what is happening in the biology of our samples. This is the case whether we choose to align directly to the genome or use transcriptome-level quantification tools.

# 1-iii. Obtaining transcriptome quantification data

So far we have taken a look specifically at the human reference information. None of this is specific to the experiments we want to analyze. Much like the reference, our data is present in some sequence databases. In contrast to the reference genome though, our data is in a more specialized database intended to provide a warehouse for large amounts of publicly-accessible high-throughput sequence data called the Short Read Archive (SRA). The SRA is hosted by the NCBI, part of the National Institute of Health in Bethesda, Maryland. It is part of an interconnected series of databases which include a variety of annotation about about experimental data. 

## Exploring the data in the NCBI databases

Some of the relevant NCBI databases for obtaining and understanding our data will be `BioProject`, `BioSample`, `SRA`, and `PubMed`. The first of these, `BioProject`, (https://www.ncbi.nlm.nih.gov/bioproject/) is designed to contain, categorize, and describe whole experimental projects. This database will contain links not only to the raw data hosted in databases like `SRA` or `GEO`, but also to any publications that may be associated with this data as well as metadata associated with the raw data itself. 

Some of this metadata is present in the `BioSample` database. This is intended to be a higher-level description of the raw data, each entry is, much like the name, a biological sample that some data derived from. This could be a particular batch of cells grown in a test tube, an environmental sample of water, a biopsy of human tissue, or essentially any other biological unit of physical *stuff* we might sample to generate data. The `BioSample` entries themselves may then be associated with raw data that was derived from that particular biological sample. The quality of this annotation may vary from experiment to experiment because it is up to the individual scientist submitting this to the database to decide which information to include in the submission to the database.

The `BioProject` entries themselves often contain only a brief summary of the experimental goals and data content. In many cases, to really understand what all of the samples represent and how they were generated, you may need to read the associated publications to get enough information. 

We will be working with data that largely resides in two `BioProject` entries. The accession number for the *Leishmania* datasets has the `BioProject` accession PRJNA290995. The data for the *Trypanosoma cruzi* is associated with `BioProject` PRJNA251582. Searching for these in the link provided above should bring up the summary page for each project. Clicking on the number of SRA experiments brings you to a result list in the SRA database. A link at the top of the page should bring you to the SRA Run Selector which lets you interact and download a table containing more information about the data itself.
 
![Run selector](Images/SRASelector.png)
 
Q7) How many RNA-seq datasets are present in each of the two projects? For the *Leishmania major* experiment, how many of these correspond to samples containing parasite infections and how many to non-infected control samples? What other potentially import biological variables can you find in the annotation for these experiments?
 
Q8) Who is the lead author of the journal article associated with the *Leishmania* project? What publication was this published it?
 
## Downloading raw data
 
There are several ways to download the raw sequence data from the SRA database. These can be directly downloaded from the website, from several cloud services, or through command line tools called `sra-tools`. Because together these datasets take up almost an entire terabyte even compressed, these have already been downloaded and are present in a folder on tod-compute at `/mnt/storage/data/` in the `leishmania` and `tcruzi` folders respectively. If you do need to download any additional datasets, instructions for using the `fasterq-dump` program to directly download SRA datasets is found at https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
 
You can list the downloaded datasets:
 
```bash
 
ls -lh /mnt/storage/data/leishmania
ls -lh /mnt/storage/data/tcruzi
 
```
 
Since these are very large files, we would like to keep them compressed as much as possible. Forunately, through the magic of piping, in many cases it is not necessary to run gzip before you run additional commands, the `gzip` program can be chained just like other unix programs. Adding the `-c` flag to the gzip command forces it to read input or send output to the console, allowing us to either directly print or pipe the data to other programs. Each dataset consists of two files, ending in either `_1.fastq.gz` or `_2.fastq.gz`. This is a consequence of the method of sequencing our samples used, which generates sequence from both ends of each DNA fragment, in a process called "paired-end sequencing". Each set of four lines in each file represent a single read. 

```bash

gzip -dc /mnt/storage/data/tcruzi/SRR1346026_1.fastq.gz | head
gzip -dc /mnt/storage/data/tcruzi/SRR1346026_2.fastq.gz | head

```

Q9) Given the output of the previous command, how can you tell which read from file 1 corresponds to its paired read in file 2?

## Initial quality control

There are several programs designed to generally analyze the output from high-throughput sequencing to screen for a variety of common problems or contamination issues that arise. Perhaps the most popular is a program called FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). You can use FastQC in two modes. One mode involves running FastQC like a normal GUI program. Download it onto your computer, unzip, and run the `run_fastqc.bat` file. You can then open up `fastq.gz` files you have downloaded to your computer and analyze them directly. The other mode runs via the command line and instead generates HTML output that you can download and view in a web browser. For the final part of the assignment, you have a choice to follow either method. Both will ultimately involve transferring files from the server to your own computer with WinSCP, Fetch, or other SFTP program.

### Run FastQC graphically on your own computer

If you choose to run FastQC on your own computer, ensure that you have at least 15 GB of free space on your hard drive first. This method is also only advisable if you are on a fast on-campus internet connection. Finally, you will need Java installed on your computer. If none of these apply, please don't use this method and skip to method 2. Download FastQC from the link above, choosing either "FastQC v0.11.8 (Win/Linux zip file)" if you run Windows or "FastQC v0.11.8 (Mac DMG image)" if you have a Mac. Extract the file you downloaded and run FastQC following the instructions for either version.

Next, you will pick one of the SRA datasets from *Leishmania* folder. **Do not use the SRR1346026 dataset in the examples.** Connect to the server with your SFTP program of choice and navigate to `/mnt/storage/data/`. Enter either the `tcruzi` or `leishmania` folders and download one dataset (both the `_1.fastq.gz` and `_2.fastq.gz` files). Finally, in the FastQC program go to the open file menu, navigate to, and select the two `.fastq.gz` files you downloaded from the server. This will take a few minutes to analyze both files and should show you the progress as it works.

### Run FastQC remotely on tod-compute and download the results
 
As an alternative to running FastQC on your computer and directly viewing the results, you can use the FastQC program installed on the tod-compute server. Pick one of the SRA datasets from the *Leishmania* folder. **Do not use the SRR1346026 dataset in the examples.** Next, you should run the `fastqc` program and give it the full path to the files you chose, as well as the location you want to store the output. On the command line the `.` character as a file represents the current directory you are in. Either of the following commands will analyze both files at the same time. This will take a few minutes to analyze both files and should show you the progress as it works.
 
```bash

fastqc /mnt/storage/data/tcruzi/[dataset]_1.fastq.gz /mnt/storage/data/tcruzi/[dataset]_2.fastq.gz --outdir=. 
fastqc /mnt/storage/data/tcruzi/[dataset]_*.fastq.gz --outdir=. 

#Examples
fastqc /mnt/storage/data/tcruzi/SRR1346026_1.fastq.gz /mnt/storage/data/tcruzi/SRR1346026_2.fastq.gz --outdir=. 
fastqc /mnt/storage/data/tcruzi/SRR1346026_*.fastq.gz --outdir=. 
 
```

After this completes, running `ls` in your current (output) directory should show several additional files. There should be two additional `.html` files and two `.zip` files. You will need to download these to your own computer to view them. Connect to the server with your SFTP program of choice and navigate to the folder your are currently in. If you followed this exactly this will probably be `ASN2` in your home directory. Download the `.html` files to a folder on your computer, and open these in your preferred web browser.

### In either case...

Once you have the results for both files, take a look through the different results sections. For the following questions you may want to refer to this guide from Michigan State University for help interpreting the results: https://rtsf.natsci.msu.edu/genomics/tech-notes/fastqc-tutorial-and-faq/

Q10) Compare the html files results from both of the read 1 and read 2 files you analyzed. Do you notice any differences between them? How many reads were present in the dataset you chose? (Make sure to tell me which dataset you chose)

Q11) Did any of the FastQC analysis modules flag potential problems? For each problem, refer to the MSU guide. Do you think we should worry about these problems? If not, why not? If so, what might we do to address the problem?

