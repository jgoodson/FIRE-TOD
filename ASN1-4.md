 # 3-i. Feature counting with HTSeq-count

## How do we know which read corresponds to which feature?

We have effectively generated a list of genomic positions where we found RNA molecules in our experiment. Remembering last week, we looked at genome annotation files that similarly were a list of genomic positions of different types of biologically-relevant features. If we was very patient and thorough, we could simply go through these lists and cross reference the positions, find a read at position 82539982 of chromosome 15, look up this position in a genome browser, and see that this corresponds to an exon of the human ribosomal protein RPS17. This would take a *very* long time to do for a billion sequence reads. Fortunately, this process can be automated. Use a simple algorithm to run through each list and find which genomic regions overlap. 

![Overlap types](Images/count_modes.png)

There are some complications to this. What do you do when a read only partially overlaps a feature you care about? What if a read overlaps two genes? Is this because the read spans them or because the two genes themselves overlap? The above chart shows a number of these possibilities. The tool we will be using, HTSeq-count has three modes that deal with these in different ways. 

**Q7) Both intersection_strict and intersection_nonempty have largely similar results. In what situations do these assign reads differently? In what situations might the less-strict intersection_nonempty keep more information?**

## Running htseq-count

To start counting our features, we only need two things. The feature annotations and the aligned reads that will correspond to them. We just generated the SAM alignment file in the previous section. The annotation we will use is similar to what we looked at last week, but we will use the full genome this time.

```bash

wget ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.gtf.gz
gzip -d Homo_sapiens.GRCh38.86.gtf.gz

```

HTSeq-count has been pre-installed for you and is ready to go. We will run it with a command something like this:

```bash

htseq-count --mode intersection-strict --order name --type exon --idattr gene_id [SRA Accession].sam Homo_sapiens.GRCh38.86.gtf > [SRA Accession].counts

head [SRA Accession].counts
tail [SRA Accession].counts

```

We can see that HTSeq-count outputs two tab-delimited columns, with a label on the left and a number on the right. The top of the file looks pretty consistent, while the end has some unusual labels starting with two underscores `_` corresponding to some categories that could not be mapped to specific features.

Going back to a technique we covered last week, we can extract single columns from these kind of files using either `cut` or `awk`.

```bash

cut -f 1 [SRA Accession].counts | head
cut -f 2 [SRA Accession].counts | head

awk '{print $1}' [SRA Accession].counts | head
awk '{print $2}' [SRA Accession].counts | head

```

There is a companion to the `cut` command, aptly named `paste`. Instead of separating out a single column out of one file, `paste` can combine multiple inputs with columns into one columns. Make sure you have done the complete process of alignment and gene counting for both of our datasets to this point.

```bash 

cut -f 1 [SRA Accession 1].counts > gene_labels.txt
cut -f 2 [SRA Accession 1].counts > [SRA Accession 1].values
cut -f 2 [SRA Accession 2].counts > [SRA Accession 2].values
paste -d ',' gene_labels.txt [SRA Accession 1].values [SRA Accession 2].values > total_counts.csv

```

**Q8) Download the total_counts.csv filet computer and open it in Excel, Google Docs other spreadsheet program. How many columns are there in this file? Find the sum of each numerical column. What do these columns add up to?

Now that we have this sort of tabular summary of each gene in the genome, with values corresponding to the amount of RNA in each sample, we are ready to begin our statistical analysis! 

Next week, we will have count tables like this prepared for all of the samples in both experiments. You will do a DataCamp assignment (ASN4) to get an introduction to R, a statistical scripting language that will allow us to use some very powerful analysis tools to follow up on these tables.