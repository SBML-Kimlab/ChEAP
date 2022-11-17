"""ChIP-exo_analysis_pipeline ChEAP

# **1. Environment setting**

Please run the cell without modification
"""
from os import system

system('pip install pandas')
system('pip install pysam')
system('pip install numpy')
system('pip install biopython')
system('pip install bowtie')
system('pip install samtools')

from numpy import zeros
from Bio import SeqIO
import pandas as pd
import os, math, pysam, shutil


#Function
def make_gff_from_ncbigb( file_gb, file_output, types = [ "CDS", "tRNA", "rRNA", "tmRNA", "misc_RNA" ] ):
    handle = open( file_gb, "r" )
    genomes = SeqIO.parse( handle, "gb" )    
    fout = open( file_output, "w" )

    for genome in genomes:  
        genome_id = genome.name
        print ( "genome id=" + genome_id, "length=" + str( len( genome ) ) )
        
        for f in genome.features:            
            if f.type in types :
                start = str( f.location.start.position )
                end = str( f.location.end.position )
                strand = f.strand
                try:
                    attr_gene = "/".join( [ q for q in f.qualifiers[ "gene" ] ] )            
                except:
                    attr_gene = ""
                try:
                    attr_locus_tag = f.qualifiers[ "locus_tag" ][ 0 ]
                except:
                    attr_locus_tag = ""

                if strand == 1:
                    strand = '+'
                else:
                    strand = '-'

                try:
                    attr_product = f.qualifiers[ "product" ][0]
                except:
                    attr_product = ""

                inner = "gene=" +  attr_gene + ";locus_tag=" + attr_locus_tag + ";type=" + f.type + ";product=" + attr_product 
                out = [ genome_id, "from_genbank", "gene_ncbi", start, end, ".", strand, ".", inner ]
                fout.write( "\t".join( out ) + "\n" )

    handle.close()
    fout.close()

def gb2fasta(file_gb, file_output) :
  handle = open( file_gb, "r" )
  genomes = SeqIO.parse( handle, "gb" )    
  fout = open( file_output, "w" )

  for genome in genomes:  
    fout.write( ">%s\n%s\n"%( genome.name, genome.seq ) )
  
  handle.close()
  fout.close()

def count_coverage(samfile, chromosome_size=6000000, flip=False):
    """counts coverage per base in a strand-specific manner

    For paired-end reads, the insert betwen the mapped reads is
    also counted.

    flip: Whether or not the strands should be flipped.
    This should be true for RNA-seq, and false for ChIP-exo

    chromsome_size: This value should be larger than the largest chromosome"""
    if samfile.has_index():
        return count_coverage_indexed(samfile, chromosome_size=chromosome_size, flip=flip)
    all_counts = {}
    plus_strands = []
    minus_strands = []
    if "SQ" in samfile.header:
        chromosome_sizes = {}
        for entry in samfile.header["SQ"]:
            chromosome_sizes[entry["SN"]] = int(entry["LN"]) + 1
    else:
        for reference in samfile.references:
            chromosome_sizes[reference] = chromosome_size
    for reference in samfile.references:  # create an array for each reference
        plus_strands.append(zeros((chromosome_sizes[reference],)))
        minus_strands.append(zeros((chromosome_sizes[reference],)))
    # iterate through each mapped read
    for i, read in enumerate(samfile):
        if read.is_unmapped:
            continue
        if not read.is_proper_pair:
            if read.is_reverse:
                minus_strands[read.tid][read.pos:read.aend] += 1
            else:
                plus_strands[read.tid][read.pos:read.aend] += 1
        # for paired-end data, get entire insert from only read1
        elif read.is_read1:
            if read.is_reverse:
                minus_strands[read.tid][read.pnext:read.aend] += 1
            else:
                plus_strands[read.tid][read.pos:read.pos + read.isize] += 1
    # store the results per reference
    for i, reference in enumerate(samfile.references):
        all_counts[reference] = {}
        if flip:
            all_counts[reference]["-"] = plus_strands[i]
            all_counts[reference]["+"] = minus_strands[i]
        else:
            all_counts[reference]["+"] = plus_strands[i]
            all_counts[reference]["-"] = minus_strands[i]
    return all_counts

def count_coverage_indexed(samfile, chromosome_size=6000000, flip=False):
    """counts coverage per base in a strand-specific manner

    For paired-end reads, the insert betwen the mapped reads is
    also counted.

    flip: Whether or not the strands should be flipped.
    This should be true for RNA-seq, and false for ChIP-exo

    chromsome_size: This value should be larger than the largest chromosome"""
    if "SQ" in samfile.header:
        chromosome_sizes = {}
        for entry in samfile.header["SQ"]:
            chromosome_sizes[entry["SN"]] = int(entry["LN"]) + 1
    else:
        for reference in samfile.references:
            chromosome_sizes[reference] = chromosome_size
    all_counts = {}
    for reference in samfile.references:  # go through each chromosome
        plus_strand = zeros((chromosome_sizes[reference],))
        minus_strand = zeros((chromosome_sizes[reference],))
        # iterate through each mapped read
        for i, read in enumerate(samfile.fetch(reference=reference)):
            if not read.is_proper_pair:
                if read.is_reverse:
                    minus_strand[read.pos:read.aend] += 1
                else:
                    plus_strand[read.pos:read.aend] += 1
            # for paired-end data, get entire insert from only read1
            elif read.is_read1:
                if read.is_reverse:
                    minus_strand[read.pnext:read.aend] += 1
                else:
                    plus_strand[read.pos:read.pos + read.isize] += 1
            all_counts[reference] = {}
            if flip:
                all_counts[reference]["-"] = plus_strand
                all_counts[reference]["+"] = minus_strand
            else:
                all_counts[reference]["+"] = plus_strand
                all_counts[reference]["-"] = minus_strand
    return all_counts


def write_samfile_to_gff(samfile, output, chromosome_size=6000000,
                         separate_strand=True, flip=False, log_scale=False):
    """write samfile object to an output object in a gff format

    flip: Whether or not the strands should be flipped.
    This should be true for RNA-seq, and false for ChIP-exo

    chromsome_size: This value should be larger than the largest chromosome

    separate_strand: Whether the forward and reverse strands should be made
    into separate tracks (True) or the negative strand should be rendered
    as negative values (False)
    """
    all_counts = count_coverage(samfile, chromosome_size=chromosome_size,
                                flip=flip)

    name = os.path.split(samfile.filename)[1].decode()
    for reference in all_counts:
        for strand in all_counts[reference]:
            counts = all_counts[reference][strand]			
            for i in counts.nonzero()[0]:			
                if log_scale:
                    count = math.log(float(counts[i]), 2)
                else:
                  count = counts[i]
                if separate_strand:
                    output.write("%s\t\t%s\t%d\t%d\t%.2f\t%s\t.\t.\n" %
                        (reference, "%s_(%s)" % (name, strand), i, i,count, strand))
                else:
                    output.write("%s\t\t%s\t%d\t%d\t%.2f\t%s\t.\t.\n" %
                        (reference, name, i, i, count,strand))
						

def convert_samfile_to_gff(sam_filename, out_filename, chromosome_size=6000000,
                           separate_strand=True, flip=False, log_scale=False):
    """read in the a samfile from a path, and write it out to a gff filepath

    flip: Whether or not the strands should be flipped.
    This should be true for RNA-seq, and false for ChIP-exo.

    chromsome_size: This value should be larger than the largest chromosome.

    separate_strand: Whether the forward and reverse strands should be made
    into separate tracks (True) or the negative strand should be rendered
    as negative values (False)
    """
    samfile = pysam.Samfile(sam_filename)
    with open(out_filename, "w") as outfile:
        write_samfile_to_gff(samfile, outfile, chromosome_size=chromosome_size,
                             separate_strand=separate_strand, flip=flip, 
							 log_scale=log_scale)
    samfile.close()

"""# **2. Directory Setting**

If you want to run with the example file, download this file and upload that.

[Example download link](https://drive.google.com/drive/folders/1Kg8DZGWFcUPoR-n65QofTZdatk55jX54?usp=sharing)

---
**List of files**

>**1.   reference genome genbank file**
>> Upload the genbank file (format:.gb) downloaded from NCBI

>**2. sequencing data file**
>> Upload compressed fastq file (format:.gz)

>**3. datasheet.csv file**
>> Upload the csv file created according to the following format

---
About **datasheet.csv** file

>Each row should contain information about one sample.

>Enter files into a Tab-separeted file with 4columns:

>> **sample_id**: Make these easy to understand to you. The replicates should end with _1,_2, etc.

>> **R1, R2** : Sequencing file name (In case of sequencing by single-end, R2 is not required)

>> **organism** : reference genome ganbank file name

"""# **3. Running the pipeline**"""
system('mkdir 0_raw 1_reference 2_alignment 3_gff')
system('mv ./*.gb ./1_reference/')
system('mv ./*.gz ./0_raw/')

#Check the datasheet
datasheet = pd.read_csv('datasheet.csv')
datasheet

for i,row in datasheet.iterrows():
  print ( "%s processing..."%(row['sample_id']) )
  dir_gb = "./1_reference/%s.gb"%(row['organism'])
  dir_fasta = "./1_reference/%s.fasta"%(row['organism'])
  dir_ref = "./1_reference/%s"%(row['organism'])
  dir_refgff = "./3_gff/%s.gff"%(row['organism'])
  dir_gff_output = "./3_gff/%s.gff"%(row['sample_id'])

  #Unzip gz file
  dir_R1 = './0_raw/%s'%(row['R1'])
  system('gzip -d %s'%dir_R1)

  #build reference gff
  if os.path.isfile( dir_refgff ) :
    print ( 'Reference for bowtie already exist : %s'%(row['organism']))
  else:
    make_gff_from_ncbigb( dir_gb, dir_refgff )
    print ( 'Making reference gff : %s'%(row['organism']))

  #build reference file
  if os.path.isfile( dir_fasta ) :
    print ( 'Reference %s.gff already exist'%(row['organism']))
  else:
    gb2fasta( dir_gb, dir_fasta )
    system('bowtie-build %s %s'%(dir_fasta, dir_ref))
  
  #bowtie alignment
  bowtie_option = '-S'
  bowtie_input = "./0_raw/%s"%(row['R1'].split(".gz")[0])
  bowtie_output = "./2_alignment/%s.sam"%(row['sample_id'])
  bowtie_unaligned = bowtie_output + "_unaligned.fastq"
  
  print ("Mapping result:")
  system( 'bowtie %s %s %s --un %s > %s'%(bowtie_option, dir_ref, bowtie_input, bowtie_unaligned, bowtie_output))
  print ("Bowtie running is done")

  #samtools
  sam_bam_unsorted = bowtie_output + ".unsorted.bam"
  sam_bam_sorted = "./2_alignment/%s.bam"%(row['sample_id']) 


  system('samtools view -bS %s -o %s'%(bowtie_output, sam_bam_unsorted))

  system('samtools sort %s -o %s'%(sam_bam_unsorted, sam_bam_sorted ))
  system('rm %s'%sam_bam_unsorted)
  system('samtools index %s'%sam_bam_sorted)
  print ("Samtools running is done")
  
  convert_samfile_to_gff( bowtie_output, dir_gff_output )
  print ("Making the gff file is done")

  
  print ( "==================================================")

"""#**4. Visualization of gff files**

Visualize the files created in the '3_Gff' folder using metascope software

Download the metascope: [Metascope download link](https://sites.google.com/view/systemskimlab/software)

"""