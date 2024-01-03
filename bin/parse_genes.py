#!/apps/python/3.8.3/bin/python3
from distutils import filelist
import sys, os, re
from Bio import SeqIO


# working_dir = /srv/scratch/z5206348/

working_dir = os.getcwd() + "/"
busco_genes = []

os.mkdir(working_dir + "parsed_genes/")



# store sample names
samples = filter(os.path.isfile, os.listdir( os.curdir ) )
print(samples)

for sample in samples:
    sample_dir = working_dir + sample
    with open(sample_dir, "r") as buscomp_output:
        for record in SeqIO.parse(buscomp_output, "fasta"):
            output_filename = working_dir + "parsed_genes/" + record.id + ".fasta"
            record.id = sample
            mode = 'a' if os.path.exists(output_filename) else 'w'
            with open(output_filename, mode) as output_file:
                SeqIO.write(record, output_file, "fasta")


#Only keep BUSCOs shared among all samples
parsed_genes = os.listdir(working_dir + "parsed_genes/")

for gene in parsed_genes:
    parsed_dir = working_dir + "parsed_genes/" + gene
    if sum(1 for _ in SeqIO.parse(parsed_dir, "fasta")) < 10:
        os.remove(parsed_dir)
        print("removing" + parsed_dir)
