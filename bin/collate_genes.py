from distutils import filelist
import sys, os, re
from Bio import SeqIO

working_dir = os.getcwd() + "/"
seq_path = sys.argv[1]
sample_id = sys.argv[2] + ".collate.faa"

print(working_dir)
print(seq_path)
print(sample_id)

genes = [s for s in os.listdir( seq_path ) if "faa" in s]


with open(sample_id, "a") as collated_file:
    for gene in genes:
        gene_path = seq_path + '/' + gene
        gene_id = gene.split('.')[0]
        with open(gene_path) as sequence:
            for record in SeqIO.parse(sequence, "fasta"):
                record.id = gene_id
                record.description = ">" + record.description
                print(record.description)
                SeqIO.write(record, collated_file, "fasta")