from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from distutils import filelist
import sys, os, re
import time
import subprocess




#argv[1] = the main run directory that contains full table
#argv[2] = fasta file of assembly


work_dir = os.getcwd()
run_dir = work_dir + "/" + sys.argv[1]
genome = sys.argv[2]

protien = run_dir + 'gene_marker.fasta'

full_table = open(run_dir+'full_table.tsv', 'r')

sequence_dict = {}

### Index protiens
for record in SeqIO.parse(protien, "fasta"):
	sequence_dict[record.id] = record.seq

print(len(sequence_dict))

best_genes = []

### Get best genes from full table
for line in full_table:
	if not line.startswith('#'):
		parts = line.split('\t')
		if parts[1].strip() == 'Single':
			best = parts[11]+':' + parts[2] + ':' + parts[3] + '-' + parts[4]
			best_genes.append(best)



### get corresponding ID's of best genes from gff file
### write these to a file and also map them to the best genes identified above

with open('remove_stops.gff', 'r') as file:
    lines = file.readlines()

idmatch = {}
with open('best_genes.txt', 'w') as output_file:
	for i, line in enumerate(lines):
		if line.startswith('##PAF'):
			parts = line.split('\t')
			gene =  parts[1] + ':' + parts[6] + ':' + parts[8] + '-' + parts[9]
			if gene in best_genes:
				ID = re.search(r'ID=([^;\s]+)', lines[i+1])[1]
				idmatch[ID] = gene
				output_file.write(ID+'\n')


### Remove PAF lines from gff file (gffread doesn't like them)
### Run gffread, with -x to get spliced CDS, --ids to only extract the best genes

f = open("new.gff", "w")
subprocess.call(['grep', '-vE', '^##PAF', 'remove_stops.gff'], stdout=f)
subprocess.run(['gffread', '-x', 'seq.fa', '--ids', 'best_genes.txt', '-g', genome, 'new.gff'])


### Write directory of single copy fna sequences with correct IDs

with open('seq.fa', "r") as single_copy_fnas:
	for record in SeqIO.parse(single_copy_fnas, "fasta"):
		name = idmatch[record.id].split(':')
		filename = name[0]
		if '_' in filename:
			filename = filename.split('_')[0]
		#Write fnas
		new_record = SeqRecord(
    		record.seq,
			id=name[1] + ':' + name[2],
			description="",
			name=""
		)
		outfile = run_dir + 'busco_sequences/single_copy_busco_sequences/' + filename + '.fna'
		with open(outfile, 'w') as output_file:
			SeqIO.write(new_record, output_file, "fasta")

		#Write proteins
		prot_seq = sequence_dict[name[0]]
		prot_record = SeqRecord(
    		prot_seq,
			id=name[1] + ':' + name[2],
			description="",
			name=""
		)
		outfile_prot = run_dir + 'busco_sequences/single_copy_busco_sequences/' + filename + '.faa'
		with open(outfile_prot, 'w') as output_file_prot:
			SeqIO.write(prot_record, output_file_prot, "fasta")