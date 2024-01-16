from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from distutils import filelist
import sys, os, re
import time

#argv[1] = the main run directory that contains full table
#argv[2] = fasta file of assembly

start = time.time()
working_dir = os.getcwd() + "/" + sys.argv[1]


table_file = working_dir+'full_table_busco_format.tsv'
full_table = open(table_file, 'r')
protien = working_dir + 'gene_marker.fasta'
assembly = SeqIO.index(sys.argv[2], 'fasta')

os.rename(table_file, working_dir+"full_table.tsv")

sequence_dict = {}

for record in SeqIO.parse(protien, "fasta"):
	gene = record.id.split('_')[0]
	sequence_dict[gene] = str(record.seq)


for line in full_table:
	if not line.startswith('#'):
		parts = line.split('\t')
		if parts[1].strip() != 'Missing':
			print(parts[1])
			record_id = parts[2] + ':' + parts[3] + '-' + parts[4]
			if parts[1] == 'Complete':
				#Write single copy protiens
				out_dir = 'single_copy_busco_sequences/'
				prot_record = SeqRecord(
    				Seq(sequence_dict[parts[0]]),
    				id=record_id,
					description="",
					name=""
				)
				prot_filename = working_dir + "/busco_sequences/" + out_dir + parts[0] + '.faa'
				with open(prot_filename, 'w') as output_file:
					SeqIO.write(prot_record, output_file, "fasta")

			elif parts[1] == 'Duplicated':
			 	out_dir = 'multi_copy_busco_sequences/'

			elif parts[1] == 'Fragmented':
				out_dir = 'fragmented_busco_sequences/'

			
			#Write fna sequences
			contig_record = assembly[parts[2]]
			contig_record.seq = contig_record.seq[int(parts[3])-1:int(parts[4])]
			if parts[5] == '-':
				contig_record.seq = contig_record.seq.reverse_complement()
			nt_record = SeqRecord(
    			contig_record.seq,
    			id=record_id,
				description="",
				name=""
			)
			nt_filename = working_dir + "busco_sequences/" + out_dir + parts[0] + '.fna'
			with open(nt_filename, 'a') as output_file:
				SeqIO.write(nt_record, output_file, "fasta")

print(start-time.time())