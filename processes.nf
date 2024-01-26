process DOWNLOAD_LINEAGES {
	output: 
    val true

    """
    compleasm download -L ${params.lineage_folder} ${params.busco_lineage}
    """
}

process BUSCO {
	publishDir params.busco_out
	errorStrategy 'retry'

	input:
	val ready
	tuple val(sample), path(assembly)

	output:
	tuple val(sample), path("run_*")

	"""
	busco --in ${assembly} -l ${params.busco_lineage} --mode genome --cpu ${task.cpus * 2} -o run_${assembly.baseName}
	"""
}


process COMPLEASM {
	publishDir params.compleasm_out
	errorStrategy 'retry'
	
	input:
	val ready
	tuple val(sample), path(assembly)

	output:
	tuple val(sample), path("run_*"), path(assembly)
	
	"""
	compleasm run -m 'busco' -a ${assembly} -l ${params.busco_lineage} -L ${params.lineage_folder} -t ${task.cpus * 2} -o run_${assembly.baseName}
	"""

}


process PARSE_TABLE {
	stageInMode "copy"

	input:
	tuple val(sample), path(run), path(assembly)
	output:
	tuple val(sample), path(run)

	"""
	# Change run directory to have run_ prefix
	mv ${run}/*/ ${run}/run_${params.busco_lineage}/

	# Create directories for busco sequences
	mkdir ${run}/run_${params.busco_lineage}/busco_sequences/
	mkdir ${run}/run_${params.busco_lineage}/busco_sequences/single_copy_busco_sequences/
	mkdir ${run}/run_${params.busco_lineage}/busco_sequences/multi_copy_busco_sequences/
	mkdir ${run}/run_${params.busco_lineage}/busco_sequences/fragmented_busco_sequences/
	
	grep -vE '^##STA' run_*/*/miniprot_output.gff > remove_PAF.gff
	awk -F'\t' '\$3 != "stop_codon"' remove_PAF.gff > remove_stops.gff

	python $projectDir/bin/parse_table.py ${run}/run_${params.busco_lineage}/ ${assembly}

	echo "# BUSCO version is: 5.4.7" | cat - ${run}/run_${params.busco_lineage}/full_table_busco_format.tsv > ${run}/run_${params.busco_lineage}/full_table.tsv
	"""
}

process BUSCOMP {
	publishDir params.buscomp_out
	input:
	tuple val(sample), path(data)

	output:
	path "${sample}.buscomp.faa"

	"""
	python /buscomp/code/buscomp.py runs="run_*" i=-1 basefile=${sample} buscompseq=F rmdreport=F 
	"""
}


process PARSE_GENES {
	publishDir params.parse_genes_out

	input:
	path buscomp_runs

	output:
	path "parsed_genes/*"

	"""
	python $projectDir/bin/parse_genes.py
	"""
}


process MAFFT {
	publishDir params.mafft_out

	input:
	path gene

    output:
	path "${gene}.fasta"

    """
    mafft --quiet ${gene} > ${gene}.fasta
    """
}


process IQTREE2 {
	publishDir params.iqtree_out

	input:
	path msa

	output:
	path "${msa}.treefile"

    """
    iqtree2 -s ${msa} -T ${task.cpus}
    """
}


process ASTRAL {
	publishDir params.astral_out
	cache false

	input:
	path trees

	output:
	path "species.treefile"

	"""
	astral -i ${trees} -o species.treefile
	"""
}


process IQTREE2_CONCAT_CONSENSUS {
	publishDir params.concat_out
	input:
	path msa

	output:
	path "species.treefile"

	"""
	mkdir Fasta
	mv *.fasta Fasta
	iqtree2 -p Fasta --prefix species -B 1000 -T 128
	"""
}

process GCF {
	publishDir params.final_out
	input:
	path supertree
	path locitree

	output:
	path "concord.cf.tree"

	"""
	iqtree2 -t ${supertree} --gcf ${locitree} --prefix concord
	"""
}