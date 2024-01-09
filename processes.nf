process BUSCO {
	publishDir params.busco_out
	errorStrategy 'ignore'

	input:
	tuple val(sample), path(assembly)

	output:
	tuple val(sample), path("run_*")

	"""
	busco --in ${assembly} -l ${params.busco_lineage} --mode genome --cpu ${task.cpus * 2} -o run_${assembly.baseName}
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
	python $projectDir/bin/parse_genes.py buscomp
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
	input:
	path msa

	output:
	path "concat.treefile"

	"""
	mkdir Fasta
	mv *.fasta Fasta
	iqtree2 -p Fasta --prefix concat -B 1000 -T 128
	"""
}


process IQTREE2_COALESCENT_CONSENSUS {
	input:
	path trees

	output:
	path "loci.treefile.contree"

	"""
	iqtree2 -t ${trees} -con -T 128
	"""
}