

process BUSCO {
	cpus 32
	input:
	tuple val(sample), path(assembly)
	output:
	tuple val(sample), path("run_*")

	"""
	busco --in ${assembly} -l ${params.busco_lineage} --mode genome --cpu ${task.cpus * 2} -o run_${assembly.baseName}
	"""
}

process BUSCOMP {
	input:
	tuple val(sample), path(data)
	output:
	path "buscomp.faa"

	"""
	python /buscomp/code/buscomp.py runs="run_*" fastadir=${task.workDir} i=-1 buscompseq=F basefile=buscomp rmdreport=F basefile=${sample}
	"""
}

process PARSE_GENES {
	input:
	path "buscomp_runs"
	output:
	path "parsed_genes/*"

	"""
	python $projectDir/parse_genes.py ${buscomp_runs}
	"""
}

process MAFFT {
	cpus 1
	input:
	path "gene"
    output:
	path "MSA.fasta"

    """
    mafft --quiet ${gene} > MSA.fasta
    """
}

process IQTREE2 {
	cpus 1
	input:
	path "msa"
	output:
	path "msa.treefile"

    """
    iqtree2 -s ${msa}
    """
}

process ASTRAL {
	input:
	path "trees"
	output:
	path "species.treefile"

	"""
	astral -i ${trees} -o species.treefile
	"""
}

workflow {
	def assemblies = []
	params.samplePaths.each { key, value ->
		def files = file(value)
		files.each { fasta ->
		assemblies.add([key, fasta])
		}
	}

	genomes = Channel.from(assemblies)
	busco_runs = BUSCO(genomes).mix(genomes).groupTuple()
	buscomp_runs = BUSCOMP(busco_runs)

	// def buscomp = Channel.fromPath(params.buscomp_runs)
	// genes = PARSE_GENES(buscomp).flatten()
	// msa = MAFFT(genes)
    // trees = IQTREE2(msa).collectFile(name: 'loci.treefile')
	// tree = ASTRAL(trees)
	// tree.view()

}