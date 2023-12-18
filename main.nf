params.assemblies = "/scratch/pawsey0812/pmisiun/assemblies"
params.busco_lineage = "vertebrata_odb10"

process BUSCO {
	cpus 32
	input:
	path assembly
	output:
	path "run_*"

	"""
	mkdir run_${assembly.baseName}
	"""
}

process BUSCOMP {
	input:
	tuple path "assembles"
	output
	path "buscomp.faa"

	"""
	python /buscomp/code/buscomp.py runs=$busco_runs fastadir=$assemblies i=-1 buscompseq=F basefile=buscomp rmdreport=F
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
	genomes = Channel.fromPath("${params.assemblies}/**")
	println params.samplePaths
	sampleChannels = Channel.fromPaths(params.samplePaths)
	
	/*busco_runs = BUSCO(genomes)
	def buscomp = Channel.fromPath(params.buscomp_runs)
	*genes = PARSE_GENES(buscomp).flatten()
	*msa = MAFFT(genes)
    trees = IQTREE2(msa).collectFile(name: 'loci.treefile')
	tree = ASTRAL(trees)
	tree.view()
	*/
}