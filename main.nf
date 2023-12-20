
process BUSCO {
	errorStrategy 'ignore'
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
	path "${sample}.buscomp.faa"

	"""
	python /buscomp/code/buscomp.py runs="run_*" i=-1 buscompseq=F basefile=buscomp rmdreport=F basefile=${sample}
	"""
}

process PARSE_GENES {
	input:
	path buscomp_runs
	output:
	path "parsed_genes/*"

	"""
	python $projectDir/bin/parse_genes.py buscomp
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
	def buscos = []
	params.sampleBuscoPaths.each { key, value ->
		def runs = files(value, type:'dir')
		runs.each {run ->
		buscos.add([key, run])
		}
	}

	def busco_base = buscos.collect {sample, fasta -> return "$fasta.name"}
	println(busco_base)

	def assemblies = []
	params.samplePaths.each { key, value ->
		def files = files(value)
		files.each { fasta ->
			if("run_$fasta.baseName" !in busco_base){
				assemblies.add([key, fasta])
			}
		}
	}
	// fasta_base = busco.collect {it[1]}
	

	genomes = Channel.from(assemblies)
	busco_input = Channel.from(buscos)
	busco_runs = BUSCO(genomes).mix(genomes, busco_input).groupTuple()

	
	buscomp_runs = BUSCOMP(busco_runs).collect()
	genes = PARSE_GENES(buscomp_runs).flatten()
	msa = MAFFT(genes)
    trees = IQTREE2(msa).collectFile(name: 'loci.treefile')
	tree = ASTRAL(trees)
	tree.view()

}