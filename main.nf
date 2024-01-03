
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
	python /buscomp/code/buscomp.py runs="run_*" i=-1 buscompseq=F basefile=buscomp rmdreport=F basefile=${sample}
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
	input:
	path gene
    output:
	path "${gene}.fasta"

    """
    mafft --quiet ${gene} > ${gene}.fasta
    """
}

process IQTREE2 {
	input:
	path msa
	output:
	path "${msa}.treefile"

    """
    iqtree2 -s ${msa}
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

process ASTRAL {
	input:
	path trees
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
	def assemblies_for_busco = []
	params.samplePaths.each { key, value ->
		def files = files(value)
		files.each { fasta ->
			assemblies.add([key, fasta])
			if("run_$fasta.baseName" !in busco_base){
				assemblies_for_busco.add([key, fasta])
			}
		}
	}
	// fasta_base = busco.collect {it[1]}
	

	for_busco = Channel.from(assemblies_for_busco)
	busco_input = Channel.from(buscos)
	asm = Channel.from(assemblies)
	busco_runs = BUSCO(for_busco).mix(asm, busco_input).groupTuple()

	
	buscomp_runs = BUSCOMP(busco_runs).collect()
	genes = PARSE_GENES(buscomp_runs).flatten()
	msa = MAFFT(genes)
	// tree = IQTREE2_CONCAT_CONSENSUS(msa)
    trees = IQTREE2(msa).collectFile(name: 'loci.treefile')
	tree = IQTREE2_COALESCENT_CONSENSUS(trees)
	// tree = ASTRAL(trees)
	// tree.view()

}