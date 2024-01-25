include { DOWNLOAD_LINEAGES; COMPLEASM; BUSCO; PARSE_TABLE; BUSCOMP; PARSE_GENES; MAFFT; IQTREE2; IQTREE2_CONCAT_CONSENSUS; ASTRAL; GCF } from './processes.nf'

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

	//determines which assemblies need BUSCO run on them (BUSCO runs have not been given)
	params.samplePaths.each { key, value ->
		def files = files(value)
		files.each { fasta ->
			assemblies.add([key, fasta])
			if("run_$fasta.baseName" !in busco_base){
				assemblies_for_busco.add([key, fasta])
			}
		}
	}


	//Channel of assemblies for BUSCO
	for_busco = Channel.from(assemblies_for_busco)
	//Channel of given BUSCOs
	busco_input = Channel.from(buscos)
	//Channel of given assemblies
	asm = Channel.from(assemblies)

	
	DOWNLOAD_LINEAGES()

	if (params.withBusco == true) {
		run = BUSCO(DOWNLOAD_LINEAGES.out, for_busco).mix(asm, busco_input).groupTuple()
	}
	else {
		run = PARSE_TABLE(COMPLEASM(DOWNLOAD_LINEAGES.out, for_busco)).mix(asm, busco_input).groupTuple()
	}

	buscomp_runs = BUSCOMP(run).collect()
	genes = PARSE_GENES(buscomp_runs).flatten()
	msa = MAFFT(genes)
    trees = IQTREE2(msa).collectFile(name: 'loci.treefile')
	if (params.withConcat == true) {
		tree = IQTREE2_CONCAT_CONSENSUS(msa.collect())
	}
	else {
		tree = ASTRAL(trees)
	}
	final_tree = GCF(tree, trees)
	final_tree.view()

}