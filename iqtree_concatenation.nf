include { BUSCO; BUSCOMP; PARSE_GENES; MAFFT; IQTREE2_CONCAT_CONSENSUS } from './processes.nf'


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
	

	for_busco = Channel.from(assemblies_for_busco)
	busco_input = Channel.from(buscos)
	asm = Channel.from(assemblies)
	busco_runs = BUSCO(for_busco).mix(asm, busco_input).groupTuple()

	
	buscomp_runs = BUSCOMP(busco_runs).collect()
	genes = PARSE_GENES(buscomp_runs).flatten()
	msa = MAFFT(genes).collect()
	tree = IQTREE2_CONCAT_CONSENSUS(msa)
	tree.view()
}