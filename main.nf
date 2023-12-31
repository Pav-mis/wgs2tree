
include { BUSCO; BUSCOMP; PARSE_GENES; MAFFT; IQTREE2; ASTRAL; IQTREE2_COALESCENT_CONSENSUS } from './processes.nf'



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
	//run BUSCO on assemblies with missing runs, combine assembly and busco channels, group them by sample name for input into buscomp
	busco_runs = BUSCO(for_busco).mix(asm, busco_input).groupTuple()

	
	buscomp_runs = BUSCOMP(busco_runs).collect()
	genes = PARSE_GENES(buscomp_runs).flatten()
	msa = MAFFT(genes)
    trees = IQTREE2(msa).collectFile(name: 'loci.treefile')
	tree = ASTRAL(trees)
	tree.view()

}