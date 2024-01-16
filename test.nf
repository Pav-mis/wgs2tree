include { COMPLEASM; PARSE_TABLE; BUSCOMP; PARSE_GENES; MAFFT; IQTREE2; ASTRAL } from './processes.nf'

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

	run = PARSE_TABLE(COMPLEASM(for_busco)).mix(asm, busco_input).groupTuple()
	//def tuple = new Tuple('OG88', "/scratch/pawsey0812/pmisiun/busco2/run_OG88G_v230728.hic1.0.hifiasm.hap2.p_ctg", "/scratch/pawsey0812/pmisiun/assemblies/OG88/OG88G_v230728.hic1.0.hifiasm.hap2.p_ctg.fasta")
	//run = PARSE_TABLE(tuple)
	buscomp_runs = BUSCOMP(run).collect()
	genes = PARSE_GENES(buscomp_runs).flatten()
	msa = MAFFT(genes)
    trees = IQTREE2(msa).collectFile(name: 'loci.treefile')
	tree = ASTRAL(trees)
	tree.view()

}