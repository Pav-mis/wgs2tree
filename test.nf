include { DOWNLOAD_LINEAGES; COMPLEASM; BUSCO; PARSE_TABLE; BUSCOMP; PARSE_GENES; MAFFT; IQTREE2; IQTREE2_CONCAT_CONSENSUS; ASTRAL; GCF } from './processes.nf'

workflow {

    def assemblies = []
    def assemblies_for_busco = []
    def buscos = []

    file(params.inPath).eachDir {it ->
        def key = it.getBaseName()
        def busco_bases = []
        it.eachFile {busco_run ->
            if( busco_run.isDirectory()) {
                busco_bases.add("$busco_run.Name")
                buscos.add([key, busco_run])
            }
        }
        it.eachFile {assembly ->
            if( assembly.isFile() ) {
                assemblies.add([key, assembly])
                if("run_$assembly.baseName" !in busco_bases) {
                    assemblies_for_busco.add([key,assembly])
                    }
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