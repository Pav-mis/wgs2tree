// include modules
include { DOWNLOAD_LINEAGES_COMPLEASM; DOWNLOAD_LINEAGES_BUSCO; COMPLEASM; BUSCO; PARSE_TABLE; BUSCOMP; COLLATE_GENES; PARSE_GENES; MAFFT; IQTREE2; IQTREE2_CONCAT_CONSENSUS; ASTRAL; GCF } from './processes.nf'


log.info """
###============###
    wgs-2-tree         
###============###

> inPath:           ${params.inPath}
> outPath:          ${params.outPath}

> withBusco:        ${params.withBusco}         
> withConcat:       ${params.withConcat} 
> geneThres:        ${params.geneThres}
> lineage:          ${params.lineage}

> buscoDir:         ${params.buscoDir}
> compleasmDir:     ${params.compleasmDir}
> updateLineage:    ${params.updateLineage}

"""


workflow {

    // Determine which assemblies need busco/compleasm, based on what inputs are given
    def assemblies = []
    def assemblies_for_busco = []
    def buscos = []
    def busco_bases = []

    file(params.inPath).eachDir {it ->
        def key = it.getBaseName()
        if( key.contains("run_")) {
            busco_bases.add("$it.Name")
            buscos.add([key.replaceAll("run_", ""), it])
        }
        else {
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
    }

    file(params.inPath).eachFile {it ->
        if( it.isFile() ) {
            assemblies.add([it.getBaseName(), it])
            if("run_$it.baseName" !in busco_bases) {
                assemblies_for_busco.add([it.getBaseName(),it])
            }
        }
    }
    

	//Channel of assemblies for BUSCO
	for_busco = Channel.from(assemblies_for_busco)
	//Channel of given BUSCOs
	busco_input = Channel.from(buscos)
	//Channel of given assemblies
	asm = Channel.from(assemblies)

	
    //Get single copy orthologues with either busco or compleasm
	if (params.withBusco == true) {
	 	DOWNLOAD_LINEAGES_BUSCO()
	 	run = BUSCO(DOWNLOAD_LINEAGES_BUSCO.out, for_busco).mix(asm, busco_input).groupTuple()
	}
	else {
	 	DOWNLOAD_LINEAGES_COMPLEASM()
	 	run = PARSE_TABLE(COMPLEASM(DOWNLOAD_LINEAGES_COMPLEASM.out, for_busco)).mix(asm, busco_input).groupTuple()
	}


    //Split sample tuples into those that consist of a single assemblies and those with more than one assembly
    run.branch {
        single: it[1].size() == 2
        multi: it[1].size() > 2
    }.set{runs}


    singles = COLLATE_GENES(runs.single) //Collate genes into the needed format for single assembly samples 
    multis = BUSCOMP(runs.multi) //Get best non-redundant set of buscos from samples with multiple assemblies/runs
    
	best_genes = singles.mix(multis).collect()

	genes = PARSE_GENES(best_genes).flatten() //Group samples by gene
	msa = MAFFT(genes) //Align each gene
    trees = IQTREE2(msa).collectFile(name: 'loci.treefile') //Generate tree for each alignment
	if (params.withConcat == true) {
	    tree = IQTREE2_CONCAT_CONSENSUS(msa.collect()) //Generate supertree with concatenation if specified
    }
	else {
	    tree = ASTRAL(trees) //Generate supertree with coalescence 
	}
    final_tree = GCF(tree, trees) //Calculate GCF values
	final_tree.view()

}