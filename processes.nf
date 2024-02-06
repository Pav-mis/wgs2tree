process DOWNLOAD_LINEAGES_COMPLEASM {
	cache false

	output: 
    val true

    """
	if [ ! -d "${params.compleasmDir}/${params.lineage}" ] || [ "${params.updateLineage}" = "true" ]; then
    	compleasm download -L ${params.compleasmDir} ${params.lineage}
	fi
    """
}

process DOWNLOAD_LINEAGES_BUSCO {
	cache false

	output: 
    val true

    """
	if [ ! -d "${params.buscoDir}/lineages/${params.lineage}" ] || [ "${params.updateLineage}" = "true" ]; then
		mkdir -p ${params.buscoDir}
    	busco --download ${params.lineage}
		cp -r busco_downloads/* ${params.buscoDir}
	fi
    """
}

process BUSCO {
	publishDir "${params.outPath}/nf_busco", mode: 'copy', overwrite: true

	input:
	val ready
	tuple val(sample), path(assembly)

	output:
	tuple val(sample), path("run_*")

	"""
	busco --in ${assembly} -l ${params.buscoDir}/lineages/${params.lineage} --mode genome --cpu ${task.cpus * 2} -o run_${assembly.baseName}
	"""
}


process COMPLEASM {
	input:
	val ready
	tuple val(sample), path(assembly)

	output:
	tuple val(sample), path("run_*"), path(assembly)
	
	"""
	compleasm run -m 'busco' -a ${assembly} -l ${params.lineage} -L ${params.compleasmDir} -t ${task.cpus * 2} -o run_${assembly.baseName}
	"""

}


process PARSE_TABLE {
	publishDir "${params.outPath}/nf_compleasm", mode: 'copy', overwrite: true
	stageInMode "copy"

	input:
	tuple val(sample), path(run), path(assembly)
	output:
	tuple val(sample), path(run)

	"""
	# Change run directory to have run_ prefix
	mv ${run}/*/ ${run}/run_${params.lineage}/

	# Create directories for busco sequences
	mkdir ${run}/run_${params.lineage}/busco_sequences/
	mkdir ${run}/run_${params.lineage}/busco_sequences/single_copy_busco_sequences/
	mkdir ${run}/run_${params.lineage}/busco_sequences/multi_copy_busco_sequences/
	mkdir ${run}/run_${params.lineage}/busco_sequences/fragmented_busco_sequences/
	
	grep -vE '^##STA' run_*/*/miniprot_output.gff > remove_PAF.gff
	awk -F'\t' '\$3 != "stop_codon"' remove_PAF.gff > remove_stops.gff

	python $projectDir/bin/parse_table.py ${run}/run_${params.lineage}/ ${assembly}

	echo "# BUSCO version is: 5.4.7" | cat - ${run}/run_${params.lineage}/full_table_busco_format.tsv > ${run}/run_${params.lineage}/full_table.tsv
	"""
}

process BUSCOMP {
	publishDir "${params.outPath}/nf_buscomp", mode: 'copy', overwrite: true
	input:
	tuple val(sample), path(data)

	output:
	path "${sample}.buscomp.faa"

	"""
	python /buscomp/code/buscomp.py runs="run_*" i=-1 basefile=${sample} buscompseq=F rmdreport=F 
	"""
}

process COLLATE_GENES {
	publishDir "${params.outPath}/nf_buscomp", mode: 'copy', overwrite: true
	input:
	tuple val(sample), path(run)

	output:
	path "${sample}.collate.faa"

	"""
	python $projectDir/bin/collate_genes.py run_*/run_${params.lineage}/busco_sequences/single_copy_busco_sequences/ ${sample}
	"""

}


process PARSE_GENES {
	publishDir "${params.outPath}/nf_parsed_genes", mode: 'copy', overwrite: true

	input:
	path buscomp_runs

	output:
	path "parsed_genes/*"

	"""
	python $projectDir/bin/parse_genes.py ${params.geneThres}
	"""
}


process MAFFT {
	publishDir "${params.outPath}/nf_mafft", mode: 'copy', overwrite: true

	input:
	path gene

    output:
	path "${gene}.fasta"

    """
    mafft --quiet ${gene} > ${gene}.fasta
    """
}


process IQTREE2 {
	publishDir "${params.outPath}/nf_iqtree", mode: 'copy', overwrite: true

	input:
	path msa

	output:
	path "${msa}.treefile"

    """
    iqtree2 -s ${msa} -T ${task.cpus * 2}
    """
}


process ASTRAL {
	publishDir "${params.outPath}/nf_astral", mode: 'copy', overwrite: true
	cache false

	input:
	path trees

	output:
	path "species.treefile"

	"""
	astral -i ${trees} -o species.treefile
	"""
}


process IQTREE2_CONCAT_CONSENSUS {
	publishDir "${params.outPath}/nf_concat", mode: 'copy', overwrite: true
	input:
	path msa

	output:
	path "species.treefile"

	"""
	mkdir Fasta
	mv *.fasta Fasta
	iqtree2 -p Fasta --prefix species -B 1000 -T ${task.cpus * 2}
	"""
}

process GCF {
	publishDir params.outPath, mode: 'copy', overwrite: true
	input:
	path supertree
	path locitree

	output:
	path "concord.cf.tree"

	"""
	iqtree2 -t ${supertree} --gcf ${locitree} --prefix concord
	"""
}