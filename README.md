# wgs2tree
This is a nextflow implementation of a workflow for the rapid generation of consensus phylogenomic trees from whole genome sequence data.

The workflow is as follows:

1. Samples, either consisting of a single assembly, or multiple assemblies (ie, hap1 and hap2) are used as input
2. For each assembly, single-copy orthologues are identified with BUSCO
3. For samples consisting of more than one assembly, the maximal set of Busco orthologues is retrieved with Buscomp
4. For each retrieved BUSCO gene, a multiple sequence alignment is performed with MAFFT
5. A phylogenetic tree is generated for each gene with IQTree2
6. A consensus phylogenomic tree is generated with ASTRAL

#Config
The nextflow config file should be modified to specify relevant parameters, using the provided 'template.config'. Within 'params', various specifications can be made. 

1. samplePaths:
This parameter is a groovy map, used to specify assembly inputs. The keys are sample identifiers (used as labels on leaves on the final phylogenomic tree), and the mapped value is a file path or a list of paths (glob operators allowed).

2. sampleBuscoPaths:
This parameter is a groovy map, used to specify busco inputs. When running the workflow, provided busco run directories will be paired to a respective assembly in samplePaths. If an assembly does not have a corresponding busco run as input, busco will be executed on the assembly to attain it. In order for pairing of assembly and busco inputs, consistent naming is needed (see below).

	Assembly = /path/to/assembly/sample1.fasta
	Busco = /path/to/busco/run_sample1/
